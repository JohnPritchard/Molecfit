#!/usr/bin/env python3

import os
import shutil
import glob
import sys
import tarfile
import math
import time
import numpy as np
import threading
import dill
#from multiprocessing import multiprocessing.Process
import multiprocessing
from datetime import datetime

import logging
logger=logging.getLogger()
logger.handlers = []
# see https://stackoverflow.com/questions/11581794/how-do-i-change-the-format-of-a-python-log-message-on-a-per-logger-basis
logFMT = "%(asctime)s %(module)20s[%(process)5d] [%(levelname)s] %(message)s"
ch = logging.StreamHandler()
formatter = logging.Formatter(logFMT, "%Y-%m-%dT%H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)


# GLOBAL VARIABLES:
ODA_OPTION = os.getenv('ODA_OPTION')
if (ODA_OPTION == None):
    ODA_OPTION="STD"
LBLRTM_BIN      = 'lblrtm'
LNFL_BIN        ="lnfl"
MOLPRFSDIR      = 'INDIVIDUAL_MOLECULE_PROFS'
TAPE5_FILENAME  = 'TAPE5'
TAPE1_FILENAME  = 'TAPE1'
TAPE3_FILENAME  = 'TAPE3'
TAPE6_FILENAME  = 'TAPE6'
TAPE7_FILENAME  = 'TAPE7'
TAPE10_FILENAME = 'TAPE10'
TAPE28_FILENAME = 'TAPE28'
TAPE28_GLOBSTR  = 'TAPE28_*'
TAPE3_DIRNAME   = 'TAPE3_DIR'
TAPE28STD_FILENAME='TAPE28STD'
TAPE28ODA_FILENAME='TAPE28Oda'
OPTICALDEPTHS_BINFILE='OPTICAL_DEPTHS.npy'
OPTICALWVNUMS_BINFILE='OPTICAL_WVNUMS.npy'
ALLMOLS_DIRNAME = "ALL"
N_MOL_TYPES=47

MOL_STR="\
( 1)  H2O  ( 2)  CO2  ( 3)    O3 ( 4)   N2O ( 5)    CO ( 6)   CH4 ( 7)    O2 \
( 8)   NO  ( 9)  SO2  (10)   NO2 (11)   NH3 (12)  HNO3 (13)    OH (14)    HF \
(15)  HCL  (16)  HBR  (17)    HI (18)   CLO (19)   OCS (20)  H2CO (21)  HOCL \
(22)   N2  (23)  HCN  (24) CH3CL (25)  H2O2 (26)  C2H2 (27)  C2H6 (28)   PH3 \
(29) COF2  (30)  SF6  (31)   H2S (32) HCOOH (33)   HO2 (34)     O (35)CLONO2 \
(36)  NO+  (37) HOBR  (38)  C2H4 (39) C3HOH (40) CH3Br (41) CH3CN (42)   CF4 \
(43) C4H2  (44) HC3N  (45)    H2 (46)    CS (47)   SO3                       \
"
# Parse the MOLECULAR set string into a workable list
def ParseMOLSTR(MOL_STR):

    ret_lst=[]
    for i in range(1,N_MOL_TYPES+1):
        if (i<10):
            substr="( "+str(i)+")"
        else:
            substr="("+str(i)+")"
        idx =MOL_STR.find(substr)
        mols=MOL_STR[idx+4:idx+11].strip()
        ret_lst.append(mols)
    # Add an all moles dirname to the list
    ret_lst.append(ALLMOLS_DIRNAME)
    return ret_lst
MOLECULES_FULL_LST=ParseMOLSTR(MOL_STR)
# Store the ALLMLS idx (which will be at the end of the list)
ALLMOLS_IDX=len(MOLECULES_FULL_LST)-1

def Mol2Idx(mol_str):
    idx=-1
    n_molecules=len(MOLECULES_FULL_LST)
    for i in range(n_molecules):
        if(MOLECULES_FULL_LST[i]==mol_str):
            idx=i+1
            return idx
    return idx

# =========================================================
#                    DEBUG UTIL
# =========================================================

def dump4plot(fname,xV,yV,n):
    fid=open(fname,'w')
    for i in range(n):
        #lst = []
        #lst.append(xV[i])
        #lst.append(yV[i])
        #xstr='%.12g' % xV[i]
        #ystr='%.12g' % yV[i]
        #fid.write(xstr + ' ' + ystr + ' ' + str(yV[i]) + '\n')
        fid.write(str(xV[i])+' '+str(yV[i])+'\n')

    fid.close()

# =========================================================
#                    READ TAPE 28
# =========================================================

def GetTAPE28Contents(dirname):

    # --------
    # COMMENTS:
    # --------
    # Read all lines in ascis file TAPE28 in dir dirname and return as list
    # Have to add additional smartness as required TAPE28 may have been moved
    # by molecfit to an upper directory or by ODA scripts to TAPE28STD

    # Assume that TAPE28 exists in the working directory as TAPE28
    # or TAPE28STD
    filename    =os.path.join(dirname,TAPE28_FILENAME)
    std_filename=os.path.join(dirname,TAPE28STD_FILENAME)

    # Check if either are present if not assume it has been moved to ../TAPE28_?
    # and make a soft link
    if (not os.path.exists(filename) and not os.path.exists(std_filename) ):
        upperdir=os.path.dirname(os.path.realpath(dirname))
        globstr=os.path.join(upperdir,TAPE28_GLOBSTR)
        lst=glob.glob(globstr)
        target=lst[0]
        os.symlink(target,filename)
    elif (os.path.exists(std_filename)):
        # There is no TAPE28 in the working directory but there is a TAPE28STD
        # so we use this
        filename=std_filename

    # Having got the rquired TAPE28 filename open and read
    fid=open(filename,'r')
    lines=fid.readlines()
    fid.close()

    # Return the file lines as a list
    return lines

def GetTAPE28Header(dirname):

    # Return the Header information of a TAPE28 file

    # Get the lines contained in this file
    lines=GetTAPE28Contents(dirname)

    # Iterate through lines up to the "WAVENUMBER" "TANSMISSION" line
    ret_lst=[]
    for line in lines:
        ret_lst.append(line)
        # If this line contains "WAVENUMBER" and "TRANSMISSION" then start extracting
        # on the next iteration
        if (line.find("WAVENUMBER")!=-1 and line.find("TRANSMISSION")!=-1):
                break
    return ret_lst

def ReadTAPE28(dirname):

    # Get the lines contained in this file
    lines=GetTAPE28Contents(dirname)

    # Iterate through all lines and extract Wavenumbers and Transmission values
    cnt=0;            # No of lines with values extracted
    extracting=False; # Boolean if in extracting mode
    wave_lst=[];      # The wavenumber extracted values list
    trans_lst=[]       # The transmission extracted values list
    for line in lines:

        # Skip empty lines
        if (len(line.strip())==0):
            continue

        # If this line contains "WAVENUMBER" and "TRANSMISSION" then start extracting
        # on the next iteration
        if (line.find("WAVENUMBER")!=-1 and line.find("TRANSMISSION")!=-1):
                extracting=True
                continue

        # Skip if we are not extracting
        if (not extracting):
            continue

        # Real Values are to be extracted from this line
        cnt+=1
        lst=line.split()
        wave_lst.append (lst[0])
        trans_lst.append(lst[1])

    # Convert to numpy vectors
    waveV=np.zeros(cnt)
    transV=np.zeros(cnt)
    for i in range (cnt):
        waveV [i] = wave_lst[i]
        transV[i] =trans_lst[i]

    return waveV, transV, cnt

def GenCombTrans(mol_lst,r_lst,npts):

    n=len(r_lst)
    print("Combining ", n, " molecule transitivity profiles")
    tV=np.zeros(npts)
    for idx in range(n):
        sV=mol_lst[idx]
        r =r_lst[idx]
        #r=1.0/r
        for j in range(npts):
            trans=max(sV[j],1.0e-08)
            tau=-1.0*math.log(trans)
            tV[j]=tV[j]+r*tau

    for j in range(npts):
        tau=tV[j]
        tV[j]=math.exp(-tau)

    return tV


# =========================================================
#         END READ TAPE 28
# =========================================================

# =========================================================
#         READ TAPE 5
# =========================================================

MAX_NMOLS=50

def READ_TAPE5(dirname):

    verbose=False
    filename=os.path.join(dirname,'TAPE5')

    # Read TAPE5 file and dump int a data list
    fid=open(filename,'r')
    data=fid.readlines()
    fid.close()


    # Iterate through each line and parse
    cnt=0

    max_nmols=0
    record_section=False
    new_record    =False
    for line in data:

        cnt=cnt+1
        if (cnt<10 and verbose) :
            print(cnt, ") ", line)

        if (record_section and new_record):

            # Parse the line as the header record in the form "val val val AA AAAA??"
            lst=line.split()
            hgth=lst[0]
            pres=lst[1]
            temp=lst[2]
            mols=lst[4]         # A string of 'A' characters one for each molecules
            nmols=len(mols)     # No of molecules = no of 'A' characters
            max_nmols=max(max_nmols,nmols)
            nmol_lines=(nmols-1)//8+1 # Record is broken up into groups of 8 molecules
            hvec[layer]=hgth
            pvec[layer]=pres
            tvec[layer]=temp
            new_record=False
            record_part=0
            if (verbose):
                print (layer+1, hgth, pres, temp, mols,nmol_lines)

        elif (record_section and layer<nlayers) :
            record_part=record_part+1
            # Parse the line as a list of abundencies
            lst=line.split()
            for i in range(len(lst)):
                j = 8*(record_part-1) + i
                abun[layer,j]=lst[i]
            if (record_part==nmol_lines):
                new_record=True
                layer=layer+1
            if (layer==nlayers):
                new_record=False
                record_section=False

        if (cnt==8) :
            # This is a datum line containing the no of layers
            nlayers=int(line)
            # Now allocate the data vactors and arrays
            abun=np.zeros((nlayers,MAX_NMOLS))
            hvec=np.zeros((nlayers))
            tvec=np.zeros((nlayers))
            pvec=np.zeros((nlayers))

            # Decalre we are now in the record section and we are reading a new record
            new_record    =True
            record_section=True
            layer=0     # Start the layer index at 0

    if (verbose) :
        print('N lines=',cnt)
        print('n_layers=',nlayers)

    return hvec,pvec,tvec,abun,nlayers,max_nmols

# end def READ_TAPE5(dirname)



def CompareTAPE5(dirname1,dirname2):

    hvec1,pvec1,tvec1,abunA1,nlayers1,nmols1=READ_TAPE5(dirname1)
    hvec2,pvec2,tvec2,abunA2,nlayers2,nmols2=READ_TAPE5(dirname2)
    different_nlayers=False
    different_nmols  =False
    if (nlayers1!=nlayers2):
        different_nlayers=True
    if (nmols1!=nmols2):
        different_nmols  =True
    print("Check layers: ",nlayers1,nlayers2, "Diff=", different_nlayers)
    print("Check nmols: ", nmols1,nmols2, "Diff=", different_nmols)
    if (different_nlayers or different_nmols):
        return
    nlayers=nlayers1
    nmols=nmols1
    sA=np.zeros((nlayers,nmols))
    for layer in range(nlayers):
        lst=[]
        for mol_idx in range(nmols):
            a1=abunA1[layer][mol_idx]
            a2=abunA2[layer][mol_idx]
            if (a2==0.0):
                scale=0.0
            else:
                scale=a1/a2
            lst.append(scale)
            sA[layer][mol_idx]=scale
    lst=[]
    avgV=np.median(sA,axis=0)
    for mol_idx in range(nmols):
        avg=avgV[mol_idx]
        lst.append(avg)
    return lst


# =========================================================
#         END READ TAPE 5
# =========================================================

# =========================================================
#       TAPE3_DIR PATH and Molecule Profile Data Path
# =========================================================

def GetTAPE3_DIRPath(dirname):
    tape3=os.path.join(dirname,TAPE3_FILENAME)
    realpath=os.path.realpath(tape3)
    realdir=os.path.dirname(realpath)
    tape3dir=os.path.join(realdir,TAPE3_DIRNAME)
    return tape3dir

def GetMolPrfDir(dirname):
    tape3_dir=GetTAPE3_DIRPath(dirname)
    hdir=os.path.dirname(tape3_dir)
    mprfdir=os.path.join(hdir,MOLPRFSDIR)
    return mprfdir

def MolPrfDirExists(dirname):
    path=GetMolPrfDir(dirname)
    if os.path.exists(path):
        return True
    return False

def MakeMOLPRFSDir(dirname,tape5,tape28,tape3dir):
    molprfs_dir=os.path.join(dirname,MOLPRFSDIR)
    print ("From dirname=",dirname, "-> molprfs=",molprfs_dir)
    if (os.path.exists(molprfs_dir)):
        print("Molecule Profile Directory allready exists. Will delete.")
        shutil.rmtree(molprfs_dir)
    os.mkdir(molprfs_dir)
    target=os.path.join(molprfs_dir,TAPE5_FILENAME)
    shutil.copy(tape5,target)
    target=os.path.join(molprfs_dir,TAPE28_FILENAME)
    shutil.copy(tape28,target)
    target=os.path.join(molprfs_dir,TAPE3_DIRNAME)
    os.symlink(tape3dir,target)
    return molprfs_dir

# =========================================================
#    END TAPE3_DIR PATH and Molecule Profile Data Path
# =========================================================

# ----------------------------------------------------------------------------
def run_dill_encoded(payload):
    i, fun, args = dill.loads(payload)
    r=fun( args )
    return i, r
# ----------------------------------------------------------------------------
def apply_async(pool, i, fun, args):
    payload = dill.dumps((i, fun, args))
    return pool.apply_async(run_dill_encoded, (payload,))
# ----------------------------------------------------------------------------
def doit_SysCall( mod, cmds ) :
    sdt=datetime.now()
    the_pool = multiprocessing.Pool(int(os.getenv('MF_LNFL_NUM_THREADS',8)))
    jobs = []
    for i,cmd in enumerate(cmds) :
        #the_queue.put([i,[mod, a1, a2])
        job = apply_async(the_pool, i, SysCall, (cmd,))
        jobs.append(job)

    retvals={}
    for job in jobs:
        retval = job.get()
        print('\t', retval[0], retval[1])
        i=retval[0]
        r=retval[1]
        retvals[i]=r

    the_pool.terminate()
    ddt=datetime.now()-sdt
    t_elapsed_time=datetime.now()
    t_elapsed_time=t_elapsed_time-t_elapsed_time

    logger.info(f'Total equivalent elapsed time [{str(t_elapsed_time)}]')
    logger.info(f'          Actual elapsed time [{str(ddt)}]')
# ----------------------------------------------------------------------------
