def MergeTransData(mol_lst,r_lst,npts):
    n=len(r_lst)
    print("Combining ", n, " molecule transitivity profiles")
    tV=np.zeros(npts)
    for idx in range(n):
        sV=mol_lst[idx]
        r =r_lst[idx]
        #r=1.0/r
        if (r==0.0):
            continue
        for j in range(npts):
            trans=max(sV[j],1.0e-08)
            tau=-1.0*math.log(trans)
            tV[j]=tV[j]+r*tau

    for j in range(npts):
        tau=tV[j]
        tV[j]=math.exp(-tau)

    return tV

def MergeOpDepData(mol_lst,r_lst,npts):

    n=len(r_lst)
    print("Combining ", n, " molecule transitivity profiles")
    dfile="OD_DEBUG.txt"
    if (os.path.exists(dfile)):
        debug=False
    else:
        debug=True
        fid=open(dfile,"w")
        fid.write("DEBUG INFO n combs = " + " " + str(n) + "\n")

    if (debug==True):
        fid.write("DEBUG==TRUE\n")
    else:
        fid.write("DEBUG==FALSE\n")

    tV=np.zeros(npts)
    st=time.time()
    for idx in range(n):
        tauV=mol_lst[idx]
        r =r_lst[idx]
        #if (r==0.0):
        #    continue
        if (debug==True):
            fid.write(str(idx))
            fid.write(" ")
            fid.write(str(r))
            fid.write(" ")
            fid.write(str(tauV[0]))
            fid.write("\n")
        tV=tV+r*tauV
    if (debug==True):
        fid.write(str(tV[0]))
        fid.write("\n")
    tV=np.exp(-tV)
    if (debug==True):
        fid.write(str(tV[0]))
        fid.close()
    et=time.time()
    print("Merge Op Dep Data Time=", 1000*(et-st), "ms")
    return tV

def GetMolTv(mprfdir,molname):
    # Get the TAPE28 file from the relevant directory
    homedir=os.path.join(mprfdir,molname)
    print("PATH=",homedir)
    waveV, transV, npts=ReadTAPE28(homedir)
    #tV=np.zeros(n)
    return waveV,transV,npts

def GetMolOptDep(mprfdir,molname):
    # Store this as a numpy binary file
    outfile=os.path.join(mprfdir,molname,OPTICALDEPTHS_BINFILE)
    tauv=np.load(outfile)
    outfile=os.path.join(mprfdir,molname,OPTICALWVNUMS_BINFILE)
    wavev=np.load(outfile)
    n=np.size(wavev)
    return wavev,tauv,n


def GenerateTAPE28FMT (filename,header,waveV,mergedtransV,npts):

    fid=open(filename,'w')
    for line in header:
        fid.write(line)
    fid.write('\n')

    for i in range(npts):
        fval=waveV[i]
        wavestr=f'{fval:.8f}'
        fval=mergedtransV[i]
        transtr=f'{fval:.8f}'
        line="     " + wavestr + "        " + transtr +"\n"
        fid.write(line)
    fid.close()


if __name__ == "__main__":
    mprfdir=GetMolPrfDir('./')
    print("Inididual molecular profile components is in : ", mprfdir)

    lst=CompareTAPE5('./',mprfdir)
    print("lst=",lst)
    n=len(lst)
    #molTV_lst=[]
    molODV_lst=[]
    r_lst=[]
    for idx in range(n):
        r=lst[idx]
        if (r==0.0):
            continue
        molname=MOLECULES_FULL_LST[idx]
    #    print(idx,lst[idx],r,molname)
    #    waveV,tV,npts=GetMolTv    (mprfdir,molname)
        waveV,odV,npts=GetMolOptDep(mprfdir,molname)
    #    molTV_lst.append(tV)
        molODV_lst.append(odV)
        r_lst.append(r)

    #mergedtransV0   = MergeTransData(molTV_lst, r_lst,npts)
    mergedtransV  = MergeOpDepData(molODV_lst,r_lst,npts)
    #dump4plot("mrgd.dat",waveV,mergedtransV,npts)
    # Locate the TAPE28 file for ALL
    alldir=os.path.join(mprfdir,'ALL')

    header=GetTAPE28Header(alldir)
    #for line in header:
    #    print(line)
    GenerateTAPE28FMT (TAPE28ODA_FILENAME,header,waveV,mergedtransV,npts)

    #waveV,tV,npts=ReadTAPE28('./')
    #dump4plot("TAPE28.dat",waveV,tV,npts)
