
# WIP RELATED HACKS
TOPL_DIR=os.getcwd()
ALL_FLAG=True
#ALL_FLAG=False



# ==========
# Functions:
# ==========

def CleanOldMolecularDirectories(dir_names):

    if (os.path.exists(TAPE3_FILENAME)):
            os.remove(TAPE3_FILENAME)

    filename_lst=[TAPE1_FILENAME,TAPE5_FILENAME,TAPE10_FILENAME,TAPE3_FILENAME,TAPE6_FILENAME,TAPE7_FILENAME]

    for dirname in dir_names:
        if (os.path.exists(dirname)):
            print("Directory ", dirname, " exists")
            for filename in filename_lst:
                pathname=os.path.join(dirname,filename)
                if (os.path.exists(pathname)):
                    print("File ", pathname, " exists")
                    os.remove(pathname)
                if (os.path.islink(pathname)):
                    print("File ", pathname, " exists as link")
                    os.unlink(pathname)

            os.rmdir(dirname)


def ForceMakeDir(dir):

    # Forecably make a directory, ie if it already exists then delete first.

    if os.path.exists(dir):
            shutil.rmtree(dir)
    os.mkdir(dir)

def SysCall(cmd_str):
    print(cmd_str)
    start = time.time()
    os.system(cmd_str)
    duration = time.time() -start
    print ("Finished " + cmd_str + " Runtime = " + str(duration))

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

if __name__ == "__main__":
    # -------------
    # MAIN ROUTINE:
    # -------------
    # Clean out any leftover Molecular Specific Directories
    #CleanOldMolecularDirectories(MOLECULES_FULL_LST)
    print("lnflOda ODA_OPTION=", ODA_OPTION)

    # Check if option is to run the default method
    if (ODA_OPTION=="BOTH"):
        cmd_str=LNFL_BIN
        print("Running Standard:", cmd_str)
        os.system(cmd_str)
    #    shutil.copyfile(TAPE3_FILENAME,'TAPE3_BACKUP')

    # Check if option is to run the Oda method
    #if (ODA_OPTION!="ODA" and ODA_OPTION!="BOTH" and ODA_OPTION!="BOTH2"):
    #    # No request to run ODA so exist
    #    sys.exit()


    # If the output directory already exists then remove it
    if (os.path.exists(TAPE3_DIRNAME)):
        shutil.rmtree(TAPE3_DIRNAME)

    # Open TAPE 5 and read the 5 records that it contains
    fid=open(TAPE5_FILENAME,'r')
    record1=fid.readline()
    record2=fid.readline()
    record3=fid.readline()
    record4=fid.readline()
    record5=fid.readline()
    fid.close()

    # Parse record3 of the TAPE5 file into a list of molecules and define
    # sel_lst as the list ofindicies of the molecules flagged in record3
    lst=record3.split()
    field=lst[0]
    #print ("TAPE5 RECORD 3          =",record3)
    #print ("MOLECULE SELECTION FIELD=",field)
    sel_lst=[]
    for idx in range(0,N_MOL_TYPES):
        flag=field[idx]
        if (flag=="1"):
            sel_lst.append(idx)
            print("Selected:", MOLECULES_FULL_LST[idx])


    # Create a subdirectory for each selected molecule and populate with
    # a copy of TAPE5 but with record3 only having flag for that molecule
    for idx in sel_lst:
        mol_name=MOLECULES_FULL_LST[idx]
        field=""
        for idx2 in range(0,N_MOL_TYPES):
            if (idx2==idx):
                field=field+"1"
            else:
                field=field+"0"
        new_record=field+record3[N_MOL_TYPES:len(record3)]

        ForceMakeDir(mol_name)
        newfile=os.path.join(mol_name,TAPE5_FILENAME)
        fid=open(newfile,'w')
        fid.write(record1)
        fid.write(record2)
        fid.write(new_record)
        fid.write(record4)
        fid.write(record5)
        fid.close()

        # Now make a softlink for TAPE1
        newfile=os.path.join(mol_name,TAPE1_FILENAME)
        srcfile=os.path.join("../",TAPE1_FILENAME)
        os.symlink(srcfile,newfile)

    # Now execute lnfl in each molecule directory
    process_lst=[]
    for idx in sel_lst:
        mol_name=MOLECULES_FULL_LST[idx]
        cmd_str="cd " + mol_name + " ; " + LNFL_BIN
        process=multiprocessing.Process(target=SysCall,args=(cmd_str,))
        process_lst.append(process)
        process.start()

    # Join the processes
    for process in process_lst:
        process.join()


    '''
    cmd_lst=[]
    for idx in sel_lst:
        mol_name=MOLECULES_FULL_LST[idx]
        cmd_lst+=["cd " + mol_name + " ; " + LNFL_BIN,]

    doit_SysCall(mod=SysCall, cmds=cmd_lst)
    '''
