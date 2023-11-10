# GLOBAL VARIABLES:
LBLRTM_MRG_BIN ="ODintMerge"
OD_MERGE_FILE  ="OD_merged.dat"
TOPL_DIR=os.getcwd()

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# ==========
# Functions:
# ==========


# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

if __name__ == "__main__":
    # -------------
    # MAIN ROUTINE:
    # -------------
    print("lblrtmOda ODA_OPTION=", ODA_OPTION)

    mbmOdaHasBeenInvoked=MolPrfDirExists(".")

    # Check if option is to run the default method
    if (ODA_OPTION=="NONE" or ODA_OPTION=="STD" or ODA_OPTION=="BOTH" or ODA_OPTION=="BOTH2" or not mbmOdaHasBeenInvoked):

        cmd_str=LBLRTM_BIN
        print(cmd_str)
        st=time.time()
        os.system(cmd_str)
        et=time.time()
        print("LBLRTM exec Time=", 1000*(et-st), "ms")
    # Check if option is to run the Oda method. If not then exit
    #if (ODA_OPTION!="ODA" and ODA_OPTION!="BOTH" and ODA_OPTION!="BOTH2"):
    #    # No request to run ODA so exist
    #    sys.exit()

    # Check if the Molecular Profile Exists and if not then  create it with mbmOda
    if (not mbmOdaHasBeenInvoked):
        cmd_str="mbmOda"
        print(cmd_str)
        os.system(cmd_str)

    # Create the Oda version of TAPE28 via mrgmOda
    cmd_str="mrgOda"
    print(cmd_str)
    os.system(cmd_str)

    # If this is "ODA" option or BOTH2 option then copy TAPE28Oda
    # to TAPE28 so that molecfit will use the Oda version
    if (ODA_OPTION=="ODA" or ODA_OPTION=="BOTH2"):

        # First if there is already a TAPE28 file then this is the standard one
        # so rename it to TAPE28STD for debug purposes
        if os.path.exists(TAPE28_FILENAME):
            os.rename(TAPE28_FILENAME,   TAPE28STD_FILENAME)

        # Now copy TAPE28Oda to TAPE28
        shutil.copyfile(TAPE28ODA_FILENAME,TAPE28_FILENAME)
