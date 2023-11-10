import matplotlib.pyplot as plt
TMPDIR="plot28_tmpdir"

def ReadTAPE28Oda():

    # A hack as parsers work only on files called TAPE28
    # So we make a temp directory and copy TAPE28Oda
    # into it as TAPE28, parse that TAPE28 then remove
    # the tmp directtory

    # Safety catch
    if(os.path.exists(TMPDIR)):
        shutil.rmtree(TMPDIR)

    # Make tmp directory
    os.mkdir(TMPDIR)

    # Copy file`
    target=os.path.join(TMPDIR,TAPE28_FILENAME)
    shutil.copyfile(TAPE28ODA_FILENAME,target)

    # Parse copied file
    xV, yV, n = ReadTAPE28(TMPDIR)

    # Remove temporray directory
    shutil.rmtree(TMPDIR)

    # Return required data
    return xV,yV,n


if __name__ == "__main__":
    if (not os.path.exists(TAPE28_FILENAME)):
            print("No TAPE28 files in the current directory")
        # sys.exit()

    xV,  yV,  n  = ReadTAPE28("./")

    if os.path.exists(TAPE28ODA_FILENAME):

        xV2, yV2, n2 = ReadTAPE28Oda()
        plt.subplot(1, 2, 1) # row 1, col 2 index 1
        plt.plot(xV, yV)
        plt.title(TAPE28_FILENAME)
        plt.xlabel('Wavenumber')
        plt.ylabel('Transmission')

        plt.subplot(1, 2, 2) # index 2
        plt.plot(xV2, yV2)
        plt.title(TAPE28ODA_FILENAME)
        plt.xlabel('Wavenumber')
        plt.ylabel('Transmission')

    else:
        plt.plot(xV, yV)
        plt.title(TAPE28_FILENAME)
        plt.xlabel('Wavenumber')
        plt.ylabel('Transmission')

    plt.show()

    #plt.plot(xV,yV)
    #plt.plot(xV2,yV2)
    #plt.show()
