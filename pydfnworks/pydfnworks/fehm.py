"""
Functions for using FEHM in dfnWorks
"""

def correct_perm_for_fehm():
    """ FEHM wants an empty line at the end of the perm file
    This functions adds that line return
    
    Parameters
    ----------
        None

    Returns
    ---------
        None

    Notes
    ------------
        Only adds a new line if the last line is not empty
    """
    fp = open("perm.dat")
    lines = fp.readlines()
    fp.close()
    # Check if the last line of file is just a new line
    # If it is not, then add a new line at the end of the file
    if len(lines[-1].split()) != 0:
        print("--> Adding line to perm.dat")
        fp = open("perm.dat","a")
        fp.write("\n")
        fp.close()

def fehm(self):
    """Run FEHM 

    Parameters
    ----------
        self : object 
            DFN Class
   
    Returns
    -------
        None

    Notes
    -----
    See https://fehm.lanl.gov/ for details about FEHM

    """
    print("--> Running FEHM")
    if self.flow_solver != "FEHM":
        sys.exit("ERROR! Wrong flow solver requested")
    try: 
        shutil.copy(self.dfnFlow_file, os.getcwd())
    except:
        sys.error("-->ERROR copying FEHM run file: %s"%self.dfnFlow_file)
    path = self.dfnFlow_file.strip(self.local_dfnFlow_file)
    fp = open(self.local_dfnFlow_file)
    line = fp.readline()
    fehm_input = line.split()[-1]
    fp.close()
    try: 
        shutil.copy(path+fehm_input, os.getcwd())
    except:
        sys.exit("-->ERROR copying FEHM input file:"%fehm_input)
    
    correct_perm_for_fehm()
    tic = time()
    subprocess.call(os.environ["FEHM_EXE"]+" "+self.local_dfnFlow_file, shell = True)
    print('='*80)
    print("FEHM Complete")
    print("Time Required %0.2f Seconds"%(time()-tic))
    print('='*80)
