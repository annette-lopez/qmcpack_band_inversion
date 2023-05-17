#!/usr/bin/env python3
import sys
import numpy as np


#Access within runs/nscf/pwscf_output/pwscf.save/atomic_proj.xml
def collectValuesFromAtomicProj(xmlfile):

    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()
    
    header = root.find('.//HEADER')
    
    # Find number of bands
    nBands = int(header.attrib['NUMBER_OF_BANDS'])
    # Find number of kpoints
    nKpoints = int(header.attrib['NUMBER_OF_K-POINTS'])
    # Find number of atomic wave functions
    nAtomicWFC = int(header.attrib['NUMBER_OF_ATOMIC_WFC'])
    # Find number of spin components
    nSpin = int(header.attrib['NUMBER_OF_SPIN_COMPONENTS'])

    kWeights = np.empty((nKpoints),dtype=float)
    
    ###Create a list of the k-points
    kValues = []

    for i in root.findall('EIGENSTATES/K-POINT'):
        k_tuple = list(map(float, i.text.split()))
        kValues.append(k_tuple) 
    #print("In collectValuesFromAtomicProj: ", kValues)


    atomicProjections = np.empty((nKpoints,nSpin,nAtomicWFC,nBands),dtype=complex)
    # Find atomic projections
    for k in range(nKpoints):
        kWeights[k] = float(root.findall('EIGENSTATES/K-POINT')[k].attrib['Weight'])
        for s in range(nSpin):
            for awfc in range(nAtomicWFC):
                if nSpin==1:
                    for b, text in enumerate(root.findall('EIGENSTATES/PROJS')[k][awfc].text.strip().splitlines()):
                        proj = float(text.split()[0])
                        proj = proj+complex(0,float(text.split()[1]))
                        # zeroth element below is for spin-type. In this case there is only one
                        atomicProjections[k][0][awfc][b]=proj
                    #end for
                else:
                    for b, text in enumerate(root.findall('EIGENSTATES/PROJS')[s*nKpoints+k][awfc].text.strip().splitlines()):
                        proj = float(text.split()[0])
                        proj = proj+complex(0,float(text.split()[1]))
                        atomicProjections[k][s][awfc][b]=proj
                    #end for
                    #for b, text in enumerate(root.find('EIGENSTATES/PROJS')[k][s][awfc].text.strip().splitlines()):
                    #    proj = float(text.split()[0])
                    #    proj = proj+complex(0,float(text.split()[1]))
                    #    atomicProjections[k][s][awfc][b]=proj
                    ##end for
                #end if
            #end for
        #end for
    #end for

    atomicOverlaps = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)

    # Find atomic overlaps
    for k in range(nKpoints):
        for s in range(nSpin):
            if nSpin==1:
                for o, text in enumerate(root.findall('OVERLAPS/OVPS')[k].text.strip().splitlines()):
                    ovlp = float(text.split()[0])
                    ovlp = ovlp+complex(0,float(text.split()[1]))
                    atomicOverlaps[k][0][o//nAtomicWFC][o%nAtomicWFC]=ovlp
                #end for
            else:
                for o, text in enumerate(root.findall('OVERLAPS/OVPS')[s*nKpoints+k].text.strip().splitlines()):
                    ovlp = float(text.split()[0])
                    ovlp = ovlp+complex(0,float(text.split()[1]))
                    atomicOverlaps[k][s][o//nAtomicWFC][o%nAtomicWFC]=ovlp
                #end for
            #end if
        #end for
    #end for

    invAtomicOverlaps = np.copy(atomicOverlaps)
    tmp = np.copy(atomicOverlaps)
    # Store inverse of atomic overlaps
    for k in range(nKpoints):
        for s in range(nSpin):
            invAtomicOverlaps[k][s] = np.linalg.inv(tmp[k][s])
        #end for
    #end for

    ###Added a return for kValues list
    return nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps,kValues

#end def

#Access within runs/nscf/pwscf_output/pwscf.xml
def collectValuesFromXML(xmlfile):

    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()

    ###Collect atom names
    atoms = []
    for i in root.iter('atom'):
        atoms.append(str(i.get('name'))+str(i.get('index')))
    #print(atoms)
    
    orbitals = ['s', 'pz', 'px', 'py', 'dz2', 'dzx','dzy','dx2-y2','dxy', 'f1', 'f2','f3','f4','f5','f6','f7']
    #print(orbitals)


###Collect total magnetization, if not present take absolute magnetization
    totmag = 0.0 # default value
    assert root.find('.//magnetization') is not None
    s = root.find('.//magnetization/total')
    if s is not None:
        totmag = int(float(s.text))
    else:
        s = root.find('.//magnetization/absolute')
        if s is not None:
            assert float(s.text)<1e-3
        #end if
    #end if

    #totmag = int(float(root.find('.//magnetization/total').text))
    nElec = int(float(root.find('.//nelec').text))
    nAtom = int(float(root.find('.//atomic_structure').attrib['nat']))

    #Take only unique names
    del atoms[nAtom:]
    #print(atoms)

    ###Added a return for orbitals and atoms list
    return nAtom,nElec,int((nElec+totmag)/2),int((nElec-totmag)/2),orbitals,atoms

#end def

###Not used
def matprint(m):
    for row in m:
        for element in row:
            print("%0.5f" % element),
        #end for
        print("\n")
    #end for
#end def

if __name__ == '__main__':

    from developer import ci
    from qmcpack_analyzer import QmcpackAnalyzer
    from uncertainties import ufloat,unumpy

    # Exit if atomic_proj.xml and outdir locations not given
    if(len(sys.argv)<5):
        print("Usage: lowdin.py <pw_prefix> <pw_outdir> <qmc_directory> <qmc_identifier> <spin>")
        quit()
    #end if

    pw_prefix = sys.argv[1]
    pw_outdir = sys.argv[2]

    qmc_directory = sys.argv[3]
    qmc_identifier = sys.argv[4]

    # spin (up=0,down=1)
    sp = int(sys.argv[5])

    if not sp in (0,1):
        print('Invalid spin specfied: {}'.format(sp))
        print('Must be either 0 (up) or 1 (down)')
        quit()
    #end if

    # Collect parameters from atomic_proj.xml.
    # Note: if atomic_proj.xml was not generated from projwfc.x or if the <OVERLAPS> element is not present in atomic_proj.xml then
    #       you should try re-running projwfc.x on a single core and single thread with "-ndiag 1" -- this can sometimes help
    ###Added kValues
    nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps,kValues = collectValuesFromAtomicProj(pw_outdir+"/"+pw_prefix+".save/atomic_proj.xml")

    # Collect parameters from <prefix>.xml
    ###Added orbitals and atoms
    nAtom,nElec,nOccUp,nOccDown,orbitals,atoms = collectValuesFromXML(pw_outdir+"/"+pw_prefix+".xml")

    print('\nNumber of up electrons: {}'.format(nOccUp))
    print('Number of down electrons: {}'.format(nOccDown))


    ###Modification for supercell qmc runs
    import os
    import pandas as pd

    #Defintition to enter into each einspline file to map from supercell twists to primitive cell twists
    #vmcdirect =r"/gpfs/alpine/mat151/proj-shared/annette/bi2te3/runs_supercell_beststats/vmc_1rdm_noJ"
    
    cwd = os.getcwd()
    #print("The current working directory is: ", cwd)
    
    os.chdir(qmc_directory)   
    #nwd = os.getcwd() #new working directory (qmc directory)
    #print("The current working directory is: ", nwd)
    
    def read_files(file_path):
        with open(file_path, 'r') as f:
            #print(f.readline())
            df = pd.read_fwf(f, delim_whitespace = True)
            #Must separate single column into two
            df[['BandIndex', 'Energy']] = df['BandIndex Energy  '].apply(lambda x: pd.Series(str(x).split()))
            #print("Column names ", df.columns.tolist())
            #Remove unncessary columns
            df = df.drop(columns = ['#','BandIndex Energy  ','   Kx  ','   Ky  ', '   Kz  ','KmK', 'Unnamed: 12' ],axis =1)
            #print(df)

            #Access the relevant columns in the einspline.dat file
            twistseries = df['TwistIndex']
            k1series = df['   K1  ']
            k2series = df['   K2  ']
            k3series = df.iloc[:,5]
            stateseries = df['State ']

            #For each different TwistIndex create a list of k-point tuples
            twistIndex = pd.unique(twistseries) #these are our future dictionary keys
            #print(twistseries.values)
            #print("Unique twist indices: ", twistIndex)
            #print(type(twistIndex))
            #Identify the unique k-points
            kvals = []
            for i in range(len(twistIndex)):
                for j in range(len(twistseries)):
                    if twistseries[j] in twistIndex:
                        mytuple = (k1series[j], k2series[j],k3series[j])
                        kvals.append(mytuple)
            kpoints = pd.unique(kvals)
            #print("The associated k-points for each twist index: ", kpoints)

            #Construct the set of states associated with each twistIndex
            state_collection = []
            for i in range(len(twistIndex)):
                templist = []
                for j in range(len(twistseries)):
                    if twistseries[j] == twistIndex[i]:
                        templist.append(stateseries[j])
                        #print(templist)
                state_collection.append(templist)
            #print("The set of states per k-point and twist index: ", state_collection)

            #Next, create a dictionary, where for each supercell twist we can access the k-point
            k_dict = dict(zip(twistIndex, kpoints))
            #print("The associated k-points for each twist index: ", k_dict)
            #Next, create a dictionary, where for each supercell twist we can access the set of states (orbitals)
            s_dict = dict(zip(twistIndex, state_collection))
            #print("The set of states per k-point and twist index: ", s_dict)
            
            #Check the number of orbitals is correct
            nStates = len(s_dict[twistIndex[0]])
            #print(nStates)
            #if nStates == 18: print("Good job!")
            return (k_dict, s_dict)
    #end def 
    
    ###Original (nonsupercell) transformation
    if False: 

        # Analyze QMC data
        qa = [] # qmcpack_analyzer instance
        nm = [] # number matrix
        for tn in range(nKpoints):
            qa_tmp = QmcpackAnalyzer('{}/{}.g{:03d}.twistnum_{}.in.xml'.format(qmc_directory,qmc_identifier,tn,tn),verbose=False)
            qa_tmp.analyze()
            qa.append(qa_tmp)

            # get the density matrix (called the number matrix here)
            nm_tmp = []

            if sp==0:
                nm_tmp.append(qa[tn].qmc[0].DensityMatrices.number_matrix.u.data)
            else:
                nm_tmp.append(qa[tn].qmc[0].DensityMatrices.number_matrix.d.data)
            #end if

            nm.append(nm_tmp)
        #end for

        nm = np.array(nm)
        
        # Obtain dimensions of number matrices
    
        nblocks,nstates,nstates = nm[0][0].shape
    
    
    ###New supercell modification
    count = 0
    # Analyze QMC data
    qa = []  # qmcpack_analyzer instance
    nm = []  # number matrix
    
    #print("The current working directory is: ", cwd)
    #Reading all einspline files
    for f in os.listdir():
        #print(f)
        if f.endswith('bandinfo.dat'):
            file_path = "{}/{}".format(cwd + "/" + qmc_directory,f)
            k_dict, s_dict = read_files(file_path)
            tn = count
            print("\n Supercell Twist #", count)
            print("Unique twist indices: ", list(k_dict.keys())) 
            print("K-point per twist: \n",k_dict)
            print("Set of primitive cell states per twist: \n", s_dict)
            #Make a Qmcpack Analyzer instance
            qa_tmp = QmcpackAnalyzer('{}.g{:03d}.twistnum_{}.in.xml'.format(qmc_identifier,tn,tn),verbose=False)
            qa_tmp.analyze()
            qa.append(qa_tmp)

            # get the density matrix (called the number matrix here)
            nm_tmp = []

            if sp==0:
                nm_tmp.append(qa[tn].qmc[0].DensityMatrices.number_matrix.u.data)
            else:
                nm_tmp.append(qa[tn].qmc[0].DensityMatrices.number_matrix.d.data)

            #nm.append(nm_tmp)
            nmt = nm_tmp
            nmt = np.array(nmt)
            #print(nmt.shape)
            for key,value in k_dict.items():
                #assign a list the set of states for each key in the dictionary
                s = s_dict[key]
                s = np.array(s)
                #print(s)
                #print(nmt[:,:,s,s])
                nmk = nmt[:,:,s][:,:,:,s]
                #print(nmk.shape)
                nm.append(nmk)
            #end for
            count += 1

    nm = np.array(nm)
    #print("Final number matrix dimensions: ", nm.shape)
    nblocks,nstates,nstates = nm[0][0].shape
    
    #quit()
 
    # Store stats of number matrix corresponding to single determinant with no jastrow, projected
    # on MO basis

    from numerics import simstats

    ###Not used below: commented out m_mo
    #m_mo,v_mo,e_mo,k_mo = simstats(nm,dim=2) # stats over blocks
 
    # Perform "unitary" transform on each block's number matrix individually
    # and store in nmqmcu (i.e., up component of number matrix prime)
    # After the transformation, number matrix has been transformed from
    # the MO basis to the AO basis

    
    ###Original (nonsupercell) transformation; same as current implementation
    if False:
        s = sp


        nmqmc = np.empty((nKpoints,nSpin,nblocks,nAtomicWFC,nAtomicWFC),dtype=complex)
        for k in range(nKpoints):
            for b in range(nblocks):
                nmqmc[k][s][b] = kWeights[k]*np.matmul(atomicProjections[k][s][:,:],np.matmul(nm[k][0][b][:,:],np.conj(atomicProjections[k][s][:,:].T)))
            #end for
        #end for

    ###New supercell modification
    ###Spin here should be generalized later
    nmqmc = np.empty((nKpoints,nSpin,nblocks,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):#Problem here is that nKpoints its 144 from nscf, but 36 twists
        for b in range(nblocks):
            nmqmc[k][sp][b] = kWeights[k]*np.matmul(atomicProjections[k][sp][:,:],np.matmul(nm[k][sp][b][:,:],np.conj(atomicProjections[k][sp][:,:].T)))
        #end for
    #end for

    os.chdir(cwd) #we were in path: vmcdirect
    #here = os.getcwd()
    #print("The current working directory is: ", here)
    #quit()
    
    m_ao,v_ao,e_ao,k_ao = simstats(nmqmc,dim=2)
    
    ###K-point resolved m_ao
    copy_means = m_ao
    #print(copy_means)
    #print("The dimensions of matrix are " + str(copy_means.shape))
    
    ###Not used: commented out m_mo_avg
    #m_mo_avg = np.sum(unumpy.uarray(m_mo.real,e_mo.real),axis=0)
    
    m_ao_avg = np.sum(unumpy.uarray(m_ao.real,e_ao.real),axis=0)
    #print(m_ao_avg)
    #print("The dimensions of matrix are " + str(m_ao_avg.shape))
    
    ###K-point resolved m_ao with uncertainty
    m_ao = unumpy.uarray(m_ao.real,e_ao.real)
    #print("The dimensions of array are " + str(m_ao.shape))
    #print(m_ao)


    # Obtain exact number matrix corresponding to single determinant with no jastrow, projected
    # on AO basis.

    exct_nmqmc = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        exct_nmqmc[k][sp] = kWeights[k]*np.matmul(atomicProjections[k][sp][:,:nOccUp],np.conj(atomicProjections[k][sp][:,:nOccUp].T))
    #end for
    exavg = np.sum(exct_nmqmc,axis=0)

    cutoff = nAtomicWFC//nAtom
    #print("cutoff calc: nAtomicWFC (" + str(nAtomicWFC)+ ") + nAtom ("+ str(nAtom)+ ") -> cutoff = "+ str(cutoff))
    
    # Print real part of mean of number matrix in MO basis
    print('nElec',nElec)

    print("\n     Total Charge of system (QMCPACK): " + str(np.trace(m_ao_avg[sp])) +"\n")
    for a in range(nAtomicWFC):
        print("          charge on AO "+str(a)+" = "+str(m_ao_avg[sp][a][a]))
    #end for

    print("\n     Total Charge of system (QE): " + str(np.trace(exavg[sp].real)) +"\n")
    for a in range(nAtomicWFC):
        print("          charge on AO "+str(a)+" = "+str(exavg[sp][a][a].real))
    #end for

    #quit()
       
    ###Print the above for each atom and its orbital per k-point
    print("QMC" + "\n")
    for k in range(nKpoints):
        print("k point: " + str(kValues[k]) + "\n")
        for a in range(nAtom):
            if atoms[a][0:3] == "Bi1":
                cutoff = 16
                for i in range(0,cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i]) +"  = " + str(m_ao[k][sp][i][i]))
            elif atoms[a][0:3] == "Bi2":
                cutoff = 16
                for i in range(16,2*cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i - cutoff]) +"  = " + str(m_ao[k][sp][i][i]))
            elif atoms[a][0:3] == "Te3":
                cutoff = 9
                for i in range(32,32+cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i-32]) +"  = " + str(m_ao[k][sp][i][i]))
            elif atoms[a][0:3] == "Te4":
                cutoff = 9
                for i in range(41,32+2*cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i-41]) +"  = " + str(m_ao[k][sp][i][i]))
            elif atoms[a][0:3] == "Te5":
                cutoff = 9
                for i in range(50,32+3*cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i-50]) +"  = " + str(m_ao[k][sp][i][i]))

    print("Quantum Espresso"+ "\n")
    for k in range(nKpoints):
        print("k point: " + str(kValues[k]) + "\n")
        for a in range(nAtom):
            if atoms[a][0:3] == "Bi1":
                cutoff = 16
                for i in range(0,cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i]) +"  = " + str(exct_nmqmc[k][sp][i][i].real))
            elif atoms[a][0:3] == "Bi2":
                cutoff = 16
                for i in range(16,2*cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i - cutoff]) +"  = " + str(exct_nmqmc[k][sp][i][i].real))
            elif atoms[a][0:3] == "Te3":
                cutoff = 9
                for i in range(32,32+cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i-32]) +"  = " + str(exct_nmqmc[k][sp][i][i].real))
            elif atoms[a][0:3] == "Te4":
                cutoff = 9
                for i in range(41,32+2*cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i-41]) +"  = " + str(exct_nmqmc[k][sp][i][i].real))
            elif atoms[a][0:3] == "Te5":
                cutoff = 9
                for i in range(50,32+3*cutoff):
                    print("          charge on "+ str(atoms[a]) + " , " + str(orbitals[i-50]) +"  = " + str(exct_nmqmc[k][sp][i][i].real))


    ###Write the resolved data to a file
    import csv

    print("printing QMC resolved data")
    for a in range(nAtom):
        with open("{name}_QMC_resolved.dat".format(name = atoms[a]), "w") as f:
            for k in range(nKpoints):
                data = []
                occupations = []
                if atoms[a][0:3] == "Bi1":
                    cutoff = 16
                    for i in range(0,cutoff):
                        occupations.append(copy_means[k][sp][i][i].real)
                elif atoms[a][0:3] == "Bi2":
                    cutoff = 16
                    for i in range(16,2*cutoff):
                        occupations.append(copy_means[k][sp][i][i].real)
                elif atoms[a][0:3] == "Te3":
                    cutoff = 9
                    for i in range(32,32+cutoff):
                        occupations.append(copy_means[k][sp][i][i].real)
                elif atoms[a][0:3] == "Te4":
                    cutoff = 9
                    for i in range(41,32+2*cutoff):
                        occupations.append(copy_means[k][sp][i][i].real)
                elif atoms[a][0:3] == "Te5":
                    cutoff = 9
                    for i in range(50,32+3*cutoff):
                        occupations.append(copy_means[k][sp][i][i].real)
                writer = csv.writer(f, delimiter = " ")
                kValues[k].extend(occupations)
                writer.writerow(kValues[k])
            f.close()
            for j in range(nKpoints):
                del kValues[j][3:]

    print("printing QE resolved data")
    for a in range(nAtom):
        with open("{name}_QE_resolved.dat".format(name = atoms[a]), "w") as f:
            for k in range(nKpoints):
                data = []
                occupations = []
                if atoms[a][0:3] == "Bi1":
                    cutoff = 16
                    for i in range(0,cutoff):
                        occupations.append(exct_nmqmc[k][sp][i][i].real)
                elif atoms[a][0:3] == "Bi2":
                    cutoff = 16
                    for i in range(16,2*cutoff):
                        occupations.append(exct_nmqmc[k][sp][i][i].real)
                elif atoms[a][0:3] == "Te3":
                    cutoff = 9
                    for i in range(32,32+cutoff):
                        occupations.append(exct_nmqmc[k][sp][i][i].real)
                elif atoms[a][0:3] == "Te4":
                    cutoff = 9
                    for i in range(41,32+2*cutoff):
                        occupations.append(exct_nmqmc[k][sp][i][i].real)
                elif atoms[a][0:3] == "Te5":
                    cutoff = 9
                    for i in range(50,32+3*cutoff):
                        occupations.append(exct_nmqmc[k][sp][i][i].real)
                writer = csv.writer(f, delimiter = " ")
                kValues[k].extend(occupations)
                writer.writerow(kValues[k])
            f.close()
            for j in range(nKpoints):
                del kValues[j][3:]
#end if
