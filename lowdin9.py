#!/usr/bin/env python3
import sys
import numpy as np

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
    
    ###ANNETTE
    kValues = []
    #for k in range (nKpoints):
    for i in root.findall('EIGENSTATES/K-POINT'):
       # print(i.attrib)
       # print(i.text)        
        my_result = tuple(map(float, i.text.split()))
       # print(my_result)

    for i in root.findall('EIGENSTATES/K-POINT'):
        k_tuple = list(map(float, i.text.split()))
        kValues.append(k_tuple) 
    #print(kValues)


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

    return nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps,kValues

#end def

def collectValuesFromXML(xmlfile):

    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()

    ###ANNETTE
    atoms = []
    for i in root.iter('atom'):
        atoms.append(str(i.get('name'))+str(i.get('index')))
   # print(atoms)
   # print('\n')
    
    orbitals = ['s', 'pz', 'px', 'py', 'dz2', 'dzx','dzy','dx2-y2','dxy', 'f1', 'f2','f3','f4','f5','f6','f7']
    #print(orbitals)
    #print(orbitals[0])    

   # kValues = []
   # for i in root.iter('k_point'):
       # print(i.attrib)
       # print(i.text)
        #my_result = tuple(map(float, i.text.split()))
        #print(my_result)


###ADDED FROM JARON
    totmag = 0.0
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

    ###Annette
    del atoms[nAtom:]
    #print(atoms)


    return nAtom,nElec,int((nElec+totmag)/2),int((nElec-totmag)/2),orbitals,atoms

#end def

def matprint(m):
    for row in m:
        for element in row:
            print("%0.5f" % element),
        #end for
        print("\n")
    #end for
#end def

if __name__ == '__main__':

#    from developer import ci
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
    nBands,nKpoints,kWeights,nAtomicWFC,nSpin,atomicProjections,atomicOverlaps,invAtomicOverlaps,kValues = collectValuesFromAtomicProj(pw_outdir+"/"+pw_prefix+".save/atomic_proj.xml")

    # Collect parameters from <prefix>.xml
    nAtom,nElec,nOccUp,nOccDown, orbitals,atoms = collectValuesFromXML(pw_outdir+"/"+pw_prefix+".xml")

    print('\nNumber of up electrons: {}'.format(nOccUp))
    print('Number of down electrons: {}'.format(nOccDown))


    ###Annette modifying for supercell qmc runs
    import os
    import pandas as pd

    #Defintition to enter into each einspline file for manipulation
    vmcdirect =r"/gpfs/alpine/mat151/proj-shared/annette/bi2te3/runs_supercell_beststats/vmc_1rdm_noJ"
    os.chdir(vmcdirect)    
    
    def read_files(file_path):
        with open(file_path, 'r') as f:
            #print(f.readline())
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
            print("Unique twist indices: ", twistIndex)
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
            print("The associated k-points for each twist index: ", k_dict)
            #Next, create a dictionary, where for each supercell twist we can access the set of states (orbitals)
            s_dict = dict(zip(twistIndex, state_collection))
            print("The set of states per k-point and twist index: ", s_dict)
            
            #Check the number of orbitals is correct
            nStates = len(s_dict[twistIndex[0]])
            #print(nStates)
            #if nStates == 18: print("Good job!")
    
    #Reading all einspline files
    count = 0
    for f in os.listdir():
        if f.endswith('bandinfo.dat'):
            file_path = "{}/{}".format(vmcdirect,f)
            read_files(file_path)
            count += 1
            print(count)
    #print(count)

    quit()

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
    print(nm.shape)
    # Obtain dimensions of number matrices
    
    nblocks,nstates,nstates = nm[0][0].shape
    # Store stats of number matrix corresponding to single determinant with no jastrow, projected
    # on MO basis

    from numerics import simstats

    m_mo,v_mo,e_mo,k_mo = simstats(nm,dim=2) # stats over blocks
 
    # Perform "unitary" transform on each block's number matrix individually
    # and store in nmqmcu (i.e., up component of number matrix prime)
    # After the transformation, number matrix has been transformed from
    # the MO basis to the AO basis

    s=sp
    ###Annette
    print(atomicProjections.shape)
    nm2 =[] #new number matrix
    for tn in range(n_suptwists):
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

        nm2.append(nm_tmp)
    #end for
    nm2 = np.array(nm2)
    print(nm2.shape)  #shape is (36,1,400,72,72)
    for b in range(nblocks):
        nm_temp =[]
        for s in range(nStates):
            print(s)
            nm_temp.append(nm[0][0][b][s_dict[0][s],s_dict[0][s]]) #this is referencing the twist index dictionary for supertwist 0 to the primitive states
            print(s_dict[0][s])
    x = np.sum(nm_temp,axis = 0)
    print(x)

    quit()
    nmqmc = np.empty((nKpoints,nSpin,nblocks,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        for b in range(nblocks):
            nmqmc[k][s][b] = kWeights[k]*np.matmul(atomicProjections[k][s][:,:],np.matmul(nm[k][0][b][:,:],np.conj(atomicProjections[k][s][:,:].T)))
        #end for
    #end for
    
    #quit()
    m_ao,v_ao,e_ao,k_ao = simstats(nmqmc,dim=2)
    ###Annette added
    copy_means = m_ao
    #print(copy_means)
    #print("The dimensions of matrix are " + str(copy_means.shape))

    m_mo_avg = np.sum(unumpy.uarray(m_mo.real,e_mo.real),axis=0)
    m_ao_avg = np.sum(unumpy.uarray(m_ao.real,e_ao.real),axis=0)
    #print(m_ao_avg)
    #print("The dimensions of matrix are " + str(m_ao_avg.shape))
    
    ###Annette
    m_ao = unumpy.uarray(m_ao.real,e_ao.real)
    #print("The dimensions of array are " + str(m_ao.shape))
    #print(m_ao)


    # Obtain exact number matrix corresponding to single determinant with no jastrow, projected
    # on AO basis.

    exct_nmqmc = np.empty((nKpoints,nSpin,nAtomicWFC,nAtomicWFC),dtype=complex)
    for k in range(nKpoints):
        exct_nmqmc[k][s] = kWeights[k]*np.matmul(atomicProjections[k][s][:,:nOccUp],np.conj(atomicProjections[k][s][:,:nOccUp].T))
    #end for
    ###ANNETTE
    #for k in range(nKpoints):
    #    print(exct_nmqmc[k])
    #    print("/n")
    #print(exct_nmqmc)
    #print("The dimensions of matrix are " + str(exct_nmqmc.shape))
    exavg = np.sum(exct_nmqmc,axis=0)
    #print("The dimensions of the printed array are " + str(exavg.shape))
    #print(exavg)

    #for k in range(nKpoints):
    #    print(str(k) + " \n")
    #    for a in range(nAtomicWFC):
    #        print(exct_nmqmc[k][sp][a][a].real)
    #charges = np.zeros(nAtomicWFC,float)
    #for k in range(nKpoints):
    #    for a in range(nAtomicWFC):
    #        charges[a] = charges[a]+(exct_nmqmc[k][sp][a][a].real)
    #print(charges)
    #print("checking math"+"\n")
    #charges = np.zeros(nAtomicWFC,float)
    #for k in range(nKpoints):
    #    for a in range(nAtomicWFC):
    #        charges[a] = charges[a]+(m_ao[k][sp][a][a].real)
    #print(charges)

    cutoff = nAtomicWFC//nAtom
    #print("cutoff calc: nAtomicWFC (" + str(nAtomicWFC)+ ") + nAtom ("+ str(nAtom)+ ") -> cutoff = "+ str(cutoff))
    # Print real part of mean of number matrix in MO basis
    #print('nElec',nElec)

    quit()

    print("\n     Total Charge of system (QMCPACK): " + str(np.trace(m_ao_avg[s])) +"\n")
    for a in range(nAtomicWFC):
        print("          charge on AO "+str(a)+" = "+str(m_ao_avg[sp][a][a]))
    #end for

    print("\n     Total Charge of system (QE): " + str(np.trace(exavg[s].real)) +"\n")
    for a in range(nAtomicWFC):
        print("          charge on AO "+str(a)+" = "+str(exavg[sp][a][a].real))
    #end for

    
    quit()
    
    ###Annette's updated printing v3.0
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
     
    
    #print()
    ###Annette's newest additions
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
