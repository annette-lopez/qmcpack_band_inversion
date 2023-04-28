#!/usr/bin/env python3

from nexus import settings, job, run_project
from nexus import generate_physical_system, read_structure
from nexus import generate_pwscf
from nexus import generate_projwfc
from nexus import generate_pw2qmcpack
from nexus import generate_qmcpack
from nexus import vmc

from qmcpack_input import dm1b
from qmcpack_input import sposet

settings(
    pseudo_dir = './pseudopotentials',
    runs = 'runs',
    results = '',
    sleep = 10,
    machine = 'andes',
    account = 'MAT151',
    status_only = 0,
    generate_only = 0,	#1 means only generate files do not submit run
)

qe_exe1 = '/ccs/proj/mat151/opt/qe_7.0_cmake/build_mpi/bin/pw.x'
qe_exe1 += ' -pd .true.'

qe_exe2 = '/ccs/proj/mat151/opt/qe_7.0_cmake/build_mpi/bin/projwfc.x'
qe_exe2 += ' -pd .true.'

qe_exe3  = '/ccs/proj/mat151/opt/qe_7.0_cmake/build_mpi/bin/pw2qmcpack.x'
qe_exe3 += ' -pd .true.'


qe_modules = '''
module load intel
module load hdf5/1.10.7
module load fftw/3.3.8
module load netlib-lapack
module load python  #should this be python3?
module list
'''

scf_job = job(nodes = 50, hours = 10, presub = qe_modules, app = qe_exe1, processes_per_node = 8)
projwfc_job = job(nodes = 50, hours = 3, presub = qe_modules, app = qe_exe2, processes_per_node = 8)
nscf_job = job(nodes = 50, hours = 10, presub = qe_modules, app = qe_exe1, processes_per_node = 8)
conv_job = job(nodes = 50, hours = 3, presub = qe_modules, app = qe_exe3, processes_per_node = 8)

#Now define material
structure = read_structure('./structures/Bi2Te3_primitive.xsf')
structure.change_units('B')		#units should be in bohr

system = generate_physical_system(
   structure = structure,
   kgrid = (6,6,1),
   kshift = (0,0,0),
   tiling = (2,2,1),
   Bi = 5, # Zeff/valence charge for Bi ECP
   Te = 6,
)

scf = generate_pwscf(
    identifier = 'scf',	
    path = 'scf',	    
    job = scf_job,
    calculation = 'scf',
    system = system,
    kgrid = (12,12,1),
    input_type = 'generic',
    input_dft = 'pbe',
    ecutwfc = 300,	    
    conv_thr = 1e-7,
    occupations  = 'smearing',
    smearing  = 'gaussian', 	
    electron_maxstep = 500,
    pseudos = 'Bi.ccECP.AREP.upf Te.ccECP.AREP.upf'.split(),
    kshift = (0,0,0),
    degauss = 0.01,
    mixing_beta = 0.6,
    nosym        = False,
    wf_collect   = False,
    nspin = 1,
    nbnd = 18
)

nscf = generate_pwscf(
   identifier  = 'nscf',
   path = 'nscf',
   job = nscf_job,
   verbosity = 'high',
   system = system,
   calculation = 'nscf',
   input_type = 'generic',
   input_dft = 'pbe',
   ecutwfc = 300,
   conv_thr = 1e-7,
   nbnd = 18, 
   occupations = 'smearing',
   smearing = 'gaussian',
   electron_maxstep = 500,
   degauss = 0.01,
   mixing_beta = 0.6,
   pseudos = 'Bi.ccECP.AREP.upf Te.ccECP.AREP.upf'.split(),
   nspin = 1,
   nosym = True,
   wf_collect = True,
   dependencies = (scf, 'charge_density'),
)

pwf = generate_projwfc(
   identifier = 'pwf',
   path = 'nscf',
   job = projwfc_job,
   lwrite_overlaps = True,
   lsym = False,
   dependencies = (nscf,'other')
)

conv = generate_pw2qmcpack(
   identifier   = 'conv',
   path         = 'nscf',
   job          = conv_job,
   write_psir   = False,
   dependencies = (nscf,'orbitals'),
)

qmc_exe = '/ccs/proj/mat151/jtkrogel/apps/andes/qmcpack/qmcpack-3.11.0/bin_andes/qmcpack_complex_cpu'
qmc_exe += ' -pd .true.'

qmc_modules = '''
module load gcc/9.3.0
#module load intel/19.0.3
module load openmpi/4.0.4
#module load essl
module load openblas/0.3.12
module load netlib-lapack
#module load netlib-scalapack
module load hdf5
module load fftw
export FFTW_ROOT=$OLCF_FFTW_ROOT
module load cmake/3.18.4
module load boost/1.74.0
#module load cuda
module load python/3.7-anaconda3
'''

qmc_job = job(nodes=50, processes_per_node=8, presub=qmc_modules, app=qmc_exe, hours=10)


dm_estimator = dm1b(
   energy_matrix = False,
   integrator = 'uniform_grid',
   points = 6,
   scale = 1.0,
   basis = sposet(type = 'bspline',size = 18*4, spindataset = 0),
   evaluator = 'matrix',
   center = (0,0,0),
   check_overlap = False,
)

qmc = generate_qmcpack(
   identifier = 'vmc_1rdm_noJ',
   path = 'vmc_1rdm_noJ',
   job = qmc_job,
   input_type = 'basic',
   system = system,
   pseudos = 'Bi.ccECP.AREP.xml Te.ccECP.AREP.xml'.split(),
   estimators = [dm_estimator],
   jastrows = [],
   calculations = [vmc(
       walkers = 1,
       warmupsteps = 20,
       blocks = 400,
       steps = 40,
       substeps = 2,
       timestep = 0.4
       )],
   dependencies = (conv, 'orbitals')
)

run_project()

