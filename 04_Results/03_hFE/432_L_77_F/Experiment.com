from driverConstants import *
from driverStandardMPI import StandardMPIAnalysis
import driverUtils, sys
options = {
    'SIMExt':'.sim',
    'ams':OFF,
    'analysisType':STANDARD,
    'applicationName':'analysis',
    'aqua':OFF,
    'ask_delete':OFF,
    'beamSectGen':OFF,
    'biorid':OFF,
    'cavityTypes':[],
    'cavparallel':OFF,
    'complexFrequency':OFF,
    'contact':OFF,
    'cosimulation':OFF,
    'coupledProcedure':OFF,
    'cpus':42,
    'cse':OFF,
    'cyclicSymmetryModel':OFF,
    'directCyclic':OFF,
    'direct_solver':DMP,
    'dsa':OFF,
    'dynStepSenseAdj':OFF,
    'dynamic':OFF,
    'excite':OFF,
    'externalField':OFF,
    'externalFieldCSEAux':OFF,
    'externalFieldExtList':['.sim', '.SMAManifest'],
    'externalFieldFiles':[],
    'externalFieldSimReader':None,
    'fieldImport':OFF,
    'filPrt':[],
    'fils':[],
    'finitesliding':OFF,
    'flexiblebody':OFF,
    'foundation':OFF,
    'geostatic':OFF,
    'geotech':OFF,
    'heatTransfer':OFF,
    'impJobExpVars':{},
    'importJobList':[],
    'importer':OFF,
    'importerParts':OFF,
    'includes':['/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step01.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step02.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step03.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step04.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step05.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step06.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step07.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step08.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step09.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step10.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step11.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step12.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step13.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step14.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step15.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step16.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step17.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step18.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step19.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step20.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step21.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step22.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step23.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step24.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step25.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step26.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step27.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step28.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step29.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step30.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step31.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step32.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step33.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step34.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step35.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step36.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step37.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step38.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step39.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step40.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step41.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step42.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step43.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step44.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step45.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step46.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step47.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step48.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step49.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/432_L_77_F/Steps/Step50.inp'],
    'initialConditionsFile':OFF,
    'input':'Simulation',
    'inputFormat':INP,
    'interactive':None,
    'interpolExtList':['.odb', '.sim', '.SMAManifest'],
    'job':'Experiment',
    'keyword_licenses':[],
    'lanczos':OFF,
    'libs':[],
    'magnetostatic':OFF,
    'massDiffusion':OFF,
    'materialresponse':OFF,
    'modifiedTet':OFF,
    'moldflowFiles':[],
    'moldflowMaterial':OFF,
    'mp_file_system':(DETECT, DETECT),
    'mp_head_node':('artorg-cortex.campus.unibe.ch', 'artorg-cortex', '130.92.125.37'),
    'mp_host_list':(('artorg-cortex.campus.unibe.ch', 42),),
    'mp_mode':MPI,
    'mp_mode_requested':MPI,
    'mp_mpi_validate':OFF,
    'mp_mpirun_path':'/usr/SIMULIA/EstProducts/2021/linux_a64/code/bin/SMAExternal/pmpi/bin/mpirun',
    'mp_rsh_command':'ssh -n -l ms20s284 %H %C',
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'noMultiHostElemLoop':[],
    'no_domain_check':1,
    'outputKeywords':ON,
    'parameterized':OFF,
    'partsAndAssemblies':OFF,
    'parval':OFF,
    'pgdHeatTransfer':OFF,
    'postOutput':OFF,
    'preDecomposition':ON,
    'restart':OFF,
    'restartEndStep':OFF,
    'restartIncrement':0,
    'restartStep':0,
    'restartWrite':OFF,
    'rezone':OFF,
    'runCalculator':ON,
    'simPack':OFF,
    'soils':OFF,
    'soliter':OFF,
    'solverTypes':['DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT'],
    'standard_parallel':ALL,
    'staticNonlinear':ON,
    'steadyStateTransport':OFF,
    'step':ON,
    'stepSenseAdj':OFF,
    'stressExtList':['.odb', '.sim', '.SMAManifest'],
    'subGen':OFF,
    'subGenLibs':[],
    'subGenTypes':[],
    'submodel':OFF,
    'substrLibDefs':OFF,
    'substructure':OFF,
    'symmetricModelGeneration':OFF,
    'tempNoInterpolExtList':['.fil', '.odb', '.sim', '.SMAManifest'],
    'thermal':OFF,
    'tmpdir':'/tmp',
    'tracer':OFF,
    'transientSensitivity':OFF,
    'unfold_param':OFF,
    'unsymm':ON,
    'user':'/home/ms20s284/FRACTIB/03_Scripts/UMAT.f',
    'visco':OFF,
    'xplSelect':OFF,
}
analysis = StandardMPIAnalysis(options)
status = analysis.run()
sys.exit(status)
