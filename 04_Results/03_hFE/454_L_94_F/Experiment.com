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
    'cpus':24,
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
    'freqSimReq':OFF,
    'geostatic':OFF,
    'geotech':OFF,
    'heatTransfer':OFF,
    'impJobExpVars':{},
    'importJobList':[],
    'importer':OFF,
    'importerParts':OFF,
    'includes':['/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step01.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step02.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step03.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step04.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step05.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step06.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step07.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step08.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step09.inp', '/home/ms20s284/FRACTIB/04_Results/03_hFE/454_L_94_F/Steps/Step10.inp'],
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
    'mp_head_node':('artorg-mtrosa.campus.unibe.ch', 'artorg-mtrosa', '130.92.125.21'),
    'mp_host_list':(('artorg-mtrosa.campus.unibe.ch', 24),),
    'mp_mode':MPI,
    'mp_mode_requested':MPI,
    'mp_mpi_validate':OFF,
    'mp_rsh_command':'ssh -n -l ms20s284 %H %C',
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'noMultiHostElemLoop':[],
    'no_domain_check':1,
    'onestepinverse':OFF,
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
    'solverTypes':['DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT'],
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
