#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """

    Read .dat simulation results file and store force-displacement
    curve together with experimental data. Then, use force result
    to assess which increment to extract deformation gradient. Write
    ODB reader from template and execute it.

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import yaml
import argparse
from Utils import *
from numba import njit

Show.ShowPlot = False

#%% Functions
# Define functions
def ReadConfigFile(Filename):

    """ Read configuration file and store to dictionary """

    print('\n\nReading configuration file', Filename)
    with open(Filename, 'r') as File:
        Configuration = yaml.load(File, Loader=yaml.FullLoader)

    return Configuration
def AdjustImageSize(Image, CoarseFactor, CropZ='Crop'):

    """
    Adapted from Denis's utils_SA.py
    Images are adjusted according to CropType:
    0 = CropType.expand     (Expand image by copying layers)
    1 = CropType.crop       (Crop image)
    2 = CropType.variable   (Either crop or expand, depending on what includes less layers)
    """

    # Get array
    Array = sitk.GetArrayFromImage(Image)
    Array = Array.transpose(2, 1, 0)

    # Measure image shape
    IMDimX = np.shape(Array)[0]
    IMDimY = np.shape(Array)[1]
    IMDimZ = np.shape(Array)[2]

    AddDimX = CoarseFactor - (IMDimX % CoarseFactor)
    AddDimY = CoarseFactor - (IMDimY % CoarseFactor)

    # adjust in x and y direction
    Shape_Diff = [AddDimX, AddDimY]
    IMG_XY_Adjusted = np.lib.pad(Array,
                                 ((0, Shape_Diff[0]), (0, Shape_Diff[1]), (0, 0)),
                                 'constant', constant_values=(0),)

    if CropZ == 'Crop':
        Image_Adjusted = IMG_XY_Adjusted

    if CropZ == 'Expand':
        AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
        Shape_Diff = [AddDimX, AddDimY, AddDimZ]
        Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                    ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                    'edge')

    if CropZ == 'Variable':
        Limit = CoarseFactor / 2.0
        if IMDimZ % CoarseFactor > Limit:
            AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
            Shape_Diff = [AddDimX, AddDimY, AddDimZ]
            Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                        ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                        'edge')
        if IMDimZ % CoarseFactor < Limit:
            Image_Adjusted = IMG_XY_Adjusted

    Image_Adjusted = sitk.GetImageFromArray(Image_Adjusted.transpose(2, 1, 0))
    Image_Adjusted.SetSpacing(Image.GetSpacing())
    Image_Adjusted.SetOrigin(Image.GetOrigin())
    Image_Adjusted.SetDirection (Image.GetDirection())

    return Image_Adjusted

@njit
def Decomposition(JacobianArray):

    # ProcessTiming(1, 'Decompose Jacobian of deformation')

    Terms = JacobianArray.shape[-1]
    ArrayShape = JacobianArray.shape[:-1]

    if Terms == 4:

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(ArrayShape)
        # VonMises_Strain = np.zeros(ArrayShape)
        # MaxShear = np.zeros(ArrayShape)

        # ProcessLength = ArrayShape[0] * ArrayShape[1]
        # Progress = 0
        for j in range(0, ArrayShape[0]):
            for i in range(0, ArrayShape[1]):
                F_d = JacobianArray[j, i, :].reshape((2,2))

                ## Unimodular decomposition of F
                J = np.linalg.det(F_d)
                SphericalCompression[j, i] = J

                if J > 0:
                    F_tilde = J ** (-1 / 3) * F_d
                    Norm_F_tilde = np.linalg.norm(F_tilde)
                else:
                    Norm_F_tilde = 0.0

                IsovolumicDeformation[j, i] = Norm_F_tilde

                # ## Optional: decomposition of F_tilde
                # R_tilde, U_tilde = polar(F_tilde)
                # Norm_U_tilde = np.sqrt(np.sum(U_tilde ** 2))

                ## Hydrostatic and deviatoric strain
                # I_d = np.matrix(np.eye(F_d.shape[0]))
                # E = 1/2 * (F_d.T * F_d - I_d)
                # Hydrostatic_E = -1/3 * np.trace(E) * I_d
                # Deviatoric_E = E - Hydrostatic_E
                #
                # HydrostaticStrain[k,j,i] = Hydrostatic_E[0,0]
                # MaxShear[k,j,i] = E.diagonal().max() - E.diagonal().min()
                #
                # VM_Strain = np.sqrt(3/2) * np.linalg.norm(Deviatoric_E)
                # VonMises_Strain[k,j,i] = VM_Strain

                # Progress += 1
                # ProgressNext(Progress/ProcessLength*20)

    elif Terms == 9:

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(ArrayShape)
        # VonMises_Strain = np.zeros(ArrayShape)
        # MaxShear = np.zeros(ArrayShape)

        # ProcessLength = ArrayShape[0] * ArrayShape[1] * ArrayShape[2]
        # Progress = 0
        for k in range(0, ArrayShape[0]):
            for j in range(0, ArrayShape[1]):
                for i in range(0, ArrayShape[2]):

                    F_d = JacobianArray[k, j, i, :].reshape((3, 3))

                    ## Unimodular decomposition of F
                    J = np.linalg.det(F_d)
                    SphericalCompression[k, j, i] = J

                    if J > 0:
                        F_tilde = J ** (-1 / 3) * F_d
                        Norm_F_tilde = np.linalg.norm(F_tilde)
                    else:
                        Norm_F_tilde = 0.0

                    IsovolumicDeformation[k, j, i] = Norm_F_tilde

                    # ## Optional: decomposition of F_tilde
                    # R_tilde, U_tilde = polar(F_tilde)
                    # Norm_U_tilde = np.sqrt(np.sum(U_tilde ** 2))

                    ## Hydrostatic and deviatoric strain
                    # I_d = np.matrix(np.eye(F_d.shape[0]))
                    # E = 1/2 * (F_d.T * F_d - I_d)
                    # Hydrostatic_E = -1/3 * np.trace(E) * I_d
                    # Deviatoric_E = E - Hydrostatic_E
                    #
                    # HydrostaticStrain[k,j,i] = Hydrostatic_E[0,0]
                    # MaxShear[k,j,i] = E.diagonal().max() - E.diagonal().min()
                    #
                    # VM_Strain = np.sqrt(3/2) * np.linalg.norm(Deviatoric_E)
                    # VonMises_Strain[k,j,i] = VM_Strain
                    
                    # Progress += 1
                    # ProgressNext(Progress/ProcessLength*20)

    return SphericalCompression, IsovolumicDeformation

def DecomposeJacobian(JacobianArray):
  
    Time.Process(1,'Decompose Jacobian')
    SC, ID = Decomposition(JacobianArray)
    Time.Process(0)

    return SC, ID

#%% For testing purpose
class Arguments:

    def __init__(self):
        self.Folder = 'FRACTIB'
        self.Sample = '441_R_64_M' # or 448_L_80_M

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    print('\nRead hFE results for sample ' + Arguments.Sample)

    # Set directories
    CWD, DD, SD, RD = SetDirectories('FRACTIB')
    FEADir = RD / '03_hFE' / Arguments.Sample
    ExpDir = RD / '02_Experiment' / Arguments.Sample

    # Read files
    print('\tPlot force-displacement hFE vs experiment')
    FEAData = Abaqus.ReadDAT(str(FEADir / 'Experiment.dat'))
    ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))

    # Truncate experiment to monotonic loading part and set to 0
    Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=1)
    MaxForce = ExpData['FZ'].idxmin()
    MaxDisp = ExpData['Z'].idxmax()
    DeltaTime = 10
    DeltaIndex = np.argmin(np.abs(ExpData['T']-DeltaTime))
    Start = Peaks[Peaks < MaxForce - DeltaIndex][-1]
    Stop = Peaks[Peaks > MaxDisp][0]
    
    ExpData = ExpData[Start:Stop].reset_index(drop=True)
    ExpData -= ExpData.loc[0]

    # Sum up fea steps results
    Indices = FEAData.groupby('Step')['Increment'].idxmax()
    FEA = FEAData.loc[Indices].cumsum()

    Show.Signal([FEA['Z']],[FEA['FZ']])

    # Truncate again experiment to match hFE
    Show.Signal([ExpData['Z'], FEA['Z']],
                [ExpData['Theta'], FEA['Theta']*180/np.pi],
                Labels=['Experiment', 'hFE'])

    # Compute min force location
    MinForceIdx = FEAData['FZ'][1:].abs().idxmin()
    Frame = FEAData.loc[MinForceIdx,'Increment']
    Step = FEAData.loc[MinForceIdx,'Step']

    # Write and execute ODB reader
    os.chdir(str(FEADir))
    with open(str(SD / 'Template_ODBReader.txt')) as Temp:
        Script = Temp.read()

    Context = {'Sample':Arguments.Sample,
               'Step':Step,
               'Frame':Frame}

    with open('ReadODB.py', 'w') as File:
        File.write(Script.format(**Context))

    print('\tRead odb')
    os.system('abaqus python ReadODB.py')

    # Read resulting files
    ElementsPositions = pd.read_csv('ElementsPositions.csv',names=['X','Y','Z'])
    DeformationGradients = pd.read_csv('DeformationGradients.csv',names=['F11','F12','F13','F21','F22','F23','F31','F32','F33'])

    # Build arrays
    X = np.unique(ElementsPositions['X'].values)
    Y = np.unique(ElementsPositions['Y'].values)
    Z = np.unique(ElementsPositions['Z'].values)

    F = np.zeros((len(Z),len(Y),len(X),9))
    for Index in DeformationGradients.index:
        
        Position = ElementsPositions.loc[Index]
        X_Index = list(X).index(Position['X'])
        Y_Index = list(Y).index(Position['Y'])
        Z_Index = list(Z).index(Position['Z'])
        
        F[Z_Index,Y_Index,X_Index] = DeformationGradients.loc[Index].values

    # Decompose deformation
    print('\tDecompose Jacobian and write to disk')
    SphericalCompression, IsovolumicDeformation = DecomposeJacobian(F)

    # Write MHDs
    # Compute metadata
    Spacing = np.array([X[1]-X[0],Y[1]-Y[0],Z[1]-Z[0]])
    Origin = np.array([X.min(), Y.min(), Z.min()])

    SC = sitk.GetImageFromArray(SphericalCompression)
    SC.SetSpacing(Spacing)
    SC.SetOrigin(Origin)

    ID = sitk.GetImageFromArray(IsovolumicDeformation)
    ID.SetSpacing(Spacing)
    ID.SetOrigin(Origin)

    Write.MHD(SC, 'J', PixelType='float')
    Write.MHD(ID, 'F_Tilde', PixelType='float')

    return

#%% Execution part
# Execution as main
if __name__ == '__main__':

    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)

    # Add long and short argument
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-V', '--Version', help='Show script version', action='version', version=SV)
    Parser.add_argument('-F', '--Folder', help='Root folder of the project', default='FRACTIB', type=str)
    Parser.add_argument('Sample', help='Sample to analyze (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)