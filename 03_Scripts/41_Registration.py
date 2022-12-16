#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used to perform pre/post-test registration in 3 steps:
        1. Masks rigid registration
        2. Masks affine registration
        3. Gray scale B-Spline registration

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
import numba
import argparse
from Utils import *
from numba import njit
Read.Echo = False

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

def DecomposeJacobian(JacobianImage):

    # Determine 2D of 3D jacobian array
    JacobianArray = sitk.GetArrayFromImage(JacobianImage)
    
    print('\nDecompose Jacobian')
    Tic = time.time()
    SC, ID = Decomposition(JacobianArray)
    Toc = time.time()
    PrintTime(Tic, Toc)

    SphericalCompression = sitk.GetImageFromArray(SC)
    IsovolumicDeformation = sitk.GetImageFromArray(ID)

    for Image in [SphericalCompression, IsovolumicDeformation]:
        Image.SetSpacing(JacobianImage.GetSpacing())
        Image.SetDirection(JacobianImage.GetDirection())
        Image.SetOrigin(JacobianImage.GetOrigin())

    return SphericalCompression, IsovolumicDeformation



#%% Classes
# Define classes

class Arguments(): # for testing purpose

    def __init__(self):
        self.Sample = '432_L_77_F'
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    Show.ShowPlot = Arguments.Show

    print('\nRegister Pre / Post-test scans')
    TIC = time.time()

    # Set directories
    WD, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    DataDir = Data / '02_uCT' / Arguments.Sample
    ResultsDir = Results / '04_Registration' / Arguments.Sample
    os.makedirs(ResultsDir, exist_ok=True)


    # Read hFE config file
    ConfigFile = str(Scripts / '3_hFE' / 'ConfigFile.yaml')
    Config = ReadConfigFile(ConfigFile)


    # Read AIMs 
    ProcessTiming(1, '\nRead AIMs')

    Files = [File for File in os.listdir(DataDir) if File.endswith('DOWNSCALED.AIM')]
    Files.sort()

    for iFile, File in enumerate(Files):

        Image = Read.AIM(str(DataDir / File))[0]
        Spacing = Image.GetSpacing()
        ProgressNext(1+iFile*3 / 6 * 6)

        Cort = Read.AIM(str(DataDir / (File[:-4] + '_CORT_MASK.AIM')))[0]
        ProgressNext(2+iFile*3 / 6 * 6)

        Trab = Read.AIM(str(DataDir / (File[:-4] + '_TRAB_MASK.AIM')))[0]
        ProgressNext(3+iFile*3 / 6 * 6)

        Mask = Cort + Trab
        Mask.SetSpacing(Spacing)

        del Cort
        del Trab

        if iFile == 0: # Adjust size as in hFE analysis for size matching
            CoarseFactor = int(round(Config['ElementSize'] / Spacing[0]))
            PreI = AdjustImageSize(Image, CoarseFactor)
            PreM = AdjustImageSize(Mask, CoarseFactor)
        else:
            PostI = Image
            PostM = Mask

    ProcessTiming(0)


    # Pad for transformations    
    Pad = 3 * CoarseFactor
    PreI = sitk.ConstantPad(PreI, (Pad, Pad, Pad), (Pad, Pad, Pad))
    PreM = sitk.ConstantPad(PreM, (Pad, Pad, Pad), (Pad, Pad, Pad))
    PostI = sitk.ConstantPad(PostI, (Pad, Pad, 0), (Pad, Pad, 0))
    PostM = sitk.ConstantPad(PostM, (Pad, Pad, 0), (Pad, Pad, 0))


    # Align centers of gravity
    print('\nAlign centers of gravity')
    Tic = time.time()
    CenterType = sitk.CenteredTransformInitializerFilter.MOMENTS
    IniT = sitk.CenteredTransformInitializer(PreM, PostM, sitk.Euler3DTransform(), CenterType)
    Ini = sitk.Resample(PostM, PreM, IniT, sitk.sitkNearestNeighbor, PostM.GetPixelID())
    Toc = time.time()
    PrintTime(Tic, Toc)


    # Perform initial rotation
    print('\nPerform rotations for registration starting point')
    Tic = time.time()


    ## Extract slices for quick registration
    PreS = GetSlice(PreM, int(PreM.GetSize()[2]*0.8))
    PostS = GetSlice(PostM, int(PostM.GetSize()[2]*0.8))

    ## Set rotations variables
    NRotations = 8
    Angle = 2*sp.pi/NRotations
    Rotation2D = sitk.Euler2DTransform()
    PhysicalSize = np.array(PostM.GetSize()) * np.array(PostM.GetSpacing())
    Center = (PhysicalSize + np.array(PostM.GetOrigin())) / 2
    Rotation2D.SetCenter(Center[:2])

    ## Find best image initial position with successive rotations
    Measure = sitk.LabelOverlapMeasuresImageFilter()
    Dices = pd.DataFrame()
    for i in range(NRotations):

        # Set initial rotation
        M = RotationMatrix(Alpha=0, Beta=0, Gamma=i*Angle)
        Rotation2D.SetMatrix([v for v in M[:2,:2].flatten()])
        PostR = sitk.Resample(PostS, Rotation2D)

        # Register images
        Dict = {'MaximumNumberOfIterations': [256]}
        Result, TPM = Registration.Register(PreS, PostR, 'rigid', Dictionary=Dict)
        Result = sitk.Cast(Result, PreS.GetPixelID())

        # Compute dice coefficient
        Measure.Execute(PreS, Result)
        Dice = Measure.GetDiceCoefficient()
        NewData = pd.DataFrame({'Angle':float(i*Angle/sp.pi*180), 'DSC':Dice}, index=[i])
        Dices = pd.concat([Dices, NewData])

        if Dice == Dices['DSC'].max():
            # Show.Slice(Moving_Bin)
            Show.Overlay(PreS, Result)
            BestAngle = float(i*Angle)
            Parameters = np.array(TPM[0]['TransformParameters'], 'float')

    ## Apply best rotation
    T = sitk.Euler3DTransform()
    R = RotationMatrix(Gamma=BestAngle + Parameters[0])
    T.SetMatrix([Value for Value in R.flatten()])
    T.SetTranslation((Parameters[0], Parameters[1], 0))
    T.SetCenter(Center)

    PostI = sitk.Resample(PostI, PreI, IniT, sitk.sitkNearestNeighbor, PostI.GetPixelID())
    PostI = sitk.Resample(PostI, T)

    Toc = time.time()
    PrintTime(Tic, Toc)


    # Perform rigid registration and transform mask
    print('\nPerform rigid registration')
    Tic = time.time()
    RigidI, TPM = Registration.Register(PreI, PostI, 'rigid', PreM,  Path=str(ResultsDir))
    RigidM = Registration.Apply(PostM, TPM) 
    Toc = time.time()
    PrintTime(Tic, Toc)

    if Arguments.Type == 'Rigid':
        Show.FName = str(ResultsDir / 'RigidRegistration')
        Show.Overlay(PreI, RigidI, AsBinary=True)


    # Perform bspline registration
    if Arguments.Type == 'BSpline':
        print('\nPerform b-spline registration')
        Tic = time.time()

        ## Specific parameters
        Schedule = np.repeat([64, 32, 16, 8, 4, 2, 1],3)
        Dictionary = {'FixedImagePyramidSchedule':Schedule,
                    'MovingImagePyramidSchedule':Schedule,
                    'NewSamplesEveryIteration':['true'],
                    'SP_a':['0.1']}

        ## Match b-spline interpolation with elements size
        JFile = sitk.ReadImage(str(Results / '03_hFE' / Arguments.Sample / 'J.mhd'))
        Dictionary['FinalGridSpacingInPhysicalUnits'] = JFile.GetSpacing()

        ## Perform b-spline registration
        BSplineI, TPM = Registration.Register(PreI, RigidI, 'bspline', RigidM, str(ResultDir), Dictionary)
        Show.FName = str(ResultsDir / 'BSplineRegistration')
        Show.Overlay(PreI, BSplineI, AsBinary=True)

    # Compute deformation jacobian
    if Arguments.Jac == True:

        ## Resample mask to match hFE element size and apply transform to compute jacobian
        RigidR = Resample(RigidM, Factor=CoarseFactor)
        BSplineM = Registration.Apply(RigidR, TPM, str(ResultsDir), Jacobian=True)
        Toc = time.time()
        PrintTime(Tic, Toc)

        ## Read Jacobian
        JacobianFile = str(ResultsDir / 'fullSpatialJacobian.nii')
        JacobianImage = sitk.ReadImage(JacobianFile)
        JacobianImage.SetSpacing(PreM.GetSpacing())

        ## Resample Jacobian image
        NewSpacing = JFile.GetSpacing()
        ResampledJacobian = Resample(JacobianImage, Spacing=NewSpacing)

        ## Perform jacobian unimodular decomposition
        SphericalCompression, IsovolumicDeformation = DecomposeJacobian(JacobianImage)

        ## Write results
        JFile = str(ResultsDir / 'J')
        FFile = str(ResultsDir / 'F_Tilde')
        Write.MHD(SphericalCompression, JFile, PixelType='float')
        Write.MHD(IsovolumicDeformation, FFile, PixelType='float')

    TOC = time.time()
    print('Images registered')
    PrintTime(TIC,TOC)

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
    Parser.add_argument('Sample', help='Sample number (required)', type=str)

    # Add defaults arguments
    Parser.add_argument('-F', '--Folder', help='Root folder name', type=str, default='FRACTIB')
    Parser.add_argument('-T','--Type', help='Registration type', type=str, default='Rigid')
    Parser.add_argument('-S','--Show', help='Show plots', type=bool, default=False)
    Parser.add_argument('-J','--Jac', help='Compute deformation Jacobian', type=bool, default=False)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)