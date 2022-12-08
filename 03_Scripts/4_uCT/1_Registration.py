#%% #!/usr/bin/env python3
# Initialization

import yaml
import pandas as pd
from numba import njit

from Utils import *
Show = Show()
Read = Read()
Write = Write()

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Functions
# Define functions
def ReadConfigFile(Filename):

    """ Read configuration file and store to dictionary """

    print('\n\nReading initialization file', Filename)
    with open(Filename, 'r') as File:
        Configuration = yaml.load(File, Loader=yaml.FullLoader)

    return Configuration
def Adjust_Image_Size(Image, CoarseFactor, CropZ='Crop'):

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
    
    SC, ID = Decomposition(JacobianArray)

    SphericalCompression = sitk.GetImageFromArray(SC)
    IsovolumicDeformation = sitk.GetImageFromArray(ID)

    for Image in [SphericalCompression, IsovolumicDeformation]:
        Image.SetSpacing(JacobianImage.GetSpacing())
        Image.SetDirection(JacobianImage.GetDirection())
        Image.SetOrigin(JacobianImage.GetOrigin())

    return SphericalCompression, IsovolumicDeformation


#%% Loading
# Set paths and load data

CW, Data, Scripts, Results = SetDirectories('FRACTIB')

SampleList = pd.read_csv(str(Data / '02_uCT' / 'BMD_Values.csv'))
LogFile = open(str(Results / '04_Registration' / 'Registration.log'),'w+')
LogFile.write('Registration Log File\n\n')

ResultsData = pd.DataFrame()

#%% Set index
for Index in range(len(SampleList)):

    SampleTime = time.time()
#%% uCT files loading
# Load uCT files
Index = 0
print('\nLoad uCT files')
Tic = time.time()

Sample = SampleList.loc[Index,'Sample']
ResultsDirectory = str(Results / '04_Registration' / Sample)
os.makedirs(ResultsDirectory,exist_ok=True)

SampleData = {'Sample': Sample}
LogFile.write('Sample: ' + Sample + '\n')

SampleDirectory = str(Data / '02_uCT' / Sample) + '/'
Files = [File for File in os.listdir(SampleDirectory) if File.endswith('DOWNSCALED.AIM')]
Files.sort()

FixedImage = Read.AIM(SampleDirectory + Files[0])[0]
MovingImage = Read.AIM(SampleDirectory + Files[1])[0]
FixedCort = Read.AIM(SampleDirectory + Files[0][:-4] + '_CORT_MASK.AIM')[0]
FixedTrab = Read.AIM(SampleDirectory + Files[0][:-4] + '_TRAB_MASK.AIM')[0]
FixedMask = FixedCort + FixedTrab
FixedMask.SetSpacing(FixedImage.GetSpacing())

Toc = time.time()
PrintTime(Tic, Toc)
LogFile.write('Files loaded in %.3f s'%(Toc-Tic) + '\n')

#%% Preprocessing
# Mean moving image value
Array = sitk.GetArrayFromImage(MovingImage)
MeanValue = np.mean(Array)
del Array

# Cast fixed images to float
FixedImage = sitk.Cast(FixedImage, 8)
FloatMask = sitk.Cast(FixedMask, 8)
MovingImage = sitk.Cast(MovingImage, 8)

#%% Adapt image size to hFE meshing
ConfigFile = str(Scripts / '3_hFE' / 'ConfigFile.yaml')

# Read config and store to dictionary
Config = ReadConfigFile(ConfigFile)

# coarsening factor = FE element size / CT voxel size
Spacing = FixedImage.GetSpacing()
CoarseFactor = int(round(Config['ElementSize'] / Spacing[0]))
FixedImage = Adjust_Image_Size(FixedImage, CoarseFactor)
FixedMask = Adjust_Image_Size(FixedMask, CoarseFactor)
FloatMask = Adjust_Image_Size(FloatMask, CoarseFactor)


#%% Cog alignment
# Align centers of gravity
print('\nAlign centers of gravity')
Tic = time.time()
CenterType = sitk.CenteredTransformInitializerFilter.MOMENTS

IniTransform = sitk.CenteredTransformInitializer(FixedImage, MovingImage, sitk.Euler3DTransform(), CenterType)
IniMove = sitk.Resample(MovingImage, FixedImage, IniTransform, sitk.sitkLinear, 0.0, MovingImage.GetPixelID())
Toc = time.time()
PrintTime(Tic, Toc)
LogFile.write('Align centers of gravity in %.3f s' % (Toc - Tic) + '\n')

#%% Initial rotation
# Perform initial rotation
###!! If needed compute ellipse and align directions using cog !!
print('\nPerform rotations for registration starting point')
Tic = time.time()

# Pad for rotations
Pad = 100
MovingPad = sitk.ConstantPad(IniMove, (Pad, Pad, 0), (Pad, Pad, 0))

# Extract slices for rapid estimation
Fixed_Slice = GetSlice(FixedImage, int(FixedImage.GetSize()[2]*0.8))
Moving_Slice = GetSlice(MovingPad, int(MovingPad.GetSize()[2]*0.8))

# Initialize Otsu for results binarization
Otsu = sitk.OtsuThresholdImageFilter()
Otsu.SetInsideValue(0)
Otsu.SetOutsideValue(1)
Fixed_Bin = Otsu.Execute(Fixed_Slice)

# Set variables
NRotations = 8
Angle = 2*sp.pi/NRotations
Rotation2D = sitk.Euler2DTransform()
PhysicalSize = np.array(MovingPad.GetSize()) * np.array(MovingPad.GetSpacing())
Center = (PhysicalSize + np.array(MovingPad.GetOrigin())) / 2
Rotation2D.SetCenter(Center[:2])

# Find best image initial position with successive rotations
Measure = sitk.LabelOverlapMeasuresImageFilter()
Dices = pd.DataFrame()
for i in range(NRotations):

    # Set initial rotation
    print('\nRotate sample of %i degrees' % (i*Angle/sp.pi*180))
    M = RotationMatrix(Alpha=0, Beta=0, Gamma=i*Angle)
    Rotation2D.SetMatrix([v for v in M[:2,:2].flatten()])
    Moving_Rotated = sitk.Resample(Moving_Slice, Rotation2D)

    # Register images
    Dict = {'MaximumNumberOfIterations': [256]}
    RI, TPM = Register.Rigid(Fixed_Slice, Moving_Rotated, Path=ResultsDirectory, Dictionary=Dict)
    Moving_Bin = Otsu.Execute(RI)

    # Compute dice coefficient
    Measure.Execute(Fixed_Bin, Moving_Bin)
    Dice = Measure.GetDiceCoefficient()
    print('Registration dice coefficient %.3f' % (Dice))
    NewData = pd.DataFrame({'Angle':float(i*Angle/sp.pi*180), 'DSC':Dice}, index=[i])
    Dices = pd.concat([Dices, NewData])

    if Dice == Dices['DSC'].max():
        # Show.Slice(Moving_Bin)
        # Show.Registration(Fixed_Bin, Moving_Bin)
        BestAngle = float(i*Angle)
        Parameters = np.array(TPM[0]['TransformParameters'], 'float')



#%% Perform initial rotation

# Rotation

T = sitk.Euler3DTransform()
R = RotationMatrix(Gamma=BestAngle + Parameters[0])
T.SetMatrix([Value for Value in R.flatten()])
T.SetTranslation((Parameters[0], Parameters[1], 0))
T.SetCenter(Center)

Moving_R = sitk.Resample(IniMove, T)


#%% Rigid registration
# Perform rigid registration
RigidResult, TPM = Register.Rigid(FixedImage, Moving_R, FixedMask)

NegativeMask = (1-FloatMask) * MeanValue
Rigid_Bin = Otsu.Execute(RigidResult * FloatMask + NegativeMask)

#%% Non-rigid registration
# Perform non-rigid registration
Schedule = np.repeat([64, 32, 16, 8, 4, 2, 1],3)
Schedule = np.repeat([50, 20, 10, 5],3)
Dictionary = {'FixedImagePyramidSchedule':Schedule,
              'MovingImagePyramidSchedule':Schedule,
              'NewSamplesEveryIteration':['true']}

# Match b-spline interpolation with elements size
JFile = sitk.ReadImage(str(Results / '03_hFE' / Sample / 'J.mhd'))
Dictionary['FinalGridSpacingInPhysicalUnits'] = JFile.GetSpacing()

ResultImage, TPM = Register.NonRigid(FixedImage, RigidResult, FixedMask, ResultsDirectory, Dictionary)

BSpline_Bin = Otsu.Execute(ResultImage * FloatMask + NegativeMask)
Fixed_Bin = Otsu.Execute(FixedImage * FloatMask)

#%% Registration results
Show.Registration(FixedImage, RigidResult, Axis='X')
Show.Registration(FixedImage, ResultImage, Axis='X')

# Registration Dice coefficient
Measure.Execute(Fixed_Bin, BSpline_Bin)
Dice = Measure.GetDiceCoefficient()
print('\nDice coefficient of the full image: %.3f' % (Dice))

#%% Inverse
# Compute the inverse transform
InitialTransform = str(Path(ResultsDirectory, 'TransformParameters.0.txt'))
InverseTPM = Register.ComputeInverse(FixedImage, InitialTransform, FixedMask, Path=ResultsDirectory)

#%% Transformix
## Use transformix to compute spatial jacobian
ResultImage = Register.Apply(RigidResult, TPM, ResultsDirectory, Jacobian=True)


#%% Jacobian resampling
# Resample jacobian to match with hFE
JacobianFile = str(Path(ResultsDirectory, 'fullSpatialJacobian.nii'))
JacobianImage = sitk.ReadImage(JacobianFile)
JacobianImage.SetSpacing(FixedImage.GetSpacing())

## Resample Jacobian image
NewSpacing = JFile.GetSpacing()
ResampledJacobian = Resample(JacobianImage, Spacing=NewSpacing)

#%% Jacobian decomposition
## Perform jacobian unimodular decomposition
SphericalCompression, IsovolumicDeformation = DecomposeJacobian(JacobianImage)

#%% Write results
## Write results
JFile = str(Path(ResultsDirectory, 'J'))
FFile = str(Path(ResultsDirectory, 'F_Tilde'))
Write.MHD(SphericalCompression, JFile, PixelType='float')
Write.MHD(IsovolumicDeformation, FFile, PixelType='float')
# os.remove(os.path.join(ResultsDirectory, 'fullSpatialJacobian.nii'))


# %%
