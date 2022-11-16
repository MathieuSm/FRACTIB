#%% Initialization
#!/usr/bin/env python3

import pandas as pd
from Utils import *

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Functions
# Define functions

def Build3DTransform(FI, MI, TPM, Path=None, Alpha=0, Tx=0, Ty=0):

    """
    Write 3D transform parameter file from 2D estimation
    
    :param FI: Fixed image
    :param MI: Moving image
    :param TPM: 2D resulting transform parameter map
    :param Alpha: Additional intial angle
    :param Tx: Additional initial x translation
    :param Ty: Additional initial y translation
    """

    TP = TPM[0]['TransformParameters']
    TP = ('0','0',str(float(TP[0]) + Alpha),str(float(TP[1]) + Tx),str(float(TP[2]) + Ty),'0')
    CR = TPM[0]['CenterOfRotationPoint']

    IT = TPM[0]
    IT['CenterOfRotationPoint'] = (CR[0], CR[1], '0.0')
    IT['Direction'] = [str(d) for d in FI.GetDirection()]
    IT['FixedImageDimension'] = str(FI.GetDimension())
    IT['Index'] = ('0', '0', '0')
    IT['MovingImageDimension'] = str(MI.GetDimension())
    IT['NumberOfParameters'] = '6'
    IT['Origin'] = [str(o) for o in FI.GetOrigin()]
    IT['Size'] = [str(s) for s in FI.GetSize()]
    IT['Spacing'] = [str(s) for s in FI.GetSpacing()]
    IT['TransformParameters'] = TP

    if Path:
        TransformFile = open(str(Path / 'TransformParameters.0.txt'))
    else:
        TransformFile = open('TransformParameters.0.txt')

    Text = TransformFile.read()
    TransformFile.close()
    for K in IT.keys():
        Start = Text.find(K)
        Stop = Text[Start:].find(')')
        OldText = Text[Start:Start+Stop]
        NewText = Text[Start:Start+len(K)]
        for i in IT[K]:
            NewText += ' ' + str(i)
        Text = Text.replace(OldText,NewText)
    
    if Path:
        TFileName = str(Path / 'InitialTransform.txt')
    else:
        TFileName = 'InitialTransform.txt'

    TransformFile = open(TFileName, 'w')
    TransformFile.write(Text)
    TransformFile.close()

    return TFileName


#%% Paths
# Set path variables

WD, Data, Scripts, Results = SetDirectories('FRACTIB')
hFE_Data = Data / '01_HRpQCT'
uCT_Data = Data / '02_uCT'
ResultsPath = Results / '05_FractureLinePrediction'
SampleList = pd.read_csv(str(Data / 'SampleList.csv'))

#%% File Loading
# Load files

iSample = 0
Sample = SampleList.loc[iSample, 'Internal ID']

uCT_Number = SampleList.loc[iSample, 'MicroCT pretest file number']
uCT_File = 'C000' + str(uCT_Number) + '_reso_0.098_DOWNSCALED_FULLMASK.mhd'
uCT_Mask = sitk.ReadImage(str(uCT_Data / Sample / uCT_File))

HRpQCT_Number = SampleList.loc[iSample, 'HRpQCT File 2 number']
HRpQCT_File = 'C000' + str(HRpQCT_Number)
Cort_File = str(hFE_Data / Sample / (HRpQCT_File + '_CORT_MASK_UNCOMP.AIM'))
Trab_File = str(hFE_Data / Sample / (HRpQCT_File + '_TRAB_MASK_UNCOMP.AIM'))
HRpQCT_Cort_Mask, AdditionalData = Read.AIM(Cort_File)
HRpQCT_Trab_Mask, AdditionalData = Read.AIM(Trab_File)
HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask


#%% Preprocessing
# Rotate image to 180°
Rotation = sitk.VersorRigid3DTransform()
M = RotationMatrix(Alpha=0, Beta=sp.pi, Gamma=0)
M = [v for v in M.flatten()]

# Symetry (have to find out why)
M[4] *= -1
Rotation.SetMatrix(M)
Center = np.array(HRpQCT_Mask.GetSize()) / 2 * np.array(HRpQCT_Mask.GetSpacing())
CO = Center + np.array(HRpQCT_Mask.GetOrigin())
Rotation.SetCenter(CO)
HRpQCT_Rotated = sitk.Resample(HRpQCT_Mask, Rotation)


#%% Resampling
# Perform resampling for similar resolution
Direction = HRpQCT_Rotated.GetDirection()
Orig_Size = np.array(HRpQCT_Rotated.GetSize())
Orig_Spacing = np.array(HRpQCT_Rotated.GetSpacing())

New_Spacing = np.array(uCT_Mask.GetSpacing())
New_Size = np.array(uCT_Mask.GetSize())
# New_Size = np.ceil(Orig_Size * Orig_Spacing / New_Spacing)

Offset = (Orig_Size * Orig_Spacing - New_Size * New_Spacing) / 2

Resample = sitk.ResampleImageFilter()
Resample.SetInterpolator = sitk.sitkLinear
Resample.SetOutputDirection(Direction)
Resample.SetOutputOrigin(Offset)
Resample.SetOutputSpacing(New_Spacing)
# Resample.SetSize([int(S) for S in New_Size])
Resample.SetSize(uCT_Mask.GetSize())

HRpQCT_Resampled = Resample.Execute(HRpQCT_Rotated)


#%% Starting position estimation
# Find best image initial position with successive rotations
Measure = sitk.LabelOverlapMeasuresImageFilter()
Dices = pd.DataFrame()

# Extract slices for rapid estimation
uCT_Slice = GetSlice(uCT_Mask, int(uCT_Mask.GetSize()[2]*0.9))
HRpQCT_Slice = GetSlice(HRpQCT_Resampled, int(HRpQCT_Resampled.GetSize()[2]*0.9))

# Binarize masks
Otsu = sitk.OtsuThresholdImageFilter()
Otsu.SetInsideValue(0)
Otsu.SetOutsideValue(1)
uCT_Bin = Otsu.Execute(uCT_Slice)

# Set variables
NRotations = 8
Angle = 2*sp.pi/NRotations
Rotation2D = sitk.Euler2DTransform()
Rotation2D.SetCenter(CO[:2])

for i in range(NRotations):

    # Set initial rotation
    print('\nRotate sample of %i degrees' % (i*Angle/sp.pi*180))
    M = RotationMatrix(Alpha=0, Beta=0, Gamma=i*Angle)
    Rotation2D.SetMatrix([v for v in M[:2,:2].flatten()])
    HRpQCT_Rotated = sitk.Resample(HRpQCT_Slice, Rotation2D)

    # Register images
    Dict = {'MaximumNumberOfIterations': [256]}
    RI, TPM = Register.Rigid(uCT_Slice, HRpQCT_Rotated, None, str(ResultsPath / Sample), Dict)
    HRpQCT_Bin = Otsu.Execute(RI)

    # Compute dice coefficient
    Measure.Execute(uCT_Bin, HRpQCT_Bin)
    Dice = Measure.GetDiceCoefficient()
    print('Registration dice coefficient %.3f' % (Dice))
    NewData = pd.DataFrame({'Angle':float(i*Angle/sp.pi*180), 'DSC':Dice}, index=[i])
    Dices = pd.concat([Dices, NewData])

    if Dice == Dices['DSC'].max():
        BestAngle = float(i*Angle)
        TransformParameterMap = TPM


# %% Full registration
# Perform full 3D registration with initial transform

ITFile = Build3DTransform(uCT_Mask,
                          HRpQCT_Resampled,
                          TransformParameterMap,
                          ResultsPath / Sample,
                          Alpha = BestAngle)
RI, TPM = Register.Rigid(uCT_Mask, HRpQCT_Resampled, ITFile, str(ResultsPath / Sample))
HRpQCT_Bin = Otsu.Execute(RI)
uCT_Bin = Otsu.Execute(uCT_Mask)


# %% Common region
# Determine common region with parallel surfaces
Common = HRpQCT_Bin * uCT_Bin
CommonArray = sitk.GetArrayFromImage(Common) == 1
Sum = np.sum(CommonArray, axis=(1,2))

Start = np.argwhere(Sum)[0][0]
Stop = np.argmax(Sum)
Length = Stop - Start
Shift = round((Length * 0.05) / 2, 0).astype('int')
Start += Shift
Stop -= Shift

Figure, Axis = plt.subplots(1,1)
Axis.plot(Sum, color=(0,0,0))
Axis.plot([Start, Stop], [Sum[Start], Sum[Stop]], linestyle='none', marker='o', color=(1,0,0), fillstyle='none')
plt.show()

Ones = np.zeros(CommonArray.shape)
Ones[Start:Stop] = 1
Final = CommonArray * Ones
FinalImage = sitk.GetImageFromArray(Final)

Show.Registration(HRpQCT_Bin, FinalImage, Axis='X')
Show.Registration(uCT_Bin, FinalImage, Axis='X')


#%% Write images and masks
# Write images and masks

Cort_Mask = Register.Apply(HRpQCT_Cort_Mask, TPM)
Trab_Mask = Register.Apply(HRpQCT_Trab_Mask, TPM)

Cort_Bin = Otsu.Execute(Cort_Mask)
Trab_Bin = Otsu.Execute(Trab_Mask)

Final_Cort = Cort_Bin * FinalImage
Final_Trab = Trab_Bin * FinalImage

Final_HRpQCT = 9

Write.MHD(Final_Cort, str(ResultsPath / Sample / 'Cort_Mask'))
Write.MHD(Final_Trab, str(ResultsPath / Sample / 'Trab_Mask'))

### Transform original mask, multiply with common region and write

