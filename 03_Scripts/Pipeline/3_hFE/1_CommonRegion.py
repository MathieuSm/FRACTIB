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

def GetTopAndBot(Image):

    """
    Return mean height of top and bottom surface
    """

    Array = sitk.GetArrayFromImage(Image).astype('bool')
    MidXSlice = Array[:,:,Array.shape[2] // 2]
    Sum = np.sum(MidXSlice, axis=0)
    Counts = np.bincount(Sum[Sum > 0])
    MeanSampleHeigth = np.argmax(Counts)
    MeanHeightPositions = np.where(Sum == MeanSampleHeigth)[0]

    TopNodes = []
    BotNodes = []
    for Position in MeanHeightPositions:
        Nodes = np.argwhere(MidXSlice[:,Position])
        TopNodes.append(Nodes.min())
        BotNodes.append(Nodes.max())

    MeanTop = int(np.mean(TopNodes).round(0))
    MeanBot = int(np.mean(BotNodes).round(0))
    
    return MeanTop, MeanBot


#%% Paths
# Set path variables

WD, Data, Scripts, Results = SetDirectories('FRACTIB')
hFE_Data = Data / '01_HRpQCT'
uCT_Data = Data / '02_uCT'
ResultsPath = Results / '05_Localizations'
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
Reader = Read()
HRpQCT_Cort_Mask, AdditionalData = Reader.AIM(Cort_File)
HRpQCT_Trab_Mask, AdditionalData = Reader.AIM(Trab_File)
HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask


#%% Preprocessing
# Rotate image to 180°
Rotation = sitk.VersorRigid3DTransform()
M = RotationMatrix(Alpha=0, Beta=sp.pi, Gamma=0)
M = [v for v in M.flatten()]

Rotation.SetMatrix(M)
Center = np.array(HRpQCT_Mask.GetSize()) / 2 * np.array(HRpQCT_Mask.GetSpacing())
CO = Center + np.array(HRpQCT_Mask.GetOrigin())
Rotation.SetCenter(CO)
HRpQCT_Resampled = sitk.Resample(HRpQCT_Mask, Rotation)


# Pad for rotations
Pad = 100
HRpQCT_Pad = sitk.ConstantPad(HRpQCT_Resampled, (Pad, Pad, 0), (Pad, Pad, 0))
#%% Starting position estimation
# Find best image initial position with successive rotations
Measure = sitk.LabelOverlapMeasuresImageFilter()
Dices = pd.DataFrame()

# Extract slices for rapid estimation
uCT_Slice = GetSlice(uCT_Mask, int(uCT_Mask.GetSize()[2]*0.8))
HRpQCT_Slice = GetSlice(HRpQCT_Pad, int(HRpQCT_Pad.GetSize()[2]*0.8))

# Binarize masks
Otsu = sitk.OtsuThresholdImageFilter()
Otsu.SetInsideValue(0)
Otsu.SetOutsideValue(1)
uCT_Bin = Otsu.Execute(uCT_Slice)

# Set variables
NRotations = 8
Angle = 2*sp.pi/NRotations
Rotation2D = sitk.Euler2DTransform()
PhysicalSize = np.array(HRpQCT_Pad.GetSize()) * np.array(HRpQCT_Pad.GetSpacing())
Center = (PhysicalSize + np.array(HRpQCT_Pad.GetOrigin())) / 2
Rotation2D.SetCenter(Center[:2])

for i in range(NRotations):

    # Set initial rotation
    print('\nRotate sample of %i degrees' % (i*Angle/sp.pi*180))
    M = RotationMatrix(Alpha=0, Beta=0, Gamma=i*Angle)
    Rotation2D.SetMatrix([v for v in M[:2,:2].flatten()])
    HRpQCT_Rotated = sitk.Resample(HRpQCT_Slice, Rotation2D)

    # Register images
    Dict = {'MaximumNumberOfIterations': [256]}
    RI, TPM = Register.Rigid(uCT_Slice, HRpQCT_Rotated, str(ResultsPath / Sample), Dict)
    HRpQCT_Bin = Otsu.Execute(RI)

    # Compute dice coefficient
    Measure.Execute(uCT_Bin, HRpQCT_Bin)
    Dice = Measure.GetDiceCoefficient()
    print('Registration dice coefficient %.3f' % (Dice))
    NewData = pd.DataFrame({'Angle':float(i*Angle/sp.pi*180), 'DSC':Dice}, index=[i])
    Dices = pd.concat([Dices, NewData])

    if Dice == Dices['DSC'].max():
        BestAngle = float(i*Angle)
        Parameters = np.array(TPM[0]['TransformParameters'], 'float')


#%% Full registration
# Perform full 3D registration with initial transform

T = sitk.Euler3DTransform()
R = RotationMatrix(Gamma=Angle + Parameters[0])
T.SetMatrix([Value for Value in R.flatten()])
T.SetTranslation((Parameters[0], Parameters[1], 0))
T.SetCenter(Center)
sitk.WriteTransform(T, str(ResultsPath / Sample / 'InitialTransform.txt'))

HRpQCT_I = sitk.Resample(HRpQCT_Pad, T)
HRpQCT_I = Otsu.Execute(HRpQCT_I)
RI, TPM = Register.Rigid(uCT_Mask, HRpQCT_I, str(ResultsPath / Sample))
HRpQCT_Bin = Otsu.Execute(RI)
uCT_Bin = Otsu.Execute(uCT_Mask)

Show.Registration(HRpQCT_Bin, uCT_Bin, Axis='X')
Show.Registration(HRpQCT_Bin, uCT_Bin)

#%% Common region
# Determine common region with parallel surfaces
Common_Raw = HRpQCT_Bin * uCT_Bin
Common_Array = sitk.GetArrayFromImage(Common_Raw)

# Get uCT top and bottom surface
uCT_Top, uCT_Bot = GetTopAndBot(uCT_Bin)

# Get HRpQCT top and bottom surface
S = HRpQCT_Resampled.GetSize()
HRpQCT_Top = sitk.Slice(HRpQCT_Resampled, (0, 0, 1), (S[0], S[1], 2))
HRpQCT_Bot = sitk.Slice(HRpQCT_Resampled, (0, 0, S[2]-2), (S[0], S[1], S[2]-1))

# Transform them into uCT space
HRpQCT_TopR = Register.Apply(HRpQCT_Top, TPM)
HRpQCT_BotR = Register.Apply(HRpQCT_Bot, TPM)

# Get lower top surface point and higher bottom surface point
HRpQCT_TopA = sitk.GetArrayFromImage(HRpQCT_TopR)
HRpQCT_BotA = sitk.GetArrayFromImage(HRpQCT_BotR)

if HRpQCT_TopA.sum() == 0: # If top surface is outside the image
    HRpQCT_TopL = 0
else:
    HRpQCT_TopL = np.max(np.argwhere(HRpQCT_TopA)[:,0])

if HRpQCT_BotA.sum() == 0: # If bottom surface is ouside the image
    HRpQCT_BotH = S[2]
else:
    HRpQCT_BotH = np.min(np.argwhere(HRpQCT_BotA)[:,0])

# Keep the most restrictive height
Common_Top = max([uCT_Top, HRpQCT_TopL])
Common_Bot = min([uCT_Bot, HRpQCT_BotH])
Common_Array[:Common_Top] = 0
Common_Array[Common_Bot:] = 0

Common = sitk.GetImageFromArray(Common_Array)
Common.SetDirection(Common_Raw.GetDirection())
Common.SetOrigin(Common_Raw.GetOrigin())
Common.SetSpacing(Common_Raw.GetSpacing())

Show.Registration(Common_Raw, Common, Axis='X')
Writer = Write()
Writer.MHD(Common, str(ResultsPath / Sample / 'CommonMask'))


#%% Test for inverse transform

Show.Registration(HRpQCT_Bin, Common, Axis='X')
# HRpQCT_A = sitk.GetArrayFromImage(HRpQCT_Bin)
# Fill = np.sum(Common_Array * HRpQCT_A) / Common_Array.sum()

# Test = Register.Apply(HRpQCT_I, TPM)
# Test = Otsu.Execute(Test)

RI, TPM_Inv = Register.Rigid(HRpQCT_I, uCT_Mask)
uCT_Inv = Otsu.Execute(RI)
Show.Registration(HRpQCT_I, uCT_Inv, Axis='X')

Common_Inv = Register.Apply(Common, TPM_Inv)
Common_Inv = Otsu.Execute(Common_Inv)
Show.Registration(HRpQCT_I, Common_Inv, Axis='X')
Writer.MHD(Common_Inv, 'Rotated')

Common_T = sitk.Resample(Common_Inv, T.GetInverse())
HRpQCT_R = Otsu.Execute(HRpQCT_Pad)
Show.Registration(HRpQCT_R, Common_T, Axis='X')

S = Common_T.GetSize()
Common_S = sitk.Slice(Common_T, (Pad, Pad, 0), (S[0]-Pad, S[1]-Pad, S[2]))
Common_S.SetOrigin((0, 0, 0))

Common_F = sitk.Resample(Common_S, Rotation)
Writer.MHD(Common_F, str(Results / '03_hFE' / Sample / 'CommonMask'))

# Tot test
# uCT_R = Resample(uCT_Mask, HRpQCT_I)
# T_Inv = Register.ComputeInverse(uCT_R, str(ResultsPath / Sample / 'TransformParameters.0.txt'))
# uCT_Inv = Otsu.Execute(RI)
# Show.Registration(HRpQCT_I, uCT_Inv, Axis='X')


# %%
