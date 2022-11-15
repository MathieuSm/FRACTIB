#%%
#!/usr/bin/env python3

# 00 Initialization
import os
import time
import numpy as np
import sympy as sp
import pandas as pd
from scipy import ndimage
import SimpleITK as sitk
from pathlib import Path
import matplotlib.pyplot as plt
from Reader import ReadAIM, ShowSlice, PrintTime, ShowBinaryRegistration, GetSlice

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Functions
## Define functions
def RotationMatrix(Alpha=0, Beta=0, Gamma=0):

    Rx = sp.Matrix([[1,             0,              0],
                    [0, sp.cos(Alpha), -sp.sin(Alpha)],
                    [0, sp.sin(Alpha),  sp.cos(Alpha)]])

    Ry = sp.Matrix([[ sp.cos(Beta), 0, sp.sin(Beta)],
                    [0,             1,              0],
                    [-sp.sin(Beta), 0, sp.cos(Beta)]])

    Rz = sp.Matrix([[sp.cos(Gamma), -sp.sin(Gamma), 0],
                    [sp.sin(Gamma),  sp.cos(Gamma), 0],
                    [0,             0,              1]])

    R = Rz * Ry * Rx

    return np.array(R, dtype='float')
def WriteRaw(ImageArray, OutputFileName, PixelType):

    if PixelType == 'uint':

        MinValue = ImageArray.min()

        if MinValue < 0:
            ShiftedImageArray = ImageArray + abs(MinValue)
            MaxValue = ShiftedImageArray.max()
        elif MinValue > 0:
            ShiftedImageArray = ImageArray - MinValue
            MaxValue = ShiftedImageArray.max()
        else :
            ShiftedImageArray = ImageArray
            MaxValue = ShiftedImageArray.max()

        ScaledImageArray = ShiftedImageArray / MaxValue * 255

        CastedImageArray = ScaledImageArray.astype(np.uint8)

    elif PixelType == 'short':
        CastedImageArray = ImageArray.astype(np.short)
    elif PixelType == 'float':
        CastedImageArray = ImageArray.astype('float32')

    File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape)
    File[:] = CastedImageArray[:]
    del File

    return
def WriteMHD(ImageArray, Spacing, Offset, Path, FileName, PixelType='uint'):

    if PixelType == 'short' or PixelType == 'float':
        if len(ImageArray.shape) == 2:

            Array_3D = np.zeros((1,ImageArray.shape[0],ImageArray.shape[1]))

            for j in range(ImageArray.shape[0]):
                for i in range(ImageArray.shape[1]):
                    Array_3D[0,j,i] = ImageArray[j,i]

            ImageArray = Array_3D

    nz, ny, nx = np.shape(ImageArray)

    lx = float(Spacing[0])
    ly = float(Spacing[1])
    lz = float(Spacing[2])

    TransformMatrix = '1 0 0 0 1 0 0 0 1'
    X_o, Y_o, Z_o = float(Offset[0]), float(Offset[1]), float(Offset[2])
    CenterOfRotation = '0 0 0'
    AnatomicalOrientation = 'LPS'

    outs = open(os.path.join(Path, FileName) + '.mhd', 'w')
    outs.write('ObjectType = Image\n')
    outs.write('NDims = 3\n')
    outs.write('BinaryData = True\n')
    outs.write('BinaryDataByteOrderMSB = False\n')
    outs.write('CompressedData = False\n')
    outs.write('TransformMatrix = %s \n' % TransformMatrix)
    outs.write('Offset = %g %g %g\n' % (X_o, Y_o, Z_o))
    outs.write('CenterOfRotation = %s \n' % CenterOfRotation)
    outs.write('AnatomicalOrientation = %s \n' % AnatomicalOrientation)
    outs.write('ElementSpacing = %g %g %g\n' % (lx, ly, lz))
    outs.write('DimSize = %i %i %i\n' % (nx, ny, nz))

    if PixelType == 'uint':
        outs.write('ElementType = %s\n' % 'MET_UCHAR')
    elif PixelType == 'short':
        outs.write('ElementType = %s\n' % 'MET_SHORT')
    elif PixelType == 'float':
        outs.write('ElementType = %s\n' % 'MET_FLOAT')

    outs.write('ElementDataFile = %s\n' % (FileName + '.raw'))
    outs.close()

    WriteRaw(ImageArray, os.path.join(Path, FileName) + '.raw', PixelType)

    return
def TransformixTransformations(MovingImage,TransformParameterMap,ResultsDirectory):

    ## Compute jacobian of deformation field using transformix
    TransformixImageFilter = sitk.TransformixImageFilter()
    TransformixImageFilter.ComputeDeformationFieldOff()
    TransformixImageFilter.ComputeSpatialJacobianOff()
    TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOff()
    TransformixImageFilter.SetMovingImage(MovingImage)
    TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
    TransformixImageFilter.SetOutputDirectory(ResultsDirectory)

    TransformixImageFilter.Execute()

    ResultImage = TransformixImageFilter.GetResultImage()

    return ResultImage

# 01 Set variables
HRpQCT_Data = Path.cwd() / '../../02_Data/01_HRpQCT'
uCT_Data = Path.cwd() / '../../02_Data/02_uCT'
ResultsPath = Path.cwd() / '../../04_Results/05_FractureLinePrediction'
SampleList = pd.read_csv(str(Path.cwd() / '../../02_Data/SampleList.csv'))

#%%
for iSample, Sample in enumerate(SampleList['Internal ID']):

    iSample = 0
    Sample = SampleList.loc[iSample, 'Internal ID']

    # 03 Load Masks
    uCT_Number = SampleList.loc[iSample, 'MicroCT pretest file number']
    uCT_File = 'C000' + str(uCT_Number) + '_reso_0.098_DOWNSCALED_FULLMASK.mhd'
    uCT_Mask = sitk.ReadImage(str(uCT_Data / Sample / uCT_File))

    HRpQCT_Number = SampleList.loc[iSample, 'HRpQCT File 2 number']
    HRpQCT_File = 'C000' + str(HRpQCT_Number)
    Cort_File = str(HRpQCT_Data / Sample / (HRpQCT_File + '_CORT_MASK_UNCOMP.AIM'))
    Trab_File = str(HRpQCT_Data / Sample / (HRpQCT_File + '_TRAB_MASK_UNCOMP.AIM'))
    HRpQCT_Cort_Mask, AdditionalData = ReadAIM(Cort_File)
    HRpQCT_Trab_Mask, AdditionalData = ReadAIM(Trab_File)
    HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask

#%%
    # 04 Resample HR-pQCT image
    Direction = HRpQCT_Mask.GetDirection()
    Orig_Size = np.array(HRpQCT_Mask.GetSize())
    Orig_Spacing = np.array(HRpQCT_Mask.GetSpacing())

    New_Size = np.array(uCT_Mask.GetSize())
    New_Spacing = np.array(uCT_Mask.GetSpacing())

    Offset = (Orig_Size * Orig_Spacing - New_Size * New_Spacing) / 2

    Resample = sitk.ResampleImageFilter()
    Resample.SetInterpolator = sitk.sitkLinear
    Resample.SetOutputDirection(Direction)
    Resample.SetOutputOrigin(Offset)
    Resample.SetOutputSpacing(New_Spacing)
    Resample.SetSize(uCT_Mask.GetSize())

    HRpQCT_Resampled = Resample.Execute(HRpQCT_Mask)

#%%   
    # Rotate HR-pQCT image to 180°
    Rotation = sitk.VersorRigid3DTransform()
    M = RotationMatrix(Alpha=0, Beta=sp.pi, Gamma=0)
    M = [v for v in M.flatten()]

    # Symetry
    M[4] *= -1
    Rotation.SetMatrix(M)
    Center = np.array(HRpQCT_Resampled.GetSize()) / 2 * np.array(HRpQCT_Resampled.GetSpacing())
    CO = Center + np.array(HRpQCT_Resampled.GetOrigin())
    Rotation.SetCenter(CO)
    HRpQCT_Rotated = sitk.Resample(HRpQCT_Resampled, Rotation)

    # Extract slices for quick testing
    HRpQCT_Slice = GetSlice(HRpQCT_Rotated, int(HRpQCT_Rotated.GetSize()[2]*0.9))
    uCT_Slice = GetSlice(uCT_Mask, int(uCT_Mask.GetSize()[2]*0.9))

    # Binarize masks
    Otsu = sitk.OtsuThresholdImageFilter()
    Otsu.SetInsideValue(0)
    Otsu.SetOutsideValue(1)
    uCT_Bin = Otsu.Execute(uCT_Slice)

    # 07 Find best HR-pQCT image initial rotation with successive 90° rotations
    Measure = sitk.LabelOverlapMeasuresImageFilter()
    DSCs = pd.DataFrame()
    NRotations = 8
    Angle = 2*sp.pi/NRotations
    Rotation2D = sitk.Euler2DTransform()
    Rotation2D.SetCenter(CO[:2])
    for i in range(NRotations):

        print('\nRotate slice of %i degrees and register it' % (i*Angle/sp.pi*180))
        Tic = time.time()

        M = RotationMatrix(Alpha=0, Beta=0, Gamma=i*Angle)
        Rotation2D.SetMatrix([v for v in M[:2,:2].flatten()])
        Rotated_Slice = sitk.Resample(HRpQCT_Slice, Rotation2D)

        ## Register images
        ParameterMap = sitk.GetDefaultParameterMap('rigid')
        # ParameterMap['MaximumNumberOfIterations'] = '2000'
        ParameterMap['FixedImagePyramidSchedule'] = ['50', '20', '10']
        ParameterMap['MovingImagePyramidSchedule'] = ['50', '20', '10']
        ParameterMap['SP_alpha'] = ['0.6']
        ParameterMap['SP_A'] = ['1000']

        ## Set Elastix and perform registration
        ElastixImageFilter = sitk.ElastixImageFilter()
        ElastixImageFilter.SetParameterMap(ParameterMap)
        ElastixImageFilter.SetFixedImage(uCT_Slice)
        ElastixImageFilter.SetMovingImage(Rotated_Slice)
        ElastixImageFilter.SetOutputDirectory(str(ResultsPath / Sample))
        ElastixImageFilter.LogToConsoleOff()
        ElastixImageFilter.LogToFileOn()
        ElastixImageFilter.Execute()

        ## Get results
        Result_Image = ElastixImageFilter.GetResultImage()  # How moving image is deformed
        Result_Bin = Otsu.Execute(Result_Image)

        Toc = time.time()
        PrintTime(Tic, Toc)

        # ShowBinaryRegistration(uCT_Bin, Result_Bin)

        ## Compute dice coefficient
        Measure.Execute(uCT_Bin, Result_Bin)
        DSC = Measure.GetDiceCoefficient()
        NewData = pd.DataFrame({'Angle':float(i*Angle/sp.pi*180), 'DSC':DSC}, index=[i])
        DSCs = pd.concat([DSCs, NewData])

        if DSC == DSCs['DSC'].max():
            BestAngle = float(i*Angle)
            TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    #%% Build 3D initial transform
    # Build initial transform
    TransformParameters = TransformParameterMap[0]['TransformParameters']
    TransformParameters = ('0',
                           '0',
                           str(float(TransformParameters[0]) + BestAngle),
                           str(TransformParameters[1]),
                           str(TransformParameters[2]),
                           '0')
    CR = TransformParameterMap[0]['CenterOfRotationPoint']

    InitialTransform = TransformParameterMap[0]
    InitialTransform['CenterOfRotationPoint'] = (CR[0], CR[1], '0.0')
    InitialTransform['Direction'] = [str(d) for d in uCT_Mask.GetDirection()]
    InitialTransform['FixedImageDimension'] = str(uCT_Mask.GetDimension())
    InitialTransform['Index'] = ('0', '0', '0')
    InitialTransform['MovingImageDimension'] = str(HRpQCT_Mask.GetDimension())
    InitialTransform['NumberOfParameters'] = '6'
    InitialTransform['Origin'] = [str(o) for o in uCT_Mask.GetOrigin()]
    InitialTransform['Size'] = [str(s) for s in uCT_Mask.GetSize()]
    InitialTransform['Spacing'] = [str(s) for s in uCT_Mask.GetSpacing()]
    InitialTransform['TransformParameters'] = TransformParameters

    TransformFile = open(str(ResultsPath / Sample) + '/TransformParameters.0.txt')
    Text = TransformFile.read()
    TransformFile.close()
    for K in InitialTransform.keys():
        Start = Text.find(K)
        Stop = Text[Start:].find(')')
        OldText = Text[Start:Start+Stop]
        NewText = Text[Start:Start+len(K)]
        for i in InitialTransform[K]:
            NewText += ' ' + str(i)
        Text = Text.replace(OldText,NewText)
    TFileName = str(ResultsPath / Sample) + '/InitialTransform.txt'
    TransformFile = open(TFileName, 'w')
    TransformFile.write(Text)
    TransformFile.close()


#%%
    # Full registration using estimated parameters
    ParameterMap = sitk.GetDefaultParameterMap('rigid')
    ParameterMap['MaximumNumberOfIterations'] = ['2000']
    ParameterMap['FixedImagePyramidSchedule'] = ['50', '20', '10']
    ParameterMap['MovingImagePyramidSchedule'] = ['50', '20', '10']
    ParameterMap['SP_alpha'] = ['0.6']
    ParameterMap['SP_A'] = ['1000']

    ## Set Elastix and perform registration
    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetParameterMap(ParameterMap)
    ElastixImageFilter.SetFixedImage(uCT_Mask)
    ElastixImageFilter.SetMovingImage(HRpQCT_Rotated)
    ElastixImageFilter.SetInitialTransformParameterFileName(TFileName)
    ElastixImageFilter.SetOutputDirectory(str(ResultsPath / Sample))
    ElastixImageFilter.LogToConsoleOff()
    ElastixImageFilter.LogToFileOn()
    ElastixImageFilter.Execute()

    Result_Image = ElastixImageFilter.GetResultImage()

    ShowBinaryRegistration(uCT_Mask, Result_Image)

#%%
    
    Result_Mask = TransformixTransformations(Rotated_Mask, TransformParameterMap, ResultsDirectory)
    WriteMHD(Result_Mask, New_Spacing, Origin, ResultsDirectory, 'HR-pQCT_Mask', PixelType='float')

    # 09 Load uCT scan and write it
    uCT_Scan = sitk.ReadImage(uCT_SamplePath + '/' + uCT_Scan[0][:-13] + '.mhd')
    uCT_Scan = sitk.GetArrayFromImage(uCT_Scan)
    WriteMHD(uCT_Scan, New_Spacing, Origin, ResultsDirectory, 'uCT', PixelType='float')
    WriteMHD(uCT_Mask, New_Spacing, Origin, ResultsDirectory, 'uCT_Mask', PixelType='float')




# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
# Axes.imshow(Rotated_Scan[:,:,int(Rotated_Scan.shape[2]/2)], cmap='bone', alpha=0.5)
# Axes.imshow(uCT_Scan[:,:,int(uCT_Scan.shape[2]/2)], cmap='bone', alpha=0.5)
# plt.show()
# plt.close(Figure)
# %%
