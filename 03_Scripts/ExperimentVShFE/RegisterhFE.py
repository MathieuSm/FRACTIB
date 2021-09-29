#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from scipy import ndimage
import matplotlib.pyplot as plt


desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


## Define functions
def TransformixTransformations(MovingImage,TransformParameterMap,ResultsDirectory=False):

    ## Compute jacobian of deformation field using transformix
    TransformixImageFilter = sitk.TransformixImageFilter()
    TransformixImageFilter.ComputeDeformationFieldOff()
    TransformixImageFilter.ComputeSpatialJacobianOff()
    TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOff()
    TransformixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage))
    TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
    TransformixImageFilter.SetOutputDirectory(ResultsDirectory)

    TransformixImageFilter.Execute()

    ResultImage = TransformixImageFilter.GetResultImage()

    return sitk.GetArrayFromImage(ResultImage)
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


# 01 Set variables
WorkingDirectory = os.getcwd()
Data_Directory = os.path.join(WorkingDirectory,'02_Data/05_hFE/02_FEA/')
Results_Directory = os.path.join(WorkingDirectory,'04_Results/05_FractureLinePrediction/')

SampleList = [Dir for Dir in os.listdir(Data_Directory) if os.path.isdir(Data_Directory+Dir)]
SampleList.sort()

for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    # 02 Set Paths and Scans
    SamplePath = os.path.join(Data_Directory,Sample)
    ResultsPath = os.path.join(Results_Directory,Sample)

    # 03 Load images
    for hFE in ['J','F_Tilde']:
        hFE_Image = sitk.ReadImage(SamplePath + '/' + hFE +'.mhd')
        HRpQCT = sitk.ReadImage(ResultsPath + '/HR-pQCT.mhd')

        # 04 Rotate image according to pre-registration rotation
        hFE_Array = sitk.GetArrayFromImage(hFE_Image)
        hFE_Array[np.isnan(hFE_Array)] = 0
        hFE_Array = np.rot90(hFE_Array, 2, (0, 1))

        # 05 Rotate image according to registration estimation
        BestAngle = np.loadtxt(ResultsPath + '/HR-pQCT_RigidRotationAngle.txt').astype('int')
        RhFE_Array = ndimage.rotate(hFE_Array, BestAngle, (1,2), reshape=True)

        RhFE_Image = sitk.GetImageFromArray(RhFE_Array)
        RhFE_Image.SetSpacing(hFE_Image.GetSpacing())
        RhFE_Image.SetOrigin(hFE_Image.GetOrigin())
        RhFE_Image.SetDirection(hFE_Image.GetDirection())

        # 06 Resample hFE to get same sampling as HR-pQCT
        Offset = RhFE_Image.GetOrigin()
        Direction = RhFE_Image.GetDirection()
        Orig_Size = np.array(RhFE_Image.GetSize(), dtype=np.int)
        Orig_Spacing = RhFE_Image.GetSpacing()

        New_Spacing = HRpQCT.GetSpacing()

        Resample = sitk.ResampleImageFilter()
        Resample.SetInterpolator = sitk.sitkLinear
        Resample.SetOutputDirection(Direction)
        Resample.SetOutputOrigin(Offset)
        Resample.SetOutputSpacing(New_Spacing)
        Resample.SetSize(RhFE_Image.GetSize())

        New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(New_Spacing))
        New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
        New_Size = [int(s) for s in New_Size]
        Resample.SetSize(New_Size)

        ResampledImage = Resample.Execute(RhFE_Image)
        ResampledArray = sitk.GetArrayFromImage(ResampledImage)


        # 07 Transform hFE with same parameters as HR-pQCT
        TransformParameterMap = sitk.ReadParameterFile(ResultsPath + '/TransformParameters.0.txt')
        Transformed_Array = TransformixTransformations(ResampledArray, TransformParameterMap, ResultsDirectory=ResultsPath)

        Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
        Axes.imshow(Transformed_Array[:, :, int(Transformed_Array.shape[2] / 2)], cmap='bone')
        plt.show()
        plt.close(Figure)

        Spacing = np.array(New_Spacing)
        Origin = np.array([0, 0, 0])
        WriteMHD(Transformed_Array, Spacing, Origin, ResultsPath, hFE, PixelType='float')
