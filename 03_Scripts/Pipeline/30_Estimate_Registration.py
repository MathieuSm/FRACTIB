#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


## Define functions
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
def TransformixTransformation(MovingImage,TransformParameterMap,ResultsDirectory):

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
WorkingDirectory = os.getcwd()
DataFolder = os.path.join(WorkingDirectory, '02_Data/02_HRpQCT/')
uCTFolder = os.path.join(WorkingDirectory, '02_Data/03_uCT/')
ResultsFolder = os.path.join(WorkingDirectory, '04_Results/06_FractureLinePrediction/')

SampleList = os.listdir(DataFolder)
SampleList.sort()

Index = 0
for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    # 02 Set Paths and Scans
    SamplePath = os.path.join(DataFolder, Sample)
    Scans = [Scan for Scan in os.listdir(SamplePath) if Scan.endswith('.mhd')]
    Scans.sort()
    Scan = Scans[0][:8]

    uCT_SamplePath = os.path.join(uCTFolder, Sample)
    uCT_Scan = [Scan for Scan in os.listdir(uCT_SamplePath) if Scan.endswith('FULLMASK.mhd')]

    Results_Path = os.path.join(ResultsFolder, Sample)
    os.makedirs(Results_Path, exist_ok=True)

    # 03 Load uCT mask, segment it, and extract distal slice for quick registration
    uCT_Mask = sitk.ReadImage(uCT_SamplePath + '/' + uCT_Scan[0])
    Size = uCT_Mask.GetSize()
    uCT_Slice = sitk.Crop(uCT_Mask,(0,0,int(Size[2]*0.75)),(0,0,int(Size[2]*0.05)))
    Spacing = uCT_Slice.GetSpacing()
    Origin = uCT_Slice.GetOrigin()
    Direction = uCT_Slice.GetDirection()
    uCT_Slice_Array = sitk.GetArrayFromImage(uCT_Slice)
    uCT_Slice_Array[uCT_Slice_Array <= 0] = 0
    uCT_Slice_Array[uCT_Slice_Array > 0] = 1
    uCT_Slice = sitk.GetImageFromArray(uCT_Slice_Array)
    uCT_Slice.SetSpacing(Spacing)
    uCT_Slice.SetOrigin(Origin)
    uCT_Slice.SetDirection(Direction)



    # 04 Load HRpQCT mask, rotate it, and extract distal slice for quick registration
    HRpQCT_Cort_Mask = sitk.ReadImage(SamplePath + '/' + Scan + '_CORT_MASK_UNCOMP.mhd')
    HRpQCT_Trab_Mask = sitk.ReadImage(SamplePath + '/' + Scan + '_TRAB_MASK_UNCOMP.mhd')
    HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask
    Size = HRpQCT_Mask.GetSize()
    Spacing = HRpQCT_Mask.GetSpacing()
    Origin = HRpQCT_Mask.GetOrigin()
    Direction = HRpQCT_Mask.GetDirection()

    HRpQCT_Mask_Array = sitk.GetArrayFromImage(HRpQCT_Mask)
    HRpQCT_Mask_Array = np.rot90(HRpQCT_Mask_Array, 2, (0, 1))
    HRpQCT_Mask = sitk.GetImageFromArray(HRpQCT_Mask_Array)
    HRpQCT_Mask.SetSpacing(Spacing)
    HRpQCT_Mask.SetOrigin(Origin)
    HRpQCT_Mask.SetDirection(Direction)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(HRpQCT_Mask_Array[:, :, int(HRpQCT_Mask_Array.shape[2] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)

    HRpQCT_Slice = sitk.Crop(HRpQCT_Mask,(0,0,int(Size[2]*0.75)),(0,0,int(Size[2]*0.05)))
    Slice_Origin = HRpQCT_Slice.GetOrigin()
    HRpQCT_Slice_Array = sitk.GetArrayFromImage(HRpQCT_Slice)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(HRpQCT_Slice_Array[:, :, int(HRpQCT_Slice_Array.shape[2] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)


    # 07 Find best HR-pQCT image initial rotation
    SpacingRatio = np.array(uCT_Mask.GetSpacing())/(np.array(HRpQCT_Mask.GetSpacing()))
    Points = np.array(HRpQCT_Mask.GetSize()) / SpacingRatio
    Points_X = np.linspace(0,HRpQCT_Mask.GetSize()[0]-1,int(Points[0]))
    Points_Y = np.linspace(0,HRpQCT_Mask.GetSize()[1]-1,int(Points[1]))

    Rotations = range(0, 360, 15)
    DSCs = pd.DataFrame()
    TranslationParameters = {}
    for Rotation in Rotations:

        Rotated_Slice_Array = ndimage.rotate(HRpQCT_Slice_Array, Rotation, (1, 2), reshape=True)
        Rotated_Slice = sitk.GetImageFromArray(Rotated_Slice_Array)
        Rotated_Slice.SetSpacing(Spacing)
        Rotated_Slice.SetOrigin(Slice_Origin)
        Rotated_Slice.SetDirection(Direction)

        # # Resample rotated slice for plot
        # Points = Rotated_Slice_Array.shape / SpacingRatio
        # Points_X = np.linspace(0, Rotated_Slice_Array.shape[2] - 1, int(Points[2]))
        # Points_Y = np.linspace(0, Rotated_Slice_Array.shape[1] - 1, int(Points[1]))
        #
        # RotatedImage = Rotated_Slice_Array[int(Rotated_Slice_Array.shape[0] / 2), :, :]
        # RotatedImage = RotatedImage[Points_Y.astype('int'), :]
        # RotatedImage = RotatedImage[:, Points_X.astype('int')]
        #
        # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
        # Axes.imshow(uCT_Slice_Array[int(uCT_Slice_Array.shape[0]/2),:,:], cmap='bone', alpha=0.5)
        # Axes.imshow(RotatedImage, cmap='bone', alpha=0.5)
        # plt.show()
        # plt.close(Figure)

        ## Register images
        ParameterMap = sitk.GetDefaultParameterMap('translation')

        ## Set Elastix and perform registration
        ElastixImageFilter = sitk.ElastixImageFilter()
        ElastixImageFilter.SetParameterMap(ParameterMap)
        ElastixImageFilter.SetFixedImage(uCT_Slice)
        ElastixImageFilter.SetMovingImage(Rotated_Slice)
        ElastixImageFilter.SetOutputDirectory(Results_Path)
        ElastixImageFilter.LogToConsoleOn()
        ElastixImageFilter.Execute()

        ## Get results
        ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
        ResultImage_Array = sitk.GetArrayFromImage(ResultImage)

        ## Store translation
        TranslationFile = open(Results_Path + '/TransformParameters.0.txt', 'r')
        TranslationParameters[Rotation] = TranslationFile.read()

        ## Segment result image
        Otsu_Filter = sitk.OtsuThresholdImageFilter()
        Otsu_Filter.SetInsideValue(0)
        Otsu_Filter.SetOutsideValue(1)
        Segmentation = Otsu_Filter.Execute(ResultImage)
        R_Threshold = Otsu_Filter.GetThreshold()
        ResultImage_Array[ResultImage_Array <= R_Threshold] = 0
        ResultImage_Array[ResultImage_Array > 0] = 1

        # Plot results
        Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
        Axes.imshow(uCT_Slice_Array[int(uCT_Slice_Array.shape[0] / 2), :, :], cmap='bone', alpha=0.5)
        Axes.imshow(ResultImage_Array[int(ResultImage_Array.shape[0] / 2), :, :], cmap='bone', alpha=0.5)
        plt.show()
        plt.close(Figure)

        DSC = 2 * np.sum(ResultImage_Array * uCT_Slice_Array) / np.sum(ResultImage_Array + uCT_Slice_Array)
        DSCs = DSCs.append({'Angle': Rotation, 'DSC': DSC}, ignore_index=True)
        DSCs.to_csv(Results_Path + '/DSCs.csv',index=False)

    ## Verify best angle
    DSCs = pd.read_csv(Results_Path + '/DSCs.csv')
    BestAngle = DSCs.loc[DSCs['DSC'].idxmax(), 'Angle']
    Rotated_Slice_Array = ndimage.rotate(HRpQCT_Slice_Array, BestAngle, (1, 2), reshape=True)
    Rotated_Slice = sitk.GetImageFromArray(Rotated_Slice_Array)
    Rotated_Slice.SetSpacing(Spacing)
    Rotated_Slice.SetOrigin(Slice_Origin)
    Rotated_Slice.SetDirection(Direction)

    TranslationParameter = TranslationParameters[BestAngle]
    File = open(Results_Path + '/InitialTranslation.txt', 'w')
    File.write(TranslationParameter)
    File.close()
    ParameterMap = sitk.ReadParameterFile(Results_Path + '/InitialTranslation.txt')
    ResultImage = TransformixTransformation(Rotated_Slice, ParameterMap, Results_Path)
    ResultImage_Array = sitk.GetArrayFromImage(ResultImage)

    ## Segment result image
    Otsu_Filter = sitk.OtsuThresholdImageFilter()
    Otsu_Filter.SetInsideValue(0)
    Otsu_Filter.SetOutsideValue(1)
    Segmentation = Otsu_Filter.Execute(ResultImage)
    R_Threshold = Otsu_Filter.GetThreshold()
    ResultImage_Array[ResultImage_Array <= R_Threshold] = 0
    ResultImage_Array[ResultImage_Array > 0] = 1

    ## Plot results
    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(uCT_Slice_Array[int(uCT_Slice_Array.shape[0] / 2), :, :], cmap='bone', alpha=0.5)
    Axes.imshow(ResultImage_Array[int(ResultImage_Array.shape[0] / 2), :, :], cmap='bone', alpha=0.5)
    Axes.set_title('Best Results: ' + str(int(BestAngle)) + ' degrees')
    plt.show()
    plt.close(Figure)
