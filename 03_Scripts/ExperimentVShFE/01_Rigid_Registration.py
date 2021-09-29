#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import time

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


## Define functions
def ElastixRegistration(Dictionary):

    # Get dictionary parameters
    Transformations = Dictionary['Transformations']
    if type(Transformations) is not list:
        Transformations = list([Transformations])
    FixedImage = Dictionary['FixedImage']
    MovingImage = Dictionary['MovingImage']
    FixedMask = Dictionary['FixedMask']
    MovingMask = Dictionary['MovingMask']
    PyramidSchedule = Dictionary['PyramidSchedule']
    NIterations = Dictionary['NIterations']
    Alpha = Dictionary['Alpha']
    A = Dictionary['A']
    ResultsDirectory = Dictionary['ResultsDirectory']
    os.makedirs(ResultsDirectory,exist_ok=True)

    ## Set parameter map
    ParameterMapVector = sitk.VectorOfParameterMap()
    Dimension = len(FixedImage.shape)
    ImagePyramidSchedule = np.repeat(PyramidSchedule,Dimension)

    for Transformation in Transformations:

        ParameterMap = sitk.GetDefaultParameterMap(Transformation)
        ParameterMap['ResultImageFormat'] = ['mhd']
        ParameterMap['NewSamplesEveryIteration'] = ['true']
        ParameterMap['FixedImagePyramidSchedule'] = [str(ImagePyramidSchedule)[1:-1]]
        ParameterMap['MovingImagePyramidSchedule'] = [str(ImagePyramidSchedule)[1:-1]]
        ParameterMap['MaximumNumberOfIterations'] = [str(NIterations)]
        ParameterMap['SP_alpha'] = [str(Alpha)]
        ParameterMap['SP_A'] = [str(A)]

        ParameterMapVector.append(ParameterMap)


    ## Set Elastix and perform registration
    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetParameterMap(ParameterMapVector)
    ElastixImageFilter.SetFixedImage(sitk.GetImageFromArray(FixedImage))
    ElastixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage))
    ElastixImageFilter.SetFixedMask(sitk.GetImageFromArray(FixedMask))
    ElastixImageFilter.SetMovingMask(sitk.GetImageFromArray(MovingMask))
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
    ElastixImageFilter.LogToConsoleOn()
    ElastixImageFilter.Execute()

    ## Get results
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    return sitk.GetArrayFromImage(ResultImage), TransformParameterMap
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
Registration_Directory = os.path.join(WorkingDirectory,'04_Results/05_FractureLinePrediction/')
Mask_Directory = os.path.join(WorkingDirectory,'02_Data/03_uCT/')

## Set registration parameters
Transformation = 'rigid'
PyramidSchedules = [[50,20,10]]
NIterations = [2000]
Alphas = [0.6]
As = [10]


SampleList = [Dir for Dir in os.listdir(Registration_Directory) if os.path.isdir(Registration_Directory+Dir)]
SampleList.sort()

SampleList = ['443_L_73_F']

for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    # 02 Set Paths and Scans
    SamplePath = os.path.join(Registration_Directory,Sample)

    # 03 Load Scans
    uCT_Scan = sitk.ReadImage(SamplePath + '/uCT.mhd')
    uCT_Mask = sitk.ReadImage(SamplePath + '/uCT_Mask.mhd')
    HRpQCT_Scan = sitk.ReadImage(SamplePath + '/HR-pQCT.mhd')
    HRpQCT_Mask = sitk.ReadImage(SamplePath + '/HR-pQCT_Mask.mhd')
    Spacing = HRpQCT_Scan.GetSpacing()

    # 04 Get Array from Images
    uCT_Scan = sitk.GetArrayFromImage(uCT_Scan)
    uCT_Mask = sitk.GetArrayFromImage(uCT_Mask)
    HRpQCT_Scan = sitk.GetArrayFromImage(HRpQCT_Scan)
    HRpQCT_Mask = sitk.GetArrayFromImage(HRpQCT_Mask)


    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(uCT_Scan[:, :, int(uCT_Scan.shape[-1] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)

    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    # Axes.imshow(uCT_Mask[:, :, int(uCT_Mask.shape[-1] / 2)], cmap='bone')
    # plt.show()
    # plt.close(Figure)
    #
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    # Axes.imshow(HRpQCT_Scan[:, :, int(HRpQCT_Scan.shape[-1] / 2)], cmap='bone')
    # plt.show()
    # plt.close(Figure)
    #
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    # Axes.imshow(HRpQCT_Mask[:, :, int(HRpQCT_Mask.shape[-1] / 2)], cmap='bone')
    # plt.show()
    # plt.close(Figure)


    Dictionary = {'Transformations':Transformation,
                  'FixedImage':uCT_Scan,
                  'MovingImage':HRpQCT_Scan,
                  'FixedMask': uCT_Mask.astype('uint8'),
                  'MovingMask':HRpQCT_Mask.astype('uint8'),
                  'PyramidSchedule':PyramidSchedules[0],
                  'NIterations':NIterations[0],
                  'Alpha': Alphas[0],
                  'A': As[0],
                  'ResultsDirectory':SamplePath}
    ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(ResultImage[:, :, int(ResultImage.shape[-1] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)

    ResultImage = TransformixTransformations(HRpQCT_Scan, TransformParameterMap, ResultsDirectory=SamplePath)

    Origin = np.array([0, 0, 0])
    WriteMHD(ResultImage,Spacing,Origin,SamplePath, 'HR-pQCT_Registered2', PixelType='float')
