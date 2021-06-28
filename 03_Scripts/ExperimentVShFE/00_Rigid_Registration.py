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
DataDirectory = os.path.join(WorkingDirectory,'02_Data/05_hFE/01_AIMs/')
uCT_Directory = os.path.join(WorkingDirectory,'02_Data/03_uCT/')
Registration_Directory = os.path.join(WorkingDirectory,'04_Results/04_Registration/')


## Set registration parameters
Transformation = 'rigid'
PyramidSchedules = [[50,20,10]]
NIterations = [2000]
Alphas = [0.6]
As = [10]


SampleList = os.listdir(DataDirectory)
SampleList.sort()


for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    SamplePath = os.path.join(DataDirectory,Sample)
    Scans = [Scan for Scan in os.listdir(SamplePath) if Scan.endswith('.mhd')]
    Scans.sort()

    # 03 Load HR-pQCT files
    CortMask = sitk.ReadImage(SamplePath + '/' + Scans[0][:8] + '_CORT_MASK_UNCOMP.mhd')
    TrabMask = sitk.ReadImage(SamplePath + '/' + Scans[0][:8] + '_TRAB_MASK_UNCOMP.mhd')
    FixedImage = sitk.ReadImage(SamplePath + '/' + Scans[0][:8] + '_UNCOMP.mhd')
    Offset = FixedImage.GetOrigin()
    Direction = FixedImage.GetDirection()

    ## Build mask
    CortMask = sitk.GetArrayFromImage(CortMask).astype('uint8')
    TrabMask = sitk.GetArrayFromImage(TrabMask).astype('uint8')
    FullMask = CortMask + TrabMask
    FullMask[FullMask > 0.1] = 1
    FullMask = sitk.GetImageFromArray(FullMask)
    FullMask.SetSpacing(FixedImage.GetSpacing())
    FullMask.SetOrigin(Offset)
    FullMask.SetDirection(Direction)

    # 04 Load uCT file
    uCT_SamplePath = os.path.join(uCT_Directory, Sample)
    uCT_Scans = [Scan for Scan in os.listdir(uCT_SamplePath) if Scan.endswith('DOWNSCALED.mhd')]
    uCT_Scans.sort()

    MovingImage = sitk.ReadImage(uCT_SamplePath + '/' + uCT_Scans[0])
    Spacing = MovingImage.GetSpacing()

    ## Resample image
    Resample = sitk.ResampleImageFilter()
    Resample.SetInterpolator = sitk.sitkLinear
    Resample.SetOutputDirection(Direction)
    Resample.SetOutputOrigin(Offset)
    Resample.SetOutputSpacing(Spacing)

    Orig_Size = np.array(FixedImage.GetSize(), dtype=np.int)
    Orig_Spacing = FixedImage.GetSpacing()
    New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(Spacing))
    New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
    New_Size = [int(s) for s in New_Size]
    Resample.SetSize(New_Size)

    ResampledImage = Resample.Execute(FixedImage)
    ResampledMask = Resample.Execute(FullMask)


    ## Get arrays from images
    MovingImage = sitk.GetArrayFromImage(MovingImage)
    FixedImage = sitk.GetArrayFromImage(ResampledImage)
    FixedMask = sitk.GetArrayFromImage(ResampledMask)


    ## Binarize mask
    FixedMask[FixedMask < np.max(FixedMask)/2] = 0
    FixedMask[FixedMask > np.max(FixedMask)/2] = 1

    ## Add layers to fixed image
    ZerosArray = np.zeros(MovingImage.shape)
    ZerosMask = np.zeros(MovingImage.shape).astype('uint8')
    Diff_Z, Diff_Y, Diff_X = np.array(MovingImage.shape) - np.array(FixedImage.shape)
    for Z in range(min(FixedImage.shape[0],MovingImage.shape[0])):
        for Y in range(min(FixedImage.shape[1],MovingImage.shape[1])):
            for X in range(min(FixedImage.shape[2],MovingImage.shape[2])):
                ZerosArray[Z + int(Diff_Z/2), Y + int(Diff_Y/2), X + int(Diff_X/2)] = FixedImage[Z,Y,X]
                ZerosMask[Z + int(Diff_Z/2), Y + int(Diff_Y/2), X + int(Diff_X/2)] = FixedMask[Z,Y,X]
    FixedImage = ZerosArray
    FixedMask = ZerosMask

    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    # Axes.imshow(MovingImage[:,:,int(MovingImage.shape[-1]/2)],cmap='bone')
    # plt.show()
    # plt.close(Figure)
    #
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    # Axes.imshow(FixedImage[:, :, int(FixedImage.shape[-1] / 2)], cmap='bone')
    # plt.show()
    # plt.close(Figure)

    #Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    #Axes.imshow(FixedMask[:, :, int(FixedMask.shape[-1] / 2)], cmap='bone')
    #plt.show()
    #plt.close(Figure)

    Data = pd.DataFrame()

    ResultsDirectory = os.path.join(WorkingDirectory, '04_Results/05_FractureLinePrediction', Sample)
    os.makedirs(ResultsDirectory,exist_ok=True)
    SampleData = {'Sample': Sample}

    ## Rotate image by 180 degree around x
    MovingImage = np.rot90(MovingImage, 2, (0, 1))
    InitialParameterMap = sitk.ReadParameterFile(ResultsDirectory + '/GuessedParameters.txt')
    MovingImageGuess = TransformixTransformations(MovingImage, InitialParameterMap, ResultsDirectory)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(FixedImage[:, :, int(FixedImage.shape[-1] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(MovingImageGuess[:, :, int(MovingImageGuess.shape[-1] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)


    Dictionary = {'Transformations':Transformation,
                  'FixedImage':FixedImage,
                  'MovingImage':MovingImageGuess,
                  'FixedMask': FixedMask,
                  'PyramidSchedule':PyramidSchedules[0],
                  'NIterations':NIterations[0],
                  'Alpha': Alphas[0],
                  'A': As[0],
                  'ResultsDirectory':ResultsDirectory}
    Tic = time.clock_gettime(0)
    ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)
    Tac = time.clock_gettime(0)
    SampleData['Rigid'] = round(Tac-Tic)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(ResultImage[:, :, int(ResultImage.shape[-1] / 2)], cmap='bone')
    plt.show()
    plt.close(Figure)

    Origin = Offset - np.array([int(Diff_X/2), int(Diff_Y/2), int(Diff_Z/2)]) * Spacing
    WriteMHD(ResultImage,Spacing,Origin,ResultsDirectory, 'uCT_Registered', PixelType='float')
    
    
    ## Compute registered J and F_tilde according to the original cropped image
    Crop = [[[0, 342], [25, 500], [65, 485]],
            [[0, 326], [75, 540], [75, 500]],
            [[0, 330], [10, 500], [55, 545]],
            [[0, 330], [75, 485], [35, 525]],
            [[0, 311], [140, 605], [100, 570]],
            [[0, 325], [70, 465], [35, 505]],
            [[0, 320], [25, 549], [55, 530]],
            [[0, 325], [35, 549], [40, 525]],
            [[0, 310], [100, 630], [40, 610]],
            [[0, 325], [65, 535], [0, 549]],
            [[10, 330], [50, 475], [50, 549]],
            [[0, 337], [100, 570], [150, 600]],
            [[0, 311], [125, 585], [130, 610]],
            [[0, 311], [100, 540], [105, 580]],
            [[10, 335], [95, 515], [40, 515]],
            [[0, 330], [20, 549], [0, 530]],
            [[0, 311], [20, 630], [100, 640]],
            [[0, 330], [170, 580], [80, 570]],
            [[0, 323], [15, 510], [5, 525]],
            [[0, 325], [50, 510], [80, 549]],
            [[0, 325], [60, 460], [30, 530]],
            [[0, 325], [45, 520], [0, 549]],
            [[0, 325], [105, 500], [50, 515]],
            [[0, 325], [0, 535], [85, 520]],
            [[0, 330], [95, 500], [5, 505]]]
    Crop_Index = Crop[Index]
    CroppedImage = MovingImage[Crop_Index[0][0]:Crop_Index[0][1],
                               Crop_Index[1][0]:Crop_Index[1][1],
                               Crop_Index[2][0]:Crop_Index[2][1]]

    RegisteredSample = os.path.join(Registration_Directory,Sample)
    Image2Move = sitk.ReadImage(RegisteredSample + '/J.mhd')
    Image2Move = sitk.GetArrayFromImage(Image2Move)

    # Rotate image
    Image2Move = np.rot90(Image2Move, 2, (0, 1))

    ## Resample image
    Repeat = np.array([10,10,10])
    Image2Move = np.repeat(Image2Move, Repeat[2],axis=2)
    Image2Move = np.repeat(Image2Move, Repeat[1], axis=1)
    Image2Move = np.repeat(Image2Move, Repeat[0], axis=0)
    Image2Move = Image2Move[:CroppedImage.shape[0],:CroppedImage.shape[1],:CroppedImage.shape[2]]

    J_Array = np.ones(MovingImage.shape)
    z = 0
    for Z in range(Crop_Index[0][0],min(Crop_Index[0][1],Image2Move.shape[0])):
        y = 0
        for Y in range(Crop_Index[1][0],min(Crop_Index[1][1],Image2Move.shape[1])):
            x = 0
            for X in range(Crop_Index[2][0],min(Crop_Index[2][1],Image2Move.shape[2])):
                J_Array[Z,Y,X] = Image2Move[z,y,x]
                x += 1
            y += 1
        z += 1
    Image2Move = J_Array

    Image2Move = TransformixTransformations(Image2Move, InitialParameterMap, ResultsDirectory)
    ImageMoved = TransformixTransformations(Image2Move, TransformParameterMap, ResultsDirectory)
    MaskedImage = ImageMoved * FixedMask
    MaskedImage[MaskedImage==0] = np.nan
    WriteMHD(ImageMoved, Spacing, Origin, ResultsDirectory, 'J_Registered', PixelType='float')

    Image2Move = sitk.ReadImage(RegisteredSample + '/F_Tilde.mhd')
    Image2Move = sitk.GetArrayFromImage(Image2Move)

    # Rotate image
    Image2Move = np.rot90(Image2Move, 2, (0, 1))
    
    ## Resample image
    Repeat = np.array([10,10,10])
    Image2Move = np.repeat(Image2Move, Repeat[2],axis=2)
    Image2Move = np.repeat(Image2Move, Repeat[1], axis=1)
    Image2Move = np.repeat(Image2Move, Repeat[0], axis=0)
    Image2Move = Image2Move[:CroppedImage.shape[0],:CroppedImage.shape[1],:CroppedImage.shape[2]]

    J_Array = np.ones(MovingImage.shape)
    z = 0
    for Z in range(Crop_Index[0][0],min(Crop_Index[0][1],Image2Move.shape[0])):
        y = 0
        for Y in range(Crop_Index[1][0],min(Crop_Index[1][1],Image2Move.shape[1])):
            x = 0
            for X in range(Crop_Index[2][0],min(Crop_Index[2][1],Image2Move.shape[2])):
                J_Array[Z,Y,X] = Image2Move[z,y,x]
                x += 1
            y += 1
        z += 1
    Image2Move = J_Array
    
    Image2Move = TransformixTransformations(Image2Move, InitialParameterMap, ResultsDirectory)
    ImageMoved = TransformixTransformations(Image2Move,TransformParameterMap,ResultsDirectory)
    MaskedImage = ImageMoved * FixedMask
    MaskedImage[MaskedImage==0] = np.nan
    WriteMHD(ImageMoved,Spacing,Origin,ResultsDirectory,'F_Tilde_Registered', PixelType='float')

