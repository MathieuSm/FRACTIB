#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

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
    Dimension = FixedImage.GetDimension()
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
    ElastixImageFilter.SetFixedImage(FixedImage)
    ElastixImageFilter.SetMovingImage(MovingImage)
    ElastixImageFilter.SetFixedMask(FixedMask)
    ElastixImageFilter.SetMovingMask(MovingMask)
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
    ElastixImageFilter.SetInitialTransformParameterFileName(ResultsDirectory + 'InitialTranslation.txt')
    ElastixImageFilter.LogToConsoleOn()
    ElastixImageFilter.Execute()

    ## Get results
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    ## Rename parameter file
    os.rename(Results_Path + '/TransformParameters.0.txt', ResultsDirectory + 'TransformParameters.txt')

    return ResultImage, TransformParameterMap
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
uCTFolder = os.path.join(WorkingDirectory, '02_Data/03_uCT/')
HRpQCTFolder = os.path.join(WorkingDirectory, '02_Data/02_HRpQCT/')
ResultsFolder = os.path.join(WorkingDirectory, '04_Results/06_FractureLinePrediction/')


SampleList = [Dir for Dir in os.listdir(HRpQCTFolder) if os.path.isdir(HRpQCTFolder+Dir)]
SampleList.sort()

Index = 0
for Index in range(len(SampleList)):

    ## Set sample paths
    Sample = SampleList[Index]
    uCT_Path = os.path.join(uCTFolder + Sample + '/')
    HRpQCT_Path = os.path.join(HRpQCTFolder + Sample + '/')
    Results_Path = os.path.join(ResultsFolder + Sample + '/')

    ## List files
    uCT_Files = os.listdir(uCT_Path)
    uCT_Mask_File = [File for File in uCT_Files if File.endswith('MASK.mhd')][0]
    uCT_Scan_File = uCT_Mask_File[:-13] + '.mhd'
    HRpQCT_Files = [File for File in os.listdir(HRpQCT_Path) if File.endswith('.mhd')]
    HRpQCT_Files.sort()

    ## Load Scans
    uCT_Scan = sitk.ReadImage(uCT_Path + uCT_Scan_File)
    uCT_Mask = sitk.ReadImage(uCT_Path + uCT_Mask_File)
    uCT_Scan_Array = sitk.GetArrayFromImage(uCT_Scan)

    HRpQCT_Scan = sitk.ReadImage(HRpQCT_Path + HRpQCT_Files[3])
    HRpQCT_Cort = sitk.ReadImage(HRpQCT_Path + HRpQCT_Files[0])
    HRpQCT_Trab = sitk.ReadImage(HRpQCT_Path + HRpQCT_Files[2])
    HRpQCT_Mask = HRpQCT_Cort + HRpQCT_Trab

    ## Apply initial rotations on HRpQCT scans
    DSCs = pd.read_csv(Results_Path+'DSCs.csv')
    BestAngle = int(DSCs.loc[DSCs['DSC'].idxmax(), 'Angle'])
    for Index, Image in enumerate([HRpQCT_Scan,HRpQCT_Mask]):
        print('\n Image number:' + str(Index))
        print(Image)
        Spacing = Image.GetSpacing()
        Origin = Image.GetOrigin()
        Direction = Image.GetDirection()
        Array = sitk.GetArrayFromImage(Image)
        Array = np.rot90(Array, 2, (0, 1))
        Array = ndimage.rotate(Array, BestAngle, (1, 2), reshape=True)

        if Index == 1:
            Otsu_Filter = sitk.OtsuThresholdImageFilter()
            Otsu_Filter.SetInsideValue(0)
            Otsu_Filter.SetOutsideValue(1)
            Segmentation = Otsu_Filter.Execute(Image)
            R_Threshold = Otsu_Filter.GetThreshold()
            Array[Array <= R_Threshold] = 0
            Array[Array > 0] = 1

        Image = sitk.GetImageFromArray(Array)
        Image.SetSpacing(Spacing)
        Image.SetOrigin(Origin)
        Image.SetDirection(Direction)

        if Index == 0:
            HRpQCT_Scan = Image
        elif Index == 1:
            HRpQCT_Mask = Image

    ## Crop HRpQCT scan for lower registration memory
    HRpQCT_Scan_Cropped = sitk.Crop(HRpQCT_Scan,(0,0,322),(0,0,0))
    HRpQCT_Mask_Cropped = sitk.Crop(HRpQCT_Mask,(0,0,322),(0,0,0))

    ## Perform registration
    Dictionary = {'Transformations': 'rigid',
                  'FixedImage': uCT_Scan,
                  'MovingImage': HRpQCT_Scan_Cropped,
                  'FixedMask': sitk.Cast(uCT_Mask, sitk.sitkUInt8),
                  'MovingMask': sitk.Cast(HRpQCT_Mask_Cropped, sitk.sitkUInt8),
                  'PyramidSchedule': [50,20,10],
                  'NIterations': 2000,
                  'Alpha': 0.6,
                  'A': 10,
                  'ResultsDirectory': Results_Path}
    ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)
    ResultImage_Array = sitk.GetArrayFromImage(ResultImage)

    ## Plot registration results
    Figure, Axes = plt.subplots(1, 1, figsize=(8.25, 6.75), dpi=100)
    Mid_Position = int(round(uCT_Scan_Array.shape[2]/2))
    uCT_Show = Axes.imshow(uCT_Scan_Array[:, :, Mid_Position], cmap='bone', alpha=1)
    HRpQCT_Show = Axes.imshow(ResultImage_Array[:, :, Mid_Position], cmap='jet', alpha=0.5)

    SliderAxis = plt.axes([0.25, 0.15, 0.65, 0.03])
    Mid_Position_Slider = Slider(SliderAxis, 'Y Position', 0, uCT_Scan_Array.shape[1], valinit=Mid_Position)

    def Update(Value):
        Position = Mid_Position_Slider.Value
        uCT_Show.set_data(uCT_Scan_Array[:, :, int(Position)])
        HRpQCT_Show.set_data(ResultImage_Array[:, :, int(Position)])
        Figure.canvas.draw_idle()

    Mid_Position_Slider.on_changed(Update)

    Axes.set_title(Sample + ' Registration Results')
    plt.show()
    plt.close(Figure)
