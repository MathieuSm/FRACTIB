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


# 01 Set variables
WorkingDirectory = os.getcwd()
DataDirectory = os.path.join(WorkingDirectory,'02_Data/05_hFE/01_AIMs/')
uCT_Directory = os.path.join(WorkingDirectory,'02_Data/03_uCT/')
Registration_Directory = os.path.join(WorkingDirectory,'04_Results/04_Registration/')


SampleList = os.listdir(DataDirectory)
SampleList.sort()

for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    # 02 Set Paths and Scans
    SamplePath = os.path.join(DataDirectory,Sample)
    Scans = [Scan for Scan in os.listdir(SamplePath) if Scan.endswith('.mhd')]
    Scans.sort()
    Scan = Scans[0][:8]

    uCT_SamplePath = os.path.join(uCT_Directory, Sample)
    uCT_Scan = [Scan for Scan in os.listdir(uCT_SamplePath) if Scan.endswith('FULLMASK.mhd')]

    Registration_Path = os.path.join(Registration_Directory,Sample)

    ResultsDirectory = os.path.join(WorkingDirectory, '04_Results/05_FractureLinePrediction', Sample)
    os.makedirs(ResultsDirectory, exist_ok=True)

    # 03 Load Masks
    uCT_Mask = sitk.ReadImage(uCT_SamplePath + '/' + uCT_Scan[0])

    HRpQCT_Cort_Mask = sitk.ReadImage(SamplePath + '/' + Scan + '_CORT_MASK_UNCOMP.mhd')
    HRpQCT_Trab_Mask = sitk.ReadImage(SamplePath + '/' + Scan + '_TRAB_MASK_UNCOMP.mhd')
    HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask


    # 04 Resample HR-pQCT image
    Offset = HRpQCT_Mask.GetOrigin()
    Direction = HRpQCT_Mask.GetDirection()
    Orig_Size = np.array(HRpQCT_Mask.GetSize(), dtype=np.int)
    Orig_Spacing = HRpQCT_Mask.GetSpacing()

    New_Spacing = uCT_Mask.GetSpacing()

    Resample = sitk.ResampleImageFilter()
    Resample.SetInterpolator = sitk.sitkLinear
    Resample.SetOutputDirection(Direction)
    Resample.SetOutputOrigin(Offset)
    Resample.SetOutputSpacing(New_Spacing)

    New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(New_Spacing))
    New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
    New_Size = [int(s) for s in New_Size]
    Resample.SetSize(New_Size)

    ResampledImage = Resample.Execute(HRpQCT_Mask)

    # 05 Get arrays from images and binarize masks
    uCT_Mask = sitk.GetArrayFromImage(uCT_Mask)
    HRpQCT_Mask = sitk.GetArrayFromImage(ResampledImage)
    uCT_Mask[uCT_Mask > 0] = 1
    HRpQCT_Mask[HRpQCT_Mask > 0] = 1

    ## Rotate HR-pQCT image to 180°
    HRpQCT_Mask = np.rot90(HRpQCT_Mask, 2, (0, 1))

    HRpQCT_Slice = HRpQCT_Mask[int(HRpQCT_Mask.shape[0] / 2), :, :]
    uCT_Slice = uCT_Mask[int(uCT_Mask.shape[0] / 2), :, :]


    # 07 Find best HR-pQCT image initial rotation
    Rotations = range(0,360,15)
    DSCs = pd.DataFrame()
    for Rotation in Rotations:
        Rotated_Slice = ndimage.rotate(HRpQCT_Slice,Rotation,reshape=False)

        ## Register images
        ParameterMap = sitk.GetDefaultParameterMap('translation')

        ## Set Elastix and perform registration
        ElastixImageFilter = sitk.ElastixImageFilter()
        ElastixImageFilter.SetParameterMap(ParameterMap)
        ElastixImageFilter.SetFixedImage(sitk.GetImageFromArray(uCT_Slice))
        ElastixImageFilter.SetMovingImage(sitk.GetImageFromArray(Rotated_Slice))
        ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
        ElastixImageFilter.LogToConsoleOn()
        ElastixImageFilter.Execute()

        ## Get results
        ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
        ResultImage = sitk.GetArrayFromImage(ResultImage)

        DSC = 2 * np.sum(ResultImage * uCT_Slice) / np.sum(ResultImage + uCT_Slice)
        DSCs = DSCs.append({'Angle':Rotation,'DSC':DSC},ignore_index=True)

        # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
        # Axes.imshow(uCT_Slice, cmap='bone', alpha=0.5)
        # Axes.imshow(ResultImage, cmap='bone', alpha=0.5)
        # Axes.set_title('DSC: ' + str(DSC.round(2)) + ' Angle: ' + str(int(DSC)))
        # plt.show()
        # plt.close(Figure)

    # 08 Load HR-pQCT scan, rotate it and write it
    HRpQCT_Scan = sitk.ReadImage(SamplePath + '/' + Scan + '_UNCOMP.mhd')
    ResampledImage = Resample.Execute(HRpQCT_Scan)
    HRpQCT_Scan = sitk.GetArrayFromImage(ResampledImage)
    HRpQCT_Scan = np.rot90(HRpQCT_Scan, 2, (0, 1))

    # ## Add Layers to HR-pQCT Scan
    # ZerosArray = np.zeros(uCT_Mask.shape)
    # Diff_Z, Diff_Y, Diff_X = np.array(uCT_Mask.shape) - np.array(HRpQCT_Scan.shape)
    # for Z in range(min(HRpQCT_Scan.shape[0],uCT_Mask.shape[0])):
    #     for Y in range(min(HRpQCT_Scan.shape[1],uCT_Mask.shape[1])):
    #         for X in range(min(HRpQCT_Scan.shape[2],uCT_Mask.shape[2])):
    #             ZerosArray[Z + int(Diff_Z/2), Y + int(Diff_Y/2), X + int(Diff_X/2)] = HRpQCT_Scan[Z,Y,X]
    # HRpQCT_Scan = ZerosArray

    # Verify best angle
    BestAngle = DSCs.loc[DSCs['DSC'].idxmax(), 'Angle']
    Rotated_Slice = ndimage.rotate(HRpQCT_Slice, BestAngle, reshape=False)
    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetMovingImage(sitk.GetImageFromArray(Rotated_Slice))
    ElastixImageFilter.Execute()
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    ResultImage = sitk.GetArrayFromImage(ResultImage)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(uCT_Slice, cmap='bone', alpha=0.5)
    Axes.imshow(ResultImage, cmap='bone', alpha=0.5)
    Axes.set_title('DSC: ' + str(DSC.round(2)) + ' Angle: ' + str(int(DSC)))
    plt.show()
    plt.close(Figure)


    File = open(ResultsDirectory + '/HR-pQCT_RigidRotationAngle.txt','+w')
    File.write(str(int(BestAngle)))
    File.close()
    Rotated_Scan = ndimage.rotate(HRpQCT_Scan, int(BestAngle), (1,2), reshape=False)
    Origin = np.array([0, 0, 0])
    WriteMHD(Rotated_Scan, New_Spacing, Origin, ResultsDirectory, 'HR-pQCT', PixelType='float')
    Rotated_Mask = ndimage.rotate(HRpQCT_Mask, int(BestAngle), (1,2), reshape=False)
    WriteMHD(Rotated_Mask, New_Spacing, Origin, ResultsDirectory, 'HR-pQCT_Mask', PixelType='float')

    # 09 Load uCT scan, crop it and write it
    uCT_Scan = sitk.ReadImage(uCT_SamplePath + '/' + uCT_Scan[0][:-13] + '.mhd')
    uCT_Scan = sitk.GetArrayFromImage(uCT_Scan)
    WriteMHD(uCT_Scan, New_Spacing, Origin, ResultsDirectory, 'uCT', PixelType='float')
    WriteMHD(uCT_Mask, New_Spacing, Origin, ResultsDirectory, 'uCT_Mask', PixelType='float')




# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
# Axes.imshow(uCT_Mask[:,:,int(uCT_Mask.shape[2]/2)], cmap='bone', alpha=0.5)
# Axes.imshow(uCT_Scan[:,:,int(uCT_Scan.shape[2]/2)], cmap='bone', alpha=0.5)
# plt.show()
# plt.close(Figure)