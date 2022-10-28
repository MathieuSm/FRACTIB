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
Data_Directory = os.path.join(WorkingDirectory,'04_Results/06_FractureLinePrediction/')

SampleList = [Dir for Dir in os.listdir(Data_Directory) if os.path.isdir(Data_Directory+Dir)]
SampleList.sort()


for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    # 02 Set Paths and Scans
    SamplePath = os.path.join(Data_Directory,Sample)

    # 03 Load Scans
    uCT_Scan = sitk.ReadImage(SamplePath + '/uCT.mhd')
    uCT_Mask = sitk.ReadImage(SamplePath + '/uCT_Mask.mhd')
    HRpQCT_Scan = sitk.ReadImage(SamplePath + '/HR-pQCT_Registered.mhd')
    HRpQCT_Mask_Image = sitk.ReadImage(SamplePath + '/HR-pQCT_Mask.mhd')
    J = sitk.ReadImage(SamplePath + '/J.mhd')
    F_Tilde = sitk.ReadImage(SamplePath + '/F_Tilde.mhd')


    # Get arrays from images
    uCT_Mask_Array = sitk.GetArrayFromImage(uCT_Mask)
    HRpQCT_Mask_Array = sitk.GetArrayFromImage(HRpQCT_Mask_Image)

    Spacing = uCT_Scan.GetSpacing()

    # 04 Transform HR-pQCT mask
    TransformParameterMap = sitk.ReadParameterFile(SamplePath + '/TransformParameters.0.txt')
    HRpQCT_Mask = TransformixTransformations(HRpQCT_Mask_Array, TransformParameterMap, ResultsDirectory=SamplePath)
    ## Re-binarize mask
    HRpQCT_Mask[HRpQCT_Mask < 0.5] = 0
    HRpQCT_Mask[HRpQCT_Mask >= 0.5] = 1
    HRpQCT_Mask = sitk.GetImageFromArray(HRpQCT_Mask)
    HRpQCT_Mask.SetSpacing(HRpQCT_Mask_Image.GetSpacing())

    # 06 Resample masks to get same sampling as HR-pQCT
    Offset = HRpQCT_Mask_Image.GetOrigin()
    Direction = HRpQCT_Mask_Image.GetDirection()
    Orig_Size = np.array(HRpQCT_Mask.GetSize(), dtype=np.int)
    Orig_Spacing = HRpQCT_Mask_Image.GetSpacing()

    New_Spacing = (0.9712, 0.9712, 0.9712)

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
    HRpQCT_Mask = sitk.GetArrayFromImage(ResampledImage)

    # Do the same fot other images
    Orig_Size = np.array(uCT_Mask.GetSize(), dtype=np.int)

    New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(New_Spacing))
    New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
    New_Size = [int(s) for s in New_Size]
    Resample.SetSize(New_Size)

    ResampledImage = Resample.Execute(uCT_Mask)
    uCT_Mask = sitk.GetArrayFromImage(ResampledImage)

    ResampledImage = Resample.Execute(uCT_Scan)
    uCT_Scan = sitk.GetArrayFromImage(ResampledImage)

    ResampledImage = Resample.Execute(HRpQCT_Scan)
    HRpQCT_Scan = sitk.GetArrayFromImage(ResampledImage)

    ResampledImage = Resample.Execute(J)
    J = sitk.GetArrayFromImage(ResampledImage)

    ResampledImage = Resample.Execute(F_Tilde)
    F_Tilde = sitk.GetArrayFromImage(ResampledImage)




    # 05 Compute first slice to keep
    DSCs = pd.DataFrame()
    for Slice in range(10):

        if np.sum(uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :]) > 0:
            DSC = 2 * np.sum(uCT_Mask[Slice, :, :] * HRpQCT_Mask[Slice, :, :]) / np.sum(uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :])
            DSCs = DSCs.append({'Slice':Slice,'DSC':DSC},ignore_index=True)

    DSCs['Slope'] = 0
    for i in DSCs.index:
        if i > 0:
            DSCs.loc[i,'Slope'] = DSCs.loc[i,'DSC'] - DSCs.loc[i-1,'DSC']

    Tolerance = 1E-2
    FirstSlice = DSCs.loc[DSCs[DSCs['Slope'] > Tolerance]['Slice'].idxmax(),'Slice'].astype('int')

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(uCT_Mask[FirstSlice-1, :, :], cmap='bone',alpha=0.5)
    Axes.imshow(HRpQCT_Mask[FirstSlice-1, :, :], cmap='bone',alpha=0.5)
    plt.show()
    plt.close(Figure)

    # 06 Compute last slice to keep
    DSCs = pd.DataFrame()
    for SliceIndex in range(10):

        Slice = uCT_Mask.shape[0]-SliceIndex-1

        if np.sum(uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :]) > 0:
            DSC = 2 * np.sum(uCT_Mask[Slice, :, :] * HRpQCT_Mask[Slice, :, :]) / np.sum(
                uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :])
            DSCs = DSCs.append({'Slice': Slice, 'DSC': DSC}, ignore_index=True)

    DSCs['Slope'] = 0
    for i in DSCs.index:
        if i > 0:
            DSCs.loc[i, 'Slope'] = DSCs.loc[i, 'DSC'] - DSCs.loc[i - 1, 'DSC']

    Tolerance = 1E-2
    LastSlice = DSCs.loc[DSCs[DSCs['Slope'] > Tolerance]['Slice'].idxmin(),'Slice'].astype('int')

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(uCT_Mask[LastSlice+1, :, :], cmap='bone', alpha=0.5)
    Axes.imshow(HRpQCT_Mask[LastSlice+1, :, :], cmap='bone', alpha=0.5)
    plt.show()
    plt.close(Figure)

    # Crop images
    uCT_Cropped = uCT_Scan[FirstSlice:LastSlice, :, :]
    HRpQCT_Cropped = HRpQCT_Scan[FirstSlice:LastSlice, :, :]

    uMask_Cropped = uCT_Mask[FirstSlice:LastSlice, :, :]
    HRMask_Cropped = HRpQCT_Mask[FirstSlice:LastSlice, :, :]

    J_Cropped = J[FirstSlice:LastSlice, :, :]
    F_Tilde_Cropped = F_Tilde[FirstSlice:LastSlice, :, :]


    Origin = np.array([0,0,New_Spacing[2]]) * FirstSlice
    WriteMHD(uCT_Cropped,New_Spacing,Origin,SamplePath, 'uCT_Cropped', PixelType='float')
    WriteMHD(HRpQCT_Cropped,New_Spacing,Origin,SamplePath, 'HRpQCT_Cropped', PixelType='float')
    WriteMHD(uMask_Cropped,New_Spacing,Origin,SamplePath, 'uCT_Mask_Cropped', PixelType='float')
    WriteMHD(HRMask_Cropped,New_Spacing,Origin,SamplePath, 'HRpQCT_Mask_Cropped', PixelType='float')
    WriteMHD(J_Cropped,New_Spacing,Origin,SamplePath, 'J_Cropped', PixelType='float')
    WriteMHD(F_Tilde_Cropped,New_Spacing,Origin,SamplePath, 'F_Tilde_Cropped', PixelType='float')



    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(uCT_Cropped[:, :, int(uCT_Mask.shape[2]/2)], cmap='bone',alpha=0.5)
    Axes.imshow(HRpQCT_Scan[FirstSlice:LastSlice, :, int(HRpQCT_Mask.shape[2]/2)], cmap='bone',alpha=0.5)
    plt.show()
    plt.close(Figure)


    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.plot(DSCs['Slope'], linestyle='--', marker='o', color=(1,0,0))
    plt.show()
    plt.close(Figure)



