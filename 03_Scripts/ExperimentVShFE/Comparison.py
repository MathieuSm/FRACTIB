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
Data_Directory = os.path.join(WorkingDirectory,'04_Results/05_FractureLinePrediction/')

SampleList = [Dir for Dir in os.listdir(Data_Directory) if os.path.isdir(Data_Directory+Dir)]
SampleList.sort()


for Index in range(len(SampleList)):

    Sample = SampleList[Index]

    # 02 Set Paths and Scans
    SamplePath = os.path.join(Data_Directory,Sample)


## Do some tests
uCT_Mask = sitk.ReadImage((SamplePath + '/uCT_Mask_Cropped.mhd'))
J_Cropped = sitk.ReadImage(SamplePath + '/J_Cropped.mhd')
J_Registration = sitk.ReadImage(WorkingDirectory + '/04_Results/04_Registration/432_L_77_F/J.mhd')


Test = J_Cropped*uCT_Mask
TestArray = sitk.GetArrayFromImage(Test)

# Crop J registration with transverse plane
J_Array = sitk.GetArrayFromImage(J_Registration)
J_Zeros = np.zeros(uCT_Mask.GetSize()[::-1])

Shift = int(J_Cropped.GetOrigin()[2] / J_Cropped.GetSpacing()[2])
for i in range(J_Cropped.GetSize()[2]):
    for j in range(J_Array.shape[1]):
        for k in range(J_Array.shape[2]):
            J_Zeros[i,j,k] += J_Array[Shift + i, j, k]

TestArray2 = J_Zeros * sitk.GetArrayFromImage(uCT_Mask)

# Plot
Vmin = min(TestArray.min(),TestArray2.min())
Vmax = max(TestArray.max(),TestArray2.max())

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Cmap = Axes.imshow(TestArray[:, :, int(TestArray.shape[2]/2)], cmap='jet', vmin=Vmin, vmax=Vmax)
plt.colorbar(Cmap)
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Cmap = Axes.imshow(TestArray2[:, :, int(TestArray2.shape[2]/2)], cmap='jet', vmin=Vmin, vmax=Vmax)
plt.colorbar(Cmap)
plt.show()
plt.close(Figure)

Origin = np.array(uCT_Mask.GetOrigin()[::-1])
Spacing = np.array(uCT_Mask.GetSpacing()[::-1])

WriteMHD(TestArray, Spacing, Origin, SamplePath, 'J_hFE', PixelType='float')
WriteMHD(TestArray2, Spacing, Origin, SamplePath, 'J_Registration', PixelType='float')
