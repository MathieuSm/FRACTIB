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


def PlotsImages(Image1, Image2, C1=np.array([[0, 0, 0], [255, 0, 0]]), C2=np.array([[0, 0, 0], [0, 0, 255]]), Cbar=False):
    # Create custom color map
    import matplotlib as mpl  # in python
    ColorMap1 = mpl.colors.ListedColormap(C1 / 255.0)
    ColorMap2 = mpl.colors.ListedColormap(C2 / 255.0)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(Image1[:, :, int(Image1.shape[2] / 2)], cmap=ColorMap1, alpha=0.5)
    ColorMap = Axes.imshow(Image2[:, :, int(Image2.shape[2] / 2)], cmap=ColorMap2, alpha=0.5)
    if Cbar:
        plt.colorbar(ColorMap)
    Axes.axis('off')
    plt.show()
    plt.close(Figure)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(Image1[int(Image1.shape[0] / 2), :, :], cmap=ColorMap1, alpha=0.5)
    ColorMap = Axes.imshow(Image2[int(Image2.shape[0] / 2), :, :], cmap=ColorMap2, alpha=0.5)
    if Cbar:
        plt.colorbar(ColorMap)
    Axes.axis('off')
    plt.show()
    plt.close(Figure)

    return


# 01 Set variables
WorkingDirectory = os.getcwd()
Data_Directory = os.path.join(WorkingDirectory,'04_Results/05_FractureLinePrediction/')

SampleList = [Dir for Dir in os.listdir(Data_Directory) if os.path.isdir(Data_Directory+Dir)]
SampleList.sort()

Index = 0
# for Index in range(len(SampleList)):

Sample = SampleList[Index]
SamplePath = os.path.join(Data_Directory,Sample)

# Read images and get arrays
uCT = sitk.ReadImage((SamplePath + '/uCT.mhd'))
uCT_Array = sitk.GetArrayFromImage(uCT)
HRpQCT = sitk.ReadImage((SamplePath + '/HR-pQCT.mhd'))
HRpQCT_Array = sitk.GetArrayFromImage(HRpQCT)
HRpQCT_R = sitk.ReadImage((SamplePath + '/HR-pQCT_Registered.mhd'))
HRpQCT_R_Array = sitk.GetArrayFromImage(HRpQCT_R)
uCT_C = sitk.ReadImage((SamplePath + '/uCT_Cropped.mhd'))
uCT_C_Array = sitk.GetArrayFromImage(uCT_C)
HRpQCT_C = sitk.ReadImage((SamplePath + '/HRpQCT_Cropped.mhd'))
HRpQCT_C_Array = sitk.GetArrayFromImage(HRpQCT_C)
J_C = sitk.ReadImage((SamplePath + '/J_Cropped.mhd'))
J_C_Array = sitk.GetArrayFromImage(J_C)
J_R = sitk.ReadImage((SamplePath + '/J_Registration.mhd'))
J_R_Array = sitk.GetArrayFromImage(J_R)

# Segment it
Otsu_Filter = sitk.OtsuThresholdImageFilter()
Otsu_Filter.SetInsideValue(0)
Otsu_Filter.SetOutsideValue(1)
Segmentation = Otsu_Filter.Execute(uCT)
uCT_Threshold = Otsu_Filter.GetThreshold()
uCT_Array[uCT_Array < uCT_Threshold] = 0
uCT_Array[uCT_Array > 0] = 1
uCT_C_Array[uCT_C_Array < uCT_Threshold] = 0
uCT_C_Array[uCT_C_Array > 0] = 1

HRpQCT_Threshold = 3000
HRpQCT_Array[HRpQCT_Array < HRpQCT_Threshold] = 0
HRpQCT_Array[HRpQCT_Array > 0] = 1
HRpQCT_R_Array[HRpQCT_R_Array < HRpQCT_Threshold] = 0
HRpQCT_R_Array[HRpQCT_R_Array > 0] = 1
HRpQCT_C_Array[HRpQCT_C_Array < HRpQCT_Threshold] = 0
HRpQCT_C_Array[HRpQCT_C_Array > 0] = 1


# Plot
PlotsImages(uCT_Array,HRpQCT_Array)
PlotsImages(uCT_Array,HRpQCT_R_Array)
PlotsImages(uCT_C_Array,HRpQCT_C_Array)

J_C_Array[J_C_Array < 0.1] = np.nan
J_R_Array[J_R_Array < 0.1] = np.nan

C = np.linspace(0,255,256)
C = np.append(C,255)
Colormap = np.zeros((len(C),3))
for Step in range(len(C)):
    Colormap[Step,0] = C[::-1][Step]
    Colormap[Step, 2] = C[Step]
Colormap = Colormap.astype('int')
PlotsImages(J_C_Array,J_R_Array,C1=Colormap,C2=Colormap,Cbar=True)
