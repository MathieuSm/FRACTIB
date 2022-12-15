#%% #!/usr/bin/env python3
import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path

from Utils import *

def Adjust_Image_Size(Image, CoarseFactor, CropZ='Crop'):

    """
    Adapted from Denis's utils_SA.py
    Images are adjusted according to CropType:
    0 = CropType.expand     (Expand image by copying layers)
    1 = CropType.crop       (Crop image)
    2 = CropType.variable   (Either crop or expand, depending on what includes less layers)
    """

    # Get array
    Array = sitk.GetArrayFromImage(Image)
    Array = Array.transpose(2, 1, 0)

    # Measure image shape
    IMDimX = np.shape(Array)[0]
    IMDimY = np.shape(Array)[1]
    IMDimZ = np.shape(Array)[2]

    AddDimX = CoarseFactor - (IMDimX % CoarseFactor)
    AddDimY = CoarseFactor - (IMDimY % CoarseFactor)

    # adjust in x and y direction
    Shape_Diff = [AddDimX, AddDimY]
    IMG_XY_Adjusted = np.lib.pad(Array,
                                 ((0, Shape_Diff[0]), (0, Shape_Diff[1]), (0, 0)),
                                 'constant', constant_values=(0),)

    if CropZ == 'Crop':
        Image_Adjusted = IMG_XY_Adjusted

    if CropZ == 'Expand':
        AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
        Shape_Diff = [AddDimX, AddDimY, AddDimZ]
        Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                    ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                    'edge')

    if CropZ == 'Variable':
        Limit = CoarseFactor / 2.0
        if IMDimZ % CoarseFactor > Limit:
            AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
            Shape_Diff = [AddDimX, AddDimY, AddDimZ]
            Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                        ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                        'edge')
        if IMDimZ % CoarseFactor < Limit:
            Image_Adjusted = IMG_XY_Adjusted

    Image_Adjusted = sitk.GetImageFromArray(Image_Adjusted.transpose(2, 1, 0))
    Image_Adjusted.SetSpacing(Image.GetSpacing())
    Image_Adjusted.SetOrigin(Image.GetOrigin())
    Image_Adjusted.SetDirection (Image.GetDirection())

    return Image_Adjusted


#%% Initialization
# Set directories and list samples

CW, Data, Script, Results = SetDirectories('FRACTIB')

ImageDir = Data / '02_uCT'
hFEDir = Results / '03_hFE'
uCTDir = Results / '04_Registration'
ResDir = Results / '05_Localizations'

Samples = pd.read_csv(str(Data / 'SampleList.csv'))

#%% Load data
# Load data

Sample = Samples.loc[0, 'Internal ID']

Variables = ['J', 'F_Tilde']
uCT, hFE = [], []

for V in Variables:

    uCT_Data = sitk.ReadImage(str(uCTDir / Sample / (V + '.mhd')))
    uCT.append(uCT_Data)

    hFE_Data = sitk.ReadImage(str(hFEDir / Sample / (V + '.mhd')))
    hFE.append(hFE_Data)

PreFile = str(ImageDir / Sample / 'C0001901_DOWNSCALED.AIM')
PostFile = str(ImageDir / Sample / 'C0001939_DOWNSCALED.AIM')
Pre = Read.AIM(PreFile)[0]
Post = Read.AIM(PostFile)[0]


#%% Padding
# Padding
Pad = tuple(int(p) for p in np.array(uCT[0].GetSize()) - np.array(hFE[0].GetSize()))
PhFE = [sitk.ConstantPad(P, (0, 0, 0), Pad) for P in hFE]

Show.Overlay(uCT[0], PhFE[0], AsBinary = False)
Show.Overlay(uCT[0], uCT[0])

Show.Intensity(Pre, uCT[0], Axis='X')
Show.Intensity(Pre, Resample(hFE[0], Size=Post.GetSize()), Axis='X')

#%% Compare values
# Results

for i in range(2):
    hFE_Array = sitk.GetArrayFromImage(PhFE[i])
    uCT_Array = sitk.GetArrayFromImage(uCT[i])

    X = uCT_Array.flatten()#[hFE_Array > 0]
    Y = hFE_Array.flatten()#[hFE_Array > 0]
    X = uCT_Array[hFE_Array > 0]
    Y = hFE_Array[hFE_Array > 0]

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(X, Y, color=(1,0,0), linestyle='none', marker='o', fillstyle='none')
    Axis.set_xlabel('uCT values')
    Axis.set_ylabel('hFE values')
    plt.show()




# %%
