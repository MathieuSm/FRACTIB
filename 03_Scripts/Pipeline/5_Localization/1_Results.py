#%% #!/usr/bin/env python3
import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path

from Utils import *

#%% Initialization
# Set directories and list samples

CW, Data, Script, Results = SetDirectories('FRACTIB')

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

Mask = sitk.ReadImage(str(ResDir / Sample / 'CommonMask.mhd'))

#%% Resampling
# Resampling

RefSize = np.array(hFE[0].GetSize())

R_Mask = Resample(Mask, Size=RefSize)

a_ = sitk.GetArrayFromImage(uCT[1])
Figure, Axis = plt.subplots(1,1)
Axis.imshow(a_[a_.shape[0]//2,:,:])
plt.show()

R_uCT = Resample(uCT[0], Size=RefSize)


#%% Compare values
# Results

hFE_Array = sitk.GetArrayFromImage(hFE[0])
uCT_Array = sitk.GetArrayFromImage(R_uCT)
Mask_Array = sitk.GetArrayFromImage(R_Mask)

X = uCT_Array * Mask_Array
Y = hFE_Array * Mask_Array

Filter_uCT = X > 0
Filter_hFE = Y > 0

Xf = X[Filter_uCT * Filter_hFE]
Yf = Y[Filter_uCT * Filter_hFE]

Figure, Axis = plt.subplots(1,1)
Axis.plot(Xf, Yf, color=(1,0,0), linestyle='none', marker='o', fillstyle='none')
plt.show()




# %%
