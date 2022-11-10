#%%
# 00 Initialization
import os
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import SimpleITK as sitk

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%%
# 01 Set variables
DataDirectory = Path.cwd() / '../../04_Results/03_Registration/'
DataFolders = [Dir for Dir in os.listdir(DataDirectory) if Path.is_dir(DataDirectory / Dir)]
Results = pd.DataFrame()

#%%
## Load file
for Sample in DataFolders:
    SamplePath = DataDirectory / Sample
    FixedFile = str(SamplePath / 'Fixed.mhd')
    RegisteredFile = str(SamplePath / 'Registered.mhd')

    FixedImage = sitk.ReadImage(FixedFile)
    RegisteredImage = sitk.ReadImage(RegisteredFile)

#%%
## Perform Otsu segmentation
    OtsuFilter = sitk.OtsuThresholdImageFilter()
    OtsuFilter.SetInsideValue(1)
    OtsuFilter.SetOutsideValue(0)

    FixedSegmented = OtsuFilter.Execute(FixedImage)
    RegisteredSegmented = OtsuFilter.Execute(RegisteredImage)

#%%
# Compute mean Dice coefficient and store
    FixedArray = sitk.GetArrayFromImage(FixedSegmented)
    RegisteredArray = sitk.GetArrayFromImage(RegisteredSegmented)

    Dice = 2 * np.sum(FixedArray * RegisteredArray) / np.sum(FixedArray + RegisteredArray)

    print('Sample ' + Sample + ' dice coefficient: ' + str(round(Dice,3)))
    Results.loc[Sample, 'Dice'] = Dice

Results.to_csv(str(DataDirectory / 'RegistrationResults.csv'))    
        
# %%
