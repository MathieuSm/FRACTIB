#%%
#!/usr/bin/env python3

# 00 Initialization
import pandas as pd
from Utils import *
from pathlib import Path

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Paths
# Set path variables
HRpQCT_Data = Path.cwd() / '../../02_Data/01_HRpQCT'
uCT_Data = Path.cwd() / '../../02_Data/02_uCT'
ResultsPath = Path.cwd() / '../../04_Results/05_FractureLinePrediction'
SampleList = pd.read_csv(str(Path.cwd() / '../../02_Data/SampleList.csv'))

#%% File Loading
# Load files
iSample = 0
Sample = SampleList.loc[iSample, 'Internal ID']

# 03 Load Masks
uCT_Number = SampleList.loc[iSample, 'MicroCT pretest file number']
uCT_File = 'C000' + str(uCT_Number) + '_reso_0.098_DOWNSCALED_FULLMASK.mhd'
uCT_Mask = sitk.ReadImage(str(uCT_Data / Sample / uCT_File))

HRpQCT_Number = SampleList.loc[iSample, 'HRpQCT File 2 number']
HRpQCT_File = 'C000' + str(HRpQCT_Number)
Cort_File = str(HRpQCT_Data / Sample / (HRpQCT_File + '_CORT_MASK_UNCOMP.AIM'))
Trab_File = str(HRpQCT_Data / Sample / (HRpQCT_File + '_TRAB_MASK_UNCOMP.AIM'))
HRpQCT_Cort_Mask, AdditionalData = Read.AIM(Cort_File)
HRpQCT_Trab_Mask, AdditionalData = Read.AIM(Trab_File)
HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask

