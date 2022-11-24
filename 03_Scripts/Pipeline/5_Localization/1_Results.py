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

    uCT_Data = sitk.ReadImage(str(uCTDir, Sample, (V + '.mhd')))
    uCT.append(uCT_Data)

    

