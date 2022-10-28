# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


# 01 Set variables
WorkingDirectory = os.getcwd()
DataDirectory = os.path.join(WorkingDirectory,'02_Data/03_uCT/')

## Load Files
SampleList = pd.read_csv(DataDirectory+'BMD_Values.csv')

MHD_Files = pd.DataFrame()
for Index in SampleList.index:
    Sample = SampleList.loc[Index,'Sample']
    SampleDirectory = os.path.join(DataDirectory,Sample+'/')
    MHD_File = [File for File in os.listdir(SampleDirectory) if File.endswith('DOWNSCALED_BMD.mhd')]
    MHD_Files = MHD_Files.append({'Scan':MHD_File[0][:-8]+'.mhd'},ignore_index=True)

SampleList['Scan'] = MHD_Files['Scan']
SampleList.to_csv(DataDirectory+'BMD1_Values.csv',index=False)


## Write batch
BatchFile = open(os.path.join(WorkingDirectory, '03_Scripts/Registration', '00_ScanBatch.bash'),'w')
PipelineFile  = open(os.path.join(WorkingDirectory, '03_Scripts/Registration', '00_Masking_MedtoolCommand.txt'), 'r')
PipelineText  = PipelineFile.read()

for Index in SampleList.index:

    Scan = SampleList.loc[Index, 'Scan']
    Threshold = SampleList.loc[Index, 'BMD Threshold (HA g/cm3)']

    Text = PipelineText.split('\n')[1]
    Text = Text.replace('Scan', Scan[:-4])
    Text = Text.replace('Threshold', str(int(Threshold.round(0))))

    BatchFile.write(Text+'\n\n')

BatchFile.close()


## Write remote to local copy file
CopyFile = open(os.path.join(WorkingDirectory, '03_Scripts/Registration', '00_CopyResults.bash'),'w')
CopyFile.write('#!/bin/bash\n\n')
FullText = 'sshpass -p "Password" scp -r RemotePath LocalPath'

Password = 'Password'
RemotePath = 'ms20s284@130.92.125.21:PostMsc2/Registration_Masking/'

for Index in SampleList.index:

    # Set variables
    Sample = SampleList.loc[Index,'Sample']
    LocalPath = os.path.join(DataDirectory,Sample+'/')
    File = SampleList.loc[Index, 'Scan'][:-4] + '_FULLMASK.mhd'

    # Write text
    Text = FullText
    Text = Text.replace('Password', Password)
    Text = Text.replace('RemotePath', RemotePath + File)
    Text = Text.replace('LocalPath', LocalPath)

    CopyFile.write(Text + '\n')

    Text = Text.replace(File, File[:-4] + '.raw')

    CopyFile.write(Text + '\n\n')

CopyFile.close()
