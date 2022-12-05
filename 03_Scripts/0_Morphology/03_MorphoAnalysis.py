# 00 Initialization
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


# 01 Set variables
WorkingDirectory = os.getcwd()
MorphoPath = os.path.join(WorkingDirectory,'04_Results/02_MorphoAnalysis/')
DataPath = os.path.join(WorkingDirectory,'02_Data/')


## Load files
SampleList = pd.read_csv(DataPath + 'SampleList.csv')
MorphoData = pd.read_csv(MorphoPath + '00_Data.csv')

## Add sample numbers
MorphoData['Sample Number'] = MorphoData['Scan']
for Index in MorphoData.index:
    uCT_SN = MorphoData['Scan'].loc[Index]
    Filter = SampleList['MicroCT pretest file number'] == int(uCT_SN[-4:])
    MorphoData.loc[Index,'Sample Number'] = int(SampleList[Filter]['Internal ID'].values[0][:3])
MorphoData = MorphoData.sort_values(by='Sample Number')


# 02 Morphometric analysis
C = MorphoData.columns
DataLabel = 'Mean Tb Sp'

Figure, Axes = plt.subplots(1, 1, figsize=(3.5, 4.5),dpi=100)
Axes.boxplot(MorphoData[DataLabel],vert=True,
             showmeans=False,
             boxprops=dict(linestyle='-',color=(0,0,0)),
             medianprops=dict(linestyle='-',color=(1,0,0)),
             whiskerprops=dict(linestyle='--',color=(0,0,0)),
             meanprops=dict(marker='x',markeredgecolor=(0,0,1)))
Axes.plot([1],MorphoData.loc[7,DataLabel],linestyle='none',marker='o',color=(0,0,1))
Axes.set_ylabel(DataLabel)
Axes.set_xticks([])
Axes.set_title('')
plt.title('')
plt.suptitle('')
plt.subplots_adjust(0.2)
plt.show()
plt.close(Figure)


BoxplotOutliers = [7,19,15,6,10]
RegressionOutliers = [440,454,446,451]


Filter1 = Data['BMC (HA mg)'] > 6300
Filter2 = Data['BMC (HA mg)'] < 6600
Data[Filter1&Filter2]

Data.loc[RegressionOutliers]
MorphoData.loc[BoxplotOutliers]