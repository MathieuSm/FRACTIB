#%%
#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
from Utils import *
from ReadDAT import Main as ReadDAT

CW, Data, Scripts, Results = SetDirectories('FRACTIB')


SampleList = pd.read_csv(str(Data / 'SampleList.csv'))

Index = 0
Sample = SampleList.loc[Index, 'Internal ID']
uCT_ID = SampleList.loc[Index, 'MicroCT pretest file number']

#%% Load data

ExpData = pd.read_csv(str(Data / '03_Experiment' / '3_Matched' / (Sample + '.csv')))
SimData5 = ReadDAT('C000' + str(uCT_ID) + '_5.dat')
SimData2 = ReadDAT('C000' + str(uCT_ID) + '_2.dat')
SimData3 = ReadDAT('C000' + str(uCT_ID) + '_3.dat')
SimData1 = ReadDAT('C000' + str(uCT_ID) + '_1.dat')

#%% Force displacement
X = 'Displacement (mm)'
Y = 'Force (N)'

Figure, Axis = plt.subplots(1,1)
Axis.plot(ExpData[X], ExpData[Y], color=(0,0,0), label='Experiment')
Axis.plot(-SimData5['U3'], -SimData5['F3'], color=(1,0,0), label='All free')
Axis.plot(-SimData2['U3'], -SimData2['F3'], color=(1,0,1), label='Translations locked')
Axis.plot(-SimData3['U3'], -SimData3['F3'], color=(0,0,1), label='Rotations locked')
Axis.plot(-SimData1['U3'], -SimData1['F3'], color=(0,1,1), label='All locked')
Axis.set_xlabel(X)
Axis.set_ylabel(Y)
plt.legend()
plt.show()

#%% Force max - displacement at force max

Exp_Fmax = ExpData[Y].min()
Exp_Dmax = ExpData.loc[ExpData[Y].idxmin(),X]

Sim5_Fmax = SimData5['F3'].max()
Sim5_Dmax = -SimData5.loc[SimData5['F3'].idxmax(), 'U3']

Sim2_Fmax = SimData2['F3'].max()
Sim2_Dmax = -SimData2.loc[SimData2['F3'].idxmax(), 'U3']

Sim3_Fmax = SimData3['F3'].max()
Sim3_Dmax = -SimData3.loc[SimData3['F3'].idxmax(), 'U3']

Sim1_Fmax = SimData1['F3'].max()
Sim1_Dmax = -SimData1.loc[SimData1['F3'].idxmax(), 'U3']

Figure, Axis = plt.subplots(1,1)
Axis.plot(Exp_Dmax, Exp_Fmax, marker='o', fillstyle='none', linestyle='none', color=(0,0,0), label='Experiment')
Axis.plot(Sim5_Dmax, -Sim5_Fmax, marker='o', fillstyle='none', linestyle='none', color=(1,0,0), label='All free')
Axis.plot(Sim2_Dmax, -Sim2_Fmax, marker='o', fillstyle='none', linestyle='none', color=(1,0,1), label='Translations locked')
Axis.plot(Sim3_Dmax, -Sim3_Fmax, marker='o', fillstyle='none', linestyle='none', color=(0,0,1), label='Rotations locked')
Axis.plot(Sim1_Dmax, -Sim1_Fmax, marker='o', fillstyle='none', linestyle='none', color=(0,1,1), label='All locked')
Axis.set_xlabel(X)
Axis.set_ylabel(Y)
plt.legend()
plt.show()

# %%
