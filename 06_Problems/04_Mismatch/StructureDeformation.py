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

ExpData = pd.read_csv(str(Results / '02_Experiment' / Sample / 'MatchedSignals.csv'))
SimData5 = ReadDAT('C000' + str(uCT_ID) + '_5.dat')
SimData2 = ReadDAT('C000' + str(uCT_ID) + '_2.dat')
SimData3 = ReadDAT('C000' + str(uCT_ID) + '_3.dat')
SimData1 = ReadDAT('C000' + str(uCT_ID) + '_1.dat')

#%% Force displacement

X, Y = ['Z', 'FZ']

Figure, Axis = plt.subplots(1,1)
Axis.plot(ExpData[X], -ExpData[Y], color=(0,0,0), label='Experiment')
Axis.plot(SimData5[X], SimData5[Y], color=(1,0,0), label='All free')
Axis.plot(SimData2[X], SimData2[Y], color=(1,0,1), label='Translations locked')
Axis.plot(SimData3[X], SimData3[Y], color=(0,0,1), label='Rotations locked')
Axis.plot(SimData1[X], SimData1[Y], color=(0,1,1), label='All locked')
Axis.set_xlabel(X)
Axis.set_ylabel(Y)
plt.legend()
plt.show()

#%% Peak detection

Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=0.02)
Start = Peaks[6]
Stop = np.argmin(np.abs(ExpData['Z'] - SimData1['Z'].max()))

Figure, Axis = plt.subplots(1,1)
Axis.plot(ExpData['Z'] / ExpData['Z'].max(), color=(1,0,0))
Axis.plot(ExpData['FZ'] / ExpData['FZ'].min(), color=(0,0,1))
Axis.plot(Start, ExpData['FZ'][Start] / ExpData['FZ'].min(),
          linestyle='none', marker='o', fillstyle='none', color=(0,0,0))
Axis.plot(Stop, ExpData['FZ'][Stop] / ExpData['FZ'].min(),
          linestyle='none', marker='o', fillstyle='none', color=(0,0,0))
plt.show()

# %%


Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=True, sharey=True)
Axis[0].plot(ExpData['X'][Start:Stop], ExpData['Z'][Start:Stop], color=(0,0,0), label='Experiment')
Axis[0].plot(SimData1['X'], SimData1['Z'], color=(0,1,1), label='All locked')
Axis[0].plot(SimData2['X'], SimData2['Z'], color=(0,0,1), label='Translations locked')
Axis[0].plot(SimData3['X'], SimData3['Z'], color=(1,0,1), label='Rotations locked')
Axis[0].plot(SimData5['X'], SimData5['Z'], color=(1,0,0), label='All Free')
Axis[1].plot(ExpData['Y'][Start:Stop], ExpData['Z'][Start:Stop], color=(0,0,0), label='Experiment')
Axis[1].plot(SimData1['Y'], SimData1['Z'], color=(0,1,1), label='All locked')
Axis[1].plot(SimData2['Y'], SimData2['Z'], color=(0,0,1), label='Translations locked')
Axis[1].plot(SimData3['Y'], SimData3['Z'], color=(1,0,1), label='Rotations locked')
Axis[1].plot(SimData5['Y'], SimData5['Z'], color=(1,0,0), label='All Free')
Axis[0].set_xlabel('X')
Axis[0].set_ylabel('Z')
Axis[1].set_xlabel('Y')
plt.legend()
plt.show()

#%% Rotate coordinate system to align them
from scipy.optimize import minimize

InterpX = np.interp(ExpData['Z'][Start:Stop], SimData5['Z'], SimData5['X'])
InterpY = np.interp(ExpData['Z'][Start:Stop], SimData5['Z'], SimData5['Y'])
Ref = np.array([InterpX, InterpY, ExpData['Z'][Start:Stop]]).T

Sys = np.array(ExpData[['X','Y','Z']][Start:Stop])


def RotateSystem(Angle, Ref, Sys):
    M = RotationMatrix(Gamma=Angle[0])
    rSys = np.dot(Sys, M)
    Cost = np.linalg.norm(Ref - rSys, axis=1).sum()
    return Cost

Result = minimize(RotateSystem, [np.pi/2], args=(Ref, Sys), bounds=([0, 2*np.pi],))
Angle = Result.x[0]

RotatedSystem = np.dot(Sys, RotationMatrix(Gamma=Angle))

Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=True, sharey=True)
Axis[0].plot(RotatedSystem[:,0], -RotatedSystem[:,2], color=(0,0,0), label='Experiment')
Axis[0].plot(SimData1['X'], -SimData1['Z'], color=(0,1,1), label='All locked')
Axis[0].plot(SimData2['X'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
Axis[0].plot(SimData3['X'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
Axis[0].plot(SimData5['X'], -SimData5['Z'], color=(1,0,0), label='All Free')
Axis[1].plot(RotatedSystem[:,1], -RotatedSystem[:,2], color=(0,0,0), label='Experiment')
Axis[1].plot(SimData1['Y'], -SimData1['Z'], color=(0,1,1), label='All locked')
Axis[1].plot(SimData2['Y'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
Axis[1].plot(SimData3['Y'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
Axis[1].plot(SimData5['Y'], -SimData5['Z'], color=(1,0,0), label='All Free')
Axis[0].set_xlabel('X')
Axis[0].set_ylabel('Z')
Axis[1].set_xlabel('Y')
plt.legend()
plt.show()


# %%

# Rotate angles values

Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=True, sharey=True)
Axis[0].plot(ExpData['Phi'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,0), label='Experiment')
Axis[0].plot(SimData1['Phi'], -SimData1['Z'], color=(0,1,1), label='All locked')
Axis[0].plot(SimData2['Phi'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
Axis[0].plot(SimData3['Phi'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
Axis[0].plot(SimData5['Phi'], -SimData5['Z'], color=(1,0,0), label='All Free')
Axis[1].plot(ExpData['Theta'][Start:Stop], -ExpData['Theta'][Start:Stop], color=(0,0,0), label='Experiment')
Axis[1].plot(SimData1['Theta'], -SimData1['Z'], color=(0,1,1), label='All locked')
Axis[1].plot(SimData2['Theta'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
Axis[1].plot(SimData3['Theta'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
Axis[1].plot(SimData5['Theta'], -SimData5['Z'], color=(1,0,0), label='All Free')
Axis[0].set_xlabel('Psi')
Axis[0].set_ylabel('Z')
Axis[1].set_xlabel('Theta')
plt.legend()
plt.show()

# %%
