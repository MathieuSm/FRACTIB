#%%
#!/usr/bin/env python3

# 00 Imports

## Import libraries
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime
from scipy.signal import find_peaks
from scipy import signal
from scipy.fft import fft
from scipy import interpolate
from scipy.stats import linregress
from pathlib import Path


desired_width = 320
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width)

CurrentWorkingDirectory = Path.cwd() / '../..'
SampleList = pd.read_csv(str(CurrentWorkingDirectory / '02_Data/SampleList.csv'))


#%%
# 01 Set paths and load files
DataPath  = '02_Data/04_Experiment'
MachineData = ['1_MTS','2_ARAMIS']
FilesPath = 'FracTib_ARAMIS_270820'
ResultsPath = str(CurrentWorkingDirectory / '04_Results/FractureZoneAssessment/01_Experiment/SamplesData')

# Set machine data path
MTSPath = os.path.join(CurrentWorkingDirectory,DataPath,MachineData[0])
ARAPath = os.path.join(CurrentWorkingDirectory,DataPath,MachineData[1])

#%%
# Open MTS data and get keys
with open(MTSPath+'/AO_data_MTS.json') as DictFile:
  MTSData = json.load(DictFile)
SamplesNames = [*MTSData]

MTSPeaks = pd.read_csv(os.path.join(ResultsPath,'MTSPeaks.csv'))

#%%
# List ARAMIS .csv files and store the samples names
CSVFiles = [File for File in os.listdir(ARAPath) if File.endswith(".csv")]
CSVSamples = []
for Name in CSVFiles:
  CSVSamples.append(Name[8:18])

#%%
# Select sample to analyze and the corresponding csv file
S = 24-24
SampleName = SamplesNames[S]
print('Sample name:')
print(SampleName)

SamplesData = pd.read_csv(os.path.join(ResultsPath,'ExperimentData.csv'))


CSVFile = CSVFiles[CSVSamples.index(SampleName)]

#%%
# Load data
MTSSampleData = pd.DataFrame.from_dict(MTSData[SampleName])
Time = 'Time UTC'
ARASampleData = pd.read_csv(os.path.join(ARAPath,CSVFile),sep=';',header=1,parse_dates=[Time])


#%%
# 02 Import MTS data
MTSTimes = MTSSampleData['time'].values
MTSDisplacements = MTSSampleData['disp'].values
MTSForces = MTSSampleData['force'].values

MTSTimeSteps = np.zeros(len(MTSTimes)-1)
for i in range(len(MTSTimeSteps)-1):
  MTSTimeSteps[i]=MTSTimes[i+1]-MTSTimes[i]
StartIndex = np.where(abs(MTSTimeSteps)==max(abs(MTSTimeSteps)))[0][0]

for Index in range(StartIndex+1,len(MTSTimes)):
  MTSTimes[Index] = MTSTimes[Index] + MTSTimes[StartIndex]

## Plot adjusted time
# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
# Axes.plot(MTSTimes,color=(0,0,0))
# plt.show()
# plt.close(Figure)

SamplingFrequency = 1/(MTSTimes[1]-MTSTimes[0])

#%%
# 03 Filter signals
P1, P2 = signal.butter(2, 0.02)
MTSFilteredDisplacements = signal.filtfilt(P1, P2, MTSDisplacements)
MTSFilteredForces = signal.filtfilt(P1, P2, MTSForces)


## Plot filtered signals
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MTSFilteredDisplacements,color=(0,0,0))
Axes.set_xlabel('Index (-)')
Axes.set_ylabel('Displacement (mm)')
plt.grid('on')
plt.xticks(np.arange(0, np.ceil(len(MTSFilteredDisplacements)/1e4)*1e4+5000, 5000),rotation='vertical')
plt.yticks(np.arange(0, np.floor(min(MTSFilteredDisplacements)), -0.2))
plt.show()
# MTSFilteredDisplacements=MTSFilteredDisplacements-MTSFilteredDisplacements[0]
plt.close(Figure)

(MTSSampleData['disp'][25000]-MTSSampleData['disp'][20000])/(MTSSampleData['time'][25000]-MTSSampleData['time'][20000])*60

#%%
# 04 Find peaks and troughs
Cycle1Range = [0,11000]
PeaksThreshold = 6e-02
TroughsThreshold = 5e-02

MTSPeaks, _ = find_peaks(-MTSFilteredDisplacements[Cycle1Range[0]:Cycle1Range[1]],PeaksThreshold)
MTSFirstCyclePeaks = np.array([MTSPeaks[0],MTSPeaks[1],MTSPeaks[2],MTSPeaks[3],MTSPeaks[4]])
MTSFirstCyclePeaks = MTSFirstCyclePeaks + Cycle1Range[0]

MTSTroughs, _ = find_peaks(MTSFilteredDisplacements[Cycle1Range[0]:Cycle1Range[1]],-TroughsThreshold)
MTSFirstCycleTroughs = np.array([MTSTroughs[-5],MTSTroughs[-4],MTSTroughs[-3],MTSTroughs[-2],MTSTroughs[-1]])
MTSFirstCycleTroughs = np.array([MTSTroughs[-6],MTSTroughs[-5],MTSTroughs[-4],MTSTroughs[-3],MTSTroughs[-2]])
# MTSFirstCycleTroughs = np.array([MTSTroughs[-35],MTSTroughs[-34],MTSTroughs[-33],MTSTroughs[-15],MTSTroughs[-14]])
MTSFirstCycleTroughs = MTSFirstCycleTroughs + Cycle1Range[0]

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MTSFilteredDisplacements[Cycle1Range[0]:Cycle1Range[1]],
          color=(0,0,0),label='Filtered signal')
Axes.plot(MTSFirstCyclePeaks,MTSFilteredDisplacements[MTSFirstCyclePeaks],
          color=(1,0,0),marker='o',linestyle='none',fillstyle='none',label='Cycle peaks')
Axes.plot(Cycle1Range,[-PeaksThreshold,-PeaksThreshold],
          color=(1,0,0),linestyle='--',label='Peaks Threshold')
Axes.plot(MTSFirstCycleTroughs,MTSFilteredDisplacements[MTSFirstCycleTroughs],
          color=(0,0,1),marker='o',linestyle='none',fillstyle='none',label='Cycle troughs')
Axes.plot(Cycle1Range,[-TroughsThreshold,-TroughsThreshold],
          color=(0,0,1),linestyle='--',label='Troughs Threshold')
Axes.set_xlabel('Index (-)')
Axes.set_ylabel('Displacement (mm)')
plt.show()
plt.close(Figure)



Cycle2Range = [38000,len(MTSFilteredDisplacements)]
PeaksThreshold = 2.7
TroughsThreshold = 2.4

MTSPeaks, _ = find_peaks(-MTSFilteredDisplacements[Cycle2Range[0]:Cycle2Range[1]],PeaksThreshold)
MTSSecondCyclePeaks = np.array([MTSPeaks[0],MTSPeaks[1],MTSPeaks[2]])
MTSSecondCyclePeaks = MTSSecondCyclePeaks + Cycle2Range[0]

MTSTroughs, _ = find_peaks(MTSFilteredDisplacements[Cycle2Range[0]:Cycle2Range[1]],-TroughsThreshold)
MTSSecondCycleTroughs = np.array([MTSTroughs[-2],MTSTroughs[-1],2*MTSTroughs[-1]-MTSTroughs[-2]])
# MTSSecondCycleTroughs = np.array([MTSTroughs[-2],MTSTroughs[-1],Cycle2Range[1]-Cycle2Range[0]-2])
MTSSecondCycleTroughs = MTSSecondCycleTroughs + Cycle2Range[0]

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(np.arange(Cycle2Range[0],Cycle2Range[1]),MTSFilteredDisplacements[Cycle2Range[0]:Cycle2Range[1]],
          color=(0,0,0),label='Filtered signal')
Axes.plot(MTSSecondCyclePeaks,MTSFilteredDisplacements[MTSSecondCyclePeaks],
          color=(1,0,0),marker='o',linestyle='none',fillstyle='none',label='Cycle peaks')
Axes.plot(Cycle2Range,[-PeaksThreshold,-PeaksThreshold],
          color=(1,0,0),linestyle='--',label='Peaks Threshold')
Axes.plot(MTSSecondCycleTroughs,MTSFilteredDisplacements[MTSSecondCycleTroughs],
          color=(0,0,1),marker='o',linestyle='none',fillstyle='none',label='Cycle troughs')
Axes.plot(Cycle2Range,[-TroughsThreshold,-TroughsThreshold],
          color=(0,0,1),linestyle='--',label='Troughs Threshold')
Axes.set_xlabel('Index (-)')
Axes.set_ylabel('Displacement (mm)')
plt.show()
plt.close(Figure)

#%%
# 05 Import ARAMIS data
TopPlateAnatomical_LZ = 'TopPlate_Csys→AnatomicalCsys.LZ [mm]'
ARADisplacements = -ARASampleData[TopPlateAnatomical_LZ]+ARASampleData[TopPlateAnatomical_LZ][0]

# # Plot displacement
# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
# Axes.plot(ARADisplacements,color=(0,0,0))
# Axes.set_xlabel('Index (-)')
# Axes.set_ylabel('Displacement (mm)')
# plt.grid('on')
# plt.xticks(np.arange(0, np.ceil(len(ARADisplacements)/1e4)*1e4+5000, 5000),rotation='vertical')
# plt.yticks(np.arange(0, np.floor(min(ARADisplacements)), -0.2))
# plt.show()
# plt.close(Figure)
#
# ## If jump at the end, cut the signal
# ARADisplacements = ARADisplacements[:23500]
# ARATimes = ARATimes[:23500]

ARATimes = (ARASampleData[Time]-ARASampleData[Time][0]).dt.total_seconds()

#%%
# 06 Interpolate displacement data
ReSampledARATimes = np.linspace(0,max(ARATimes.values),int((max(ARATimes.values)-0)*SamplingFrequency+1))
InterpARADisp = interpolate.interp1d(ARATimes.values/max(ARATimes.values),ARADisplacements)
ReSampledARADisp = InterpARADisp(ReSampledARATimes/max(ReSampledARATimes))
ReSampledARADisp = pd.DataFrame(ReSampledARADisp).fillna(method='ffill').values.ravel()

#%%
# 07 Filter signal
P1, P2 = signal.butter(2, 0.005)
ARAFilteredDisplacements = signal.filtfilt(P1, P2, ReSampledARADisp)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ARAFilteredDisplacements,color=(0,0,0))
Axes.set_xlabel('Index (-)')
Axes.set_ylabel('Displacement (mm)')
plt.grid('on')
plt.xticks(np.arange(0, np.ceil(len(ARAFilteredDisplacements)/1e4)*1e4+5000, 5000),rotation='vertical')
plt.yticks(np.arange(0, np.floor(min(ARAFilteredDisplacements)), -0.2))
plt.show()
plt.close(Figure)


#%%
# 08 Find peaks
Cycle1Range = [0,13000]
PeaksThreshold = 3.8e-02

ARAPeaks, _ = find_peaks(-ARAFilteredDisplacements[Cycle1Range[0]:Cycle1Range[1]],PeaksThreshold)
ARAFirstCyclePeaks = np.array([ARAPeaks[0],ARAPeaks[1],ARAPeaks[2],ARAPeaks[3],ARAPeaks[4]])
# ARAFirstCyclePeaks = np.array([ARAPeaks[0],ARAPeaks[2],ARAPeaks[3],ARAPeaks[4],ARAPeaks[5]])
ARAFirstCyclePeaks = ARAFirstCyclePeaks + Cycle1Range[0]

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ARAFilteredDisplacements[Cycle1Range[0]:Cycle1Range[1]],
          color=(0,0,0),label='Filtered signal')
Axes.plot(ARAFirstCyclePeaks,ARAFilteredDisplacements[ARAFirstCyclePeaks],
          color=(1,0,0),marker='o',linestyle='none',fillstyle='none',label='Cycle peaks')
Axes.plot(Cycle1Range,[-PeaksThreshold,-PeaksThreshold],
          color=(1,0,0),linestyle='--',label='Peaks Threshold')
Axes.set_xlabel('Index (-)')
Axes.set_ylabel('Displacement (mm)')
plt.show()
plt.close(Figure)



Cycle2Range = [59000,72000]
PeaksThreshold = 2.6

ARAPeaks, _ = find_peaks(-ARAFilteredDisplacements[Cycle2Range[0]:Cycle2Range[1]],PeaksThreshold)
ARASecondCyclePeaks = np.array([ARAPeaks[0],ARAPeaks[1],ARAPeaks[2]])
ARASecondCyclePeaks = ARASecondCyclePeaks + Cycle2Range[0]

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(np.arange(Cycle2Range[0],Cycle2Range[1]),ARAFilteredDisplacements[Cycle2Range[0]:Cycle2Range[1]],
          color=(0,0,0),label='Filtered signal')
Axes.plot(ARASecondCyclePeaks,ARAFilteredDisplacements[ARASecondCyclePeaks],
          color=(1,0,0),marker='o',linestyle='none',fillstyle='none',label='Cycle peaks')
Axes.plot(Cycle2Range,[-PeaksThreshold,-PeaksThreshold],
          color=(1,0,0),linestyle='--',label='Peaks Threshold')
Axes.set_xlabel('Index (-)')
Axes.set_ylabel('Displacement (mm)')
plt.show()
plt.close(Figure)


#%%
# 09 Align signals
FirstProcessARAPeaks = ReSampledARATimes[ARAFirstCyclePeaks]
FirstProcessMTSPeaks = MTSTimes[MTSFirstCyclePeaks]

SecondProcessARAPeaks = ReSampledARATimes[ARASecondCyclePeaks]
SecondProcessMTSPeaks = MTSTimes[MTSSecondCyclePeaks]

FirstProcessShift = np.mean(FirstProcessARAPeaks-FirstProcessMTSPeaks)
ShiftedARATimes = ReSampledARATimes-FirstProcessShift
FirstProcessShiftIndex = np.where(abs(ShiftedARATimes) == min(abs(ShiftedARATimes)))[0][0]
ShiftedARATimes = ShiftedARATimes[FirstProcessShiftIndex:]
ShiftedARADisplacement = ARAFilteredDisplacements[FirstProcessShiftIndex:]


## Actualize the peaks positions
ARAFirstCycleShiftedPeaks = ARAFirstCyclePeaks - FirstProcessShiftIndex
ARASecondCycleShiftedPeaks = ARASecondCyclePeaks - FirstProcessShiftIndex


## Find artificial junction in MTS data
MTSSteps = np.zeros(len(MTSDisplacements)-1)
for i in range(len(MTSSteps)-1):
  MTSSteps[i]=MTSDisplacements[i+1]-MTSDisplacements[i]

# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
# Axes.plot(MTSSteps,color=(0,0,0))
# Axes.set_xlabel('Index (-)')
# Axes.set_ylabel('Relative Displacement (mm)')
# plt.show()
# plt.close(Figure)

Junction = np.where(abs(MTSSteps)==max(abs(MTSSteps)))[0][0]
Junction = np.where(MTSSteps==max(MTSSteps))[0][0]


## Compute the second shift needed and its index
SecondProcessShift = np.mean(SecondProcessARAPeaks-SecondProcessMTSPeaks) - FirstProcessShift
SecondProcessShiftIndex = int(np.round(SecondProcessShift*SamplingFrequency))

ShiftedMTSTimes = np.zeros(len(MTSTimes))
ShiftedMTSTimes = ShiftedMTSTimes + MTSTimes
ShiftedMTSTimes[Junction:] = ShiftedMTSTimes[Junction:] + SecondProcessShift
CreatedTimes = np.linspace(ShiftedMTSTimes[Junction-1],ShiftedMTSTimes[Junction],SecondProcessShiftIndex+1)
ShiftedMTSTimes = np.insert(ShiftedMTSTimes,Junction,CreatedTimes[1:-1])

ShiftedMTSDisplacements = np.zeros(len(MTSFilteredDisplacements))
ShiftedMTSDisplacements = ShiftedMTSDisplacements + MTSFilteredDisplacements
CreatedDisplacements = np.linspace(ShiftedMTSDisplacements[Junction],ShiftedMTSDisplacements[Junction+1],SecondProcessShiftIndex+1)
ShiftedMTSDisplacements = np.insert(ShiftedMTSDisplacements,Junction+1,CreatedDisplacements[1:-1])

ShiftedMTSForces = np.zeros(len(MTSFilteredForces))
ShiftedMTSForces = ShiftedMTSForces + MTSFilteredForces
CreatedForces = np.linspace(ShiftedMTSForces[Junction],ShiftedMTSForces[Junction+1],SecondProcessShiftIndex+1)
ShiftedMTSForces = np.insert(ShiftedMTSForces,Junction+1,CreatedForces[1:-1])

MTSSecondCycleShiftedPeaks = MTSSecondCyclePeaks + SecondProcessShiftIndex
MTSSecondCycleShiftedTroughs = MTSSecondCycleTroughs + SecondProcessShiftIndex

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ShiftedMTSTimes,ShiftedMTSForces,color=(0,0,0),label='Signal')
Axes.plot(ShiftedMTSTimes[MTSFirstCyclePeaks],ShiftedMTSForces[MTSFirstCyclePeaks],color=(1,0,0),marker='o',fillstyle='none',linestyle='none',label='$1^{st}$ cycle peaks')
Axes.plot(ShiftedMTSTimes[MTSFirstCycleTroughs],ShiftedMTSForces[MTSFirstCycleTroughs],color=(1,0,0),marker='x',fillstyle='none',linestyle='none',label='$1^{st}$ cycle troughs')
Axes.plot(ShiftedMTSTimes[MTSSecondCycleShiftedPeaks],ShiftedMTSForces[MTSSecondCycleShiftedPeaks],color=(0,0,1),marker='o',fillstyle='none',linestyle='none',label='$2^{nd}$ cycle peaks')
Axes.plot(ShiftedMTSTimes[MTSSecondCycleShiftedTroughs],ShiftedMTSForces[MTSSecondCycleShiftedTroughs],color=(0,0,1),marker='x',fillstyle='none',linestyle='none',label='$2^{nd}$ cycle troughs')
Axes.set_xlabel('Time (s)')
Axes.set_ylabel('Force (N)')
plt.legend()
# plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ShiftedMTSTimes,ShiftedMTSDisplacements,color=(0.6,0.6,0.6),label='MTS')
Axes.plot(ShiftedMTSTimes,ShiftedARADisplacement[:len(ShiftedMTSTimes)],color=(0,0,0),label='ARAMIS')
Axes.plot(ShiftedMTSTimes[MTSFirstCyclePeaks],ShiftedARADisplacement[MTSFirstCyclePeaks],color=(1,0,0),marker='o',fillstyle='none',linestyle='none',label='$1^{st}$ cycle peaks')
Axes.plot(ShiftedMTSTimes[MTSFirstCycleTroughs],ShiftedARADisplacement[MTSFirstCycleTroughs],color=(1,0,0),marker='x',fillstyle='none',linestyle='none',label='$1^{st}$ cycle troughs')
Axes.plot(ShiftedMTSTimes[MTSSecondCycleShiftedPeaks],ShiftedARADisplacement[MTSSecondCycleShiftedPeaks],color=(0,0,1),marker='o',fillstyle='none',linestyle='none',label='$2^{nd}$ cycle peaks')
Axes.plot(ShiftedMTSTimes[MTSSecondCycleShiftedTroughs],ShiftedARADisplacement[MTSSecondCycleShiftedTroughs],color=(0,0,1),marker='x',fillstyle='none',linestyle='none',label='$2^{nd}$ cycle troughs')
Axes.set_xlabel('Time (s)')
Axes.set_ylabel('Displacement (mm)')
plt.legend()
# plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
plt.show()
plt.close(Figure)

#%%
# 10 Evaluate force-displacement curve
PreDamageDeltaDisplacements = ShiftedARADisplacement[MTSFirstCycleTroughs] - ShiftedARADisplacement[MTSFirstCyclePeaks]
PreDamageDeltaForces = ShiftedMTSForces[MTSFirstCycleTroughs] - ShiftedMTSForces[MTSFirstCyclePeaks]
PreDamageStiffness = PreDamageDeltaForces / PreDamageDeltaDisplacements
print('Pre-damage Stiffness:')
print(np.round(PreDamageStiffness))

PostDamageDeltaDisplacements = ShiftedARADisplacement[MTSSecondCycleShiftedTroughs] - ShiftedARADisplacement[MTSSecondCycleShiftedPeaks]
PostDamageDeltaForces = ShiftedMTSForces[MTSSecondCycleShiftedTroughs] - ShiftedMTSForces[MTSSecondCycleShiftedPeaks]
PostDamageStiffness = PostDamageDeltaForces / PostDamageDeltaDisplacements
print('Post-damage stiffness:')
print(np.round(PostDamageStiffness))

# Plot curve
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ShiftedARADisplacement[:len(ShiftedMTSForces)],ShiftedMTSForces,color=(0,0,0),label='Signal')
Axes.plot(ShiftedARADisplacement[MTSFirstCyclePeaks],ShiftedMTSForces[MTSFirstCyclePeaks],color=(1,0,0),marker='o',fillstyle='none',linestyle='none',label='$1^{st}$ cycle peaks')
Axes.plot(ShiftedARADisplacement[MTSFirstCycleTroughs],ShiftedMTSForces[MTSFirstCycleTroughs],color=(1,0,0),marker='x',fillstyle='none',linestyle='none',label='$1^{st}$ cycle troughs')
Axes.plot(ShiftedARADisplacement[MTSSecondCycleShiftedPeaks],ShiftedMTSForces[MTSSecondCycleShiftedPeaks],color=(0,0,1),marker='o',fillstyle='none',linestyle='none',label='$2^{nd}$ cycle peaks')
Axes.plot(ShiftedARADisplacement[MTSSecondCycleShiftedTroughs],ShiftedMTSForces[MTSSecondCycleShiftedTroughs],color=(0,0,1),marker='x',fillstyle='none',linestyle='none',label='$2^{nd}$ cycle troughs')
Axes.set_xlabel('Displacement (mm)')
Axes.set_ylabel('Force (N)')
plt.legend()
# plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
plt.show()
plt.close(Figure)

#%%
# Compute stiffness with moving regression
UltimateForce = min(ShiftedMTSForces)
UltimateForceIndex = np.where(ShiftedMTSForces==min(ShiftedMTSForces))[0][0]
RegressionRange = int((UltimateForceIndex-MTSFirstCycleTroughs[-1])/3)

InitialStiffness = np.zeros(UltimateForceIndex-MTSFirstCycleTroughs[-1]-RegressionRange)
for Index in range(UltimateForceIndex-MTSFirstCycleTroughs[-1]-RegressionRange):
    Index = Index + MTSFirstCycleTroughs[-1]

    RegressionDisplacementData = ShiftedARADisplacement[Index:Index+RegressionRange]
    RegressionForceData = ShiftedMTSForces[Index:Index+RegressionRange]

    InitialStiffness[Index-MTSFirstCycleTroughs[-1]], Interception, RValue, PValue, StdError = linregress(RegressionDisplacementData,RegressionForceData)

# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
# Axes.plot(InitialStiffness,color=(0,0,0),label='Regression results')
# Axes.plot([0,len(InitialStiffness)],[PreDamageStiffness[0],PreDamageStiffness[0]],color=(1,0,0),linestyle='--',label='Stiffness cycle 1')
# Axes.plot([0,len(InitialStiffness)],[PreDamageStiffness[1],PreDamageStiffness[1]],color=(1,0,1),linestyle='--',label='Stiffness cycle 2')
# Axes.plot([0,len(InitialStiffness)],[PreDamageStiffness[2],PreDamageStiffness[2]],color=(0,0,1),linestyle='--',label='Stiffness cycle 3')
# Axes.plot([0,len(InitialStiffness)],[PreDamageStiffness[3],PreDamageStiffness[3]],color=(0,1,1),linestyle='--',label='Stiffness cycle 4')
# Axes.plot([0,len(InitialStiffness)],[PreDamageStiffness[4],PreDamageStiffness[4]],color=(0,1,0),linestyle='--',label='Stiffness cycle 5')
# Axes.set_xlabel('Index (-)')
# Axes.set_ylabel('Stiffness (N/mm)')
# plt.legend()
# # plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
# plt.show()
# plt.close(Figure)

#%%
InitialStiffnessMaxIndex = np.where(InitialStiffness==max(InitialStiffness))[0][0]+MTSFirstCycleTroughs[-1]

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ShiftedARADisplacement[UltimateForceIndex:len(ShiftedMTSForces)],ShiftedMTSForces[UltimateForceIndex:],color=(0,0,0),label='Experiment')
Axes.plot(ShiftedARADisplacement[MTSFirstCycleTroughs[-1]: UltimateForceIndex],ShiftedMTSForces[MTSFirstCycleTroughs[-1]: UltimateForceIndex],color=(1,0,0),label='Stiffness Regression Range')
Axes.plot(ShiftedARADisplacement[InitialStiffnessMaxIndex:InitialStiffnessMaxIndex+RegressionRange],ShiftedMTSForces[InitialStiffnessMaxIndex:InitialStiffnessMaxIndex+RegressionRange],color=(0,0,1),label='Max Stiffness Range')
Axes.set_xlabel('Displacement (mm)')
Axes.set_ylabel('Force (N)')
plt.legend()
# plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
plt.show()
plt.close(Figure)

MaxDisplacement = min(ShiftedARADisplacement)
MaxDisplacementIndex = np.where(ShiftedARADisplacement==MaxDisplacement)[0][0]
RegressionRange = int((Junction-MaxDisplacementIndex)/10)

DamagedStiffness = np.zeros(Junction - MaxDisplacementIndex - RegressionRange)
for Index in range(Junction - MaxDisplacementIndex - RegressionRange):
    Index = Index + MaxDisplacementIndex

    RegressionDisplacementData = ShiftedARADisplacement[Index:Index + RegressionRange]
    RegressionForceData = ShiftedMTSForces[Index:Index + RegressionRange]

    DamagedStiffness[Index - MaxDisplacementIndex], Interception, RValue, PValue, StdError = linregress(RegressionDisplacementData, RegressionForceData)

# Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
# Axes.plot(DamagedStiffness,color=(0,0,0),label='Regression results')
# Axes.plot([0,len(DamagedStiffness)],[PostDamageStiffness[0],PostDamageStiffness[0]],color=(1,0,0),linestyle='--',label='Stiffness cycle 1')
# Axes.plot([0,len(DamagedStiffness)],[PostDamageStiffness[1],PostDamageStiffness[1]],color=(0,0,1),linestyle='--',label='Stiffness cycle 2')
# Axes.plot([0,len(DamagedStiffness)],[PostDamageStiffness[2],PostDamageStiffness[2]],color=(0,1,0),linestyle='--',label='Stiffness cycle 3')
# Axes.set_xlabel('Index (-)')
# Axes.set_ylabel('Stiffness (N/mm)')
# plt.legend()
# # plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
# plt.show()
# plt.close(Figure)

#%%
DamagedStiffnessMaxIndex = np.where(DamagedStiffness==max(DamagedStiffness))[0][0]+MaxDisplacementIndex

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(ShiftedARADisplacement[:MaxDisplacementIndex],ShiftedMTSForces[:MaxDisplacementIndex],color=(0,0,0),label='Experiment')
Axes.plot(ShiftedARADisplacement[MaxDisplacementIndex:Junction],ShiftedMTSForces[MaxDisplacementIndex:Junction],color=(1,0,0),label='Stiffness Regression Range')
Axes.plot(ShiftedARADisplacement[DamagedStiffnessMaxIndex:DamagedStiffnessMaxIndex+RegressionRange],ShiftedMTSForces[DamagedStiffnessMaxIndex:DamagedStiffnessMaxIndex+RegressionRange],color=(0,0,1),label='Max Stiffness Range')
Axes.set_xlabel('Displacement (mm)')
Axes.set_ylabel('Force (N)')
plt.legend()
# plt.savefig(os.path.join(ResultsPath,SampleName+'.png'))
plt.show()
plt.close(Figure)

#%%
SampleSignals = pd.DataFrame({'Displacement (mm)':ShiftedARADisplacement[:len(ShiftedMTSForces)],
                              'Force (N)':ShiftedMTSForces,
                              'Time (s)':ShiftedMTSTimes})
SampleSignals.to_csv(os.path.join(ResultsPath,SampleName+'.csv'),index=False)

Filter = SampleList['Internal ID'] == SampleName
Age = SampleList[Filter]['Age'].values[0]
Sex = SampleList[Filter]['Sex'].values[0]
Thickness1 =  SampleList[Filter]['Thickness Point 1 (mm)'].values[0]
Thickness2 =  SampleList[Filter]['Thickness Point 2 (mm)'].values[0]
Thickness3 =  SampleList[Filter]['Thickness Point 3 (mm)'].values[0]

SamplesData = SamplesData.append({'Sample':SampleName,
                                  'Age':Age,
                                  'Sex':Sex,
                                  'Thickness (mm)':np.array([Thickness1,Thickness2,Thickness3]),
                                  'Protocol 1 Peaks':MTSFirstCyclePeaks,
                                  'Protocol 1 Troughs':MTSFirstCycleTroughs,
                                  'Protocol 2 Peaks': MTSSecondCycleShiftedPeaks,
                                  'Protocol 2 Troughs': MTSSecondCycleShiftedTroughs,
                                  'Ultimate Load (N)':UltimateForce,
                                  'Ultimate Load Index':UltimateForceIndex,
                                  'Max Displacement (mm)':MaxDisplacement,
                                  'Max Displacement Index':MaxDisplacementIndex,
                                  'Protocol Junction Index':Junction,
                                  'Protocol 1 Stiffness (N/mm)':PreDamageStiffness,
                                  'Protocol 2 Stiffness (N/mm)':PostDamageStiffness,
                                  'Loading Max Stiffness (N/mm)':max(InitialStiffness),
                                  'Unloading Max Stiffness (N/mm)':max(DamagedStiffness)
                                  },ignore_index=True)
SamplesData.to_csv(os.path.join(ResultsPath,'ExperimentData.csv'),index=False)


# %%
