#!/usr/bin/env python3

# 00 Imports

## Import libraries
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import linregress


CorrelationsData = pd.DataFrame()

# 00 Load Data

ResultsPath = '/home/mathieu/Documents/MscThesis/04_Results/FractureZoneAssessment'
SamplesData = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData/ExperimentData.csv'))



Reg1Stiffness = np.zeros((len(SamplesData),5))
for Index in SamplesData.index:
    Stiffness = SamplesData.loc[Index]['Protocol 1 Regression Stiffness (N/mm)']
    Stiffness = Stiffness[1:-1]
    Reg1Stiffness[Index,:] = np.array(Stiffness.split())

Reg2Stiffness = np.zeros((len(SamplesData),3))
for Index in SamplesData.index:
    Stiffness = SamplesData.loc[Index]['Protocol 2 Regression Stiffness (N/mm)']
    Stiffness = Stiffness[1:-1]
    Reg2Stiffness[Index,:] = np.array(Stiffness.split())

Protocol1Troughs = np.zeros((len(SamplesData),5))
for Index in SamplesData.index:
    Troughs = SamplesData.loc[Index]['Protocol 1 Troughs']
    Troughs = Troughs[1:-1]
    Protocol1Troughs[Index,:] = np.array(Troughs.split()).astype('float')

Thicknesses = np.zeros((len(SamplesData), 3))
for Index in SamplesData.index:
    Thickness = SamplesData.loc[Index]['Thickness (mm)']
    Thickness = Thickness[1:-1]
    Thicknesses[Index,:] = np.array(Thickness.split()).astype('float')


SamplesData['Max Strain (-)'] = SamplesData['Max Displacement (mm)']/Thicknesses.mean(axis=1)


# 01 Typical measurement protocol
MaxForceSample = SamplesData.loc[SamplesData['Ultimate Load (N)'].idxmin()]['Sample']
MaxSampleCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',MaxForceSample+'.csv'))

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MaxSampleCurves['Time (s)'],
          MaxSampleCurves['Force (N)']/min(MaxSampleCurves['Force (N)']),
          color=(0,0,1),label='Normalized force')
Axes.plot(MaxSampleCurves['Time (s)'],
          MaxSampleCurves['Displacement (mm)']/min(MaxSampleCurves['Displacement (mm)']),
          color=(1,0,0),label='Normalized displacement')
Axes.set_xlabel('Time (s)')
Axes.set_ylabel('Normalized Signal (-)')
plt.legend(loc='center',bbox_to_anchor=(0.5, 1.1),ncol=2)
plt.show()
plt.close(Figure)

LastTroughIndex = int(Protocol1Troughs[SamplesData['Ultimate Load (N)'].idxmin(),-1])

# 01 Max and min force displacement curves

MinForceSample = SamplesData.loc[SamplesData['Ultimate Load (N)'].idxmax()]['Sample']

MinSampleCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',MinForceSample+'.csv'))

Filter = SamplesData['Sample']==MaxForceSample
MaxCurveStart = int(SamplesData[Filter]['Protocol 1 Troughs'].values[0][-6:-1])
MaxCurveEnd = int(SamplesData[Filter]['Protocol Junction Index'].values[0])

Filter = SamplesData['Sample']==MaxForceSample
MinCurveStart = int(SamplesData[Filter]['Protocol 1 Troughs'].values[0][-6:-1])
MinCurveEnd = int(SamplesData[Filter]['Protocol Junction Index'].values[0])

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MaxSampleCurves['Displacement (mm)'].iloc[MaxCurveStart:MaxCurveEnd],
          MaxSampleCurves['Force (N)'].iloc[MaxCurveStart:MaxCurveEnd]/1e3,
          color=(1,0,0),label='Max force reached')
Axes.plot(MinSampleCurves['Displacement (mm)'].iloc[MinCurveStart:MinCurveEnd],
          MinSampleCurves['Force (N)'].iloc[MinCurveStart:MinCurveEnd]/1e3,
          color=(0,0,1),label='Min force reached')
Axes.set_xlabel('Displacement (mm)')
Axes.set_ylabel('Force (kN)')
plt.legend()
plt.show()
plt.close(Figure)


# Initial stress computation
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MaxSampleCurves['Time (s)'].iloc[:LastTroughIndex],
          MaxSampleCurves['Force (N)'].iloc[:LastTroughIndex]/min(MaxSampleCurves['Force (N)']),
          color=(0,0,1),label='Normalized force')
Axes.plot(MaxSampleCurves['Time (s)'].iloc[:LastTroughIndex],
          MaxSampleCurves['Displacement (mm)'].iloc[:LastTroughIndex]/min(MaxSampleCurves['Displacement (mm)']),
          color=(1,0,0),label='Normalized displacement')
Axes.set_xlabel('Time (s)')
Axes.set_ylabel('Normalized Signal (-)')
plt.legend(loc='center',bbox_to_anchor=(0.5, 1.1),ncol=2)
plt.show()
plt.close(Figure)


# Stiffness distribution
Figure, Axes = plt.subplots(1, 1, figsize=(8.5, 4.5),dpi=100)
Axes.plot([],color=(0,0,1),marker='.',linestyle='none',label='5 cycles measurements')
Axes.plot(SamplesData['Sample'],Reg1Stiffness/1e3,color=(0,0,1),marker='.',linestyle='none')
Axes.plot(SamplesData['Sample'],Reg1Stiffness.mean(axis=1)/1e3,color=(0,0,1),marker='o',linestyle='none', label='Measurement mean')
Axes.plot(SamplesData['Sample'],np.median(Reg1Stiffness,axis=1)/1e3,color=(1,0,0),marker='o',linestyle='none', label='Measurement median')
Axes.plot(SamplesData['Sample'],SamplesData['Loading Max Stiffness (N/mm)']/1e3,color=(0,0,0),marker='x',linestyle='none', label='Max loading stiffness')
Axes.set_xlabel('Sample ID')
Axes.set_ylabel('Stiffness (kN/mm)')
plt.legend(loc='upper left')
plt.xticks(rotation='vertical')
plt.show()
plt.close(Figure)


# Correlation with hFE stiffness
Slope, Intercept, RValue, PValue, StdErr = linregress(Reg1Stiffness.mean(axis=1)/1e3,SamplesData['hFE Stiffness']/1e3)
x = np.array([min(Reg1Stiffness.mean(axis=1)/1e3),max(Reg1Stiffness.mean(axis=1)/1e3)])

CorrelationsData = CorrelationsData.append({'Independent variable':'Preconditioning Mean Stiffness',
                                            'Dependent variable':'hFE Stiffness',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot([0,170],[0,170],color=(0,0,0),linestyle='--',label='Diagonal')
Axes.plot(Reg1Stiffness.mean(axis=1)/1e3,
          SamplesData['hFE Stiffness']/1e3,
          marker='o',linestyle='none',fillstyle='none',color=(1,0,0),label='Data Points')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('Experiment Stiffness (kN/mm)')
Axes.set_ylabel('hFE Stiffness (kN/mm)')
Axes.set_xlim([0,170])
Axes.set_ylim([0,170])
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(20,110))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)


Slope, Intercept, RValue, PValue, StdErr = linregress(np.median(Reg1Stiffness,axis=1)/1e3,SamplesData['hFE Stiffness']/1e3)
x = np.array([min(np.median(Reg1Stiffness,axis=1)/1e3),max(np.median(Reg1Stiffness,axis=1)/1e3)])

CorrelationsData = CorrelationsData.append({'Independent variable':'Preconditioning Median Stiffness',
                                            'Dependent variable':'hFE Stiffness',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot([0,170],[0,170],color=(0,0,0),linestyle='--',label='Diagonal')
Axes.plot(np.median(Reg1Stiffness,axis=1)/1e3,
          SamplesData['hFE Stiffness']/1e3,
          marker='o',linestyle='none',fillstyle='none',color=(1,0,0),label='Data Points')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('Experiment Stiffness (kN/mm)')
Axes.set_ylabel('hFE Stiffness (kN/mm)')
Axes.set_xlim([0,170])
Axes.set_ylim([0,170])
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(20,110))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)


Slope, Intercept, RValue, PValue, StdErr = linregress(SamplesData['Loading Max Stiffness (N/mm)']/1e3,SamplesData['hFE Stiffness']/1e3)
x = np.array([min(SamplesData['Loading Max Stiffness (N/mm)']/1e3),max(SamplesData['Loading Max Stiffness (N/mm)']/1e3)])

CorrelationsData = CorrelationsData.append({'Independent variable':'Loading Max Stiffness',
                                            'Dependent variable':'hFE Stiffness',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot([0,170],[0,170],color=(0,0,0),linestyle='--',label='Diagonal')
Axes.plot(SamplesData['Loading Max Stiffness (N/mm)']/1e3,
          SamplesData['hFE Stiffness']/1e3,
          marker='o',linestyle='none',fillstyle='none',color=(1,0,0),label='Data Points')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('Experiment Stiffness (kN/mm)')
Axes.set_ylabel('hFE Stiffness (kN/mm)')
Axes.set_xlim([0,170])
Axes.set_ylim([0,170])
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(10,80))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)



# Post-damage stiffness
Figure, Axes = plt.subplots(1, 1, figsize=(8.5, 4.5),dpi=100)
Axes.plot([],color=(0,0,1),marker='.',linestyle='none',label='3 cycles measurements')
Axes.plot(SamplesData['Sample'],Reg2Stiffness/1e3,color=(0,0,1),marker='.',linestyle='none')
Axes.plot(SamplesData['Sample'],Reg2Stiffness.mean(axis=1)/1e3,color=(0,0,1),marker='o',linestyle='none', label='Measurement mean')
Axes.plot(SamplesData['Sample'],np.median(Reg2Stiffness,axis=1)/1e3,color=(1,0,0),marker='o',linestyle='none', label='Measurement median')
Axes.plot(SamplesData['Sample'],SamplesData['Unloading Max Stiffness (N/mm)']/1e3,color=(0,0,0),marker='x',linestyle='none', label='Max unloading stiffness')
Axes.set_xlabel('Sample ID')
Axes.set_ylabel('Stiffness (kN/mm)')
plt.legend(loc='upper left')
plt.xticks(rotation='vertical')
plt.show()
plt.close(Figure)

MeanSlope, MeanIntercept, MeanRValue, MeanPValue, MeanStdErr = linregress(SamplesData['Unloading Max Stiffness (N/mm)']/1e3,Reg2Stiffness.mean(axis=1)/1e3)
MeanX = np.array([SamplesData['Unloading Max Stiffness (N/mm)'].min()/1e3,SamplesData['Unloading Max Stiffness (N/mm)'].max()/1e3])
MedianSlope, MedianIntercept, MedianRValue, MedianPValue, MedianStdErr = linregress(SamplesData['Loading Max Stiffness (N/mm)']/1e3,np.median(Reg2Stiffness,axis=1)/1e3)
MedianX = np.array([SamplesData['Unloading Max Stiffness (N/mm)'].min()/1e3,SamplesData['Unloading Max Stiffness (N/mm)'].max()/1e3])

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=1000)
Axes.plot([0,110],[0,110],color=(0,0,0),linestyle='--',label='Diagonal')
Axes.plot(SamplesData['Unloading Max Stiffness (N/mm)']/1e3,
          Reg2Stiffness.mean(axis=1)/1e3,
          marker='o',linestyle='none',fillstyle='none',color=(1,0,0),label='Mean')
Axes.plot(MeanX, MeanX*MeanSlope+MeanIntercept,
          linestyle='--',color=(1,0,0),label='Mean Linear Regression')
Axes.annotate('Slope : ' + str(np.round(MeanSlope,2)),(10,70))
# Axes.plot(SamplesData['Unloading Max Stiffness (N/mm)']/1e3,
#           np.median(Reg2Stiffness,axis=1)/1e3,
#           marker='o',linestyle='none',fillstyle='none',color=(1,0,0),label='Median')
# Axes.plot(MedianX, MedianX*MedianSlope+MedianIntercept,
#           linestyle='--',color=(1,0,0),label='Median Linear Regression')
# Axes.annotate('Slope : ' + str(np.round(MedianSlope,2)),(10,70))
Axes.set_xlabel('Unloading Max Stiffness (kN/mm)')
Axes.set_ylabel('Unloading stiffness (kN/mm)')
Axes.set_xlim([0,110])
Axes.set_ylim([0,110])
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)


# Stiffness ratio vs strain reached and their distribution
SamplesData['Regression Stiffness Ratio (-)'] = np.median(Reg2Stiffness,axis=1) / SamplesData['Loading Max Stiffness (N/mm)']

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SamplesData['Max Strain (-)'],
          SamplesData['Regression Stiffness Ratio (-)'],color=(0,0,0),marker='o',linestyle='none',fillstyle='none')
Axes.set_xlabel('Max strain applied (-)')
Axes.set_ylabel('Stiffness ratio (-)')
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=1000)
SamplesData.boxplot(column='Max Strain (-)',ax=Axes,grid=False,vert=False,
                    color=dict(boxes='b', whiskers='k', medians='r', caps='k'),
                    boxprops=dict(linestyle='-'),
                    medianprops=dict(linestyle='-'),
                    whiskerprops=dict(linestyle='--')
                    )
Axes.set_yticks([0])
Axes.set_xlabel('Max strain applied (-)')
plt.show()
plt.close(Figure)


Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=1000)
SamplesData.boxplot(column='Regression Stiffness Ratio (-)',ax=Axes,grid=False,vert=True,
                    color=dict(boxes='b', whiskers='k', medians='r', caps='k'),
                    boxprops=dict(linestyle='-'),
                    medianprops=dict(linestyle='-'),
                    whiskerprops=dict(linestyle='--')
                    )
Axes.set_ylabel('Stiffness ratio (-)')
Axes.set_xticks([0])
plt.show()
plt.close(Figure)


# Analysis of outliers
MaxRatioSample = SamplesData.loc[SamplesData['Regression Stiffness Ratio (-)'].idxmax()]['Sample']
MaxRatioCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',MaxRatioSample+'.csv'))

Filter = SamplesData['Sample']==MaxRatioSample
MaxRatioCurveStart = int(SamplesData[Filter]['Protocol 1 Troughs'].values[0][-6:-1])
MaxRatioCurveEnd = int(SamplesData[Filter]['Protocol Junction Index'].values[0])

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MaxRatioCurves['Displacement (mm)'],
          MaxRatioCurves['Force (N)']/1e3,
          color=(1,0,0),label='Max force reached')
Axes.set_xlabel('Displacement (mm)')
Axes.set_ylabel('Force (kN)')
plt.show()
plt.close(Figure)

# Ultimate load - hFE vs Experiment
Slope, Intercept, RValue, PValue, StdErr = linregress(SamplesData['Ultimate Load (N)']/1e3,SamplesData['hFE Ultimate Load (N)']/1e3)
x = np.array([min(SamplesData['Ultimate Load (N)']/1e3),max(SamplesData['Ultimate Load (N)']/1e3)])

CorrelationsData = CorrelationsData.append({'Independent variable':'Ultimate Load',
                                            'Dependent variable':'hFE Ultimate Load',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot([0,-25],[0,-25],color=(0,0,0),linestyle='--',label='Diagonal')
Axes.plot(SamplesData['Ultimate Load (N)']/1e3,
          SamplesData['hFE Ultimate Load (N)']/1e3,
          marker='o',linestyle='none',fillstyle='none',color=(1,0,0),label='Data Points')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('Experiment Ultimate Load (kN)')
Axes.set_ylabel('hFE Ultimate Load (kN)')
Axes.set_xlim([-25,0])
Axes.set_ylim([-25,0])
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(-22,-9))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)


# hFE vs Experiment curve
SampleList = pd.read_csv('/home/mathieu/Documents/MscThesisData/FractureZoneAssessment/SampleList.csv')

SampleIndex = 6
ExampleSample = SampleList['Internal ID'].loc[SampleIndex]

HRpQCT = SampleList[SampleList['Internal ID']==ExampleSample]['HrpQCT File 2 number'].values[0]
hFE = pd.read_csv(os.path.join('/home/mathieu/Documents/MscThesis/04_Results/FractureZoneAssessment/02_StandardhFE',ExampleSample,'C000'+str(HRpQCT)+'_HFE.CSV'))
hFE.columns = ['Displacement (mm)', 'Force (N)']
hFEData = pd.read_csv(os.path.join('/home/mathieu/Documents/MscThesis/04_Results/FractureZoneAssessment/02_StandardhFE',ExampleSample,'C000'+str(HRpQCT)+'_FEARESULTS_HFE_XT2_STD.TXT'),sep='\t')
hFEData['F.Load']

SampleCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',ExampleSample+'.csv'))
CurveStart = int(Protocol1Troughs[SampleIndex,-1])
CurveEnd = int(SampleCurves['Force (N)'].idxmin())


Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SampleCurves['Displacement (mm)'].iloc[CurveStart:CurveEnd],
          SampleCurves['Force (N)'].iloc[CurveStart:CurveEnd]/1e3,
          color=(1,0,0),label='Experiment')
Axes.plot(hFE['Displacement (mm)'],hFE['Force (N)']/1e3,color=(0,0,1),label='hFE')
Axes.set_xlabel('Displacement(mm)')
Axes.set_ylabel('Force (kN)')
plt.legend()
plt.show()
plt.close(Figure)
SamplesData['Ultimate Load (N)']



# Max and min stress strain curves
SamplesData['Apparent Strength (MPa)'] = SamplesData['Ultimate Load (N)']/SamplesData['Mean Area (mm2)']
MaxForceSample = SamplesData.loc[SamplesData['Apparent Strength (MPa)'].idxmin()]['Sample']
MinForceSample = SamplesData.loc[SamplesData['Apparent Strength (MPa)'].idxmax()]['Sample']

MaxSampleCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',MaxForceSample+'.csv'))
MinSampleCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',MinForceSample+'.csv'))

Filter = SamplesData['Sample']==MaxForceSample
MaxCurveStart = int(SamplesData[Filter]['Protocol 1 Troughs'].values[0][-6:-1])
MaxCurveEnd = int(SamplesData[Filter]['Protocol Junction Index'].values[0])

Filter = SamplesData['Sample']==MinForceSample
MinCurveStart = int(SamplesData[Filter]['Protocol 1 Troughs'].values[0][-6:-1])
MinCurveEnd = int(SamplesData[Filter]['Protocol Junction Index'].values[0])

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(MaxSampleCurves['Mean Strain (-)'].iloc[MaxCurveStart:MaxCurveEnd],
          MaxSampleCurves['Apparent Stress (MPa)'].iloc[MaxCurveStart:MaxCurveEnd],
          color=(1,0,0),label='Max stress reached')
Axes.plot(MinSampleCurves['Mean Strain (-)'].iloc[MinCurveStart:MinCurveEnd],
          MinSampleCurves['Apparent Stress (MPa)'].iloc[MinCurveStart:MinCurveEnd],
          color=(0,0,1),label='Min stress reached')
Axes.set_xlabel('Apparent Strain (-)')
Axes.set_ylabel('Apparent Stress (MPa)')
plt.legend()
plt.show()
plt.close(Figure)




# Extensive properties correlation: BMC vs Stiffness
Slope, Intercept, RValue, PValue, StdErr = linregress(SamplesData['BMC (HA mg)'],SamplesData['Loading Max Stiffness (N/mm)']/1e3)
x = np.array([min(SamplesData['BMC (HA mg)']),max(SamplesData['BMC (HA mg)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMC',
                                            'Dependent variable':'Stiffness',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SamplesData['BMC (HA mg)'],SamplesData['Loading Max Stiffness (N/mm)']/1e3,
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('BMC (HA mg)')
Axes.set_ylabel('Stiffness(kN/mm)')
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(8000,20))
Axes.annotate('Standard error = ' + str(np.round(StdErr,3)),(8000,15))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)


SlopeFiltered, InterceptFiltered, RValueFiltered, PValueFiltered, StdErrFiltered = linregress(SamplesData['BMC (HA mg)'].drop(8),SamplesData['Loading Max Stiffness (N/mm)'].drop(8)/1e3)
xFiltered = np.array([min(SamplesData['BMC (HA mg)']),max(SamplesData['BMC (HA mg)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMC Filtered',
                                            'Dependent variable':'Stiffness Filtered',
                                            'RSquare':RValueFiltered**2,
                                            'Standard error':StdErrFiltered}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=1000)
Axes.plot(SamplesData['BMC (HA mg)'],SamplesData['Loading Max Stiffness (N/mm)']/1e3,
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.plot(SamplesData['BMC (HA mg)'].loc[8],SamplesData['Loading Max Stiffness (N/mm)'].loc[8]/1e3,
          linestyle='none',marker='x',color=(0,0,1),fillstyle='none',label='Outlier',markersize=10)
Axes.plot(x, x*SlopeFiltered+InterceptFiltered,linestyle='--',color=(0,0,1),label='Linear Regression Filtered')
Axes.set_xlabel('BMC (HA mg)')
Axes.set_ylabel('Stiffness(kN/mm)')
Axes.annotate('Complete :',(7000,27.5))
# Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(9000,30))
Axes.annotate('R$^2$ = 0.920',(9000,30))
Axes.annotate('Standard error = ' + str(np.round(StdErr,4)),(9000,25))
Axes.annotate('Filtered :',(7000,12.5))
Axes.annotate('R$^2$ = ' + str(np.round(RValueFiltered**2,3)),(9000,15))
Axes.annotate('Standard error = ' + str(np.round(StdErrFiltered,4)),(9000,10))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)



# Extensive properties correlation: BMC vs Ultimate load
Slope, Intercept, RValue, PValue, StdErr = linregress(SamplesData['BMC (HA mg)'],SamplesData['Ultimate Load (N)']/1e3)
x = np.array([min(SamplesData['BMC (HA mg)']),max(SamplesData['BMC (HA mg)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMC',
                                            'Dependent variable':'Ultimate Load',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SamplesData['BMC (HA mg)'],SamplesData['Ultimate Load (N)']/1e3,
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('BMC (HA mg)')
Axes.set_ylabel('Ultimate Load (kN)')
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(3500,-20))
Axes.annotate('Standard error = ' + str(np.round(StdErr,4)),(3500,-21.5))
plt.legend(loc='upper right')
plt.show()
plt.close(Figure)


SlopeFiltered, InterceptFiltered, RValueFiltered, PValueFiltered, StdErrFiltered = linregress(SamplesData['BMC (HA mg)'].drop(8),SamplesData['Ultimate Load (N)'].drop(8)/1e3)
xFiltered = np.array([min(SamplesData['BMC (HA mg)']),max(SamplesData['BMC (HA mg)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMC Filtered',
                                            'Dependent variable':'Ultimate Load Filtered',
                                            'RSquare':RValueFiltered**2,
                                            'Standard error':StdErrFiltered}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=1000)
Axes.plot(SamplesData['BMC (HA mg)'],SamplesData['Ultimate Load (N)']/1e3,
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.plot(SamplesData['BMC (HA mg)'].loc[8],SamplesData['Ultimate Load (N)'].loc[8]/1e3,
          linestyle='none',marker='x',color=(0,0,1),fillstyle='none',label='Outlier',markersize=10)
Axes.plot(x, x*SlopeFiltered+InterceptFiltered,linestyle='--',color=(0,0,1),label='Linear Regression Filtered')
Axes.set_xlabel('BMC (HA mg)')
Axes.set_ylabel('Ultimate Load (kN)')
Axes.annotate('Complete :',(3000,-17.3))
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(5000,-16.7))
Axes.annotate('Standard error = ' + str(np.round(StdErr,5)),(5000,-17.9))
Axes.annotate('Filtered :',(3000,-20.9))
Axes.annotate('R$^2$ = ' + str(np.round(RValueFiltered**2,3)),(5000,-20.3))
Axes.annotate('Standard error = ' + str(np.round(StdErrFiltered,5)),(5000,-21.5))
plt.legend(loc='upper right')
plt.show()
plt.close(Figure)




# Intensive properties correlation: BMD vs Apparent Modulus
Slope, Intercept, RValue, PValue, StdErr = linregress(SamplesData['vBMD (HA g/cm3)'],SamplesData['Apparent Modulus (MPa)']/1e3)
x = np.array([min(SamplesData['vBMD (HA g/cm3)']),max(SamplesData['vBMD (HA g/cm3)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMD',
                                            'Dependent variable':'Apparent Modulus',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SamplesData['vBMD (HA g/cm3)'],SamplesData['Apparent Modulus (MPa)']/1e3,
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('vBMD (HA g/cm3)')
Axes.set_ylabel('Apparent Modulus (GPa)')
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(110,2.05))
Axes.annotate('Standard error = ' + str(np.round(StdErr,3)),(110,1.90))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)



SlopeFiltered, InterceptFiltered, RValueFiltered, PValueFiltered, StdErrFiltered = linregress(SamplesData['vBMD (HA g/cm3)'].drop(8),SamplesData['Apparent Modulus (MPa)'].drop(8)/1e3)
xFiltered = np.array([min(SamplesData['vBMD (HA g/cm3)']),max(SamplesData['vBMD (HA g/cm3)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMD Filtered',
                                            'Dependent variable':'Apparent Modulus Filtered',
                                            'RSquare':RValueFiltered**2,
                                            'Standard error':StdErrFiltered}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SamplesData['vBMD (HA g/cm3)'],SamplesData['Apparent Modulus (MPa)']/1e3,
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.plot(SamplesData['vBMD (HA g/cm3)'].loc[8],SamplesData['Apparent Modulus (MPa)'].loc[8]/1e3,
          linestyle='none',marker='x',color=(0,0,1),fillstyle='none',label='Outlier',markersize=10)
Axes.plot(x, x*SlopeFiltered+InterceptFiltered,linestyle='--',color=(0,0,1),label='Linear Regression Filtered')
Axes.set_xlabel('vBMD (HA g/cm3)')
Axes.set_ylabel('Apparent Modulus (GPa)')
Axes.annotate('Complete :',(110,2.475))
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(160,2.55))
Axes.annotate('Standard error = ' + str(np.round(StdErr,4)),(160,2.40))
Axes.annotate('Filtered :',(110,2.025))
Axes.annotate('R$^2$ = ' + str(np.round(RValueFiltered**2,3)),(160,2.1))
Axes.annotate('Standard error = ' + str(np.round(StdErrFiltered,4)),(160,1.95))
plt.legend(loc='lower right')
plt.show()
plt.close(Figure)



# Intensive properties correlation: BMC vs Ultimate load
Slope, Intercept, RValue, PValue, StdErr = linregress(SamplesData['vBMD (HA g/cm3)'],SamplesData['Apparent Strength (MPa)'])
x = np.array([min(SamplesData['vBMD (HA g/cm3)']),max(SamplesData['vBMD (HA g/cm3)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMD',
                                            'Dependent variable':'Apparent Strength',
                                            'RSquare':RValue**2,
                                            'Standard error':StdErr}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(SamplesData['vBMD (HA g/cm3)'],SamplesData['Apparent Strength (MPa)'],
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.set_xlabel('vBMD (HA g/cm3)')
Axes.set_ylabel('Apparent Strength (MPa)')
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(110,13))
Axes.annotate('Standard error = ' + str(np.round(StdErr,3)),(110,12))
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)


SlopeFiltered, InterceptFiltered, RValueFiltered, PValueFiltered, StdErrFiltered = linregress(SamplesData['vBMD (HA g/cm3)'].drop(8),SamplesData['Apparent Strength (MPa)'].drop(8))
xFiltered = np.array([min(SamplesData['vBMD (HA g/cm3)']),max(SamplesData['vBMD (HA g/cm3)'])])

CorrelationsData = CorrelationsData.append({'Independent variable':'BMD Filtered',
                                            'Dependent variable':'Apparent Strength Filtered',
                                            'RSquare':RValueFiltered**2,
                                            'Standard error':StdErrFiltered}, ignore_index=True)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=1000)
Axes.plot(SamplesData['vBMD (HA g/cm3)'],SamplesData['Apparent Strength (MPa)'],
          linestyle='none',marker='o',color=(0,0,0),fillstyle='none',label='Data')
Axes.plot(x, x*Slope+Intercept,linestyle='--',color=(1,0,0),label='Linear Regression')
Axes.plot(SamplesData['vBMD (HA g/cm3)'].loc[8],SamplesData['Apparent Strength (MPa)'].loc[8],
          linestyle='none',marker='x',color=(0,0,1),fillstyle='none',label='Outlier',markersize=10)
Axes.plot(x, x*SlopeFiltered+InterceptFiltered,linestyle='--',color=(0,0,1),label='Linear Regression Filtered')
Axes.set_xlabel('vBMD (HA g/cm3)')
Axes.set_ylabel('Apparent Strength (MPa)')
Axes.annotate('Complete :',(110,16.5))
Axes.annotate('R$^2$ = ' + str(np.round(RValue**2,3)),(160,17))
Axes.annotate('R$^2$ = 0.790',(160,17))
Axes.annotate('Standard error = ' + str(np.round(StdErr,4)),(160,16))
Axes.annotate('Filtered :',(110,13.5))
Axes.annotate('R$^2$ = ' + str(np.round(RValueFiltered**2,3)),(160,14))
Axes.annotate('Standard error = ' + str(np.round(StdErrFiltered,4)),(160,13))
plt.legend(loc='lower right')
plt.show()
plt.close(Figure)