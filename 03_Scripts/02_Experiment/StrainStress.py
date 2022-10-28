import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

ResultsPath = '/home/mathieu/Documents/MscThesis/04_Results/FractureZoneAssessment'
SamplesData = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData/ExperimentData.csv'))


for SampleName in SamplesData['Sample']:
    SampleCurves = pd.read_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',SampleName+'.csv'))
    Filter = SamplesData['Sample']==SampleName

    SampleCurves['Apparent Stress (MPa)'] = SampleCurves['Force (N)'] / SamplesData[Filter]['Mean Area (mm2)'].values


    Thickness = SamplesData[Filter]['Thickness (mm)']
    Thickness = Thickness.values[0][1:-1]
    Thickness = np.array(Thickness.split()).astype('float')

    SampleCurves['Mean Strain (-)'] = SampleCurves['Displacement (mm)'] / Thickness.mean()

    SampleCurves.to_csv(os.path.join(ResultsPath,'01_Experiment/SamplesData',SampleName+'.csv'),index=False)
