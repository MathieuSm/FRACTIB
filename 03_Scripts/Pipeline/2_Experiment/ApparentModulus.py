import os
import pandas as pd
import matplotlib.pyplot as plt


ResultsPath = '/home/mathieu/Documents/MscThesis/04_Results/FractureZoneAssessment'
SamplesData = pd.read_csv(os.path.join(ResultsPath, '01_Experiment/SamplesData/ExperimentData.csv'))

Data = pd.DataFrame()

for Sample in SamplesData['Sample']:
    # Sample = SamplesData['Sample'][0]
    SampleCurves = pd.read_csv(os.path.join(ResultsPath, '01_Experiment/SamplesData', Sample + '.csv'))

    Filter = SamplesData['Sample'] == Sample
    CurveStart = int(SamplesData[Filter]['Protocol 1 Troughs'].values[0][-6:-1])
    CurveEnd = int(SamplesData[Filter]['Ultimate Load Index'].values[0])

    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Axes.plot(SampleCurves['Mean Strain (-)'].iloc[CurveStart:CurveEnd],
    #           SampleCurves['Apparent Stress (MPa)'].iloc[CurveStart:CurveEnd],
    #           color=(0,0,0))
    # Axes.set_xlabel('Mean Strain (-)')
    # Axes.set_ylabel('Apparent Stress (MPa)')
    # plt.show()
    # plt.close(Figure)

    RegressionRange = int((CurveEnd - CurveStart) / 3)
    ApparentModulus = np.zeros((CurveEnd - CurveStart) - RegressionRange)

    for Index in range(CurveStart, CurveEnd - RegressionRange):
        RegressionStrainData = SampleCurves['Mean Strain (-)'].iloc[Index:Index + RegressionRange]
        RegressionStressData = SampleCurves['Apparent Stress (MPa)'].iloc[Index:Index + RegressionRange]

        ApparentModulus[Index - CurveStart], Interception, RValue, PValue, StdError = linregress(RegressionStrainData,
                                                                                                 RegressionStressData)

    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Axes.plot(ApparentModulus,
    #           color=(0,0,0))
    # Axes.set_xlabel('Start index (-)')
    # Axes.set_ylabel('Apparent Modulus (MPa)')
    # plt.show()
    # plt.close(Figure)

    ApparentModulusMaxIndex = np.where(ApparentModulus == max(ApparentModulus))[0][0] + CurveStart

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.plot(SampleCurves['Mean Strain (-)'].iloc[CurveStart:CurveEnd],
              SampleCurves['Apparent Stress (MPa)'].iloc[CurveStart:CurveEnd],
              color=(0, 0, 0), label='Loading curve')
    Axes.plot(SampleCurves['Mean Strain (-)'].iloc[ApparentModulusMaxIndex:ApparentModulusMaxIndex + RegressionRange],
              SampleCurves['Apparent Stress (MPa)'].iloc[
              ApparentModulusMaxIndex:ApparentModulusMaxIndex + RegressionRange],
              color=(1, 0, 0), label='Max regression range')
    Axes.set_xlabel('Mean Strain (-)')
    Axes.set_ylabel('Apparent Stress (MPa)')
    plt.legend()
    plt.show()
    plt.close(Figure)

    Data = Data.append({'Apparent Modulus (MPa)': ApparentModulus[ApparentModulusMaxIndex - CurveStart]},
                       ignore_index=True)

SamplesData['Apparent Modulus (MPa)'] = Data['Apparent Modulus (MPa)']
SamplesData.to_csv(os.path.join(ResultsPath, '01_Experiment/SamplesData/ExperimentData.csv'), index=False)