#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script aimed to analyse the relationship of structural variables
    (force and stiffness namely) between experiment and simulation

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: January 2023
    """

#%% Imports
# Modules import

import argparse
from Utils import *
Show.ShowPlot = False

#%% Functions
# Define functions

def ComputeStiffness(Force, Displacement, RelativeRange=1/3, StepSize=None):

    Length = int(round(len(Force) * RelativeRange))
    if Length < 3:
        Length = 3

    if StepSize == None:
        StepSize = Length // 100
    
    if StepSize == 0:
        StepSize = 1

    Slopes = []
    for i in range(0, len(Force)-Length-StepSize+1, StepSize):
        Y = Force[i:i+Length]
        X = Displacement[i:i+Length]
        Show.FName = None
        Fit = Show.OLS(X, Y)
        Slopes.append(Fit.params['X'])

    Stiffness = max(Slopes)

    return Stiffness

#%% Main
# Main code

def Main(Arguments):

    # Set directories and read sample list
    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))
    ResultsDir = RD / '05_Comparison'

    # Set assessed variables
    Variables = ['Force','Displacement at Fmax','Stiffness']
    DataSets = ['hFE', 'Experiment']
    Columns = pd.MultiIndex.from_product([Variables, DataSets])
    Data = pd.DataFrame(index=SampleList.index, columns=Columns)

    # Compute variables from experiment and simulation    
    for Index in SampleList.index:

        Text = '\nGet structural results for sample ' + str(Index+1) + '/' + str(len(SampleList))
        ProcessTiming(1, Text)

        Sample = SampleList.loc[Index, 'Internal ID']
        FEADir = RD / '03_hFE' / Sample
        ExpDir = RD / '02_Experiment' / Sample

        # Read files
        FEAData = Abaqus.ReadDAT(str(FEADir / (Sample + '.dat')))
        ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))

        # Plot force-displacement and min force point
        MinForceIdx = FEAData['FZ'][1:].abs().idxmin()
        Xs = [FEAData['Z'], ExpData['Z']]
        Ys = [FEAData['FZ'], -ExpData['FZ']]
        Axes = ['Displacement (mm)', 'Force(N)']
        Labels = ['hFE', 'Experiment']

        Show.FName = str(ResultsDir / (Sample + '_FD'))
        Show.Signal(Xs, Ys, Axes=Axes, Labels=Labels, Points=[MinForceIdx,[]])
        Show.FName = None
        ProgressNext(1)

        # Max force
        Data[Variables[0],'hFE'].loc[Index] = max(FEAData['FZ'])
        Data[Variables[0],'Experiment'].loc[Index] = max(ExpData['FZ'].abs())

        # Displacement at max force
        FEAIndex = FEAData['FZ'].idxmax()
        Data[Variables[1],'hFE'].loc[Index] = FEAData.loc[FEAIndex,'Z']
        ExpIndex = ExpData['FZ'].idxmin()
        Data[Variables[1],'Experiment'].loc[Index] = ExpData.loc[ExpIndex,'Z']
        ProgressNext(2)

        # Stiffness
        F = FEAData['FZ'].values[:FEAIndex]
        D = FEAData['Z'].values[:FEAIndex]
        FEAStiff = ComputeStiffness(F, D)
        Data[Variables[2],'hFE'].loc[Index] = FEAStiff
        ProgressNext(3)
        
        Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
        F = ExpData['FZ'].values[Peaks[4]:ExpIndex]
        D = ExpData['Z'].values[Peaks[4]:ExpIndex]
        ExpStiff = ComputeStiffness(-F, D)
        Data[Variables[2],'Experiment'].loc[Index] = ExpStiff
        ProcessTiming(0)

    # Assess relationships (linear regression)
    X = Data[Variables[0]]['Experiment'].values.astype(float) / 1E3
    Y = Data[Variables[0]]['hFE'].values.astype(float) / 1E3
    Labels = ['hFE (kN)', 'Experiment (kN)']
    Show.FName = str(RD / 'UltimateForce')
    Results = Show.OLS(X, Y, Labels=Labels)

    X = Data[Variables[1]]['Experiment'].values.astype(float)
    Y = Data[Variables[1]]['hFE'].values.astype(float)
    Labels = ['hFE (mm)', 'Experiment (mm)']
    Show.FName = str(RD / 'DispUtlForce')
    Results = Show.OLS(X, Y, Labels=Labels)

    X = Data[Variables[2]]['Experiment'].values.astype(float) / 1E3
    Y = Data[Variables[2]]['hFE'].values.astype(float) / 1E3
    Labels = ['hFE (kN/mm)', 'Experiment (kN/mm)']
    Show.FName = str(RD / 'Stiffness')
    Results = Show.OLS(X, Y, Labels=Labels)
    Show.FName = None

    return

#%% Execution part
# Execution as main
if __name__ == '__main__':

    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)

    # Add long and short argument
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-V', '--Version', help='Show script version', action='version', version=SV)
    Parser.add_argument('--Folder', help='Root folder of the project (required)', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)