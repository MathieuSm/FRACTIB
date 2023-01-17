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
Stats.Echo = False

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
        Slope = Stats.LinearRegression(X, Y)[1]
        Slopes.append(Slope)

    Stiffness = max(Slopes)

    return Stiffness

#%% Main
# Main code

def Main():

    # Set directories and read sample list
    CWD, DD, SD, RD = SetDirectories('FRACTIB')
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))

    # Set assessed variables
    Variables = ['Force','Displacement at Fmax','Stiffness']
    DataSets = ['hFE', 'Experiment']
    Columns = pd.MultiIndex.from_product([Variables, DataSets])
    Data = pd.DataFrame(index=SampleList.index, columns=Columns)

    # Compute variables from experiment and simulation    
    ProcessTiming(1, 'Get experiment and simulation max force, diplacement at max force and stiffness')
    for Index in SampleList.index:

        Sample = SampleList.loc[Index, 'Internal ID']
        FEADir = RD / '03_hFE' / Sample
        ExpDir = RD / '02_Experiment' / Sample

        # Read files
        FEAData = Abaqus.ReadDAT(str(FEADir / (Sample + '.dat')))
        ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))

        # Max force
        Data[Variables[0],'hFE'].loc[Index] = max(FEAData['FZ'])
        Data[Variables[0],'Experiment'].loc[Index] = max(ExpData['FZ'].abs())

        # Displacement at max force
        FEAIndex = FEAData['FZ'].idxmax()
        Data[Variables[1],'hFE'].loc[Index] = FEAData.loc[FEAIndex,'Z']
        ExpIndex = ExpData['FZ'].idxmin()
        Data[Variables[1],'Experiment'].loc[Index] = ExpData.loc[ExpIndex,'Z']

        # Stiffness
        F = FEAData['FZ'].values[:FEAIndex]
        D = FEAData['Z'].values[:FEAIndex]
        FEAStiff = ComputeStiffness(F, D)
        Data[Variables[2],'hFE'].loc[Index] = FEAStiff
        
        Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
        F = ExpData['FZ'].values[Peaks[4]:ExpIndex]
        D = ExpData['Z'].values[Peaks[4]:ExpIndex]
        ExpStiff = ComputeStiffness(-F, D)
        Data[Variables[2],'Experiment'].loc[Index] = ExpStiff 

        ProgressNext(Index / len(SampleList) * 20)
    ProcessTiming(0)

    # Assess relationships (linear regression)
    X = Data[Variables[0]]['Experiment'].values.astype(float)
    Y = Data[Variables[0]]['hFE'].values.astype(float)
    Labels = ['hFE (N/mm)', 'Experiment (N/mm)']
    Results = Stats.LinearRegression(X, Y, Labels=Labels, Plot=True)

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
    Parser.add_argument('File', help='File to process (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments.File)