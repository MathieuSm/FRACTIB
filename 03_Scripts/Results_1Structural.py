#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """

    Read .dat simulation results file and store force-displacement
    curve together with experimental data. Then, use force result
    to assess which increment to extract deformation gradient. Write
    ODB reader from template and execute it.

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import argparse
import numpy as np
import pandas as pd
from scipy import signal as sig
from Utils import SetDirectories, Time, Show, Abaqus, Signal

Show.ShowPlot = False


#%% Main
# Main code

def Main():

    CWD, DD, SD, RD = SetDirectories('FRACTIB')
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))

    Text = 'Nodal Results'
    Time.Process(1, Text)
    Variables = ['Max Displacement (mm)', 'Max Force (N)', 'Stiffness (N/mm)']
    Columns = pd.MultiIndex.from_product([['hFE','Experiment'],Variables])
    Data = pd.DataFrame(index=SampleList.index, columns=Columns)
    for Index, Sample in enumerate(SampleList['Internal ID']):

        Time.Update((Index + 1) / len(SampleList), Sample)
        FEADir = RD / '03_hFE' / Sample
        ExpDir = RD / '02_Experiment' / Sample

        # Read files
        FEAData = Abaqus.ReadDAT(str(FEADir / 'Simulation.dat'))
        ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))

        # Truncate experiment to monotonic loading part and set to 0
        Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=1)
        MaxForce = ExpData['FZ'].idxmin()
        MaxDisp = ExpData['Z'].idxmax()
        DeltaTime = 10
        DeltaIndex = np.argmin(np.abs(ExpData['T']-DeltaTime))
        Start = Peaks[Peaks < MaxForce - DeltaIndex][-1]
        Stop = Peaks[Peaks > MaxDisp][0]
        
        ExpData = ExpData[Start:Stop].reset_index(drop=True)
        ExpData -= ExpData.loc[0]

        # Truncate FEA if force became negative
        FEAData = FEAData[FEAData['FZ'] >= 0]

        # Plot force displacement curves
        Show.FName = str(RD / '05_Comparison' / (Sample + '_Curve.png'))     
        Show.Signal([FEAData['Z'],ExpData['Z']],
                    [FEAData['FZ'] / 1E3, -ExpData['FZ'] / 1E3],
                    Axes=['Displacement (mm)', 'Force (kN)'],
                    Labels=['hFE','Experiment'])
        
        # Store stiffess, force at max(ExpForce), max displacement
        Data.loc[Index]['Experiment',Variables[0]] = ExpData['Z'].max()
        Data.loc[Index]['Experiment',Variables[1]] = abs(ExpData['FZ'].min())

        WindowWidth = ExpData['FZ'].idxmin() // 3
        X = ExpData.loc[:ExpData['FZ'].idxmin(), 'Z']
        Y = -ExpData.loc[:ExpData['FZ'].idxmin(), 'FZ']
        Data.loc[Index]['Experiment',Variables[2]] = Signal.MaxSlope(X, Y, WindowWidth)

        # Compute min force location
        Data.loc[Index]['hFE',Variables[0]] = FEAData['Z'].max()
        Data.loc[Index]['hFE',Variables[1]] = FEAData['FZ'].max()

        WindowWidth = FEAData['FZ'].idxmax() // 3
        if WindowWidth < 3:
            WindowWidth = 3
        X = FEAData.loc[:FEAData['FZ'].idxmax(), 'Z']
        Y = FEAData.loc[:FEAData['FZ'].idxmax(), 'FZ']
        Data.loc[Index]['hFE',Variables[2]] = Signal.MaxSlope(X, Y, WindowWidth)

    Show.ShowPlot = True
    Show.FName = str(RD / '05_Comparison' / ('MaxForce.png')) 
    Show.OLS(Data['Experiment',Variables[1]].astype('float') / 1E3,
             Data['hFE',Variables[1]].astype('float') / 1E3,
             Labels=['Experiment (kN)', 'hFE (kN)'])
    
    Show.FName = str(RD / '05_Comparison' / ('Stiffness.png')) 
    Show.OLS(Data['Experiment',Variables[2]].astype('float') / 1E3,
             Data['hFE',Variables[2]].astype('float') / 1E3,
             Labels=['Experiment (kN/mm)', 'hFE (kN/mm)'])

    Data.to_csv(str(RD / 'StructuralResults.csv'))

    Time.Process(0, Text)

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

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main()