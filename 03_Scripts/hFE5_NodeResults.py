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
    Columns = pd.MultiIndex.from_product([['hFE','Experiment'],['MD', 'Fm', 'Fmd', 'S']])
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

        # Plot force displacement curves
        Show.FName = str(RD / '05_Comparison' / (Sample + '_Curve.png'))     
        Show.Signal([FEAData['Z'],ExpData['Z']],
                    [FEAData['FZ'], -ExpData['FZ']],
                    Axes=['Displacement (mm)', 'Force (N)'],
                    Labels=['hFE','Experiment'])
        
        # Store stiffess, force at max(ExpForce), max displacement
        Data.loc[Index]['Experiment','MD'] = ExpData['Z'].max()
        Data.loc[Index]['Experiment','Fm'] = abs(ExpData['FZ'].min())
        Data.loc[Index]['Experiment','Fmd'] = ExpData.loc[ExpData['FZ'].idxmin(),'Z']

        WindowWidth = ExpData['FZ'].idxmin() // 3
        X = ExpData.loc[:ExpData['FZ'].idxmin(), 'Z']
        Y = -ExpData.loc[:ExpData['FZ'].idxmin(), 'FZ']
        Data.loc[Index]['Experiment','S'] = Signal.MaxSlope(X, Y, WindowWidth)

        # Compute min force location
        Data.loc[Index]['hFE','MD'] = FEAData['Z'].max()
        Data.loc[Index]['hFE','Fm'] = FEAData['FZ'].max()
        Data.loc[Index]['hFE','Fmd'] = FEAData.loc[FEAData['FZ'].idxmax(),'Z']

        WindowWidth = FEAData['FZ'].idxmax() // 3
        if WindowWidth < 3:
            WindowWidth = 3
        X = FEAData.loc[:FEAData['FZ'].idxmax(), 'Z']
        Y = FEAData.loc[:FEAData['FZ'].idxmax(), 'FZ']
        Data.loc[Index]['hFE','S'] = Signal.MaxSlope(X, Y, WindowWidth)

    # Show.ShowPlot = True
    Show.FName = str(RD / '05_Comparison' / ('MaxForce.png')) 
    Show.OLS(Data['Experiment','Fm'].astype('float'),
             Data['hFE','Fm'].astype('float'),
             Labels=['Experiment (N)', 'hFE (N)'])
    
    Show.FName = str(RD / '05_Comparison' / ('Stiffness.png')) 
    Show.OLS(Data['Experiment','S'].astype('float'),
             Data['hFE','S'].astype('float'),
             Labels=['Experiment (N/mm)', 'hFE (N/mm)'])
    
    Show.FName = str(RD / '05_Comparison' / ('DispFmax.png')) 
    Show.OLS(Data['Experiment','Fmd'].astype('float'),
             Data['hFE','Fmd'].astype('float'),
             Labels=['Experiment (mm)', 'hFE (mm)'])

    Data.to_csv(str(RD / 'StrucuralResults.csv'))

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
    Parser.add_argument('-F', '--Folder', help='Root folder of the project', default='FRACTIB', type=str)
    Parser.add_argument('Sample', help='Sample to analyze (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)