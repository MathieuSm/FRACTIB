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
from Utils import SetDirectories, Time, Show, Abaqus

Show.ShowPlot = False

#%% Functions


#%% Main
# Main code

def Main():

    CWD, DD, SD, RD = SetDirectories('FRACTIB')
    Data = pd.read_csv(str(DD / 'SampleList.csv'))

    Text = 'Nodal Results'
    Time.Process(1, Text)
    for Index, Sample in enumerate(Data['Internal ID']):

        Time.Update((Index + 1) / len(Data))
        FEADir = RD / '03_hFE' / Sample
        ExpDir = RD / '02_Experiment' / Sample

        # Read files
        FEAData = Abaqus.ReadDAT(str(FEADir / 'Experiment.dat'))
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

        # Sum up fea steps results
        Indices = FEAData.groupby('Step')['Increment'].idxmax()
        FEA = FEAData.loc[Indices].cumsum()

        # Plot force displacement curves
        Show.FName = str(RD / '05_Comparison' / (Sample + '_Curve.png'))     
        Show.Signal([FEA['Z'],ExpData['Z']],
                    [FEA['FZ'], -ExpData['FZ']],
                    Axes=['Displacement (mm)', 'Force (N)'],
                    Labels=['hFE','Experiment'],
                    Normalize=True)
        
        # Store stiffess, force at max(ExpForce), max displacement
        
        
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