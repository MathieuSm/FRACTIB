#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """

    Script aimed to investigate how the reference corrdinate system is defined
    in order to put "correct" boundary conditions for the simulation

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: Match 2023
    """

#%% Imports
# Modules import

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import scipy.signal as sig
import matplotlib.pyplot as plt

#%% Funtions
# Define functions

def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results

def ReadDAT(File):

    """
    Read .dat file from abaqus and extract reference point data
    """

    try:
        with open(File) as F:
            Text = F.read()

            Values = []
            Condition = Text.find('U1') + 1
            Start = Text.find('U1')
            Steps, Increments = [], []
            Step, Increment = 1, 1
            while Condition:

                Inc = int(Text[Start-349:Start-347])
                if Inc < Increment:
                    Step += 1
                Increment = Inc
                Increments.append(Increment)
                Steps.append(Step)

                for i in range(6):
                    iStart = Start + 125 + 15*i
                    iStop = iStart + 14
                    Values.append(float(Text[iStart : iStop]))

                    iStart += Text[Start:].find('RF1')
                    iStop += Text[Start:].find('RF1')
                    Values.append(float(Text[iStart : iStop]))

                Start = iStop + Text[iStop:].find('U1')
                Condition = Text[iStop:].find('U1') + 1

            Values = np.array(Values)
            Cols = 12
            Rows = Values.size // Cols
            Values = np.reshape(Values,(Rows,Cols))

            ColNames = []
            for i in range(3):
                for V in ['U', 'F']:
                    ColNames.append(V + str(i+1))
            for i in range(3):
                for V in ['R', 'M']:
                    ColNames.append(V + str(i+1))

            Data = pd.DataFrame()
            for iName, Name in enumerate(ColNames):
                Data[Name] = np.concatenate([[0], Values[:,iName]])
            
            Data.columns = ['X', 'FX', 'Y', 'FY', 'Z', 'FZ', 'Phi', 'MX', 'Theta', 'MY', 'Psi', 'MZ']

            Data['Step'] = np.concatenate([[0],Steps])
            Data['Increment'] = np.concatenate([[0],Increments])

        return Data

    except FileNotFoundError:
        print('File' + File + 'does not exist')

        return


#%% Main
# Main code

def Main():

    WD, DD, SD, RD = SetDirectories('FRACTIB')
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))

    Index = 0
    Sample = SampleList.loc[Index, 'Internal ID']

    ExpData = pd.read_csv(RD / '02_Experiment' / Sample / 'MatchedSignals.csv')
    SimData = ReadDAT(str(RD / '03_hFE' / Sample / 'MaxLoad.dat'))

    # Truncate experiment to monotonic loading part
    Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
    Start = Peaks[4]
    Stop = ExpData['FZ'].idxmin()
    ExpData = ExpData[Start:Stop].reset_index(drop=True)

    # Interpolate data to match hFE increments
    InterpData = pd.DataFrame()
    Coords = np.concatenate([np.array([0]),ExpData['Z']])
    for C in ExpData.columns:
        Values = np.concatenate([np.array([0]),ExpData[C]])
        Interpolated = np.interp(SimData['Z'], Coords, Values)
        InterpData[C] = Interpolated

    # Compute cost after rotation
    for A in Angles:
        R = RotationMatrix(Psi=A/180*np.pi)
        RSys = np.dot(R, InterpData[['X','Y','Z']].values.T).T
        Costs.loc[(Sample,'dX'), A] = np.array(np.abs(RSys[:,0] - SimData['X']))
        Costs.loc[(Sample,'dY'), A] = np.array(np.abs(RSys[:,1] - SimData['Y']))

    # # Plot costs
    # Show.Signal([Costs.loc[Sample,A]['dX'] for A in Angles], Labels=[str(A)for A in Angles])



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