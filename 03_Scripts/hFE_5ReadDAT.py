#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Read .dat file of simulations, rotate signals into experiment
    coordinate system and compare signals

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


#%% Functions
# Define functions

def AFunction(Argument):

    return Something


#%% Classes
# Define classes

class Arguments():

    def __init__(self):
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    # Set folders
    CW, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    hFE = Results / '03_hFE'
    Exp = Results /  '02_Experiment'
    Folders = [F for F in os.listdir(str(hFE)) if os.path.isdir(str(hFE/F))]

    # Read data
    hFEData = []
    ExpData = []
    X, Y = ['Z', 'FZ']

    for F in Folders:
        Dexp = Abaqus.ReadDAT(str(hFE / F / (F + '.dat')))
        hFEData.append(Dexp)
        Dsim = pd.read_csv(str(Exp / F / 'MatchedSignals.csv'))
        ExpData.append(Dsim)

        Xs = [Dexp[X], Dsim[X]]
        Ys = [Dexp[Y], -Dsim[Y]]
        Axes = ['Displacement (mm)', 'Force(N)']
        Labels = ['hFE', 'Experiment']

        Show.FName = str(hFE / F / 'ForceDisplacement.png')
        Show.Signal(Xs, Ys, Axes=Axes, Labels=Labels)
        del Show.FName



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
    Parser.add_argument('-F', '--Folder', help='Root folder of the project', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)