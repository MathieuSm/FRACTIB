#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Modify the boundary conditions file to set the Z displacement
    to the maximum displacement applied during experiment

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
import pandas as pd
from Utils import *


#%% Main
# Main code

def Main(Arguments):

    ProcessTiming(1, 'Set loading')

    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)

    Data = pd.read_csv(str(DD / 'SampleList.csv'))

    for Index in Data.index:

        Sample = Data.loc[Index, 'Internal ID']
        
        FileName = str(RD / '02_Experiment' / Sample / 'MatchedSignals.csv')
        Experiment = pd.read_csv(FileName)

        FileName = str(RD / '03_hFE' / Sample / 'boundary_conditions_FZ_MAX.inp')
        with open(FileName) as File:

            Text = File.read()
            Line = Text.split('\n')[7]
            
        NewLine = Line[:16] + str(round(Experiment['Z'].max(),2))
        NewText = Text.replace(Line, NewLine)

        with open(FileName, 'w') as File:
            File.write(NewText)

        ProgressNext(Index / len(Data) * 10)

    ProcessTiming(0)

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

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)