#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Add a step at the end of the .inp file with BCs to come back
    to zeros displacement according to the maximum displacement
    reached in loading simulation

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
from Utils import *


#%% For testing purpose
class Arguments:

    def __init__(self):
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    ProcessTiming(1, 'Set full simulation')

    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)

    Data = pd.read_csv(str(DD / 'SampleList.csv'))

    for Index in Data.index:

        Sample = Data.loc[Index, 'Internal ID']
        SamplePath = RD / '03_hFE' / Sample
        
        # Read max displacement reached
        FileName = str(SamplePath / (Sample + '.dat'))
        Simulation = Abaqus.ReadDAT(FileName)
        MaxDisp = Simulation.loc[Simulation['Z'].idxmax(),'Z']

        # Replace in original loading BCs file
        FileName = str(SamplePath / 'boundary_conditions_FZ_MAX.inp')
        with open(FileName) as File:

            Text = File.read()
            Line = Text.split('\n')[7]
            
        NewLine = Line[:16] + str(round(MaxDisp - 0.005,2))
        NewText = Text.replace(Line, NewLine)

        with open(FileName, 'w') as File:
            File.write(NewText)

        # Write new BCs file for unloading
        FileName = str(SamplePath / 'boundary_conditions_FZ_MIN.inp')
            
        NewLine = Line[:16] + str(-round(MaxDisp - 0.005,2))
        NewText = Text.replace(Line, NewLine)

        with open(FileName, 'w') as File:
            File.write(NewText)

        # Add unloading step to simulation input file
        uCTFile = 'C000' + str(Data.loc[Index, 'MicroCT pretest file number']) + '_DOWNSCALED_00_FZ_MAX.inp'
        SimFile = str(SamplePath / uCTFile)

        with open(SimFile) as File:
            Text = File.read()
            Start = Text.find('*STEP')
            NewText = Text[Start:]
        
        NewText = '**\n' + NewText
        NewText = NewText.replace('boundary_conditions_FZ_MAX.inp',FileName[-30:])
        with open(SimFile,'a') as File:
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