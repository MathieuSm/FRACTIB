#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    This script read a dat file from abaqus and output
    displacements, rotations, forces, and moments of
    the hFE analysis reference node

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Modules import
import argparse
import numpy as np
import pandas as pd

#%% Main function
def Main(File):

    with open(File) as F:
        Text = F.read()

        Values = []
        Condition = Text.find('U1') + 1
        Start = Text.find('U1')
        while Condition:

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

        Data = pd.DataFrame()
        ColNames = []
        for i in range(3):
            for V in ['U', 'F']:
                ColNames.append(V + str(i+1))
        for i in range(3):
            for V in ['R', 'M']:
                ColNames.append(V + str(i+1))

        for iName, Name in enumerate(ColNames):
            Data[Name] = Values[:,iName]
        
        Data.columns = ['X', 'FX', 'Y', 'FY', 'Z', 'FZ', 'Phi', 'MX', 'Theta', 'MY', 'Psi', 'MZ']

        return Data


#%% Main part
if __name__ == '__main__':

    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add long and short argument
    ScriptVersion = Parser.prog + ' version ' + Version
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('File', help='Dat file (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments.File)

