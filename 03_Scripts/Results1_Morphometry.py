#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Analysis of the morphometric measurements

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: April 2023
    """

#%% Imports
# Modules import

import os
import argparse
import pandas as pd
from pathlib import Path
from Utils import SetDirectories, Show


#%% Functions
# Define functions

def AFunction(Argument):

    return Something

#%% Main
# Main code

def Main():

    WD, DD, SD, RD = SetDirectories('FRACTIB')
    Data = pd.read_csv(str(RD / 'Morphometry.csv'))
    Check = pd.read_csv(str(RD / '01_Morphology' / '00_Data.csv'))
    Map = pd.read_csv(str(DD / 'SampleList.csv'))
    Data.loc[Map['MicroCT pretest file number'].argsort(),'BV/TV (-)']
    

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