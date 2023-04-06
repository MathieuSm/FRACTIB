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

    return 

#%% Main
# Main code

def Main():

    WD, DD, SD, RD = SetDirectories('FRACTIB')
    Data = pd.read_csv(str(RD / 'Morphometry.csv'))
    Check = pd.read_csv(str(RD / '01_Morphology' / '00_Data.csv'))
    Map = pd.read_csv(str(DD / 'SampleList.csv'))
    Sort = Map['MicroCT pretest file number'].argsort().sort_values().index
    Structural = pd.read_csv(str(RD / 'Structural.csv'), header=[0,1])

    print(Check.mean())
    print(Data.mean())
    print(Data.std())

    OLS = Show.OLS(Data['BV/TV (-)'], Structural['Experiment']['Stiffness (N/mm)']/1E3)
    OLS = Show.OLS(Check.loc[Sort, 'BV/TV'], Structural['Experiment']['Stiffness (N/mm)']/1E3)

    OLS = Show.OLS(Data['BMC (mgHA)'], Structural['Experiment']['Stiffness (N/mm)']/1E3)
    OLS = Show.OLS(Data['BMC (mgHA)']/1E3, Structural['Experiment']['Max Force (N)']/1E3)

    

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