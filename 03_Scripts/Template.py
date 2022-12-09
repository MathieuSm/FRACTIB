#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Short description of the analysis performed by this script

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: Month Year
    """

#%% Imports
# Modules import

import argparse


#%% Functions
# Define functions

def AFunction(Argument):

    return Something


#%% Classes
# Define classes

def AClass():

    def __init__(self):
        self.Attribute = 'DefaultValue'

#%% Main
# Main code

def Main(File):

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
    Parser.add_argument('File', help='File to process (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments.File)