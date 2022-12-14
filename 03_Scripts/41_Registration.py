#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used to perform pre/post-test registration in 3 steps:
        1. Masks rigid registration
        2. Masks affine registration
        3. Gray scale B-Spline registration

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

#%% Functions
# Define functions

def AFunction(Argument):

    return Something


#%% Classes
# Define classes

class Arguments(): # for testing purpose

    def __init__(self):
        self.Pre = 'C0001901'
        self.Post = 'C0001939'


#%% Main
# Main code

def Main(Pre, Post):

    print('\nRegister Pre / Post-test scans')
    TIC = time.time()

    # Set directories and load data
    WD, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    DataDir = Data / '02_uCT'
    ResultsDir = Results / '02_Experiment' / Arguments.Sample
    os.makedirs(ResultsDir, exist_ok=True)

    print('\nLoad data')
    Tic = time.time()

    SampleDirectory = str(Data / '02_uCT' / Sample) + '/'

    FixedImage = Read.AIM(SampleDirectory + Files[0])[0]
    MovingImage = Read.AIM(SampleDirectory + Files[1])[0]
    FixedCort = Read.AIM(SampleDirectory + Files[0][:-4] + '_CORT_MASK.AIM')[0]
    FixedTrab = Read.AIM(SampleDirectory + Files[0][:-4] + '_TRAB_MASK.AIM')[0]
    FixedMask = FixedCort + FixedTrab
    FixedMask.SetSpacing(FixedImage.GetSpacing())

    Toc = time.time()
    PrintTime(Tic, Toc)


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
    Parser.add_argument('Sample', help='Sample number (required)', type=str)

    # Add defaults arguments
    Parser.add_argument('-F', '--Folder', help='Root folder name', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments.Sample)