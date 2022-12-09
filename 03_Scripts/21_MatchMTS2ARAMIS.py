#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Match data recorded by MTS with records of the ARAMIS device
    using peaks alignment

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import json
import argparse

from Utils import *


#%% Functions
# Define functions

def AFunction(Argument):

    return Something


#%% Classes
# Define classes

class AClass():

    def __init__(self):
        self.Attribute = 'DefaultValue'

#%% Main
# Main code

class Arguments():

    def __init__(self):
        self.Folder = 'FRACTIB'
        self.Sample = '432_L_77_F'
Arguments = Arguments()

def Main(Arguments):

    # Set directories and load data
    WD, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    DataDir = Data / '03_Experiment'
    ResultsDir = Results / '02_Experiment' / Arguments.Sample
    os.makedirs(ResultsDir, exist_ok=True)

    print('\nLoad data')
    Tic = time.time()

    MTSFile = str(DataDir / '1_MTS' / 'AO_data_MTS.json')
    with open(MTSFile) as File:
        Data = json.load(File)
    MTSData = pd.DataFrame(Data[Arguments.Sample])
    MTSData.columns = ['T', 'D', 'F'] # rename time displacement and force

    ARAMISFile = str(DataDir / '2_ARAMIS' / (Arguments.Sample + '.csv'))
    ARAMISData = pd.read_csv(ARAMISFile, index_col=0, sep=';', header=1, parse_dates=['Time UTC'])

    Toc = time.time()
    PrintTime(Tic, Toc)

    # Preprocessing of ARAMIS data
    Z = 'TopPlate_Csys→AnatomicalCsys.LZ [mm]'
    ARAMISData['D'] = ARAMISData[Z][0] - ARAMISData[Z]
    ARAMISData['T'] = (ARAMISData['Time UTC'] - ARAMISData.loc[0,'Time UTC']).dt.total_seconds()


    # Signals filtering
    print('\nFilter signals')
    Tic = time.time()

    FilteredMTS = pd.DataFrame()
    
    Sampling = 1 / (MTSData.loc[1,'T'] - MTSData.loc[0,'T'])
    CutOff = 2.5 # Cut-off frequency in Hz

    FilteredMTS['D'] = Signal.Filter(MTSData['D'], Sampling, CutOff)
    FilteredMTS['F'] = Signal.Filter(MTSData['F'], Sampling, CutOff)

    FilteredARAMIS = pd.DataFrame()



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
    Parser.add_argument('Sample', help='Sample to process (required)', type=str)

    # Add defaults arguments
    Parser.add_argument('-F', '--Folder', help='Root folder name', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)