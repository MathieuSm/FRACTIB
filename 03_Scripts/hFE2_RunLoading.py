#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Run hFE simulation
    Adapted from Denis hFE pipeline

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import sh
import argparse
from Utils import *


#%% Main
# Main code

def Main(Arguments):

    Text = 'Run ' + Arguments.Sample
    Time.Process(1, Text)

    # Set working directory
    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)
    FEADir = RD / '03_hFE' / Arguments.Sample
    os.chdir(str(FEADir))

    # Build some variables
    if hasattr(Arguments, 'InputFile'):
        InputFile = Arguments.InputFile
    else:
        InputFile = 'Simulation.inp'

    # Absolutely necessary to start abaqus job
    try:
        os.environ.pop('PYTHONIOENCODING')
    except KeyError:
        pass

    # Write script
    Script = '''
#!/bin/bash
abaqus interactive job={Job} inp={InputFile} user={UMAT} cpus={nCPUs} ask_delete=OFF
'''
    
    Context = {'Job':'MaxLoad',
               'InputFile':InputFile,
               'UMAT':str(SD / Arguments.UMAT),
               'nCPUs':Arguments.nCPUs}

    FileName = 'Run.sh'
    with open(FileName, 'w') as File:
        File.write(Script.format(**Context))

    # Run simulation with script
    try:
        sh.bash(str(FEADir / FileName))
        Completed = True
    except:
        Completed = False

    if Completed:
        # Remove unnecessary files
        os.remove('MaxLoad.com')
        os.remove('MaxLoad.msg')
        os.remove('MaxLoad.odb')
        os.remove('MaxLoad.prt')
        os.remove('MaxLoad.sim')
        os.remove('MaxLoad.sta')

    # Print time
    Time.Process(0)

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
    Parser.add_argument('-I', '--Input', help='Input file for abaqus', type=str)
    Parser.add_argument('-U', '--UMAT', help='UMAT used for simulation', default='UMAT.f', type=str)
    Parser.add_argument('-N', '--CPUs', help='Number of CPUs to use', default=4, type=int)

    Parser.add_argument('Sample', help='Sample for which run simulation (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)