#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Read output of ReadOdb.py and plot results

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: May 2023
    """

#%% Imports
# Modules import

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#%% Main
# Main code

def Main():

    # Read data
    PP = pd.read_csv('PP.csv')
    SS = pd.read_csv('SS.csv')
    PP_NL = pd.read_csv('PP_NL.csv')
    SS_NL = pd.read_csv('SS_NL.csv')

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(PP['Uz'], PP['Fz']/1E3, color=(1,0,0), label='Perfect Plasticity')
    Axis.plot(SS['Uz'], SS['Fz']/1E3, color=(0,0,1), label='Simple Softening')
    Axis.plot(PP_NL['Uz'], PP_NL['Fz']/1E3, color=(1,0,0), linestyle='--', label='NL Perfect Plasticity')
    Axis.plot(SS_NL['Uz'], SS_NL['Fz']/1E3, color=(0,0,1), linestyle='--', label='NL Simple Softening')
    Axis.set_xlabel('Displacement (mm)')
    Axis.set_ylabel('Force (kN)')
    plt.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.2))
    plt.show()

        
    SS_0X = pd.read_csv('RefNode.csv')

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(SS['Uz'], SS['Fz']/1E3, color=(0,0,1), label='Simple Softening')
    Axis.plot(SS_NL['Uz'], SS_NL['Fz']/1E3, color=(0,0,1), linestyle='--', label='NL Simple Softening')
    Axis.plot(SS_0X['Uz'], SS_0X['Fz']/1E3, color=(1,0,0), label='NL Simple Softening 0X')
    Axis.set_xlabel('Displacement (mm)')
    Axis.set_ylabel('Force (kN)')
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.1))
    plt.show()

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