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

    Date: Mai 2023
    """

#%% Imports
# Modules import


import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#%% Functions
# Define functions

def PlotCube(Nodes):

    Faces = [[0,1,3,2], [4,5,7,6],
            [0,3,4,7], [1,2,5,6],
            [0,1,4,5], [2,3,6,7]]
    
    Figure = plt.figure()
    Axis = Figure.add_subplot(111, projection='3d')
    for Face in Faces:
        X = Nodes[Face,0].reshape((2,2))
        Y = Nodes[Face,1].reshape((2,2))
        Z = Nodes[Face,2].reshape((2,2))
        Axis.plot_surface(X, Y, Z, alpha=0.25, color=(0,0,0), edgecolor=(0,0,0))
    Axis.scatter3D(Nodes[:,0], Nodes[:,1], Nodes[:,2], color=(0,0,0))
    Axis.set_xlabel('X')
    Axis.set_ylabel('Y')
    Axis.set_zlabel('Z')

    # scaling hack
    Min, Max = -0.5, 0.5
    Axis.auto_scale_xyz([Min, Max], [Min, Max], [Min, Max])
    
    # make the panes transparent
    Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    # make the grid lines transparent
    Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    
    # modify ticks
    Axis.set_xticks([Min, 0, Max])
    Axis.set_yticks([Min, 0, Max])
    Axis.set_zticks([Min, 0, Max])
    Axis.xaxis.set_ticklabels([Min, 0, Max])
    Axis.yaxis.set_ticklabels([Min, 0, Max])
    Axis.zaxis.set_ticklabels([Min, 0, Max])

    plt.show()

    return

#%% Classes
# Define classes


#%% Main
# Main code

def Main():

    Nodes = pd.read_csv('Nodes.csv')
    Elements = pd.read_csv('Elements.csv')

    TopRef = Nodes[Nodes['Label']== 1]

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(TopRef['Ux'], color=(1,0,0))
    Axis.plot(TopRef['Uy'], color=(0,0,1))
    Axis.plot(TopRef['Uz'], color=(0,0,0))
    plt.show()

    Frames = Nodes.groupby(by=['Step','Frame'])

    for (Step, Frame), Data in Frames:

        if Frame == 502 or Frame == 0:
            P = Data[['X','Y','Z']][:8]
            D = Data[['Ux','Uy','Uz']][:8]
            N = P.values + D.values
            PlotCube(N)

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