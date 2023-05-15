#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Compute and plot the yield surface
    Adapted from the script of Gabriela Gerber, Thanks!

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: May 2023
    """

#%% Imports
# Modules import

import mpld3
import argparse
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

mpld3.enable_notebook()
%matplotlib widget   # type:ignore

#%% Functions
# Define functions

def PlotSurface(X, Y, Z):

    Figure = plt.figure(figsize=(5.5, 4))
    Axis = Figure.add_subplot(111, projection='3d')
    Axis.plot_surface(X, Y, Z, cmap='jet')

    # Scaling hack
    Bbox_min = np.nanmin([X, Y, Z])
    Bbox_max = np.nanmax([X, Y, Z])
    Axis.auto_scale_xyz([Bbox_min, Bbox_max], [Bbox_min, Bbox_max], [Bbox_min, Bbox_max])

    # Make the panes transparent
    Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Make the grid lines transparent
    Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)

    # # Modify ticks
    # MinX, MaxX = -1, 1
    # MinY, MaxY = -1, 1
    # MinZ, MaxZ = -1, 1
    # Axis.set_xticks([MinX, 0, MaxX])
    # Axis.set_yticks([MinY, 0, MaxY])
    # Axis.set_zticks([MinZ, 0, MaxZ])
    # Axis.xaxis.set_ticklabels([MinX, 0, MaxX])
    # Axis.yaxis.set_ticklabels([MinY, 0, MaxY])
    # Axis.zaxis.set_ticklabels([MinZ, 0, MaxZ])

    Axis.set_xlabel('S11')
    Axis.set_ylabel('S22')
    Axis.set_zlabel('S33')

    plt.show()


#%% Main
# Main code

def Main():

    # Initialise Variables to characterise yield surface
    SIGD0P = 0.15
    SIGD0N = 0.2
    chido = 0.3
    taudo = 0.12
    S0=(SIGD0P+SIGD0N)/2/SIGD0P/SIGD0N
    pp=0
    qq=0
    rho = 1
    MM1 = 1
    MM2 = 1
    MM3 = 1

    # Variables for surface visualisation
    r = 0.4

    # Initialise Fabric Tensors
    F = np.zeros((6,6))
    F[0,0] = S0**2/(rho**pp*MM1**(2*qq))**2
    F[1,1] = S0**2/(rho**pp*MM2**(2*qq))**2
    F[2,2] = S0**2/(rho**pp*MM3**(2*qq))**2
    F[0,1] = -(chido*(MM1/MM2)**(2*qq))*S0**2/(rho**pp*MM1**(2*qq))**2
    F[0,2] = -(chido*(MM1/MM3)**(2*qq))*S0**2/(rho**pp*MM1**(2*qq))**2
    F[1,2] = -(chido*(MM2/MM3)**(2*qq))*S0**2/(rho**pp*MM2**(2*qq))**2
    F[1,0] = F[0,1]
    F[2,0] = F[0,2]
    F[2,1] = F[1,2]
    F[3,3] = 0.5/((taudo*rho**pp)**2)
    F[4,4] = 0.5/((taudo*rho**pp)**2)
    F[5,5] = 0.5/((taudo*rho**pp)**2)

    F2 = np.zeros((6,1))
    F2[0] = -(SIGD0P-SIGD0N)/2/SIGD0P/SIGD0N/(rho^pp*MM1**(2*qq))
    F2[1] = -(SIGD0P-SIGD0N)/2/SIGD0P/SIGD0N/(rho^pp*MM2**(2*qq))
    F2[2] = -(SIGD0P-SIGD0N)/2/SIGD0P/SIGD0N/(rho^pp*MM3**(2*qq))

    # Projection in the S1 S2 S3 space
    X = np.linspace(-r, r)
    Y = np.linspace(-r, r)
    X, Y = np.meshgrid(X,Y)
    P1 = (-2*F[0,2]*X - 2*F[1,2]*Y + 2*F2[0]*F2[2]*X + 2*F2[1]*F2[2]*Y - 2*F2[2])**2 - 4*(F2[2]**2 - F[2,2]) * (-F[0,0]*X**2 - F[1,1]*Y**2 - 2*F[0,2]*X*Y + F2[0]**2*X**2 + 2*F2[0]*F2[1]*X*Y - 2*F2[0]*X + F2[1]**2*Y**2 - 2*F2[1]*Y + 1)
    P2 = 2*F[0,2]*X + 2*F[1,2]*Y - 2*F2[0]*F2[2]*X - 2*F2[2]*F2[2]*Y + 2*F2[2]
    Z1 = (-np.sqrt(P1) + P2)/(2*(F2[2]**2 - F[2,2]))
    Z2 = (np.sqrt(P1) + P2)/(2*(F2[2]**2 - F[2,2]))

    Xc = np.concatenate([X, X])
    Yc = np.concatenate([Y, Y])
    Zc = np.concatenate([Z1, Z2])

    PlotSurface(Xc, Yc, Zc)

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