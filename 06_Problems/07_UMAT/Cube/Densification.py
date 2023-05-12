#%% Imports
# Modules import

import os
import argparse
import numpy as np
import pandas as pd
from skimage import io
from scipy import signal
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D


#%% Functions
# Define functions

def PlotCube(Nodes, R, Values, FName):

    C = plt.get_cmap('jet')(round(R,3))

    Faces = [[0,1,3,2], [4,5,7,6],
            [0,3,4,7], [1,2,5,6],
            [0,1,4,5], [2,3,6,7]]
    
    Nodes = Nodes - Nodes.mean(axis=0)

    Figure = plt.figure()
    Axis = Figure.add_subplot(111, projection='3d')
    for Face in Faces:
        X = Nodes[Face,0].reshape((2,2))
        Y = Nodes[Face,1].reshape((2,2))
        Z = Nodes[Face,2].reshape((2,2))
        Axis.plot_surface(X, Y, Z, alpha=0.5, color=C[:-1], edgecolor=(0,0,0))
    Axis.scatter3D(Nodes[:,0], Nodes[:,1], Nodes[:,2], color=(0,0,0))
    Axis.set_xlabel('X')
    Axis.set_ylabel('Y')
    Axis.set_zlabel('Z')

    CBar = cm.ScalarMappable(Normalize(vmin=Values[0], vmax=Values[1]), cmap='jet')
    CBar = plt.colorbar(CBar, ax=Axis, location='top',
                        label='Stress (MPa)', ticks=Values,
                        fraction=0.05)

    # scaling hack
    Min, Max = -1.0, 1.0
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

    plt.savefig(FName)
    plt.close(Figure)

    return

#%% Main
# Main code

def Main():

    # Read data
    Nodes = pd.read_csv('Nodes.csv')
    Elements = pd.read_csv('Elements.csv')

    # Input file values
    MM1 = 1.0
    MM2 = 1.0
    MM3 = 1.0
    PBVT = 1.0

    # UMAT values
    SCA1 = 0.9419
    SCA2 = 0.7801

    # Elasticity parameters (FABRGWTB PMUBC, PHZ, 2013) new version 2017 trabecular bone
    E0  = 9759.0*(12000.0/9759.0)*SCA1
    V0  = 0.2278
    MU0 = 3117.0*(12000.0/9759.0)*SCA1
    KS  = 1.91
    LS  = 1.1

    # Strength parameters (FABRGWTB, PHZ, 2013) trabecular bone
    SIGD0N = 73.1*(12000.0/9759.0)*SCA2
    SIGD0P = 57.69*(12000.0/9759.0)*SCA2
    ZETA0  = 0.28
    TAUD0  = 29.61*(12000.0/9759.0)*SCA2
    PP     = 1.82
    QQ     = 0.98

    # Densification
    ECRIT = -0.3
    GAMMAL0 = 1100.0
    RL = 2.928
    GAMMAP0 = 1300.0
    RP = 2.77
    POWERDNS = 6.0

    # Compute densification
    TRACEEN = Elements['E1']+Elements['E2']+Elements['E3']
    CRITERION = TRACEEN-ECRIT
    SDNST = np.zeros((len(Elements), 3,3))
    XIDEN = np.eye(3)
    BVTV = 1.0

    GAMMAL = GAMMAL0*(BVTV**RL)
    GAMMAP = GAMMAP0*(BVTV**RP)
    for K1 in range(3):
        for K2 in range(3):
            SDNST[:,K1,K2] = (GAMMAL*CRITERION+GAMMAP* CRITERION**(POWERDNS-1.0))*XIDEN[K1,K2]

    # Elastic trial stress
    STR = np.zeros((6,len(Elements)))
    for i in range(6):
        STR[i] = Elements['S' + str(i+1)]

    # Quadric fourth order tensor FFFF
    S0 = (SIGD0P+SIGD0N)/2.0/SIGD0P/SIGD0N
    FFFFT = np.zeros((6,6))
    FFFFT[0,0] = S0**2/((BVTV**PP)*MM1**(2.0*QQ))**2
    FFFFT[1,1] = S0**2/((BVTV**PP)*MM2**(2.0*QQ))**2
    FFFFT[2,2] = S0**2/((BVTV**PP)*MM3**(2.0*QQ))**2
    FFFFT[0,1] = -(ZETA0*(MM1/MM2)**(2.0*QQ))*S0**2/((BVTV**PP)*MM1**(2.0*QQ))**2
    FFFFT[0,2] = -(ZETA0*(MM1/MM3)**(2.0*QQ))*S0**2/((BVTV**PP)*MM1**(2.0*QQ))**2
    FFFFT[1,2] = -(ZETA0*(MM2/MM3)**(2.0*QQ))*S0**2/((BVTV**PP)*MM2**(2.0*QQ))**2
    FFFFT[1,0] = FFFFT[0,1]
    FFFFT[2,0] = FFFFT[0,2]
    FFFFT[2,1] = FFFFT[1,2]
    FFFFT[3,3] = 0.5/((TAUD0*(BVTV**PP)*(MM1*MM2)**QQ)**2)
    FFFFT[4,4] = 0.5/((TAUD0*(BVTV**PP)*(MM1*MM3)**QQ)**2)
    FFFFT[5,5] = 0.5/((TAUD0*(BVTV**PP)*(MM2*MM3)**QQ)**2)
    FFFF = FFFFT*PBVT

    # Quadratic second order tensor FF
    FFT = np.zeros((3,3))
    FFT[0,0] = -(SIGD0P-SIGD0N)/2.0/SIGD0P/SIGD0N/((BVTV**PP)*MM1**(2.0*QQ))
    FFT[1,1] = -(SIGD0P-SIGD0N)/2.0/SIGD0P/SIGD0N/((BVTV**PP)*MM2**(2.0*QQ))
    FFT[2,2] = -(SIGD0P-SIGD0N)/2.0/SIGD0P/SIGD0N/((BVTV**PP)*MM3**(2.0*QQ))
    FFT  = FFT*PBVT

    FF = np.zeros(6)
    FF[0] = FFT[0,0]
    FF[1] = FFT[1,1]
    FF[2] = FFT[2,2]
    FF[3] = 0.0
    FF[4] = 0.0
    FF[5] = 0.0

    # Yield criterion
    RAD  = 1.0 # Perfect plasticity
    FFS  = np.matmul(FFFF,STR)
    SFFS = np.sum(STR * FFS, axis=0)
    YSTR = np.sqrt(SFFS) + np.sum(FF * STR.T, axis=1) - RAD

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(Elements['E3'], Elements['S3'], color=(1,0,0), linewidth=2, label='Abaqus Results')
    Axis.plot([ECRIT, ECRIT], [0, -1.5], color=(0,0,0), linewidth=1, linestyle='--', label='E$_c$')
    Axis.set_xlabel('Strain (-)') # Height of 1 mm -> direct conversion
    Axis.set_ylabel('Stress (MPa)') # Surface of 1 mm2 -> direct conversion
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2))
    plt.show()

    Frames = Nodes.groupby(by=['Step','Frame'])
    Stress, Strain = [], []
    i = 0
    for (Step, Frame), Data in Frames:

        for Node in [5,6,7,8]:
            if Node == 5:
                Force = Data[Data['Label']== Node]['Fz'].values[0]
                Disp = Data[Data['Label']== Node]['Uz'].values[0]
            else:
                Force += Data[Data['Label']== Node]['Fz'].values[0]
                Disp += Data[Data['Label']== Node]['Uz'].values[0]
        
        # Mean displacement of the 4 nodes
        Disp /= 4

        # Symmetric top and bottom displacement
        Disp *= 2

        # Compute and store stress and strain
        Stress.append(Force / 1.0)
        Strain.append(Disp / 1.0)

        if Frame < Elements['S3'].idxmin() * 1.5 and Step == 0.0:
            Mod = 3
        else:
            Mod = 10

        if np.mod(Frame, Mod) == 0.0:

            if Force >= 0:
                R = Force / SIGDAP
                R = R / 2 + 0.5
            else:
                R = Force / SIGDAN
                R = (R + 1) / 2
        
            P = Data[['X','Y','Z']][:8]
            D = Data[['Ux','Uy','Uz']][:8]
            N = P.values + D.values
            PlotCube(N, R, [-SIGDAN, SIGDAP], 'SS/Cube-' + '%03d'%i + '.png')

            Figure, Axis = plt.subplots(1,1)
            Axis.plot(Strain, Stress, color=(1,0,0), linewidth=2, label='Abaqus Results')
            Axis.plot(Strain[-1], Stress[-1], color=(1,0,0), linewidth=2, marker='o')
            Axis.plot([Elements['E3'].min(), Elements['E3'].max()],
                    np.array([SIGDAP, SIGDAP]) * RADK.values[-1],
                    linestyle='-.', linewidth=1, color=(0,0,0))
            Axis.plot([Elements['E3'].min(), Elements['E3'].max()],
                    [-SIGDAN, -SIGDAN],
                    linestyle='-.', linewidth=1, color=(0,0,0), label='UMAT $\sigma_0^-$ / $\sigma_0^+$')
            Axis.plot(-np.linspace(0,E0), -np.linspace(0,E0) * EAA,
                    color=(0,0,0), linestyle=':', linewidth=1, label='UMAT Modulus')
            Axis.plot([-E0, -E0], [0, -SIGDAN],
                    color=(0,0,0), linestyle='--', linewidth=1, label='Article $\epsilon_0^-$')
            Axis.plot(Elements.loc[:Elements['E3'].idxmin(),'E3'],
                    -YS.loc[:Elements['E3'].idxmin()],
                    color=(0,0,1), linewidth=1, linestyle='--', label='UMAT Softening')
            Axis.plot(np.linspace(0, dE) + Elements['E3'].min(),
                    np.linspace(0, dE) * Er - SIGDAN * np.array(RADK)[-1],
                    color=(0,0,0), linestyle='-', linewidth=1, label='UMAT dModulus')
            Axis.set_xlabel('Strain (-)') # Height of 1 mm -> direct conversion
            Axis.set_ylabel('Stress (MPa)') # Surface of 1 mm2 -> direct conversion
            plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.2))
            plt.subplots_adjust(top=0.85)
            plt.savefig('SS/Sig-' + '%03d'%i + '.png')
            plt.close(Figure)

            Image1 = io.imread('SS/Cube-' + '%03d'%i + '.png')
            Image2 = io.imread('SS/Sig-' + '%03d'%i + '.png')
            Shape = (Image1.shape[0], Image1.shape[1] + Image2.shape[1], Image1.shape[2])
            Image = np.zeros(Shape, 'uint8')
            Image[:,:Image1.shape[1]] = Image1
            Image[:,Image1.shape[1]:] = Image2
            
            io.imsave('SS/UMAT-' + '%03d'%i + '.png', Image)
            os.remove('SS/Cube-' + '%03d'%i + '.png')
            os.remove('SS/Sig-' + '%03d'%i + '.png')

            i += 1

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