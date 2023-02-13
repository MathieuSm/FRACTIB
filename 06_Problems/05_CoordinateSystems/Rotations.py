#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script aimed to blabla

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: January 2023
    """

#%% Imports
# Modules import

import mpld3
import argparse
import numpy as np
import sympy as sp
import pandas as pd
from numba import njit
from pathlib import Path
import scipy.signal as sig
import matplotlib.pyplot as plt
from scipy.optimize import minimize

mpld3.enable_notebook()
np.set_printoptions(linewidth=500,suppress=True,formatter={'float_kind':'{:3}'.format})
# %matplotlib widget

#%% Functions
# Define functions

def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results

def ReadDAT(File):

    try:
        with open(File) as F:
            Text = F.read()

            Values = []
            Condition = Text.find('U1') + 1
            Start = Text.find('U1')
            Steps, Increments = [], []
            Step, Increment = 1, 1
            while Condition:

                Inc = int(Text[Start-349:Start-347])
                if Inc < Increment:
                    Step += 1
                Increment = Inc
                Increments.append(Increment)
                Steps.append(Step)

                for i in range(6):
                    iStart = Start + 125 + 15*i
                    iStop = iStart + 14
                    Values.append(float(Text[iStart : iStop]))

                    iStart += Text[Start:].find('RF1')
                    iStop += Text[Start:].find('RF1')
                    Values.append(float(Text[iStart : iStop]))

                Start = iStop + Text[iStop:].find('U1')
                Condition = Text[iStop:].find('U1') + 1

            Values = np.array(Values)
            Cols = 12
            Rows = Values.size // Cols
            Values = np.reshape(Values,(Rows,Cols))

            ColNames = []
            for i in range(3):
                for V in ['U', 'F']:
                    ColNames.append(V + str(i+1))
            for i in range(3):
                for V in ['R', 'M']:
                    ColNames.append(V + str(i+1))

            Data = pd.DataFrame()
            for iName, Name in enumerate(ColNames):
                Data[Name] = np.concatenate([[0], Values[:,iName]])
            
            Data.columns = ['X', 'FX', 'Y', 'FY', 'Z', 'FZ', 'Phi', 'MX', 'Theta', 'MY', 'Psi', 'MZ']

            Data['Step'] = np.concatenate([[0],Steps])
            Data['Increment'] = np.concatenate([[0],Increments])

        return Data

    except FileNotFoundError:
        print('File' + File + 'does not exist')

        return

def RotationMatrix(Phi=0.0, Theta=0.0, Psi=0.0, V=np.zeros(3), A=0):

    if (V != 0).any():
        a = np.cos(A) * np.eye(3)
        b = np.sin(A) * np.array([[0, -V[2], V[1]],[V[2], 0, -V[0]],[-V[1], V[0], 0]])
        c = (1-np.cos(A)) * np.outer(V, V)
        R = np.round(a + b + c, 15)

    else:

        # if list of angles, use numpy for speed
        try:
            len(Phi)
            Phi, Theta, Psi = np.array(Phi), np.array(Theta), np.array(Psi)
            Rx = np.array([[np.ones(len(Phi)),  np.zeros(len(Phi)), np.zeros(len(Phi))],
                           [np.zeros(len(Phi)),        np.cos(Phi),       -np.sin(Phi)],
                           [np.zeros(len(Phi)),        np.sin(Phi),        np.cos(Phi)]])

            Ry = np.array([[ np.cos(Theta),      np.zeros(len(Theta)),        np.sin(Theta)],
                           [np.zeros(len(Theta)), np.ones(len(Theta)), np.zeros(len(Theta))],
                           [-np.sin(Theta),      np.zeros(len(Theta)),        np.cos(Theta)]])

            Rz = np.array([[np.cos(Psi),              -np.sin(Psi), np.zeros(len(Psi))],
                           [np.sin(Psi),               np.cos(Psi), np.zeros(len(Psi))],
                           [np.zeros(len(Psi)), np.zeros(len(Psi)),  np.ones(len(Psi))]])

            R = np.einsum('ijl,jkl->lik',Rz, np.einsum('ijl,jkl->ikl',Ry, Rx))

        # if only float angles, use sympy for more accuracy
        except:
            Rx = sp.Matrix([[1,             0,              0],
                            [0, sp.cos(Phi), -sp.sin(Phi)],
                            [0, sp.sin(Phi),  sp.cos(Phi)]])

            Ry = sp.Matrix([[ sp.cos(Theta), 0, sp.sin(Theta)],
                            [0,             1,              0],
                            [-sp.sin(Theta), 0, sp.cos(Theta)]])

            Rz = sp.Matrix([[sp.cos(Psi), -sp.sin(Psi),     0],
                            [sp.sin(Psi),  sp.cos(Psi),     0],
                            [0,             0,              1]])

            R = Rz * Ry * Rx
    
    return np.array(R, dtype='float')

def Show3DSys(Sys, Vector=np.zeros(3)):

    Figure = plt.figure(figsize=(5.5, 4))
    Axis = Figure.add_subplot(111, projection='3d')
    Axis.quiver(-1,0,0,2,0,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(0,-1,0,0,2,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(0,0,-1,0,0,2,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(0,0,0,Sys[0,0], Sys[0,1], Sys[0,2],color=(1,0,0))
    Axis.quiver(0,0,0,Sys[1,0], Sys[1,1], Sys[1,2],color=(0,0,1))
    Axis.quiver(0,0,0,Sys[2,0], Sys[2,1], Sys[2,2],color=(0,1,0))

    if (Vector != 0).any():
        Axis.quiver(0, 0, 0, Vector[0], Vector[1], Vector[2], color=(0,0,0))

    # scaling hack
    Bbox_min = np.min([-1, -1, -1])
    Bbox_max = np.max([1, 1, 1])
    Axis.auto_scale_xyz([Bbox_min, Bbox_max], [Bbox_min, Bbox_max], [Bbox_min, Bbox_max])

    # make the panes transparent
    Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # make the grid lines transparent
    Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)

    # modify ticks
    MinX, MaxX = -1, 1
    MinY, MaxY = -1, 1
    MinZ, MaxZ = -1, 1
    Axis.set_xticks([MinX, 0, MaxX])
    Axis.set_yticks([MinY, 0, MaxY])
    Axis.set_zticks([MinZ, 0, MaxZ])
    Axis.xaxis.set_ticklabels([MinX, 0, MaxX])
    Axis.yaxis.set_ticklabels([MinY, 0, MaxY])
    Axis.zaxis.set_ticklabels([MinZ, 0, MaxZ])

    Axis.set_xlabel('X')
    Axis.set_ylabel('Y')
    Axis.set_zlabel('Z')

    plt.show()

    return

def Show3DPath(X, Y, Z, R):

    Figure = plt.figure(figsize=(5.5, 4))
    Axis = Figure.add_subplot(111, projection='3d')

    # Reference frame
    Axis.quiver(X[0],Y[0],Z[0],(X[-1]-X[0]),0,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[0],Y[0],Z[0],0,(Y[-1]-Y[0]),0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[0],Y[0],Z[0],0,0,(Z[-1]-Z[0]),color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[-1],Y[-1],Z[-1],(X[0]-X[-1]),0,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[-1],Y[-1],Z[-1],0,(Y[0]-Y[-1]),0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[-1],Y[-1],Z[-1],0,0,(Z[0]-Z[-1]),color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[0],Y[0],Z[-1],0,(Y[-1]-Y[0]),0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[0],Y[-1],Z[0],(X[-1]-X[0]),0,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[0],Y[0],Z[-1],(X[-1]-X[0]),0,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[-1],Y[0],Z[0],0,(Y[-1]-Y[0]),0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[-1],Y[0],Z[0],0,0,(Z[-1]-Z[0]),color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(X[0],Y[-1],Z[0],0,0,(Z[-1]-Z[0]),color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    
    # Path
    Axis.plot(X, Y, Z, color=(0,0,0))

    # scaling hack
    XSize = X.max() - X.min()
    YSize = Y.max() - Y.min()
    ZSize = Z.max() - Z.min()
    MaxSize = max([XSize, YSize, ZSize])
    X0 = X.min() - (MaxSize - XSize) / 2
    X1 = X.max() + (MaxSize - XSize) / 2
    Y0 = Y.min() - (MaxSize - YSize) / 2
    Y1 = Y.max() + (MaxSize - YSize) / 2
    Z0 = Z.min() - (MaxSize - ZSize) / 2
    Z1 = Z.max() + (MaxSize - ZSize) / 2
    Axis.auto_scale_xyz([X0, X1], [Y0, Y1], [Z0, Z1])

    # Frame orientation
    R *= MaxSize / 2
    Axis.quiver(X[-1],Y[-1],Z[-1],MaxSize/2, 0, 0,color=(1,0,0,0.25))
    Axis.quiver(X[-1],Y[-1],Z[-1],0, MaxSize/2, 0,color=(0,0,1,0.25))
    Axis.quiver(X[-1],Y[-1],Z[-1],0, 0, MaxSize/2,color=(0,1,0,0.25))
    Axis.quiver(X[-1],Y[-1],Z[-1],R[0,0], R[0,1], R[0,2],color=(1,0,0,1))
    Axis.quiver(X[-1],Y[-1],Z[-1],R[1,0], R[1,1], R[1,2],color=(0,0,1,1))
    Axis.quiver(X[-1],Y[-1],Z[-1],R[2,0], R[2,1], R[2,2],color=(0,1,0,1))

    # make the panes transparent
    Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # make the grid lines transparent
    Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)

    # modify ticks
    if MaxSize < 1:
        D = -np.floor(np.log10(MaxSize)).astype('int')
    elif MaxSize >= 1:
        D = np.ceil(np.log10(MaxSize)).astype('int')

    Axis.set_xticks([round(X0,D), round(X1,D)])
    Axis.xaxis.set_ticklabels([round(X0,D), round(X1,D)])  
    Axis.set_yticks([round(Y0,D), round(Y1,D)])
    Axis.yaxis.set_ticklabels([round(Y0,D), round(Y1,D)])  
    Axis.set_zticks([round(Z0,D), round(Z1,D)])
    Axis.zaxis.set_ticklabels([round(Z0,D), round(Z1,D)])  

    Axis.set_xlabel('X')
    Axis.set_ylabel('Y')
    Axis.set_zlabel('Z')

    plt.show()

    return

def Show3D(Paths):

    Colors = [(1,0,0), (0,0,1), (0,0,0), (0,1,0), (0,1,1), (1,0,1)]

    Figure = plt.figure(figsize=(5.5, 4))
    Axis = Figure.add_subplot(111, projection='3d')

    # Path
    XSize, YSize, ZSize = 0, 0, 0
    X0, X1, Y0, Y1, Z0, Z1 = np.zeros(6)
    for iP, Path in enumerate(Paths):
        X = Path[:,0]
        Y = Path[:,1]
        Z = Path[:,2]
        Axis.plot(X, Y, Z, color=Colors[iP])

        if abs(X.max() - X.min()) > XSize:
            XSize = abs(X.max() - X.min())
        if abs(Y.max() - Y.min()) > YSize:
            YSize = abs(Y.max() - Y.min())
        if abs(Z.max() - Z.min()) > ZSize:
            ZSize = abs(Z.max() - Z.min())

        # scaling hack
        MaxSize = max([XSize, YSize, ZSize])
        Xx = X.min() - (MaxSize - XSize) / 2
        if Xx < X0:
            X0 = Xx
        Xx = X.max() + (MaxSize - XSize) / 2
        if Xx > X1:
            X1 = Xx
        Yy = Y.min() - (MaxSize - YSize) / 2
        if Yy < Y0:
            Y0 = Yy
        Yy = Y.max() + (MaxSize - YSize) / 2
        if Yy > Y1:
            Y1 = Yy
        Zz = Z.min() - (MaxSize - ZSize) / 2
        if Zz < Z0:
            Z0 = Zz
        Zz = Z.max() + (MaxSize - ZSize) / 2
        if Zz > Z0:
            Z1 = Zz
        Axis.auto_scale_xyz([X0, X1], [Y0, Y1], [Z0, Z1])

    # make the panes transparent
    Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # make the grid lines transparent
    Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0) 

    Axis.set_xlabel('X')
    Axis.set_ylabel('Y')
    Axis.set_zlabel('Z')

    plt.show()

    return

def GetVectorAndAngle(R):
    """
    https://en.wikipedia.org/wiki/Rotation_matrix
    """

    V = np.array([R[2,1] - R[1,2], R[0,2] - R[2,0],R[1,0] - R[0,1]])

    Angle = np.arccos((np.trace(R)- 1)/2)

    return V / np.linalg.norm(V), Angle

def GetAngles(R):

    # Compute Euler angles from rotation matrix
    # Assuming R = RxRyRz
    # Adapted from https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix

    if len(R.shape) == 2:
        # Special case
        if R[0,2] == 1 or R[0,2] == -1:
            E3 = 0 # Set arbitrarily
            dlta = np.arctan2(R[0,1],R[0,2])

            if R[0,2] == -1:
                E2 = np.pi/2;
                E1 = E3 + dlta

            else:
                E2 = -np.pi/2;
                E1 = -E3 + dlta

        else:
            E2 = - np.arcsin(R[0,2])
            E1 = np.arctan2(R[1,2]/np.cos(E2), R[2,2]/np.cos(E2))
            E3 = np.arctan2(R[0,1]/np.cos(E2), R[0,0]/np.cos(E2))
    
    else:
        E1, E2, E3 = np.zeros(len(R)), np.zeros(len(R)), np.zeros(len(R))
        
        M1 = R[:,0,2] == 1
        M2 = R[:,0,2] == -1
        if sum(M1 + M2) > 0:
            dlta = np.arctan2(R[:,0,1],R[:,0,2])

            if sum(M2) > 0:
                E2[M2] = np.pi/2
                E1[M2] = E3 + dlta

            else:
                E2[M1] = -np.pi/2
                E1[M1] = -E3 + dlta

        else:
            E2 = - np.arcsin(R[:,0,2])
            E1 = np.arctan2(R[:,1,2]/np.cos(E2), R[:,2,2]/np.cos(E2))
            E3 = np.arctan2(R[:,0,1]/np.cos(E2), R[:,0,0]/np.cos(E2))

    return np.array([-E1, -E2, -E3]).T

def ComputeCost(Angle, Paths):
    R = RotationMatrix(Psi=Angle[0])
    RSys = np.dot(R, Paths[0].T).T
    Delta = np.abs(RSys - Paths[1])
    Cost = Delta.sum()
    return Cost


#%% Class
# Define classes

class Arguments:

    def __init__(self):
        self.Folder = 'FRACTIB'
        return
    
Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    # Set directories and read sample list
    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))
    ResultsDir = RD / '05_Comparison'

    # Compute general best rotation angle
    Angles = pd.DataFrame(columns=np.arange(0, 360, 10))
    # for Index in SampleList.index:
    Index = 0

    Sample = SampleList.loc[Index, 'Internal ID']
    FEADir = RD / '03_hFE' / Sample
    ExpDir = RD / '02_Experiment' / Sample

    # Read files
    FEAData = ReadDAT(str(FEADir / (Sample + '.dat')))
    ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))
    Paths = [FEAData[['X','Y','Z']].values,
             ExpData[['X','Y','Z']].values]
    Show3D(Paths)

    # Rotate system to align coordinates
    Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
    Start = Peaks[4]
    Stop = np.argmin(np.abs(ExpData['Z'] - FEAData['Z'].max()))
    RefData = ExpData['Z'][Start:Stop]
    FEStop = FEAData['Z'].idxmax()
    Interp = FEAData['Z'][:FEStop]

    InterpData = pd.DataFrame()
    for C in FEAData.columns[:-2]:
        Interpolated = np.interp(RefData, Interp, FEAData[C][:FEStop])
        InterpData[C] = Interpolated
    RefData = ExpData[Start:Stop].reset_index(drop=True)

    V = 'Z'
    Figure, Axis = plt.subplots(1,1)
    Axis.plot(InterpData[V],color=(0,0,0),label='hFE')
    Axis.plot(RefData[V],color=(1,0,0),label='Experiment')
    plt.legend()
    plt.show()

    Time = -1
    Path = RefData[['X', 'Y', 'Z']].values[:Time,:]
    Phi, Theta, Psi = RefData[['Phi', 'Theta', 'Psi']].values[Time] / 180 * np.pi
    Sys = RotationMatrix(Phi, Theta, Psi)
    Show3DPath(Path[:,0], Path[:,1], Path[:,2], Sys)

    Paths = [RefData[['X','Y','Z']].values,
             InterpData[['X','Y','Z']].values]
    Show3D(Paths)

    Optim = minimize(ComputeCost, [6*np.pi/4], args=(Paths), bounds=([0, 2*np.pi],))

    R = RotationMatrix(Psi=Optim.x[0])
    RSys = np.dot(R, Paths[0].T).T
    Show3D([RSys,Paths[1]])

    # Try to rotate angles
    PTP = RefData[['Phi', 'Theta', 'Psi']].values[:Time] / 180 * np.pi
    Phi, Theta, Psi = PTP[:,0], PTP[:,1], PTP[:,2]
    Rs = RotationMatrix(Phi, Theta, Psi)
    Rs = np.einsum('ij,ljk->lik',R,Rs)

    # Compute angles
    Angles = GetAngles(Rs)

    # Show rotated system
    Show3DSys(RotationMatrix(Angles[0,0], Angles[0,1], Angles[0,2]))


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
    Parser.add_argument('--Folder', help='Root folder of the project (required)', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)