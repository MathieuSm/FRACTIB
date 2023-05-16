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

import os
import argparse
import numpy as np
import pandas as pd
from numba import njit
from skimage import io
import SimpleITK as sitk
import scipy.signal as sig
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

#%% Functions
# Define functions

@njit
def Decomposition(JacobianArray):

    # ProcessTiming(1, 'Decompose Jacobian of deformation')

    Terms = JacobianArray.shape[-1]
    ArrayShape = JacobianArray.shape[:-1]

    if Terms == 4:

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(ArrayShape)
        # VonMises_Strain = np.zeros(ArrayShape)
        # MaxShear = np.zeros(ArrayShape)

        # ProcessLength = ArrayShape[0] * ArrayShape[1]
        # Progress = 0
        for j in range(0, ArrayShape[0]):
            for i in range(0, ArrayShape[1]):
                F_d = JacobianArray[j, i, :].reshape((2,2))

                ## Unimodular decomposition of F
                J = np.linalg.det(F_d)
                SphericalCompression[j, i] = J

                if J > 0:
                    F_tilde = J ** (-1 / 3) * F_d
                    Norm_F_tilde = np.linalg.norm(F_tilde)
                else:
                    Norm_F_tilde = 0.0

                IsovolumicDeformation[j, i] = Norm_F_tilde

                # ## Optional: decomposition of F_tilde
                # R_tilde, U_tilde = polar(F_tilde)
                # Norm_U_tilde = np.sqrt(np.sum(U_tilde ** 2))

                ## Hydrostatic and deviatoric strain
                # I_d = np.matrix(np.eye(F_d.shape[0]))
                # E = 1/2 * (F_d.T * F_d - I_d)
                # Hydrostatic_E = -1/3 * np.trace(E) * I_d
                # Deviatoric_E = E - Hydrostatic_E
                #
                # HydrostaticStrain[k,j,i] = Hydrostatic_E[0,0]
                # MaxShear[k,j,i] = E.diagonal().max() - E.diagonal().min()
                #
                # VM_Strain = np.sqrt(3/2) * np.linalg.norm(Deviatoric_E)
                # VonMises_Strain[k,j,i] = VM_Strain

                # Progress += 1
                # ProgressNext(Progress/ProcessLength*20)

    elif Terms == 9:

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(ArrayShape)
        # VonMises_Strain = np.zeros(ArrayShape)
        # MaxShear = np.zeros(ArrayShape)

        # ProcessLength = ArrayShape[0] * ArrayShape[1] * ArrayShape[2]
        # Progress = 0
        for k in range(0, ArrayShape[0]):
            for j in range(0, ArrayShape[1]):
                for i in range(0, ArrayShape[2]):

                    F_d = JacobianArray[k, j, i, :].reshape((3, 3))

                    ## Unimodular decomposition of F
                    J = np.linalg.det(F_d)
                    SphericalCompression[k, j, i] = J

                    if J > 0:
                        F_tilde = J ** (-1 / 3) * F_d
                        Norm_F_tilde = np.linalg.norm(F_tilde)
                    else:
                        Norm_F_tilde = 0.0

                    IsovolumicDeformation[k, j, i] = Norm_F_tilde

                    # ## Optional: decomposition of F_tilde
                    # R_tilde, U_tilde = polar(F_tilde)
                    # Norm_U_tilde = np.sqrt(np.sum(U_tilde ** 2))

                    ## Hydrostatic and deviatoric strain
                    # I_d = np.matrix(np.eye(F_d.shape[0]))
                    # E = 1/2 * (F_d.T * F_d - I_d)
                    # Hydrostatic_E = -1/3 * np.trace(E) * I_d
                    # Deviatoric_E = E - Hydrostatic_E
                    #
                    # HydrostaticStrain[k,j,i] = Hydrostatic_E[0,0]
                    # MaxShear[k,j,i] = E.diagonal().max() - E.diagonal().min()
                    #
                    # VM_Strain = np.sqrt(3/2) * np.linalg.norm(Deviatoric_E)
                    # VonMises_Strain[k,j,i] = VM_Strain
                    
                    # Progress += 1
                    # ProgressNext(Progress/ProcessLength*20)

    return SphericalCompression, IsovolumicDeformation

def DecomposeJacobian(JacobianArray):

    SC, ID = Decomposition(JacobianArray)

    return SC, ID

def ShowSC(Image, IRange, Slice=None, Title=None, Axis='Z', FName=None):

    try:
        Array = sitk.GetArrayFromImage(Image)
        Dimension = Image.GetDimension()
    except:
        Array = Image
        Dimension = len(Array.shape)

    if Dimension == 3:
        
        if Axis == 'Z':
            if Slice:
                Array = Array[Slice,:,:]
            else:
                Array = Array[Array.shape[0]//2,:,:]
        if Axis == 'Y':
            if Slice:
                Array = Array[:,Slice,:]
            else:
                Array = Array[:,Array.shape[1]//2,:]
        if Axis == 'X':
            if Slice:
                Array = Array[:,:,Slice]
            else:
                Array = Array[:,:,Array.shape[2]//2]

    Figure, Axis = plt.subplots()
    CMap = Axis.imshow(Array,interpolation=None, cmap='jet', vmin=IRange[0], vmax=IRange[1])
    CBar = plt.colorbar(CMap, shrink=0.85, ticks=[IRange[0], 1, IRange[1]])
    Axis.axis('Off')
    
    if (Title):
        Axis.set_title(Title)

    if (FName):
        plt.savefig(FName, bbox_inches='tight', pad_inches=0)

    plt.close(Figure)

    return

def ShowID(Image, IRange, Slice=None, Title=None, Axis='Z', FName=None):

    N = 256
    Values = np.zeros((N, 4))
    Values[:, 0] = np.linspace(0, 1, N)
    Values[:, 2] = np.linspace(1, 0, N)
    Values[:, -1] = 1.0
    CMP = ListedColormap(Values)

    try:
        Array = sitk.GetArrayFromImage(Image)
        Dimension = Image.GetDimension()
    except:
        Array = Image
        Dimension = len(Array.shape)

    if Dimension == 3:
        
        if Axis == 'Z':
            if Slice:
                Array = Array[Slice,:,:]
            else:
                Array = Array[Array.shape[0]//2,:,:]
        if Axis == 'Y':
            if Slice:
                Array = Array[:,Slice,:]
            else:
                Array = Array[:,Array.shape[1]//2,:]
        if Axis == 'X':
            if Slice:
                Array = Array[:,:,Slice]
            else:
                Array = Array[:,:,Array.shape[2]//2]

    Figure, Axis = plt.subplots()
    CMap = Axis.imshow(Array,interpolation=None, cmap=CMP, vmin=IRange[0], vmax=IRange[1])
    CBar = plt.colorbar(CMap, shrink=0.85, ticks=[IRange[0], IRange[1]])
    Axis.axis('Off')
    
    if (Title):
        Axis.set_title(Title)

    if (FName):
        plt.savefig(FName, bbox_inches='tight', pad_inches=0)

    plt.close(Figure)

    return


#%% Main
# Main code

def Main():

    # Read experimental results
    Exp = pd.read_csv('../../../04_Results/02_Experiment/435_R_90_F/MatchedSignals.csv')

    # Truncate experiment to monotonic loading part and set to 0
    Peaks, Properties = sig.find_peaks(Exp['FZ'], prominence=1)
    MaxForce = Exp['FZ'].idxmin()
    MaxDisp = Exp['Z'].idxmax()
    DeltaTime = 10
    DeltaIndex = np.argmin(np.abs(Exp['T']-DeltaTime))
    Start = Peaks[Peaks < MaxForce - DeltaIndex][-1]
    Stop = Peaks[Peaks > MaxDisp][0]
    
    Exp = Exp[Start:Stop].reset_index(drop=True)
    Exp -= Exp.loc[0]

    # Structural results
    PS = ['PP', 'PPD', 'SS', 'SSD']
    NL = ['', 'NL']
    Colors = [(1,0,0), (1,0,1), (0,0,1), (0,1,0)]
    Style = ['-', '--']

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(Exp['Z'], -Exp['FZ']/1E3, color=(0,0,0), label='Experiment')
    for i, P in enumerate(PS):
        for j, L in enumerate(NL):
            if len(L) > 0:
                Label = P + ' ' + L
                Data = pd.read_csv(P + '_' + L + '_RefNode.csv')
            else:
                Label = P
                Data = pd.read_csv(P + '_RefNode.csv')
            Axis.plot(Data['Uz'], Data['Fz']/1E3, linestyle=Style[j], color=Colors[i], label=Label)
    Axis.set_xlabel('Displacement (mm)')
    Axis.set_ylabel('Force (kN)')
    plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.225))
    plt.show()

    # Compare deformation fields
    Fs = ['F11','F12','F13','F21','F22','F23','F31','F32','F33']
    ElementsDG = pd.read_csv('PPD_NL_DG.csv')

    Coords = ElementsDG.drop_duplicates(subset=['X', 'Y', 'Z'])
    X = np.unique(Coords['X'].values)
    Y = np.unique(Coords['Y'].values)
    Z = np.unique(Coords['Z'].values)

    # Build arrays
    Frames = ElementsDG.groupby(by=['Step','Frame'])
    SC, ID, SF = [], [], []
    i = 0
    for (Step, Frame), Data in Frames:
        F = np.zeros((len(Z),len(Y),len(X),9))
        F[:,:,:,0] = 1.0
        F[:,:,:,4] = 1.0
        F[:,:,:,8] = 1.0
        for Index in Data.index:
            
            Element = Data.loc[Index]
            X_Index = list(X).index(Element['X'])
            Y_Index = list(Y).index(Element['Y'])
            Z_Index = list(Z).index(Element['Z'])
            
            F[Z_Index,Y_Index,X_Index] = Element[Fs].values

        # Decompose deformation
        SphericalCompression, IsovolumicDeformation = DecomposeJacobian(F)
        SC.append(SphericalCompression)
        ID.append(IsovolumicDeformation)
        SF.append((Step, Frame))
    SC = np.array(SC)
    ID = np.array(ID)

    Mask = SC[-1] == 1.0
    for i, sf in enumerate(SF):
        SC[i,Mask] = np.nan
        if i == 0:
            SC[i,~Mask] = 1.0
        ShowSC(SC[i], [0, 2], Axis='X', FName='SC_%02d'%i)
        # ShowID(ID[i], [1.75,3.5], Axis='X', FName='ID_%02d'%i)

        Data = pd.read_csv('PPD_NL_RefNode.csv')
        F1 = Data['Step'] == sf[0]
        F2 = Data['Frame'] == sf[1]
        Index = Data[F1 & F2].index[0]

        Figure, Axis = plt.subplots(1,1, figsize=(4.55, 3.28))
        Axis.plot(Data.loc[:Index,'Uz'], Data.loc[:Index,'Fz']/1E3, color=(1,0,0))
        Axis.plot(Data.loc[Index,'Uz'], Data.loc[Index,'Fz']/1E3, color=(1,0,0), marker='o')
        Axis.set_xlabel('Displacement (mm)')
        Axis.set_ylabel('Force (kN)')
        Axis.set_xlim([-0.05, Data['Uz'].max()*1.05])
        Axis.set_ylim([Data['Fz'].min()*1.05/1E3, Data['Fz'].max()*1.05/1E3])
        plt.subplots_adjust(left=0.15, bottom=0.15)
        plt.savefig('RN_%02d'%i, dpi=100)
        plt.close(Figure)

        Image1 = io.imread('SC_%02d'%i + '.png')
        Image2 = io.imread('RN_%02d'%i + '.png')

        Shape = (Image1.shape[0], Image1.shape[1] + Image2.shape[1], Image1.shape[2])
        Image = np.zeros(Shape, 'uint8')
        Image[:,:Image1.shape[1]] = Image1
        Image[:,Image1.shape[1]:] = Image2
        
        io.imsave('IM_%02d'%i + '.png', Image)
        os.remove('SC_%02d'%i + '.png')
        os.remove('RN_%02d'%i + '.png')

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