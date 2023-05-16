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
import scipy.signal as sig
import matplotlib.pyplot as plt

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
    Columns = ['X','Y','Z','F11','F12','F13','F21','F22','F23','F31','F32','F33']
    ElementsDG = pd.read_csv('Elements_DG.csv', names=Columns)

    # Build arrays
    X = np.unique(ElementsDG['X'].values)
    Y = np.unique(ElementsDG['Y'].values)
    Z = np.unique(ElementsDG['Z'].values)

    F = np.zeros((len(Z),len(Y),len(X),9))
    for Index in ElementsDG.index:
        


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