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
from numba import njit
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

        
    SS_X = pd.read_csv('RefNode.csv')
    Exp = pd.read_csv('../../../04_Results/02_Experiment/455_L_97_F/MatchedSignals.csv')

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(Exp['Z'], -Exp['FZ']/1E3, color=(0,0,1), label='Experiment')
    Axis.plot(SS_X['Uz'], SS_X['Fz']/1E3, color=(1,0,0), label='Simple Softening')
    Axis.plot(PP_NL['Uz'], PP_NL['Fz']/1E3, color=(1,0,0), linestyle='--', label='Perfect Plasticity')
    Axis.set_xlabel('Displacement (mm)')
    Axis.set_ylabel('Force (kN)')
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.1))
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
        
        Element = ElementsDG.loc[Index]
        X_Index = list(X).index(Element['X'])
        Y_Index = list(Y).index(Element['Y'])
        Z_Index = list(Z).index(Element['Z'])
        
        F[Z_Index,Y_Index,X_Index] = Element[Columns[3:]].values

    # Decompose deformation
    SphericalCompression, IsovolumicDeformation = DecomposeJacobian(F)

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