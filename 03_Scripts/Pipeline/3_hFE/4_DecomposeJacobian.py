#%%
#!/usr/bin/env python3

# 00 Initialization
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from Utils import *

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Functions
# Define functions
def DecomposeJacobian(JacobianArray):

    # Determine 2D of 3D jacobian array
    JacobianTerms = JacobianArray.shape[-1]
    ArrayShape = JacobianArray.shape[:-1]

    if JacobianTerms == 4:

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(JacobianArray[ArrayShape)
        # VonMises_Strain = np.zeros(JacobianArray[ArrayShape)
        # MaxShear = np.zeros(JacobianArray[ArrayShape)

        for j in range(0, ArrayShape[0]):
            for i in range(0, ArrayShape[1]):
                F_d = np.matrix(JacobianArray[j, i, :].reshape((2,2)))

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

    elif JacobianTerms == 9:

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(JacobianArray[ArrayShape)
        # VonMises_Strain = np.zeros(JacobianArray[ArrayShape)
        # MaxShear = np.zeros(JacobianArray[ArrayShape)

        for k in range(0, ArrayShape[0]):
            for j in range(0, ArrayShape[1]):
                for i in range(0, ArrayShape[2]):

                    F_d = np.matrix(JacobianArray[k, j, i, :].reshape((3, 3)))

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

    return SphericalCompression, IsovolumicDeformation


#%% Load files
# 01 Set variables
FilePath = Path.cwd() / '..' / '..' / '..' / '04_Results' / '03_hFE' / '432_L_77_F'

# 02 Load files
ElementsPositions = pd.read_csv(str(FilePath / 'ElementsPositions.csv'),names=['X','Y','Z'])
DeformationGradients = pd.read_csv(str(FilePath / 'DeformationGradients.csv'),names=['F11','F12','F13','F21','F22','F23','F31','F32','F33'])

#%% Build arrays
# Build arrays
X = np.unique(ElementsPositions['X'].values)
Y = np.unique(ElementsPositions['Y'].values)
Z = np.unique(ElementsPositions['Z'].values)

F = np.zeros((len(Z),len(Y),len(X),9))
for Index in DeformationGradients.index:
    
    Position = ElementsPositions.loc[Index]
    X_Index = list(X).index(Position['X'])
    Y_Index = list(Y).index(Position['Y'])
    Z_Index = list(Z).index(Position['Z'])
    
    F[Z_Index,Y_Index,X_Index] = DeformationGradients.loc[Index].values

#%% Decompose deformation
# Symmetry in Y plane is necessary
SphericalCompression, IsovolumicDeformation = DecomposeJacobian(F)

#%% Write MHD
# Compute metadata
Spacing = np.array([X[1]-X[0],Y[1]-Y[0],Z[1]-Z[0]])
Origin = [X.min(), Y.min(), Z.min()]
# Origin = [0, 0, 0]

SC = sitk.GetImageFromArray(SphericalCompression)
SC.SetSpacing(Spacing)
SC.SetOrigin(Origin)
ID = sitk.GetImageFromArray(IsovolumicDeformation)
ID.SetSpacing(Spacing)
ID.SetOrigin(Origin)

Writer = Write()
Writer.MHD(SC, str(FilePath / 'J'), PixelType='float')
Writer.MHD(ID, str(FilePath / 'F_Tilde'), PixelType='float')

# %%

Figure, Axis = plt.subplots(1,1)
Axis.imshow(SphericalCompression[:,:,15])
plt.show()

# %%
