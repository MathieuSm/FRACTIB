#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

def DecomposeJacobian(JacobianArray, SubSampling=1):

    # Determine 2D of 3D jacobian array
    JacobianTerms = JacobianArray.shape[-1]

    if JacobianTerms == 4:

        ArrayShape = JacobianArray[::SubSampling, ::SubSampling, 0].shape

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(JacobianArray[ArrayShape)
        # VonMises_Strain = np.zeros(JacobianArray[ArrayShape)
        # MaxShear = np.zeros(JacobianArray[ArrayShape)

        for j in range(0, ArrayShape[0]):
            for i in range(0, ArrayShape[1]):
                F_d = np.matrix(
                    JacobianArray[int(j * SubSampling), int(i * SubSampling), :].reshape((2,2)))

                ## Unimodular decomposition of F
                J = np.linalg.det(F_d)
                SphericalCompression[j, i] = J
                F_tilde = J ** (-1 / 3) * F_d
                Norm_F_tilde = np.linalg.norm(F_tilde)
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

        ArrayShape = JacobianArray[::SubSampling, ::SubSampling, ::SubSampling, 0].shape

        SphericalCompression = np.zeros(ArrayShape)
        IsovolumicDeformation = np.zeros(ArrayShape)
        # HydrostaticStrain = np.zeros(JacobianArray[ArrayShape)
        # VonMises_Strain = np.zeros(JacobianArray[ArrayShape)
        # MaxShear = np.zeros(JacobianArray[ArrayShape)

        for k in range(0, ArrayShape[0]):
            for j in range(0, ArrayShape[1]):
                for i in range(0, ArrayShape[2]):

                    F_d = np.matrix(JacobianArray[int(k*SubSampling), int(j*SubSampling), int(i*SubSampling), :].reshape((3, 3)))

                    ## Unimodular decomposition of F
                    J = np.linalg.det(F_d)
                    SphericalCompression[k, j, i] = J
                    F_tilde = J ** (-1 / 3) * F_d
                    Norm_F_tilde = np.linalg.norm(F_tilde)
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
def WriteRaw(ImageArray, OutputFileName, PixelType):

    if PixelType == 'uint':

        MinValue = ImageArray.min()

        if MinValue < 0:
            ShiftedImageArray = ImageArray + abs(MinValue)
            MaxValue = ShiftedImageArray.max()
        elif MinValue > 0:
            ShiftedImageArray = ImageArray - MinValue
            MaxValue = ShiftedImageArray.max()
        else :
            ShiftedImageArray = ImageArray
            MaxValue = ShiftedImageArray.max()

        ScaledImageArray = ShiftedImageArray / MaxValue * 255

        CastedImageArray = ScaledImageArray.astype(np.uint8)

    elif PixelType == 'short':
        CastedImageArray = ImageArray.astype(np.short)
    elif PixelType == 'float':
        CastedImageArray = ImageArray.astype('float32')

    File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape)
    File[:] = CastedImageArray[:]
    del File

    return
def WriteMHD(ImageArray, Spacing, Offset, Path, FileName, PixelType='uint'):

    if PixelType == 'short' or PixelType == 'float':
        if len(ImageArray.shape) == 2:

            Array_3D = np.zeros((1,ImageArray.shape[0],ImageArray.shape[1]))

            for j in range(ImageArray.shape[0]):
                for i in range(ImageArray.shape[1]):
                    Array_3D[0,j,i] = ImageArray[j,i]

            ImageArray = Array_3D

    nz, ny, nx = np.shape(ImageArray)

    lx = float(Spacing[0])
    ly = float(Spacing[1])
    lz = float(Spacing[2])

    TransformMatrix = '1 0 0 0 1 0 0 0 1'
    X_o, Y_o, Z_o = float(Offset[0]), float(Offset[1]), float(Offset[2])
    CenterOfRotation = '0 0 0'
    AnatomicalOrientation = 'LPS'

    outs = open(os.path.join(Path, FileName) + '.mhd', 'w')
    outs.write('ObjectType = Image\n')
    outs.write('NDims = 3\n')
    outs.write('BinaryData = True\n')
    outs.write('BinaryDataByteOrderMSB = False\n')
    outs.write('CompressedData = False\n')
    outs.write('TransformMatrix = %s \n' % TransformMatrix)
    outs.write('Offset = %g %g %g\n' % (X_o, Y_o, Z_o))
    outs.write('CenterOfRotation = %s \n' % CenterOfRotation)
    outs.write('AnatomicalOrientation = %s \n' % AnatomicalOrientation)
    outs.write('ElementSpacing = %g %g %g\n' % (lx, ly, lz))
    outs.write('DimSize = %i %i %i\n' % (nx, ny, nz))

    if PixelType == 'uint':
        outs.write('ElementType = %s\n' % 'MET_UCHAR')
    elif PixelType == 'short':
        outs.write('ElementType = %s\n' % 'MET_SHORT')
    elif PixelType == 'float':
        outs.write('ElementType = %s\n' % 'MET_FLOAT')

    outs.write('ElementDataFile = %s\n' % (FileName + '.raw'))
    outs.close()

    WriteRaw(ImageArray, os.path.join(Path, FileName) + '.raw', PixelType)

    return


# 01 Set variables
WorkingDirectory = os.getcwd()
DataPath = os.path.join(WorkingDirectory,'06_Problems/02_hFE/CubeTest/')
ResultsPath = os.path.join(WorkingDirectory,'06_Problems/02_hFE/CubeTest/')

# 02 Load files
ElementsPositions = pd.read_csv(DataPath + 'ElementsPositions.csv',names=['X','Y','Z'])
DeformationGradients = pd.read_csv(DataPath + 'DeformationGradients.csv',names=['F11','F12','F13','F21','F22','F23','F31','F32','F33'])

X = np.unique(ElementsPositions['X'].values)
Y = np.unique(ElementsPositions['Y'].values)
Z = np.unique(ElementsPositions['Z'].values)

F = np.zeros((len(X),len(Y),len(Z),9))
for Index in DeformationGradients.index:

    Position = ElementsPositions.loc[Index]
    X_Index = list(X).index(Position['X'])
    Y_Index = list(Y).index(Position['Y'])
    Z_Index = list(Z).index(Position['Z'])

    F[X_Index,Y_Index,Z_Index] = DeformationGradients.loc[Index].values

SphericalCompression, IsovolumicDeformation = DecomposeJacobian(F)

Spacing = np.array([X[1]-X[0],Y[1]-Y[0],Z[1]-Z[0]])
Offset = np.array([min(X),min(Y),min(Z)])
WriteMHD(SphericalCompression, Spacing, Offset, ResultsPath, 'J', PixelType='float')
WriteMHD(IsovolumicDeformation, Spacing, Offset, ResultsPath, 'F_Tilde', PixelType='float')
