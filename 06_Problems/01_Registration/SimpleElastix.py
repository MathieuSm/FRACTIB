#%%
# 00 Initialization
import numpy as np
import pandas as pd
from pathlib import Path
import SimpleITK as sitk
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%%
## Define functions
def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results
def GenerateCube(Dimension, CubeSize, Scale=10, Border=1):

    HalvedCubeSize = CubeSize / 2

    if Dimension == 2:

        ## Define cube vertices
        x1 = np.array([-HalvedCubeSize, -HalvedCubeSize])
        x2 = np.array([ HalvedCubeSize, -HalvedCubeSize])
        x3 = np.array([-HalvedCubeSize,  HalvedCubeSize])
        x4 = np.array([x3[0] + CubeSize, x2[1] + CubeSize])
        CubeVertices = np.array([x1, x2, x3, x4])

        ## Scale cube
        CubeScaled = CubeVertices * np.array([CubeScale])
        CubeValues = np.random.rand(CubeScale, CubeScale)
        for i in range(CubeScale):
            for j in range(CubeScale):
                Condition1 = i < Border-1 or j < Border-1
                Condition2 = i >= CubeScale - Border or j >= CubeScale - Border
                if Condition1 or Condition2:
                    CubeValues[j, i] = 1

    elif Dimension == 3:

        ## Define cube vertices
        x1 = np.array([-HalvedCubeSize, -HalvedCubeSize, -HalvedCubeSize])
        x2 = np.array([ HalvedCubeSize, -HalvedCubeSize, -HalvedCubeSize])
        x3 = np.array([-HalvedCubeSize,  HalvedCubeSize, -HalvedCubeSize])
        x4 = np.array([x3[0] + CubeSize, x2[1] + CubeSize, -HalvedCubeSize])
        x5 = np.array([x1[0],x1[1],HalvedCubeSize])
        x6 = np.array([x2[0],x2[1],HalvedCubeSize])
        x7 = np.array([x3[0],x3[1],HalvedCubeSize])
        x8 = np.array([x4[0],x4[1],HalvedCubeSize])
        CubeVertices = np.array([x1, x2, x3, x4, x5, x6, x7, x8])

        ## Scale cube
        CubeScaled = CubeVertices * np.array([CubeScale])
        CubeValues = np.random.rand(CubeScale, CubeScale, CubeScale)

        for i in range(CubeScale):
            for j in range(CubeScale):
                for k in range(CubeScale):
                    Condition1 = i < Border-1 or j < Border-1 or k < Border-1
                    Condition2 = i >= CubeScale - Border or j >= CubeScale - Border or k >= CubeScale - Border
                    if Condition1 or Condition2:
                        CubeValues[k, j, i] = 1

    return CubeVertices, CubeValues
def GenerateImage(Dimension, ImageSize, CubeVertices, CubeScale, CubeValues, b=np.matrix([[0],[0],[0]])):

    ImageCenter = np.array(ImageSize) / 2
    Image = np.zeros(ImageSize)

    if Dimension == 2:

        ## Compute cube center and lengths
        CubeCenter = (CubeVertices[3] + CubeVertices[0]) / 2
        XCubeLength = np.sqrt((CubeVertices[1,0] - CubeVertices[0,0])**2
                            + (CubeVertices[1,1] - CubeVertices[0,1])**2) * CubeScale
        YCubeLength = np.sqrt((CubeVertices[2,0] - CubeVertices[0,0])**2
                            + (CubeVertices[2,1] - CubeVertices[0,1])**2) * CubeScale

        ## Compute xy delta
        Delta11 = (CubeVertices[1, 0] - CubeVertices[0, 0]) * CubeScale / XCubeLength
        Delta12 = (CubeVertices[1, 1] - CubeVertices[0, 1]) * CubeScale / XCubeLength
        Delta21 = (CubeVertices[2, 0] - CubeVertices[0, 0]) * CubeScale / YCubeLength
        Delta22 = (CubeVertices[2, 1] - CubeVertices[0, 1]) * CubeScale / YCubeLength

        ## Compute initial coordinates (x1 vertex position)
        X0 = ImageCenter[1] + (CubeVertices[0,0]+b[0])*CubeScale
        Y0 = ImageCenter[0] + (CubeVertices[0,1]+b[1])*CubeScale

        ## Build fixed image
        for j in range(int(round(YCubeLength))+1):
            for i in range(int(round(XCubeLength))+1):

                ## Compute xy coordinates
                X_Coordinate = int(X0 + i * Delta11 + j * Delta21)
                Y_Coordinate = int(Y0 + j * Delta22 + i * Delta12)

                ## Compute corresponding xy values
                X_Ratio = round(i / XCubeLength * CubeValues.shape[1])
                Y_Ratio = round(j / YCubeLength * CubeValues.shape[0])

                Image[Y_Coordinate, X_Coordinate] = CubeValues[int(Y_Ratio)-1, int(X_Ratio)-1]

    elif Dimension == 3:

        ## Compute cube center and lengths
        CubeCenter = (CubeVertices[7] + CubeVertices[0]) / 2

        XCubeLength = np.sqrt((CubeVertices[1,0] - CubeVertices[0,0])**2
                            + (CubeVertices[1,1] - CubeVertices[0,1])**2
                            + (CubeVertices[1,2] - CubeVertices[0,2])**2) * CubeScale
        YCubeLength = np.sqrt((CubeVertices[2,0] - CubeVertices[0,0])**2
                            + (CubeVertices[2,1] - CubeVertices[0,1])**2
                            + (CubeVertices[2,2] - CubeVertices[0,2])**2) * CubeScale
        ZCubeLength = np.sqrt((CubeVertices[4,0] - CubeVertices[0,0])**2
                            + (CubeVertices[4,1] - CubeVertices[0,1])**2
                            + (CubeVertices[4,2] - CubeVertices[0,2])**2) * CubeScale

        ## Compute xy delta
        Delta11 = (CubeVertices[1, 0] - CubeVertices[0, 0]) * CubeScale / XCubeLength
        Delta12 = (CubeVertices[1, 1] - CubeVertices[0, 1]) * CubeScale / XCubeLength
        Delta13 = (CubeVertices[1, 2] - CubeVertices[0, 2]) * CubeScale / XCubeLength
        Delta21 = (CubeVertices[2, 0] - CubeVertices[0, 0]) * CubeScale / YCubeLength
        Delta22 = (CubeVertices[2, 1] - CubeVertices[0, 1]) * CubeScale / YCubeLength
        Delta23 = (CubeVertices[2, 2] - CubeVertices[0, 2]) * CubeScale / YCubeLength
        Delta31 = (CubeVertices[4, 0] - CubeVertices[0, 0]) * CubeScale / ZCubeLength
        Delta32 = (CubeVertices[4, 1] - CubeVertices[0, 1]) * CubeScale / ZCubeLength
        Delta33 = (CubeVertices[4, 2] - CubeVertices[0, 2]) * CubeScale / ZCubeLength

        ## Compute initial coordinates (x1 vertex position)
        X0 = ImageCenter[2] + (CubeVertices[0, 0] + b[0]) * CubeScale
        Y0 = ImageCenter[1] + (CubeVertices[0, 1] + b[1]) * CubeScale
        Z0 = ImageCenter[0] + (CubeVertices[0, 2] + b[2]) * CubeScale

        ## Build fixed image
        for k in range(int(round(ZCubeLength))+1):
            for j in range(int(round(YCubeLength))+1):
                for i in range(int(round(XCubeLength))+1):

                    ## Compute xy coordinates
                    X_Coordinate = int(X0 + i * Delta11 + j * Delta21 + k * Delta31)
                    Y_Coordinate = int(Y0 + i * Delta12 + j * Delta22 + k * Delta32)
                    Z_Coordinate = int(Z0 + i * Delta13 + j * Delta23 + k * Delta33)

                    ## Compute corresponding xy values
                    X_Ratio = round(i / XCubeLength * CubeValues.shape[2])
                    Y_Ratio = round(j / YCubeLength * CubeValues.shape[1])
                    Z_Ratio = round(k / ZCubeLength * CubeValues.shape[0])

                    Image[Z_Coordinate, Y_Coordinate, X_Coordinate] = CubeValues[int(Z_Ratio)-1, int(Y_Ratio)-1, int(X_Ratio)-1]

    return Image
def Translation(Dimension, u=0, v=0, w=0):
    if Dimension == 2:
        b = np.matrix([[u], [v]])

    elif Dimension == 3:
        b = np.matrix([[u], [v], [w]])

    return b
def Rotation(Dimension, Gamma=0, Beta=0, Alpha=0):

    if Dimension == 2:
        R = np.matrix([[np.cos(Gamma), -np.sin(Gamma)],
                       [np.sin(Gamma), np.cos(Gamma)]])

    elif Dimension == 3:
        Rz = np.matrix([[np.cos(Gamma), -np.sin(Gamma), 0],
                        [np.sin(Gamma), np.cos(Gamma), 0],
                        [0, 0, 1]])
        Ry = np.matrix([[np.cos(Beta), 0, np.sin(Beta)],
                        [0, 1, 0],
                        [-np.sin(Beta), 0, np.cos(Beta)]])
        Rx = np.matrix([[1, 0, 0],
                        [0, np.cos(Alpha), -np.sin(Alpha)],
                        [0, np.sin(Alpha), np.cos(Alpha)]])
        R = Rz * Ry * Rx

    return R
def Stretch(Dimension, LambdaH=0, LambdaX=0, LambdaY=0, LambdaZ=0):

    if Dimension == 2:
        e1 = np.array([1, 0])
        e2 = np.array([0, 1])
        I = np.eye(2)
        U = np.matrix(I + LambdaH * I
                      + LambdaX * np.outer(e1, e1)
                      + LambdaY * np.outer(e2, e2))

    elif Dimension == 3:
        e1 = np.array([1, 0, 0])
        e2 = np.array([0, 1, 0])
        e3 = np.array([0, 0, 1])
        I = np.eye(3)
        U = np.matrix(I + LambdaH * I
                      + LambdaX * np.outer(e1, e1)
                      + LambdaY * np.outer(e2, e2)
                      + LambdaZ * np.outer(e3, e3))

    return U
def PureShearing(Dimension, GammaXY=0, GammaYZ=0, GammaZX=0):
    if Dimension == 2:
        e1 = np.array([1, 0])
        e2 = np.array([0, 1])
        I = np.eye(2)
        G = np.matrix(I + GammaXY / 2 * np.outer(e1, e2)
                      + GammaXY / 2 * np.outer(e2, e1))

    elif Dimension == 3:
        e1 = np.array([1, 0, 0])
        e2 = np.array([0, 1, 0])
        e3 = np.array([0, 0, 1])
        I = np.eye(3)
        G = np.matrix(I + GammaXY / 2 * np.outer(e1, e2)
                      + GammaXY / 2 * np.outer(e2, e1)
                      + GammaYZ / 2 * np.outer(e2, e3)
                      + GammaYZ / 2 * np.outer(e3, e2)
                      + GammaZX / 2 * np.outer(e1, e3)
                      + GammaZX / 2 * np.outer(e3, e1))

    return G
def SimpleShearing(Dimension, GammaXY=0, GammaYZ=0, GammaZX=0):

    if Dimension == 2:
        e1 = np.array([1, 0])
        e2 = np.array([0, 1])
        I = np.eye(2)
        G = np.matrix(I + GammaXY * np.outer(e1, e2))

    elif Dimension == 3:
        e1 = np.array([1, 0, 0])
        e2 = np.array([0, 1, 0])
        e3 = np.array([0, 0, 1])
        I = np.eye(3)
        G = np.matrix(I + GammaXY * np.outer(e1, e2)
                      + GammaYZ * np.outer(e2, e3)
                      + GammaZX * np.outer(e1, e3))

    return G
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
def WriteMHD(ImageArray, Spacing, Path, FileName, PixelType='uint'):

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
    Offset = '0 0 0'
    CenterOfRotation = '0 0 0'
    AnatomicalOrientation = 'LPS'

    outs = open(os.path.join(Path, FileName) + '.mhd', 'w')
    outs.write('ObjectType = Image\n')
    outs.write('NDims = 3\n')
    outs.write('BinaryData = True\n')
    outs.write('BinaryDataByteOrderMSB = False\n')
    outs.write('CompressedData = False\n')
    outs.write('TransformMatrix = %s \n' % TransformMatrix)
    outs.write('Offset = %s \n' % Offset)
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
def WriteVTK(VectorField,FilePath,FileName,SubSampling=1):

    # Determine 2D of 3D vector field
    Dimension = VectorField.shape[-1]
    Size = VectorField.shape[:-1]

    if Dimension == 2:

        Size = np.append(1,Size)

        # Build Grid point
        z, y, x = np.arange(0, Size[0], SubSampling), np.arange(0, Size[1], SubSampling), np.arange(0, Size[2], SubSampling)
        NumberOfElements = len(x) * len(y) * len(z)

        # Build vector arrays
        u = VectorField[:, :, 0]
        v = VectorField[:, :, 1]
        w = np.zeros(NumberOfElements).reshape(Size[1:])

    elif Dimension == 3:

        # Build Grid point
        z, y, x = np.arange(0, Size[0], SubSampling), np.arange(0, Size[1], SubSampling), np.arange(0, Size[2], SubSampling)
        NumberOfElements = len(x) * len(y) * len(z)

        # Build vector arrays
        u = VectorField[:, :, :, 0]
        v = VectorField[:, :, :, 1]
        w = VectorField[:, :, :, 2]

    File = open(os.path.join(FilePath,FileName + '.vtk'),'w')

    # ASCII file header
    File.write('# vtk DataFile Version 4.2\n')
    File.write('VTK from Python\n')
    File.write('ASCII\n\n')
    File.write('DATASET RECTILINEAR_GRID\n')
    File.write('DIMENSIONS ' + str(len(x)) + ' ' + str(len(y)) + ' ' + str(len(z)) + '\n\n')
    File.close()

    # Append ascii x,y,z
    MaxLineWidth = Size[2]*Size[1]*Size[0]
    File = open(os.path.join(FilePath,FileName + '.vtk'),'a')
    File.write('X_COORDINATES ' + str(len(x)) + ' int\n')
    File.write(np.array2string(x.astype('int'),
                               max_line_width=MaxLineWidth,
                               threshold=MaxLineWidth)[1:-1] + '\n')
    File.write('\nY_COORDINATES ' + str(len(y)) + ' int\n')
    File.write(np.array2string(y.astype('int'),
                               max_line_width=MaxLineWidth,
                               threshold=MaxLineWidth)[1:-1] + '\n')
    File.write('\nZ_COORDINATES ' + str(len(z)) + ' int\n')
    File.write(np.array2string(z.astype('int'),
                               max_line_width=MaxLineWidth,
                               threshold=MaxLineWidth)[1:-1] + '\n\n')
    File.close()

    # Append another subheader
    File = open(os.path.join(FilePath,FileName + '.vtk'),'a')
    File.write('\nPOINT_DATA ' + str(NumberOfElements) + '\n\n')
    File.write('VECTORS Deformation float\n')
    File.close()

    # Append ascii u,v,w and build deformation magnitude array
    Magnitude = np.zeros(VectorField.shape[:-1])
    File = open(os.path.join(FilePath,FileName + '.vtk'),'a')

    if Dimension == 2:
        for j in range(0,Size[1],SubSampling):
            for i in range(0,Size[2],SubSampling):
                Magnitude[j,i] = np.sqrt(u[j,i]**2 + v[j,i]**2 + w[j,i]**2)
                File.write(str(u[j,i]) + ' ' + str(v[j,i]) + ' ' + str(w[j,i]) + ' ')

    elif Dimension == 3:
        for k in range(0,Size[0],SubSampling):
            for j in range(0,Size[1],SubSampling):
                for i in range(0,Size[2],SubSampling):
                    Magnitude[k,j,i] = np.sqrt(u[k,j,i] ** 2 + v[k,j,i] ** 2 + w[k,j,i] ** 2)
                    File.write(str(u[k,j,i]) + ' ' + str(v[k,j,i]) + ' ' + str(w[k,j,i]) + ' ')

    File.close()

    # Append another subheader
    File = open(os.path.join(FilePath, FileName + '.vtk'), 'a')
    File.write('\n\nSCALARS Magnitude float\n')
    File.write('LOOKUP_TABLE default\n')

    if Dimension == 2:
        for j in range(0, Size[1], SubSampling):
            for i in range(0, Size[2], SubSampling):
                File.write(str(Magnitude[j,i]) + ' ')

    elif Dimension == 3:
        for k in range(0, Size[0], SubSampling):
            for j in range(0, Size[1], SubSampling):
                for i in range(0, Size[2], SubSampling):
                    File.write(str(Magnitude[k,j,i]) + ' ')

    File.close()

    return
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


#%%
# 01 Set variables
WD, Data, Scripts, Results = SetDirectories('FRACTIB')
ResultsDirectory = str(WD / '06_Problems/01_Registration/3D_Tests/')


# 02 Build cube values with random values and border of fix values
Dimension = 3
CubeSize = 1
CubeScale = 50
Border = 5
CubeVertices, CubeValues = GenerateCube(Dimension, CubeSize, CubeScale, Border)

#%%
## Generate fixed image
ImageSize = (CubeScale*2,CubeScale*2,CubeScale*2)
FixedImage = GenerateImage(Dimension,ImageSize,CubeVertices,CubeScale,CubeValues)

WriteMHD(FixedImage,np.array([1,1,1]),ResultsDirectory,'FixedImage', PixelType='float')

#%%
# 03 Define transformation
b = Translation(Dimension)
R = Rotation(Dimension, Gamma=np.pi/12)
U = Stretch(Dimension)
Gp = PureShearing(Dimension)
Gs = SimpleShearing(Dimension)

## Gradient of the transformation
I = np.eye(Dimension)
F = R * U * (Gp + Gs - I)

#%%
## Compute deformed cube
TransformedCubeVertices = np.zeros(CubeVertices.shape)
for VertexNumber in range(CubeVertices.shape[0]):
    Vertex = np.matrix(CubeVertices[VertexNumber]).reshape((Dimension,1))
    TransformedCubeVertices[VertexNumber] = np.array(F * Vertex + b).reshape((1,Dimension))


#%%
## Generate moving image
MovingImage = GenerateImage(Dimension,ImageSize,TransformedCubeVertices,CubeScale,CubeValues,b)
WriteMHD(MovingImage,np.array([1,1,1]),ResultsDirectory,'MovingImage', PixelType='float')

## Plot images
Figure, Axes = plt.subplots(1, 2, figsize=(5.5, 4.5),dpi=100)
Axes[0].imshow(FixedImage[FixedImage.shape[2] // 2, :, :],cmap='binary')
Axes[1].imshow(MovingImage[MovingImage.shape[2] // 2, :, :],cmap='binary')
for i in range(2):
    Axes[i].set_xlim([0,FixedImage.shape[1]])
    Axes[i].set_ylim([0,FixedImage.shape[0]])
plt.show()
plt.close(Figure)

#%%
# 02 Set the parameter map vector
ParameterMap1 = sitk.GetDefaultParameterMap('affine')
ParameterMap1['ResultImageFormat'] = ['mhd']
ParameterMap1['SP_alpha'] = ['0.6']
ParameterMap1['SP_A'] = ['100']
ParameterMap1['MaximumNumberOfIterations'] = ['2000']
ParameterMap1['FixedImagePyramidSchedule'] = ['20 20 20 10 10 10']
ParameterMap1['MovingImagePyramidSchedule'] = ['20 20 20 10 10 10']
ParameterMap1['NewSamplesEveryIteration'] = ['true']

ParameterMap2 = sitk.GetDefaultParameterMap('bspline')
ParameterMap2['ResultImageFormat'] = ['mhd']
ParameterMap2['SP_alpha'] = ['0.6']
ParameterMap2['SP_A'] = ['100']
ParameterMap2['MaximumNumberOfIterations'] = ['2000']
ParameterMap2['FixedImagePyramidSchedule'] = ['20 20 20 10 10 10']
ParameterMap2['MovingImagePyramidSchedule'] = ['20 20 20 10 10 10']
ParameterMap2['NewSamplesEveryIteration'] = ['true']

ParameterMapVector = sitk.VectorOfParameterMap()
ParameterMapVector.append(ParameterMap1)
ParameterMapVector.append(ParameterMap2)

#%%
# 03 Set Elastix and register
ElastixImageFilter = sitk.ElastixImageFilter()
ElastixImageFilter.SetFixedImage(sitk.GetImageFromArray(FixedImage))
ElastixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage))
ElastixImageFilter.SetParameterMap(ParameterMapVector)
ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
ElastixImageFilter.LogToConsoleOff()
ElastixImageFilter.LogToFileOn()
ElastixImageFilter.Execute()

#%%
## Get results
ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
WriteMHD(sitk.GetArrayFromImage(ResultImage),np.array([1,1,1]),ResultsDirectory,'ResultImage1', PixelType='float')
TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

#%%
## Compute deformation field using transformix
TransformixImageFilter = sitk.TransformixImageFilter()
TransformixImageFilter.ComputeDeformationFieldOn()
TransformixImageFilter.ComputeSpatialJacobianOn()
TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOn()

TransformixImageFilter.SetOutputDirectory(ResultsDirectory)
TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
TransformixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage))

TransformixImageFilter.Execute()

DeformationField = TransformixImageFilter.GetDeformationField()
DeformationFieldNda = sitk.GetArrayFromImage(DeformationField)
WriteVTK(DeformationFieldNda, ResultsDirectory, 'DeformationField',SubSampling=10)


#%%
## Get transformation parameters
a11, a12, a21, a22, tx, ty = np.array(TransformParameterMap[0]['TransformParameters']).astype('float')
A = np.matrix([[a11,a12],[a21,a22]])

## Polar decomposition
from scipy.linalg import polar
R_A, U_A = polar(A)


#%%
## Plot
Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5),dpi=100)
Axes[0].imshow(FixedImage + MovingImage, cmap='jet')
Axes[1].imshow(sitk.GetArrayFromImage(ResultImage)+FixedImage, cmap='jet')
for i in range(2):
    Axes[i].set_xlim([0,FixedImage.shape[1]])
    Axes[i].set_ylim([0,FixedImage.shape[0]])
plt.show()
plt.close(Figure)

y,x,v = DeformationFieldNda.shape
x1,y1 = np.arange(x), np.arange(y)
x2 = np.tile(x1,len(y1)).reshape((len(y1),len(x1)))
y2 = np.repeat(y1,len(x1)).reshape((len(y1),len(x1)))

SubSample = 21
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(sitk.GetArrayFromImage(ResultImage), cmap='binary')
Axes.quiver(x2[::SubSample,::SubSample],y2[::SubSample,::SubSample],
            DeformationFieldNda[:,:,0][::SubSample,::SubSample],
            DeformationFieldNda[:,:,1][::SubSample,::SubSample], color = (1,0,0))
Axes.set_xlim([0,DeformationFieldNda.shape[1]])
Axes.set_ylim([0,DeformationFieldNda.shape[0]])
plt.show()
plt.close(Figure)


z,y,x,v = DeformationFieldNda.shape
x1,y1,z2 = np.arange(x), np.arange(y), np.arange(z)
x2 = np.tile(x1,len(y1)).reshape((len(y1),len(x1)))
y2 = np.repeat(y1,len(x1)).reshape((len(y1),len(x1)))
x3 = np.tile(x2,len(z2)).reshape((len(z2),len(y2),len(x2)))
y3 = np.zeros(x3.shape)
for i in range(len(z2)):
    y3[i, :, :] += y2
z3 = np.repeat(np.repeat(z2,len(x1)).reshape((len(z2),len(x1))),len(y1)).reshape((len(z2),len(y1),len(x1)))

from mpl_toolkits.mplot3d import axes3d
fig = plt.figure()
# ax = fig.gca(projection='3d')
ax = fig.gca()
ax.quiver(x3[:,::55,::55],y3[:,::55,::55],z3[:,::55,::55],
          DeformationFieldNda[:,:,:,2][:,::55,::55],
          DeformationFieldNda[:,:,:,1][:,::55,::55],
          DeformationFieldNda[:,:,:,0][:,::55,::55], length=0.1, color = 'black')
plt.show()


JacobianImage = sitk.ReadImage(ResultsDirectory+'fullSpatialJacobian.mhd')
JacobianArray = sitk.GetArrayFromImage(JacobianImage)

SubSampling = 5
SphericalCompression, IsovolumicDeformation = DecomposeJacobian(JacobianArray, SubSampling=SubSampling)

WriteMHD(SphericalCompression,np.array([1,1,1])*SubSampling,ResultsDirectory,'J', PixelType='float')
WriteMHD(IsovolumicDeformation,np.array([1,1,1])*SubSampling,ResultsDirectory,'F_Tilde', PixelType='float')

JacobianImage = sitk.ReadImage(ResultsDirectory+'spatialJacobian.mhd')
JacobianArray = sitk.GetArrayFromImage(JacobianImage)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(sitk.GetArrayFromImage(ResultImage), cmap='bone')
Image = Axes.imshow(JacobianArray,cmap='jet',alpha=0.5)
ColorBar = plt.colorbar(Image)
ColorBar.set_label('Spatial jacobian (-)')
Axes.set_xlim([0, JacobianArray.shape[1]])
Axes.set_ylim([0, JacobianArray.shape[0]])
plt.subplots_adjust(bottom=0.25,top=0.75)
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Image = Axes.imshow(SphericalCompression,cmap='jet',alpha=0.5)
ColorBar = plt.colorbar(Image)
ColorBar.set_label('Spherical compression (-)')
Axes.set_xlim([0, SphericalCompression.shape[1]])
Axes.set_ylim([0, SphericalCompression.shape[0]])
plt.subplots_adjust(bottom=0.25,top=0.75)
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(sitk.GetArrayFromImage(ResultImage), cmap='bone')
Image = Axes.imshow(IsovolumicDeformation,cmap='jet',alpha=0.5)
ColorBar = plt.colorbar(Image)
ColorBar.set_label('Isovolumic deformation (-)')
Axes.set_xlim([0, IsovolumicDeformation.shape[1]])
Axes.set_ylim([0, IsovolumicDeformation.shape[0]])
plt.subplots_adjust(bottom=0.25,top=0.75)
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(sitk.GetArrayFromImage(ResultImage), cmap='bone')
Image = Axes.imshow(HydrostaticStrain,cmap='jet',alpha=0.5)
ColorBar = plt.colorbar(Image)
ColorBar.set_label('Hydrostatic strain (-)')
plt.subplots_adjust(bottom=0.25,top=0.75)
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(sitk.GetArrayFromImage(ResultImage), cmap='bone')
Image = Axes.imshow(VonMises_Strain,cmap='jet',alpha=0.5)
ColorBar = plt.colorbar(Image)
ColorBar.set_label('Von Mises strain (-)')
plt.subplots_adjust(bottom=0.25,top=0.75)
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(sitk.GetArrayFromImage(ResultImage), cmap='bone')
Image = Axes.imshow(MaxShear,cmap='jet',alpha=0.5)
ColorBar = plt.colorbar(Image)
ColorBar.set_label('Max shear (-)')
plt.subplots_adjust(bottom=0.25,top=0.75)
plt.show()
plt.close(Figure)
