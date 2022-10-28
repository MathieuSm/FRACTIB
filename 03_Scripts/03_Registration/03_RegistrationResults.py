# 00 Initialization
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import SimpleITK as sitk

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

## Define functions
def ElastixRegistration(Dictionary):

    # Get dictionary parameters
    Transformations = Dictionary['Transformations']
    if type(Transformations) is not list:
        Transformations = list([Transformations])
    FixedImage = Dictionary['FixedImage']
    MovingImage = Dictionary['MovingImage']
    FixedMask = Dictionary['FixedMask']
    PyramidSchedule = Dictionary['PyramidSchedule']
    NIterations = Dictionary['NIterations']
    Alpha = Dictionary['Alpha']
    A = Dictionary['A']
    ResultsDirectory = Dictionary['ResultsDirectory']
    os.makedirs(ResultsDirectory,exist_ok=True)

    ## Set parameter map
    ParameterMapVector = sitk.VectorOfParameterMap()
    Dimension = len(FixedImage.shape)
    ImagePyramidSchedule = np.repeat(PyramidSchedule,Dimension)

    for Transformation in Transformations:

        ParameterMap = sitk.GetDefaultParameterMap(Transformation)
        ParameterMap['ResultImageFormat'] = ['mhd']
        ParameterMap['NewSamplesEveryIteration'] = ['true']
        ParameterMap['FixedImagePyramidSchedule'] = [str(ImagePyramidSchedule)[1:-1]]
        ParameterMap['MovingImagePyramidSchedule'] = [str(ImagePyramidSchedule)[1:-1]]
        ParameterMap['MaximumNumberOfIterations'] = [str(NIterations)]
        ParameterMap['SP_alpha'] = [str(Alpha)]
        ParameterMap['SP_A'] = [str(A)]

        ParameterMapVector.append(ParameterMap)


    ## Set Elastix and perform registration
    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetParameterMap(ParameterMapVector)
    ElastixImageFilter.SetFixedImage(sitk.GetImageFromArray(FixedImage))
    ElastixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage))
    ElastixImageFilter.SetFixedMask(sitk.GetImageFromArray(FixedMask))
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
    ElastixImageFilter.LogToConsoleOn()
    ElastixImageFilter.Execute()

    ## Get results
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    return sitk.GetArrayFromImage(ResultImage), TransformParameterMap
def TransformixTransformations(MovingImage,TransformParameterMap,ResultsDirectory=False):

    ## Compute jacobian of deformation field using transformix
    TransformixImageFilter = sitk.TransformixImageFilter()
    # TransformixImageFilter.ComputeDeformationFieldOn()
    TransformixImageFilter.ComputeSpatialJacobianOn()
    # TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOn()
    TransformixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage))
    TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
    TransformixImageFilter.SetOutputDirectory(ResultsDirectory)

    TransformixImageFilter.Execute()
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


# 01 Set variables
WorkingDirectory = os.getcwd()
DataDirectory = os.path.join(WorkingDirectory,'04_Results/04_Registration/')
MaskDirectory = os.path.join(WorkingDirectory,'02_Data/03_uCT/')

## Load file
Data = pd.read_csv(DataDirectory+'RegistrationResults.csv')

## Plot results
Figure, Axes = plt.subplots(1, 1, figsize=(7.5, 4.5),dpi=100)
Axes.plot(Data['Sample'],Data['DSC'],marker='o',fillstyle='none',linestyle='none',color=(1,0,0))
Axes.set_xlabel('Sample (-)')
Axes.set_ylabel('DSC (-)')
Axes.set_ylim([0,1])
plt.xticks(rotation=90)
plt.subplots_adjust(bottom=0.32)
plt.show()
plt.close(Figure)


Samples = ['433_R_77_F','441_R_64_M','452_L_75_F']