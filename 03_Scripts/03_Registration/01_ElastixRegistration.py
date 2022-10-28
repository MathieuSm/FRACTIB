#%% #!/usr/bin/env python3
# 00 Initialization

import os
import time
import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Functions
# 01 Define functions
def PrintTime(Tic, Toc):

    """
    Print elapsed time in seconds to time in HH:MM:SS format
    :param Tic: Actual time at the beginning of the process
    :param Toc: Actual time at the end of the process
    """

    Delta = Toc - Tic

    Hours = np.floor(Delta / 60 / 60)
    Minutes = np.floor(Delta / 60) - 60 * Hours
    Seconds = Delta - 60 * Minutes - 60 * 60 * Hours

    print('Process executed in %02i:%02i:%02i (HH:MM:SS)' % (Hours, Minutes, Seconds))

    return
def ShowSlice(Image, Slice=None, Title=None, Axis='Z'):

    Array = sitk.GetArrayFromImage(Image)
    
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
    Axis.imshow(Array,interpolation=None, cmap='binary_r')
    Axis.axis('Off')
    
    if(Title):
        Axis.set_title(Title)

    plt.show(Figure)

    return
def GetSlice(Image, Slice, Axis='Z'):

    Direction = Image.GetDirection()
    Origin =  Image.GetOrigin()
    Spacing = Image.GetSpacing()
    Array = sitk.GetArrayFromImage(Image)

    if Axis == 'Z':
        Sliced = Array[Slice,:,:]
        ISliced = sitk.GetImageFromArray(Sliced)
        ISliced.SetDirection((Direction[0], Direction[1], Direction[3], Direction[4]))
        ISliced.SetOrigin((Origin[0], Origin[1]))
        ISliced.SetSpacing((Spacing[0], Spacing[1]))

    elif Axis == 'Y':
        Sliced = Array[:,Slice,:]
        ISliced = sitk.GetImageFromArray(Sliced)
        ISliced.SetDirection((Direction[0], Direction[2], Direction[6], Direction[8]))
        ISliced.SetOrigin((Origin[0], Origin[2]))
        ISliced.SetSpacing((Spacing[0], Spacing[2]))

    elif Axis == 'X':
        Sliced = Array[:,:,Slice]
        ISliced = sitk.GetImageFromArray(Sliced)
        ISliced.SetDirection((Direction[4], Direction[5], Direction[7], Direction[8]))
        ISliced.SetOrigin((Origin[1], Origin[2]))
        ISliced.SetSpacing((Spacing[1], Spacing[2]))

    return ISliced
def ElastixRotation(Dictionary):

    # Get dictionary parameters
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
    Dimension = FixedImage.GetDimension()
    ImagePyramidSchedule = np.repeat(PyramidSchedule,Dimension)

    ParameterMap = sitk.GetDefaultParameterMap('rigid')
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
    ElastixImageFilter.SetFixedImage(FixedImage)
    ElastixImageFilter.SetMovingImage(MovingImage)
    ElastixImageFilter.SetFixedMask(sitk.Cast(FixedMask, sitk.sitkUInt8))
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
    ElastixImageFilter.LogToConsoleOn()
    ElastixImageFilter.Execute()

    ## Get results
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    return ResultImage, TransformParameterMap
def ShowRegistration(FixedImage, MovingImage):

    # Perform Otsu's segmentation to compare Dice
    OtsuFilter = sitk.OtsuThresholdImageFilter()
    F_Otsu = OtsuFilter.Execute(FixedImage)
    F_Otsu = sitk.GetArrayFromImage(F_Otsu)
    M_Otsu = OtsuFilter.Execute(MovingImage)
    M_Otsu = sitk.GetArrayFromImage(M_Otsu)
    Dice = 2 * np.sum(F_Otsu * M_Otsu) / np.sum(F_Otsu + M_Otsu)

    # Plot both images
    Zeros = np.zeros((F_Otsu.shape[0], F_Otsu.shape[1], 3),'uint8')
    

    return Dice
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
def DecomposeJacobian(JacobianArray):

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

        ArrayShape = JacobianArray[:, :, :, 0].shape

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


#%% Variables
# 02 Set variables
WorkingDirectory = Path.cwd() / '../..'
DataDirectory = WorkingDirectory / '02_Data/02_uCT/'

#%% Registration parameters
# 03 Set registration parameters
Transformations = ['rigid','bspline']
PyramidSchedules = [[50,20,10]]
NIterations = [2000]
Alphas = [0.6]
As = [1000]

#%% Files loading
# 04 Load files
SampleList = pd.read_csv(str(DataDirectory /'BMD_Values.csv'))
LogFile = open(str(WorkingDirectory / '03_Scripts/03_Registration' / 'Registration.log'),'w')
LogFile.write('Registration Log File\n\n')

Fixed_Crop = [[[0,342],[25,500],[65,485]],
              [[0,326],[75,540],[75,500]],
              [[0,330],[10,500],[55,545]],
              [[0,330],[75,485],[35,525]],
              [[0,311],[140,605],[100,570]],
              [[0,325],[70,465],[35,505]],
              [[0,320],[25,549],[55,530]],
              [[0,325],[35,549],[40,525]],
              [[0,310],[100,630],[40,610]],
              [[0,325],[65,535],[0,549]],
              [[10,330],[50,475],[50,549]],
              [[0,337],[100,570],[150,600]],
              [[0,311],[125,585],[130,610]],
              [[0,311],[100,540],[105,580]],
              [[10,335],[95,515],[40,515]],
              [[0,330],[20,549],[0,530]],
              [[0,311],[20,630],[100,640]],
              [[0,330],[170,580],[80,570]],
              [[0,323],[15,510],[5,525]],
              [[0,325],[50,510],[80,549]],
              [[0,325],[60,460],[30,530]],
              [[0,325],[45,520],[0,549]],
              [[0,325],[105,500],[50,515]],
              [[0,325],[0,535],[85,520]],
              [[0,330],[95,500],[5,505]]]
Moving_Crop = [[[0,342],[70,575],[120,605]],
               [[0,326],[110,600],[120,585]],
               [[0,330],[75,560],[75,580]],
               [[0,330],[120,570],[20,520]],
               [[0,330],[110,585],[120,580]],
               [[0,320],[115,575],[60,545]],
               [[0,320],[90,625],[85,595]],
               [[0,325],[80,615],[110,590]],
               [[0,323],[95,660],[90,650]],
               [[10,325],[90,630],[70,620]],
               [[5,320],[115,545],[100,600]],
               [[0,337],[120,600],[140,670]],
               [[5,330],[100,590],[105,600]],
               [[0,335],[80,530],[120,600]],
               [[0,340],[160,580],[120,600]],
               [[15,325],[90,630],[60,635]],
               [[20,330],[30,655],[95,674]],
               [[10,330],[145,570],[80,580]],
               [[5,325],[115,630],[100,615]],
               [[5,325],[100,580],[135,600]],
               [[0,325],[110,555],[75,590]],
               [[0,325],[55,570],[30,630]],
               [[0,325],[170,580],[75,555]],
               [[0,325],[35,585],[105,595]],
               [[0,325],[140,575],[50,575]]]
Data = pd.DataFrame()
# Data = pd.read_csv(os.path.join(WorkingDirectory,'04_Results/03_Registration','RegistrationResults.csv'))

Indices = [1,9,20]
Indices = [20]

#%% Set index
# for Index in Indices:
Index = 0
SampleTime = time.time()
#%% uCT files loading
# 05 Load uCT files
print('\nLoad uCT files')
Tic = time.time()
Sample = SampleList.loc[Index,'Sample']
ResultsDirectory = os.path.join(WorkingDirectory, '04_Results/03_Registration', Sample)
os.makedirs(ResultsDirectory,exist_ok=True)
SampleData = {'Sample': Sample}
LogFile.write('Sample: ' + Sample + '\n')
SampleDirectory = os.path.join(DataDirectory,Sample+'/')
Files = [File for File in os.listdir(SampleDirectory) if File.endswith('DOWNSCALED.mhd')]
Files.sort()

FixedImage = sitk.ReadImage(SampleDirectory + Files[0])
MovingImage = sitk.ReadImage(SampleDirectory + Files[1])
FixedMask = sitk.ReadImage(SampleDirectory + Files[0][:-4] + '_FULLMASK.mhd')
Toc = time.time()
PrintTime(Tic, Toc)
LogFile.write('Files loaded in %.3f s'%(Toc-Tic) + '\n')

#%% Convert to array
# 06 Convert images to arrays
print('\nConvert images to np arrays')
Tic = time.time()
FixedImage = sitk.GetArrayFromImage(FixedImage)
MovingImage = sitk.GetArrayFromImage(MovingImage)
FixedMask = sitk.GetArrayFromImage(FixedMask).astype('uint8')
# Fixed_C = Fixed_Crop[Index]
# Moving_C = Moving_Crop[Index]
# FixedImage = FixedImage[Fixed_C[0][0]:Fixed_C[0][1],Fixed_C[1][0]:Fixed_C[1][1],Fixed_C[2][0]:Fixed_C[2][1]]
# MovingImage = MovingImage[Moving_C[0][0]:Moving_C[0][1],Moving_C[1][0]:Moving_C[1][1],Moving_C[2][0]:Moving_C[2][1]]
# FixedMask = FixedMask
Toc = time.time()
PrintTime(Tic,Toc)
LogFile.write('Image converted into arrays in %.3f s' % (Toc - Tic) + '\n')

#%% MHD writing
# 07 Write MHDs
print('\Write MHDs')
Tic = time.time()
WriteMHD(FixedImage, np.array([1, 1, 1]), ResultsDirectory, 'FixedImage', PixelType='float')
WriteMHD(MovingImage, np.array([1, 1, 1]), ResultsDirectory, 'MovingImage', PixelType='float')
Toc = time.time()
PrintTime(Tic,Toc)
LogFile.write('Write fixed and moving images in %.3f s' % (Toc - Tic) + '\n')

#%% Cog alignment
# 08 Align centers of gravity
print('\nAlign centers of gravity')
Tic = time.time()
CenterType = sitk.CenteredTransformInitializerFilter.MOMENTS
IniTransform = sitk.CenteredTransformInitializer(FixedImage, MovingImage, sitk.Euler3DTransform(), CenterType)
IniMove = sitk.Resample(MovingImage, FixedImage, IniTransform, sitk.sitkLinear, 0.0, MovingImage.GetPixelID())
Toc = time.time()
PrintTime(Tic, Toc)
LogFile.write('Align centers of gravity in %.3f s' % (Toc - Tic) + '\n')

#%% Initial rotation
print('\nPerform initial registration (2D only)')
Tic = time.time()
Slice = 60
Dictionary = {'FixedImage':GetSlice(FixedImage, Slice),
              'MovingImage':GetSlice(IniMove, Slice),
              'FixedMask': GetSlice(FixedMask, Slice),
              'PyramidSchedule':PyramidSchedules[-1],
              'NIterations':NIterations[-1],
              'Alpha': Alphas[0],
              'A': 1000,
              'ResultsDirectory':ResultsDirectory}
ResultImage, TransformParameterMap = ElastixRotation(Dictionary)
Toc = time.time()
PrintTime(Tic, Toc)



#%%
# 09 Perform rigid registration and write MHD
Tic = time.clock_gettime(0)
Dictionary = {'Transformations':Transformations[0],
                'FixedImage':FixedImage,
                'MovingImage':MovingImage,
                'FixedMask': FixedMask,
                'PyramidSchedule':PyramidSchedules[-1],
                'NIterations':NIterations[-1],
                'Alpha': Alphas[0],
                'A': 1000,
                'ResultsDirectory':ResultsDirectory}
ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)
Tac = time.clock_gettime(0)
LogFile.write('Perform rigid registration %i min %i s' % (np.floor((Tac-Tic)/60),np.mod(Tac-Tic,60)) + '\n')
SampleData['Rigid'] = round(Tac-Tic)

Tic = time.clock_gettime(0)
WriteMHD(ResultImage,np.array([1,1,1]),ResultsDirectory,'RigidResult', PixelType='float')
Tac = time.clock_gettime(0)
LogFile.write('Write rigid result image %.3f s' % (Tac - Tic) + '\n')

#%%
def Test():
    ## Build parameters dictionary for non-rigid registration
    N_k = 0
    for k in As:
        SampleData['A'] = k
        N_j = 0
        for j in NIterations:
            SampleData['N iterations'] = j
            N_i = 0
            for i in PyramidSchedules:
                SampleData['Pyramid Schedule'] = i

                Name = '_P' + str(PyramidSchedules[N_i][-1]) + '_I' + str(NIterations[N_j]) + '_A' + str(As[N_k])

                ## Perform non-rigid registration and write results
                Tic = time.clock_gettime(0)
                Dictionary = {'Transformations': Transformations[1:],
                              'FixedImage': FixedImage,
                              'MovingImage': MovingImage,
                              'FixedMask': FixedMask,
                              'PyramidSchedule': i,
                              'NIterations': j,
                              'Alpha': Alphas[0],
                              'A': k,
                              'ResultsDirectory':ResultsDirectory}
                ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)
                Tac = time.clock_gettime(0)
                LogFile.write('Perform non-rigid registration %i min %i s' % (np.floor((Tac - Tic) / 60), np.mod(Tac - Tic, 60)) + '\n')
                SampleData['Non-Rigid'] = round(Tac - Tic)

                Tic = time.clock_gettime(0)
                WriteMHD(ResultImage,np.array([1,1,1]),ResultsDirectory,'ResultImage', PixelType='float')
                Tac = time.clock_gettime(0)
                LogFile.write('Write rigid result image in %.3f s' % (Tac - Tic) + '\n')

                ## Evaluate Registration - Dice similarity coefficient (DSC)
                Tic = time.clock_gettime(0)
                SEG_R = np.zeros(ResultImage.shape)
                SEG_R += ResultImage
                SEG_R[SEG_R < 0] = 0

                SEG_F = np.zeros(FixedImage.shape)
                SEG_F += FixedImage
                SEG_F[SEG_F < 0] = 0

                Otsu_Filter = sitk.OtsuThresholdImageFilter()
                Otsu_Filter.SetInsideValue(0)
                Otsu_Filter.SetOutsideValue(1)
                Segmentation = Otsu_Filter.Execute(sitk.GetImageFromArray(SEG_R))
                R_Threshold = Otsu_Filter.GetThreshold()
                SEG_R[SEG_R < R_Threshold] = 0
                SEG_R[SEG_R > 0] = 1

                Otsu_Filter = sitk.OtsuThresholdImageFilter()
                Otsu_Filter.SetInsideValue(0)
                Otsu_Filter.SetOutsideValue(1)
                Segmentation = Otsu_Filter.Execute(sitk.GetImageFromArray(SEG_F))
                F_Threshold = Otsu_Filter.GetThreshold()
                SEG_F[SEG_F < F_Threshold] = 0
                SEG_F[SEG_F > 0] = 1

                DSC = 2*np.sum(SEG_F * SEG_R)/np.sum(SEG_F + SEG_R)
                SampleData['DSC'] = np.round(DSC,3)
                Tac = time.clock_gettime(0)
                LogFile.write('Compute DSC index in %.3f s' % (Tac - Tic) + '\n')

                ## Evaluate Registration - Plotting
                Tic = time.clock_gettime(0)
                AbsMax = 1e4
                F_Shape = np.array(FixedImage.shape)
                M_Shape = np.array(MovingImage.shape)
                Shape = np.min(np.array([F_Shape,M_Shape]),axis=0)

                F_Shaped = FixedImage[:Shape[0], :Shape[1], int(Shape[2]/2)]
                M_Shaped = MovingImage[:Shape[0], :Shape[1], int(Shape[2] / 2)]
                R_Shaped = ResultImage[:Shape[0], :Shape[1], int(Shape[2] / 2)]

                Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=100)
                Axes[0].imshow(F_Shaped - M_Shaped, cmap='jet', vmin=-AbsMax, vmax=AbsMax)
                ColorBar = Axes[1].imshow(F_Shaped - R_Shaped, cmap='jet', vmin=-AbsMax, vmax=AbsMax)
                for Axis in range(2):
                    Axes[Axis].set_xlim([0, F_Shaped.shape[1]])
                    Axes[Axis].set_ylim([0, F_Shaped.shape[0]])
                plt.subplots_adjust(bottom=0.3)
                ColorBarAxis = Figure.add_axes([0.25, 0.1, 0.5, 0.05])
                Figure.colorbar(ColorBar, cax=ColorBarAxis, orientation='horizontal')
                Axes[0].set_title('Initial')
                Axes[1].set_title('Registered')
                plt.title('DSC :' + str(round(DSC,3)))
                # plt.savefig(ResultsDirectory+'/R'+Name+'.png')
                plt.show()
                plt.close(Figure)
                Tac = time.clock_gettime(0)
                LogFile.write('Perform plotting in %.3f s' % (Tac - Tic) + '\n')

                ## Use transformix to compute spatial jacobian
                Tic = time.clock_gettime(0)
                TransformixTransformations(MovingImage, TransformParameterMap, ResultsDirectory)
                Tac = time.clock_gettime(0)
                LogFile.write('Compute transformation jacobian in %i min %i s' % (np.floor((Tac - Tic) / 60), np.mod(Tac - Tic, 60)) + '\n')

                Tic = time.clock_gettime(0)
                JacobianImage = sitk.ReadImage(ResultsDirectory + '/fullSpatialJacobian.mhd')
                JacobianImage.SetSpacing((0.098, 0.098, 0.098))

                ## Resample Jacobian image
                Offset = JacobianImage.GetOrigin()
                Direction = JacobianImage.GetDirection()
                Orig_Size = np.array(JacobianImage.GetSize(), dtype=np.int)
                Orig_Spacing = JacobianImage.GetSpacing()

                New_Spacing = (0.9712, 0.9712, 0.9712)

                Resample = sitk.ResampleImageFilter()
                Resample.SetInterpolator = sitk.sitkLinear
                Resample.SetOutputDirection(Direction)
                Resample.SetOutputOrigin(Offset)
                Resample.SetOutputSpacing(New_Spacing)

                New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(New_Spacing))
                New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
                New_Size = [int(s) for s in New_Size]
                Resample.SetSize(New_Size)

                ResampledJacobian = Resample.Execute(JacobianImage)

                ## Perform jacobian unimodular decomposition
                JacobianArray = sitk.GetArrayFromImage(ResampledJacobian)
                # SubSampling = i[-1]
                SphericalCompression, IsovolumicDeformation = DecomposeJacobian(JacobianArray)
                Tac = time.clock_gettime(0)
                LogFile.write('Read and decompose jacobian in %i min %i s' % (np.floor((Tac - Tic) / 60), np.mod(Tac - Tic, 60)) + '\n')

                ## Write results
                Tic = time.clock_gettime(0)
                WriteMHD(SphericalCompression, np.array([0.098, 0.098, 0.098]), ResultsDirectory, 'J', PixelType='float')
                WriteMHD(IsovolumicDeformation, np.array([0.098, 0.098, 0.098]), ResultsDirectory, 'F_Tilde', PixelType='float')
                Tac = time.clock_gettime(0)
                LogFile.write('Write decomposition results in %.3f s' % (Tac - Tic) + '\n\n')

                ## Add results to dataframe
                Data = Data.append(SampleData,ignore_index=True)
                Data.to_csv(os.path.join(WorkingDirectory,'04_Results/04_Registration','RegistrationResults2.csv'),index=False)

                ## Update name variable
                N_i += 1
            N_j += 1
        N_k += 1
#%%
LogFile.close()

# %%
