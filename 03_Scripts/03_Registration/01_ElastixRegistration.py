#%% #!/usr/bin/env python3
# 00 Initialization

import os
import time
import yaml
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
def GetSlice(Image, Slice=None, Axis='Z'):

    Direction = Image.GetDirection()
    Origin =  Image.GetOrigin()
    Spacing = Image.GetSpacing()
    Array = sitk.GetArrayFromImage(Image)

    if Axis == 'Z':

        if Slice:
            Sliced = Array[Slice,:,:]
        else:
            Sliced = Array[Array.shape[0] // 2, :, :]

        ISliced = sitk.GetImageFromArray(Sliced)
        ISliced.SetDirection((Direction[0], Direction[1], Direction[3], Direction[4]))
        ISliced.SetOrigin((Origin[0], Origin[1]))
        ISliced.SetSpacing((Spacing[0], Spacing[1]))

    elif Axis == 'Y':

        if Slice:
            Sliced = Array[:,Slice,:]
        else:
            Sliced = Array[:, Array.shape[1] // 2, :]

        ISliced = sitk.GetImageFromArray(Sliced)
        ISliced.SetDirection((Direction[0], Direction[2], Direction[6], Direction[8]))
        ISliced.SetOrigin((Origin[0], Origin[2]))
        ISliced.SetSpacing((Spacing[0], Spacing[2]))

    elif Axis == 'X':

        if Slice:
            Sliced = Array[:,:,Slice]
        else:
            Sliced = Array[:, :, Array.shape[2] // 2]

        ISliced = sitk.GetImageFromArray(Sliced)
        ISliced.SetDirection((Direction[4], Direction[5], Direction[7], Direction[8]))
        ISliced.SetOrigin((Origin[1], Origin[2]))
        ISliced.SetSpacing((Spacing[1], Spacing[2]))

    return ISliced
def ReadConfigFile(Filename):

    """ Read configuration file and store to dictionary """

    print('\n\nReading initialization file', Filename)
    with open(Filename, 'r') as File:
        Configuration = yaml.load(File, Loader=yaml.FullLoader)

    return Configuration
def Adjust_Image_Size(Image, CoarseFactor, CropZ='Crop'):

    """
    Adapted from Denis's utils_SA.py
    Images are adjusted according to CropType:
    0 = CropType.expand     (Expand image by copying layers)
    1 = CropType.crop       (Crop image)
    2 = CropType.variable   (Either crop or expand, depending on what includes less layers)
    """

    # Get array
    Array = sitk.GetArrayFromImage(Image)
    Array = Array.transpose(2, 1, 0)

    # Measure image shape
    IMDimX = np.shape(Array)[0]
    IMDimY = np.shape(Array)[1]
    IMDimZ = np.shape(Array)[2]

    AddDimX = CoarseFactor - (IMDimX % CoarseFactor)
    AddDimY = CoarseFactor - (IMDimY % CoarseFactor)

    # adjust in x and y direction
    Shape_Diff = [AddDimX, AddDimY]
    IMG_XY_Adjusted = np.lib.pad(Array,
                                 ((0, Shape_Diff[0]), (0, Shape_Diff[1]), (0, 0)),
                                 'constant', constant_values=(0),)

    if CropZ == 'Crop':
        Image_Adjusted = IMG_XY_Adjusted

    if CropZ == 'Expand':
        AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
        Shape_Diff = [AddDimX, AddDimY, AddDimZ]
        Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                    ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                    'edge')

    if CropZ == 'Variable':
        Limit = CoarseFactor / 2.0
        if IMDimZ % CoarseFactor > Limit:
            AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
            Shape_Diff = [AddDimX, AddDimY, AddDimZ]
            Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                        ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                        'edge')
        if IMDimZ % CoarseFactor < Limit:
            Image_Adjusted = IMG_XY_Adjusted

    Image_Adjusted = sitk.GetImageFromArray(Image_Adjusted.transpose(2, 1, 0))
    Image_Adjusted.SetSpacing(Image.GetSpacing())
    Image_Adjusted.SetOrigin(Image.GetOrigin())
    Image_Adjusted.SetDirection (Image.GetDirection())

    return Image_Adjusted
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
    Image1 = np.zeros((F_Otsu.shape[0], F_Otsu.shape[1], 4),'uint8')
    Image1[:,:,0] = (1-F_Otsu) * 255
    Image1[:,:,3] = (1-M_Otsu) * 255

    Image2 = np.zeros((M_Otsu.shape[0], M_Otsu.shape[1], 4),'uint8')
    Image2[:,:,2] = (1-M_Otsu) * 255
    Image2[:,:,3] = (1-M_Otsu) * 255

    Figure, Axis = plt.subplots(1,1)
    Axis.imshow(Image1)
    Axis.imshow(Image2, alpha=0.5)
    Axis.axis('Off')
    plt.subplots_adjust(0,0,1,1)
    plt.show()
    
    return Dice
def ElastixRegistration(Dictionary):

    # Get dictionary parameters
    Transformation = Dictionary['Transformation']
    FixedImage = Dictionary['FixedImage']
    MovingImage = Dictionary['MovingImage']
    FixedMask = Dictionary['FixedMask']
    PyramidSchedule = Dictionary['PyramidSchedule']
    NIterations = Dictionary['NIterations']
    Alpha = Dictionary['Alpha']
    A = Dictionary['A']
    ElementSize = Dictionary['ElementSize']
    ResultsDirectory = Dictionary['ResultsDirectory']
    os.makedirs(ResultsDirectory,exist_ok=True)

    ## Set parameter map
    Dimension = FixedImage.GetDimension()
    ImagePyramidSchedule = np.repeat(PyramidSchedule,Dimension)

    ParameterMapVector = sitk.VectorOfParameterMap()
    ParameterMap = sitk.GetDefaultParameterMap(Transformation)

    if Transformation == 'rigid':
        # ParameterMap['InitialTransformParametersFileName'] = 'InitialTransform.txt'
        ParameterMap['ResultImageFormat'] = ['mhd']
        ParameterMap['NewSamplesEveryIteration'] = ['true']
        ParameterMap['FixedImagePyramidSchedule'] = [str(ImagePyramidSchedule)[1:-1]]
        ParameterMap['MovingImagePyramidSchedule'] = [str(ImagePyramidSchedule)[1:-1]]
        ParameterMap['MaximumNumberOfIterations'] = [str(NIterations)]
        ParameterMap['SP_alpha'] = [str(Alpha)]
        ParameterMap['SP_A'] = [str(A)]

    if Transformation == 'bspline':
        ParameterMap = sitk.ParameterMap()
        # ParameterMap['NumberOfSpatialSamples'] = '4096'
        # ParameterMap['BSplineInterpolationOrder'] = 3

        ParameterMap['FixedInternalImagePixelType'] = ("float",)
        ParameterMap['MovingInternalImagePixelType'] = ("float",)
        ParameterMap['FixedImagePyramid'] = ("FixedRecursiveImagePyramid",)
        ParameterMap['MovingImagePyramid'] = ("MovingRecursiveImagePyramid",)
        ParameterMap['AutomaticParameterEstimation'] = ('true',)
        ParameterMap['CheckNumberOfSamples'] = ('true',)
        ParameterMap['DefaultPixelValue'] = (f'{0.0}',)
        ParameterMap['BSplineInterpolationOrder'] = ('3',)
        ParameterMap['FinalBSplineInterpolationOrder'] = ('3',)
        # ParameterMap['FinalGridSpacingInVoxels'] = ('24',)
        ParameterMap['FinalGridSpacingInPhysicalUnits'] = (str(ElementSize),str(ElementSize),str(ElementSize))
        # ParameterMap['FixedImagePyramid'] = ('FixedSmoothingImagePyramid',)
        # ParameterMap['GridSpacingSchedule'] = ('4.0', '4.0', '4.0', '2.0', '2.0', '2.0', '1.0', '1.0','1.0',)  # ('2.803221', '1.988100', '1.410000', '1.000000')

        ParameterMap['NumberOfHistogramBins'] = ('16', '32', '64')  # , '128')   ###

        ParameterMap['ImageSampler'] = ('RandomCoordinate',)
        ParameterMap['Interpolator'] = ('LinearInterpolator',)
        ParameterMap['MaximumNumberOfIterations'] = ('2000',)  # ('256',)
        ParameterMap['MaximumNumberOfSamplingAttempts'] = ('5',)  ###
        ParameterMap['Metric'] = ('AdvancedMattesMutualInformation')
        ParameterMap['Metric0Weight'] = ('1.0',)
        ParameterMap['Metric1Weight'] = ('1.0',)
        ParameterMap['MovingImagePyramid'] = ('MovingSmoothingImagePyramid',)
        ParameterMap['ImagePyramidSchedule'] = ('8', '8', '8', '4', '4', '4', '2', '2', '2', '1', '1', '1',)
        ParameterMap['NewSamplesEveryIteration'] = ('true',)
        ParameterMap['NumberOfResolutions'] = ('6',)
        ParameterMap['NumberOfSamplesForExactGradient'] = ('4096',)
        ParameterMap['NumberOfSpatialSamples'] = ('2048',)  ###
        ParameterMap['Optimizer'] = ('AdaptiveStochasticGradientDescent',)
        ParameterMap['Registration'] = ('MultiResolutionRegistration',)
        # ParameterMap['Registration'] = ("MultiResolutionRegistration",)
        ParameterMap['Interpolator'] = ("BSplineInterpolator",)
        ParameterMap['Metric'] = ('AdvancedMattesMutualInformation',)
        ParameterMap['ResampleInterpolator'] = ('FinalBSplineInterpolator',)
        ParameterMap['Resampler'] = ('DefaultResampler',)
        ParameterMap['HowToCombineTransforms'] = ('Compose',)

        ParameterMap['ResultImageFormat'] = ('mhd',)
        # ParameterMap['SP_A'] = ('100',)
        # ParameterMap['SP_alpha'] = ('0.6',)
        ParameterMap['Transform'] = ('BSplineTransform',)
        ParameterMap['WriteIterationInfo'] = ('false',)
        ParameterMap['WriteResultImage'] = ('true',)
    ParameterMapVector.append(ParameterMap)


    ## Set Elastix and perform registration
    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetParameterMap(ParameterMapVector)
    ElastixImageFilter.SetFixedImage(FixedImage)
    ElastixImageFilter.SetMovingImage(MovingImage)
    ElastixImageFilter.SetFixedMask(sitk.Cast(FixedMask, sitk.sitkUInt8))
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
    ElastixImageFilter.LogToConsoleOn()
    ElastixImageFilter.LogToFileOn()
    ElastixImageFilter.Execute()

    ## Get results
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    return ResultImage, TransformParameterMap
def TransformixTransformations(MovingImage,TransformParameterMap,ResultsDirectory=False):

    ## Compute jacobian of deformation field using transformix
    TransformixImageFilter = sitk.TransformixImageFilter()
    # TransformixImageFilter.ComputeDeformationFieldOn()
    TransformixImageFilter.ComputeSpatialJacobianOn()
    # TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOn()
    TransformixImageFilter.SetMovingImage(MovingImage)
    TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
    TransformixImageFilter.SetOutputDirectory(ResultsDirectory)

    TransformixImageFilter.Execute()
def WriteRaw(Image, OutputFileName, PixelType):

    ImageArray = sitk.GetArrayFromImage(Image)

    if PixelType == 'uint' or 'norm':

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

        ScaledImageArray = ShiftedImageArray / MaxValue

    if PixelType == 'uint':
        CastedImageArray = ScaledImageArray.astype('uint8')
    elif PixelType == 'norm':
        CastedImageArray = ScaledImageArray.astype('float32')
    elif PixelType == 'short':
        CastedImageArray = ImageArray.astype(np.short)
    elif PixelType == 'float':
        CastedImageArray = ImageArray.astype('float32')

    File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape, order='F')
    File[:] = CastedImageArray[:]
    del File

    return
def WriteMHD(Image, Path, FileName, PixelType='uint'):

    if PixelType == 'short' or PixelType == 'float':
        if Image.GetDimension() == 2:

            Array_3D = np.zeros((1,ImageArray.shape[0],ImageArray.shape[1]))

            for j in range(ImageArray.shape[0]):
                for i in range(ImageArray.shape[1]):
                    Array_3D[0,j,i] = ImageArray[j,i]

            ImageArray = Array_3D

    nz, ny, nx = Image.GetSize()

    lx, ly, lz = Image.GetSpacing()

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
    elif PixelType == 'float' or PixelType == 'norm':
        outs.write('ElementType = %s\n' % 'MET_FLOAT')

    outs.write('ElementDataFile = %s\n' % (FileName + '.raw'))
    outs.close()

    WriteRaw(Image, os.path.join(Path, FileName) + '.raw', PixelType)

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
def DecomposeJacobian(JacobianImage):

    # Determine 2D of 3D jacobian array
    JacobianArray = sitk.GetArrayFromImage(JacobianImage)
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

    SphericalCompression = sitk.GetImageFromArray(SphericalCompression)
    IsovolumicDeformation = sitk.GetImageFromArray(IsovolumicDeformation)

    for Image in [SphericalCompression, IsovolumicDeformation]:
        Image.SetSpacing(JacobianImage.GetSpacing())
        Image.SetDirection(JacobianImage.GetDirection())
        Image.SetOrigin(JacobianImage.GetOrigin())

    return SphericalCompression, IsovolumicDeformation


#%% Variables
# 02 Set variables
WorkingDirectory = Path.cwd() / '../..'
DataDirectory = WorkingDirectory / '02_Data/02_uCT/'

#%% Files loading
# 04 Load files
SampleList = pd.read_csv(str(DataDirectory /'BMD_Values.csv'))
LogFile = open(str(WorkingDirectory / '03_Scripts/03_Registration' / 'Registration.log'),'w+')
LogFile.write('Registration Log File\n\n')

Data = pd.DataFrame()

#%% Set index
for Index in range(1):
# Index = 8
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
    FixedMask.SetSpacing(FixedImage.GetSpacing())
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Files loaded in %.3f s'%(Toc-Tic) + '\n')
    
    #%% Adapt image size to hFE meshing
    ConfigFile = Path(WorkingDirectory,'03_Scripts/04_hFE/ConfigFile.yaml')

    # Read config and store to dictionary
    Config = ReadConfigFile(ConfigFile)

    # coarsening factor = FE element size / CT voxel size
    Spacing = FixedImage.GetSpacing()
    CoarseFactor = int(round(Config['ElementSize'] / Spacing[0]))
    FixedImage = Adjust_Image_Size(FixedImage, CoarseFactor)



    #%% Cog alignment
    # 08 Align centers of gravity
    print('\nAlign centers of gravity')
    Tic = time.time()
    CenterType = sitk.CenteredTransformInitializerFilter.MOMENTS
    Otsu = sitk.OtsuThresholdImageFilter()
    Otsu.SetInsideValue(1)
    Otsu.SetOutsideValue(0)
    # F_Otsu = Otsu.Execute(FixedImage)
    # M_Otsu = Otsu.Execute(MovingImage)
    IniTransform = sitk.CenteredTransformInitializer(FixedImage, MovingImage, sitk.Euler3DTransform(), CenterType)
    IniMove = sitk.Resample(MovingImage, FixedImage, IniTransform, sitk.sitkLinear, 0.0, MovingImage.GetPixelID())
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Align centers of gravity in %.3f s' % (Toc - Tic) + '\n')

    #%% Initial rotation
    # Perform initial rotation
    ###!! If needed compute ellipse and align directions using cog !!
    print('\nPerform initial registration (2D only)')
    Tic = time.time()
    Slice = 60
    Dictionary = {'FixedImage':GetSlice(FixedImage, Slice),
                'MovingImage':GetSlice(IniMove, Slice),
                'FixedMask': GetSlice(FixedMask, Slice),
                'PyramidSchedule': [50, 20, 10],
                'NIterations': 2000,
                'Alpha': 0.6,
                'A': 1000,
                'ResultsDirectory':ResultsDirectory}
    ResultImage, TransformParameterMap = ElastixRotation(Dictionary)
    Dice = ShowRegistration(GetSlice(FixedImage,Slice),ResultImage)
    print('Dice coefficient: %.3f' % (Dice))
    Toc = time.time()
    PrintTime(Tic, Toc)

    #%% Build 3D initial transform
    # Build initial transform
    TransformParameters = TransformParameterMap[0]['TransformParameters']
    TransformParameters = ('0',
                        '0',
                        str(TransformParameters[0]),
                        str(TransformParameters[1]),
                        str(TransformParameters[2]),
                        '0')
    CR = TransformParameterMap[0]['CenterOfRotationPoint']

    InitialTransform = TransformParameterMap[0]
    InitialTransform['CenterOfRotationPoint'] = (CR[0], CR[1], '0.0')
    InitialTransform['Direction'] = [str(d) for d in FixedImage.GetDirection()]
    InitialTransform['FixedImageDimension'] = str(FixedImage.GetDimension())
    InitialTransform['Index'] = ('0', '0', '0')
    InitialTransform['MovingImageDimension'] = str(MovingImage.GetDimension())
    InitialTransform['NumberOfParameters'] = '6'
    InitialTransform['Origin'] = [str(o) for o in FixedImage.GetOrigin()]
    InitialTransform['Size'] = [str(s) for s in FixedImage.GetSize()]
    InitialTransform['Spacing'] = [str(s) for s in FixedImage.GetSpacing()]
    InitialTransform['TransformParameters'] = TransformParameters

    TransformFile = open(ResultsDirectory + '/TransformParameters.0.txt')
    Text = TransformFile.read()
    TransformFile.close()
    for K in InitialTransform.keys():
        Start = Text.find(K)
        Stop = Text[Start:].find(')')
        OldText = Text[Start:Start+Stop]
        NewText = Text[Start:Start+len(K)]
        for i in InitialTransform[K]:
            NewText += ' ' + str(i)
        Text = Text.replace(OldText,NewText)
    TransformFile = open(ResultsDirectory + '/InitialTransform.txt', 'w')
    TransformFile.write(Text)
    TransformFile.close()


    #%% Rigid registration
    # 09 Perform rigid registration and write MHD
    print('\nPerform 3D rigid registration')
    Tic = time.time()
    Dictionary = {'Transformation':'rigid',
                'FixedImage':FixedImage,
                'MovingImage':IniMove,
                'FixedMask':FixedMask,
                'PyramidSchedule':[50, 20, 10],
                'NIterations':2000,
                'Alpha': 0.6,
                'A': 1000,
                'ElementSize':Config['ElementSize'],
                'ResultsDirectory':ResultsDirectory}
    ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)
    Dice = ShowRegistration(GetSlice(FixedImage, Slice), GetSlice(ResultImage, Slice))
    print('Dice coefficient: %.3f' % (Dice))
    Dice = ShowRegistration(GetSlice(FixedImage, Axis='Y'), GetSlice(ResultImage, Axis='Y'))
    print('Dice coefficient: %.3f' % (Dice))
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Perform rigid registration %i min %i s' % (np.floor((Toc-Tic)/60),np.mod(Toc-Tic,60)) + '\n')

    #%% Non-rigid registration
    # 10 Perform non-rigid registration
    print('\nPerform non-rigid registration')
    Tic = time.time()
    Dictionary = {'Transformation':'bspline',
                'FixedImage':FixedImage,
                'MovingImage':ResultImage,
                'FixedMask':FixedMask,
                'PyramidSchedule':[64, 32, 16, 8, 4, 2, 1],
                'NIterations':2000,
                'Alpha':0.6,
                'A':1000,
                'ElementSize':Config['ElementSize'],
                'ResultsDirectory':ResultsDirectory}
    DeformedImage, DeformedParameterMap = ElastixRegistration(Dictionary)
    Dice = ShowRegistration(GetSlice(FixedImage, Axis='Y'), GetSlice(DeformedImage, Axis='Y'))
    print('Dice coefficient: %.3f' % (Dice))
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Perform non-rigid registration %i min %i s' % (np.floor((Toc-Tic)/60),np.mod(Toc-Tic,60)) + '\n')


    #%% Registration Dice
    # Registration Dice coefficient

    Otsu = sitk.OtsuThresholdImageFilter()
    OtsuFixed = Otsu.Execute(FixedImage)
    OtsuDeformed = Otsu.Execute(DeformedImage)
    FixedArray = sitk.GetArrayFromImage(OtsuFixed).astype('bool') * 1
    DeformedArray = sitk.GetArrayFromImage(OtsuDeformed).astype('bool') * 1
    Dice = 2 * np.sum(FixedArray * DeformedArray) / np.sum(FixedArray + DeformedArray)
    print('\nDice coefficient of the full image: %.3f' % (Dice))
    LogFile.write('Dice coefficient of the registration: %.3f (-)' % (Dice) + '\n')


    #%% Write registration results
    # Write registration results
    print('\nWrite registration results')
    Tic = time.time()
    WriteMHD(FixedImage, ResultsDirectory, 'Fixed', PixelType='norm')
    WriteMHD(DeformedImage, ResultsDirectory, 'Registered', PixelType='norm')
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Write registration result image %.3f s' % (Toc - Tic) + '\n')

    #%% Transformix
    ## Use transformix to compute spatial jacobian
    print('\nUse transformix to compute jacobian')
    Tic = time.time()
    TransformixTransformations(MovingImage, TransformParameterMap, ResultsDirectory)
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Compute transformation jacobian in %i min %i s' % (np.floor((Toc - Tic) / 60), np.mod(Toc - Tic, 60)) + '\n')

    #%% Jacobian resampling
    # Resample jacobian to match with hFE
    print('\nResample jacobian')
    Tic = time.time()
    JacobianImage = sitk.ReadImage(ResultsDirectory + '/fullSpatialJacobian.mhd')
    JacobianImage.SetSpacing(FixedImage.GetSpacing())

    ## Resample Jacobian image
    Offset = JacobianImage.GetOrigin()
    Direction = JacobianImage.GetDirection()
    Orig_Size = np.array(JacobianImage.GetSize(), dtype='int')
    Orig_Spacing = JacobianImage.GetSpacing()

    New_Spacing = (0.9712, 0.9712, 0.9712)

    Resample = sitk.ResampleImageFilter()
    Resample.SetInterpolator = sitk.sitkLinear
    Resample.SetOutputDirection(Direction)
    Resample.SetOutputOrigin(Offset)
    Resample.SetOutputSpacing(New_Spacing)

    New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(New_Spacing))
    New_Size = np.ceil(New_Size).astype('int')
    New_Size = [int(s) for s in New_Size]
    Resample.SetSize(New_Size)

    ResampledJacobian = Resample.Execute(JacobianImage)
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Jacobian resample in %.3f s' % (Toc - Tic) + '\n')

    #%% Jacobian decomposition
    ## Perform jacobian unimodular decomposition
    print('\nUnimodular decomposition')
    Tic = time.time()
    SphericalCompression, IsovolumicDeformation = DecomposeJacobian(ResampledJacobian)
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Read and decompose jacobian in %i min %i s' % (np.floor((Toc - Tic) / 60), np.mod(Toc - Tic, 60)) + '\n')

    #%% Write results
    ## Write results
    print('\nWrite results to MHD')
    Tic = time.time()
    WriteMHD(SphericalCompression, ResultsDirectory, 'J', PixelType='float')
    WriteMHD(IsovolumicDeformation, ResultsDirectory, 'F_Tilde', PixelType='float')
    Toc = time.time()
    PrintTime(Tic, Toc)
    LogFile.write('Write decomposition results in %.3f s' % (Toc - Tic) + '\n\n')
    os.remove(os.path.join(ResultsDirectory, 'fullSpatialJacobian.mhd'))
    os.remove(os.path.join(ResultsDirectory, 'fullSpatialJacobian.raw'))

    #%% Store data
    # Store registration results
    print('\nStore registration results for sample ' + Sample)
    Data.loc[Index, 'Sample'] = Sample
    Data.loc[Index, 'Time'] = time.time() - SampleTime
    Data.loc[Index, 'Dice'] = Dice
#%%
print('\nRegistration done!')
FileName = os.path.join(WorkingDirectory, '04_Results/03_Registration', 'Results.csv')
Data.to_csv(FileName, index=False)
LogFile.close

#%% Show Dices

Figure, Axis = plt.subplots(1,1)
Axis.plot(Data['Sample'], Data['Dice'], color=(1,0,0), linestyle='none', marker='o', fillstyle='none')
Axis.set_xticklabels(Data['Sample'],rotation=90)
Axis.set_ylabel('Dice coefficient (-)')
plt.show()
# %%
