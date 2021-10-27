import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
from scipy import ndimage

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
    ElastixImageFilter.SetFixedImage(FixedImage)
    ElastixImageFilter.SetMovingImage(MovingImage)
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)
    # ElastixImageFilter.SetInitialTransformParameterFileName(ResultsDirectory + 'ReverseInitialTranslation.txt')
    ElastixImageFilter.LogToConsoleOn()
    ElastixImageFilter.Execute()

    ## Get results
    ResultImage = ElastixImageFilter.GetResultImage()  # How moving image is deformed
    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    ## Rename parameter file
    os.rename(Results_Path + '/TransformParameters.0.txt', ResultsDirectory + 'ReverseTransformParameters.txt')

    return ResultImage, TransformParameterMap
def PlotsImages(HRpQCT, uCT, C1=np.array([[0, 0, 0], [255, 0, 0]]), C2=np.array([[0, 0, 0], [0, 0, 255]]), Cbar=False, Resample=False):

    if Resample:
        # Resample HRpQCT slice for plot
        SpacingRatio = np.array([1.614505508264456, 1.614505508264456, 1.6145294470401392])
        Points = HRpQCT.shape / SpacingRatio
        Points_X = np.linspace(0, HRpQCT.shape[2] - 1, int(Points[2]))
        Points_Y = np.linspace(0, HRpQCT.shape[1] - 1, int(Points[1]))
        Points_Z = np.linspace(0, HRpQCT.shape[0] - 1, int(Points[0]))

    Resampled_HRpQCT = HRpQCT[:, :, int(HRpQCT.shape[2] / 2)]

    if Resample:
        Resampled_HRpQCT = Resampled_HRpQCT[Points_Z.astype('int'), :]
        Resampled_HRpQCT = Resampled_HRpQCT[:, Points_Y.astype('int')]

    # Create custom color map
    import matplotlib as mpl  # in python
    ColorMap1 = mpl.colors.ListedColormap(C1 / 255.0)
    ColorMap2 = mpl.colors.ListedColormap(C2 / 255.0)

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(Resampled_HRpQCT, cmap=ColorMap1, alpha=0.5)
    ColorMap = Axes.imshow(uCT[:, :, int(uCT.shape[2] / 2)], cmap=ColorMap2, alpha=0.5)
    if Cbar:
        plt.colorbar(ColorMap)
    Axes.axis('off')
    plt.show()
    plt.close(Figure)

    Resampled_HRpQCT = HRpQCT[int(HRpQCT.shape[0] / 2), :, :]

    if Resample:
        Resampled_HRpQCT = Resampled_HRpQCT[Points_Y.astype('int'), :]
        Resampled_HRpQCT = Resampled_HRpQCT[:, Points_X.astype('int')]

    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Axes.imshow(Resampled_HRpQCT, cmap=ColorMap1, alpha=0.5)
    ColorMap = Axes.imshow(uCT[int(uCT.shape[0] / 2), :, :], cmap=ColorMap2, alpha=0.5)
    if Cbar:
        plt.colorbar(ColorMap)
    Axes.axis('off')
    plt.show()
    plt.close(Figure)

    return
def TransformixTransformation(MovingImage,TransformParameterMap,ResultsDirectory):

    ## Compute jacobian of deformation field using transformix
    TransformixImageFilter = sitk.TransformixImageFilter()
    TransformixImageFilter.ComputeDeformationFieldOff()
    TransformixImageFilter.ComputeSpatialJacobianOff()
    TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOff()
    TransformixImageFilter.SetMovingImage(MovingImage)
    TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
    TransformixImageFilter.SetOutputDirectory(ResultsDirectory)

    TransformixImageFilter.Execute()

    ResultImage = TransformixImageFilter.GetResultImage()

    return ResultImage
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
uCTFolder = os.path.join(WorkingDirectory, '02_Data/03_uCT/')
HRpQCTFolder = os.path.join(WorkingDirectory, '02_Data/02_HRpQCT/')
ResultsFolder = os.path.join(WorkingDirectory, '04_Results/06_FractureLinePrediction/')

SampleList = [Dir for Dir in os.listdir(HRpQCTFolder) if os.path.isdir(HRpQCTFolder+Dir)]
SampleList.sort()

Index = 0
# for Index in range(len(SampleList)):

## Set sample paths
Sample = SampleList[Index]
uCT_Path = os.path.join(uCTFolder + Sample + '/')
HRpQCT_Path = os.path.join(HRpQCTFolder + Sample + '/')
Results_Path = os.path.join(ResultsFolder + Sample + '/')

## List files
uCT_Files = os.listdir(uCT_Path)
uCT_Mask_File = [File for File in uCT_Files if File.endswith('MASK.mhd')][0]
uCT_Scan_File = uCT_Mask_File[:-13] + '.mhd'
HRpQCT_Files = [File for File in os.listdir(HRpQCT_Path) if File.endswith('.mhd')]
HRpQCT_Files.sort()

## Load uCT mask and rotate and segment it
uCT_Mask = sitk.ReadImage(uCT_Path + uCT_Mask_File)
Spacing = uCT_Mask.GetSpacing()
Origin = uCT_Mask.GetOrigin()
Direction = uCT_Mask.GetDirection()
uCT_Mask_Array = sitk.GetArrayFromImage(uCT_Mask)

DSCs = pd.read_csv(Results_Path + 'DSCs.csv')
BestAngle = int(DSCs.loc[DSCs['DSC'].idxmax(), 'Angle'])
uCT_Mask_Array = ndimage.rotate(uCT_Mask_Array, -BestAngle, (1, 2), reshape=True)
uCT_Mask_Array = np.rot90(uCT_Mask_Array, 2, (0, 1))

Otsu_Filter = sitk.OtsuThresholdImageFilter()
Otsu_Filter.SetInsideValue(0)
Otsu_Filter.SetOutsideValue(1)
Segmentation = Otsu_Filter.Execute(uCT_Mask)
R_Threshold = Otsu_Filter.GetThreshold()
uCT_Mask_Array[uCT_Mask_Array <= R_Threshold] = 0
uCT_Mask_Array[uCT_Mask_Array > 0] = 1

## Shift uCT mask according to initial translation
ParameterMapFile = open(Results_Path + 'InitialTranslation.txt','r')
Text = ParameterMapFile.read()
Start = Text.find('TransformParameters')
Stop = Start + Text[Start:].find('\n')
Parameters = Text[Start+20:Stop-1].split()
ReverseParameters = -np.array(Parameters).astype('float')
ParameterMapFile.close()

uCT_Mask = sitk.GetImageFromArray(uCT_Mask_Array)
uCT_Mask.SetSpacing(Spacing)
uCT_Mask.SetOrigin(Origin)
uCT_Mask.SetDirection(Direction)
# uCT_Mask.SetOrigin(tuple(np.array(Origin) + ReverseParameters))

## Load HR-pQCT
HRpQCT_Cort = sitk.ReadImage(HRpQCT_Path + HRpQCT_Files[0])
HRpQCT_Trab = sitk.ReadImage(HRpQCT_Path + HRpQCT_Files[2])
HRpQCT_Mask = HRpQCT_Cort + HRpQCT_Trab
HRpQCT_Mask_Array = sitk.GetArrayFromImage(HRpQCT_Mask)

PlotsImages(HRpQCT_Mask_Array,uCT_Mask_Array,Resample=True)

## Register uCT mask on HRpQCT
Dictionary = {'Transformations': 'rigid',
              'FixedImage': HRpQCT_Mask,
              'MovingImage': uCT_Mask,
              'PyramidSchedule': [50, 20, 10],
              'NIterations': 2000,
              'Alpha': 0.6,
              'A': 10,
              'ResultsDirectory': Results_Path}
ResultImage, TransformParameterMap = ElastixRegistration(Dictionary)
ResultImage_Array = sitk.GetArrayViewFromImage(ResultImage)

Spacing = uCT_Mask.GetSpacing()
Offset = uCT_Mask.GetOrigin()
WriteMHD(uCT_Mask_Array, Spacing, Offset, Results_Path, 'uCT_Registered', PixelType='float')

PlotsImages(HRpQCT_Mask_Array,ResultImage_Array)


TransformParameterMap = sitk.ReadParameterFile(Results_Path + 'TransformParameters.txt')
uCT_Mask_Transformed = TransformixTransformation(uCT_Cropped, TransformParameterMap, ResultsDirectory=Results_Path)



## Re-binarize mask
HRpQCT_Mask[HRpQCT_Mask < 0.5] = 0
HRpQCT_Mask[HRpQCT_Mask >= 0.5] = 1
HRpQCT_Mask = sitk.GetImageFromArray(HRpQCT_Mask)
HRpQCT_Mask.SetSpacing(HRpQCT_Mask_Image.GetSpacing())

# 06 Resample masks to get same sampling as HR-pQCT
Offset = HRpQCT_Mask_Image.GetOrigin()
Direction = HRpQCT_Mask_Image.GetDirection()
Orig_Size = np.array(HRpQCT_Mask.GetSize(), dtype=np.int)
Orig_Spacing = HRpQCT_Mask_Image.GetSpacing()

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

ResampledImage = Resample.Execute(HRpQCT_Mask)
HRpQCT_Mask = sitk.GetArrayFromImage(ResampledImage)

# Do the same fot other images
Orig_Size = np.array(uCT_Mask.GetSize(), dtype=np.int)

New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(New_Spacing))
New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
New_Size = [int(s) for s in New_Size]
Resample.SetSize(New_Size)

ResampledImage = Resample.Execute(uCT_Mask)
uCT_Mask = sitk.GetArrayFromImage(ResampledImage)

ResampledImage = Resample.Execute(uCT_Scan)
uCT_Scan = sitk.GetArrayFromImage(ResampledImage)

ResampledImage = Resample.Execute(HRpQCT_Scan)
HRpQCT_Scan = sitk.GetArrayFromImage(ResampledImage)

ResampledImage = Resample.Execute(J)
J = sitk.GetArrayFromImage(ResampledImage)

ResampledImage = Resample.Execute(F_Tilde)
F_Tilde = sitk.GetArrayFromImage(ResampledImage)




# 05 Compute first slice to keep
DSCs = pd.DataFrame()
for Slice in range(10):

    if np.sum(uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :]) > 0:
        DSC = 2 * np.sum(uCT_Mask[Slice, :, :] * HRpQCT_Mask[Slice, :, :]) / np.sum(uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :])
        DSCs = DSCs.append({'Slice':Slice,'DSC':DSC},ignore_index=True)

DSCs['Slope'] = 0
for i in DSCs.index:
    if i > 0:
        DSCs.loc[i,'Slope'] = DSCs.loc[i,'DSC'] - DSCs.loc[i-1,'DSC']

Tolerance = 1E-2
FirstSlice = DSCs.loc[DSCs[DSCs['Slope'] > Tolerance]['Slice'].idxmax(),'Slice'].astype('int')

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.imshow(uCT_Mask[FirstSlice-1, :, :], cmap='bone',alpha=0.5)
Axes.imshow(HRpQCT_Mask[FirstSlice-1, :, :], cmap='bone',alpha=0.5)
plt.show()
plt.close(Figure)

# 06 Compute last slice to keep
DSCs = pd.DataFrame()
for SliceIndex in range(10):

    Slice = uCT_Mask.shape[0]-SliceIndex-1

    if np.sum(uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :]) > 0:
        DSC = 2 * np.sum(uCT_Mask[Slice, :, :] * HRpQCT_Mask[Slice, :, :]) / np.sum(
            uCT_Mask[Slice, :, :] + HRpQCT_Mask[Slice, :, :])
        DSCs = DSCs.append({'Slice': Slice, 'DSC': DSC}, ignore_index=True)

DSCs['Slope'] = 0
for i in DSCs.index:
    if i > 0:
        DSCs.loc[i, 'Slope'] = DSCs.loc[i, 'DSC'] - DSCs.loc[i - 1, 'DSC']

Tolerance = 1E-2
LastSlice = DSCs.loc[DSCs[DSCs['Slope'] > Tolerance]['Slice'].idxmin(),'Slice'].astype('int')

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.imshow(uCT_Mask[LastSlice+1, :, :], cmap='bone', alpha=0.5)
Axes.imshow(HRpQCT_Mask[LastSlice+1, :, :], cmap='bone', alpha=0.5)
plt.show()
plt.close(Figure)

# Crop images
uCT_Cropped = uCT_Scan[FirstSlice:LastSlice, :, :]
HRpQCT_Cropped = HRpQCT_Scan[FirstSlice:LastSlice, :, :]

uMask_Cropped = uCT_Mask[FirstSlice:LastSlice, :, :]
HRMask_Cropped = HRpQCT_Mask[FirstSlice:LastSlice, :, :]

J_Cropped = J[FirstSlice:LastSlice, :, :]
F_Tilde_Cropped = F_Tilde[FirstSlice:LastSlice, :, :]


Origin = np.array([0,0,New_Spacing[2]]) * FirstSlice
WriteMHD(uCT_Cropped,New_Spacing,Origin,SamplePath, 'uCT_Cropped', PixelType='float')
WriteMHD(HRpQCT_Cropped,New_Spacing,Origin,SamplePath, 'HRpQCT_Cropped', PixelType='float')
WriteMHD(uMask_Cropped,New_Spacing,Origin,SamplePath, 'uCT_Mask_Cropped', PixelType='float')
WriteMHD(HRMask_Cropped,New_Spacing,Origin,SamplePath, 'HRpQCT_Mask_Cropped', PixelType='float')
WriteMHD(J_Cropped,New_Spacing,Origin,SamplePath, 'J_Cropped', PixelType='float')
WriteMHD(F_Tilde_Cropped,New_Spacing,Origin,SamplePath, 'F_Tilde_Cropped', PixelType='float')



Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.imshow(uCT_Cropped[:, :, int(uCT_Mask.shape[2]/2)], cmap='bone',alpha=0.5)
Axes.imshow(HRpQCT_Scan[FirstSlice:LastSlice, :, int(HRpQCT_Mask.shape[2]/2)], cmap='bone',alpha=0.5)
plt.show()
plt.close(Figure)


Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.plot(DSCs['Slope'], linestyle='--', marker='o', color=(1,0,0))
plt.show()
plt.close(Figure)

