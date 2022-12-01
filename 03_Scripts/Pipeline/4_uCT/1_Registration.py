#%% #!/usr/bin/env python3
# 00 Initialization

import pandas as pd
from Utils import *

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

#%% Functions
# 01 Define functions
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
    for Image in [SphericalCompression, IsovolumicDeformation]:
        Image.SetSpacing(JacobianImage.GetSpacing())
        Image.SetDirection(JacobianImage.GetDirection())
        Image.SetOrigin(JacobianImage.GetOrigin())

    return SphericalCompression, IsovolumicDeformation


#%% Loading
# 02 Set paths and load data

CW, Data, Scripts, Results = SetDirectories('FRACTIB')

SampleList = pd.read_csv(str(Data / '02_uCT' / 'BMD_Values.csv'))
LogFile = open(str(Results / '04_Registration' / 'Registration.log'),'w+')
LogFile.write('Registration Log File\n\n')

Data = pd.DataFrame()

#%% Set index
for Index in range(len(SampleList)):

    SampleTime = time.time()
    #%% uCT files loading
    # 05 Load uCT files

    print('\nLoad uCT files')
    Tic = time.time()

    Sample = SampleList.loc[Index,'Sample']
    ResultsDirectory = str(Results / '04_Registration' / Sample)
    os.makedirs(ResultsDirectory,exist_ok=True)

    SampleData = {'Sample': Sample}
    LogFile.write('Sample: ' + Sample + '\n')

    SampleDirectory = str(Data / '02_uCT' / Sample) + '/'
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
    ConfigFile = str(Scripts / '03_hFE' / 'ConfigFile.yaml')

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

#%% Perform initial rotation
# Rotation

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


