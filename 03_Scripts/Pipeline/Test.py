#%% #!/usr/bin/env python3
# Imports
import sys
import vtk
import time
import struct
import numpy as np
import sympy as sp
import SimpleITK as sitk
from pathlib import Path
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy

#%% Functions
# Function definitions
def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results
def RotationMatrix(Alpha=0, Beta=0, Gamma=0):

    Rx = sp.Matrix([[1,             0,              0],
                    [0, sp.cos(Alpha), -sp.sin(Alpha)],
                    [0, sp.sin(Alpha),  sp.cos(Alpha)]])

    Ry = sp.Matrix([[ sp.cos(Beta), 0, sp.sin(Beta)],
                    [0,             1,              0],
                    [-sp.sin(Beta), 0, sp.cos(Beta)]])

    Rz = sp.Matrix([[sp.cos(Gamma), -sp.sin(Gamma), 0],
                    [sp.sin(Gamma),  sp.cos(Gamma), 0],
                    [0,             0,              1]])

    R = Rz * Ry * Rx

    return np.array(R, dtype='float')
def GetParameterMap(FileName):

    """
    Builds parameter map according to given file
    """

    File = open(FileName, 'r')
    Text = File.read()
    Start = Text.find('(')
    Stop = Text.find(')')

    ParameterMap = {}
    while Start-Stop+1:
        Line = Text[Start+1:Stop]
        Sep = Line.find(' ')
        Name = Line[:Sep]
        Parameter = Line[Sep+1:]

        if Line[Sep+1:].find(' ')+1:
            ParameterMap[Name] = [P for P in Parameter.split()]

        else:
            ParameterMap[Name] = [Parameter]

        Start = Stop + Text[Stop:].find('(')
        Stop = Start + Text[Start:].find(')')

    File.close()

    return ParameterMap
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
def ProgressStart(Text='Progress'):
    global CurrentProgress
    sys.stdout.write('\n' + Text + '|')
    CurrentProgress = 0
    sys.stdout.flush()
    return
def ProgressNext(Progress):
    global CurrentProgress
    if Progress > CurrentProgress:
        CurrentProgress += 1
        sys.stdout.write('=')
        sys.stdout.flush()
    return
def ProgressEnd():
    sys.stdout.write('|\n')
    sys.stdout.flush()
    return
def ProcessTiming(StartStop:bool, Process='Progress'):

    if StartStop*1 == 1:
        global Tic
        Tic = time.time()
        ProgressStart(Process)
    elif StartStop*1 == 0:
        ProgressEnd()
        Toc = time.time()
        PrintTime(Tic, Toc)
    return
def Adjust_Image_Size(Image, CoarseFactor, CropZ='Crop'):

    """
    Adapted from Denis's utils_SA.py
    Images are adjusted according to CropType:
    0 = CropType.expand     (Expand image by copying layers)
    1 = CropType.crop       (Crop image)
    2 = CropType.variable   (Either crop or expand, depending on what includes less layers)
    """

    # Measure image shape
    IMDimX = np.shape(Image)[0]
    IMDimY = np.shape(Image)[1]
    IMDimZ = np.shape(Image)[2]

    AddDimX = CoarseFactor - (IMDimX % CoarseFactor)
    AddDimY = CoarseFactor - (IMDimY % CoarseFactor)

    # adjust in x and y direction
    Shape_Diff = [AddDimX, AddDimY]
    IMG_XY_Adjusted = np.lib.pad(Image,
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

    return Image_Adjusted
def ShowRegistration(Fixed, Moving, Slice=None, Title=None, Axis='Z'):

    FixedArray = sitk.GetArrayFromImage(Fixed)
    MovingArray = sitk.GetArrayFromImage(Moving)

    if len(np.unique(FixedArray)) > 2 or len(np.unique(MovingArray)) > 2:
        Otsu = sitk.OtsuThresholdImageFilter()
        Otsu.SetInsideValue(0)
        Otsu.SetOutsideValue(1)
        Fixed_Bin = Otsu.Execute(Fixed)
        Moving_Bin = Otsu.Execute(Moving)
        FixedArray = sitk.GetArrayFromImage(Fixed_Bin)
        MovingArray = sitk.GetArrayFromImage(Moving_Bin)

    if Fixed.GetDimension() == 3:
        Array = np.zeros((Fixed.GetSize()[2], Fixed.GetSize()[1], Fixed.GetSize()[0], 3))
        Array[:,:,:,0] = FixedArray
        Array[:,:,:,1] = MovingArray
        Array[:,:,:,2] = MovingArray
        
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

    else:
        Array = np.zeros((Fixed.GetSize()[1], Fixed.GetSize()[0], 3))
        Array[:,:,0] = FixedArray
        Array[:,:,1] = MovingArray
        Array[:,:,2] = MovingArray

    Figure, Axis = plt.subplots()
    Axis.imshow(Array,interpolation=None)
    Axis.axis('Off')
    
    if(Title):
        Axis.set_title(Title)

    plt.show(Figure)

    return
def Resample(Image, Factor=None, Size=None, Spacing=None):

    Dimension = Image.GetDimension()
    OriginalSpacing = np.array(Image.GetSpacing())
    OriginalSize = np.array(Image.GetSize())
    PhysicalSize = OriginalSize * OriginalSpacing

    Origin = Image.GetOrigin()
    Direction = Image.GetDirection()
    Center = OriginalSize * OriginalSpacing / 2

    if Factor:
        NewSize = [round(Size/Factor) for Size in Image.GetSize()] 
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Size:
        NewSize = Size
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Spacing:
        NewSpacing = Spacing
        NewSize = [Size/Spacing + 1 for Size,Spacing in zip(PhysicalSize, NewSpacing)]
    
    NewImage = sitk.Image(NewSize, Image.GetPixelIDValue())
    NewImage.SetOrigin(Origin)
    NewImage.SetDirection(Direction)
    NewImage.SetSpacing(NewSpacing)
  
    Transform = sitk.TranslationTransform(Dimension)
    
    return sitk.Resample(Image, NewImage, Transform, sitk.sitkLinear, 0.0)
def ShowSlice(Image, Slice=None, Title=None, Axis='Z'):

    Array = sitk.GetArrayFromImage(Image)

    if Image.GetDimension() == 3:
        
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
def Raw(Image, OutputFileName, PixelType):

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

    File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape, order='C')
    File[:] = CastedImageArray[:]
    del File

    return
def MHD(Image, FileName, PixelType='uint'):

    print('\nWrite MHD')
    Tic = time.time()

    if PixelType == 'short' or PixelType == 'float':
        if Image.GetDimension() == 2:

            Array_3D = np.zeros((1,ImageArray.shape[0],ImageArray.shape[1]))

            for j in range(ImageArray.shape[0]):
                for i in range(ImageArray.shape[1]):
                    Array_3D[0,j,i] = ImageArray[j,i]

            ImageArray = Array_3D

    nx, ny, nz = Image.GetSize()

    lx, ly, lz = Image.GetSpacing()

    TransformMatrix = '1 0 0 0 1 0 0 0 1'
    Offset = str(np.array(Image.GetOrigin(),'int'))[1:-1]
    CenterOfRotation = '0 0 0'
    AnatomicalOrientation = 'LPS'

    outs = open(FileName + '.mhd', 'w')
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

    Raw(Image, FileName + '.raw', PixelType)

    Toc = time.time()
    PrintTime(Tic, Toc)

    return


class Read:

    def __init__(self):
        pass
    
    def Get_AIM_Ints(self, File):

        """
        Function by Glen L. Niebur, University of Notre Dame (2010)
        reads the integer data of an AIM file to find its header length
        """

        nheaderints = 32
        File.seek(0)
        binints = File.read(nheaderints * 4)
        header_int = struct.unpack("=32i", binints)

        return header_int

    def AIM(self, File):

        """
        Reads an AIM file and provides
        the corresponding itk image additional
        data (i.e. spacing, calibration data, 
        and header)

        from Denis hFE pipeline
        """

        print('\nRead AIM')
        Tic = time.time()

        # read header
        with open(File, 'rb') as f:
            AIM_Ints = self.Get_AIM_Ints(f)
            # check AIM version
            if int(AIM_Ints[5]) == 16:
                if int(AIM_Ints[10]) == 131074:
                    Format = "short"
                elif int(AIM_Ints[10]) == 65537:
                    Format = "char"
                elif int(AIM_Ints[10]) == 1376257:
                    Format = "bin compressed"
                    print("     -> format " + Format + " not supported! Exiting!")
                    exit(1)
                else:
                    Format = "unknown"
                    exit(1)
                Header = f.read(AIM_Ints[2])
                Header_Length = len(Header) + 160
                Extents = (0, AIM_Ints[14] - 1, 0, AIM_Ints[15] - 1, 0, AIM_Ints[16] - 1)
            else:
                print("     -> version 030")
                if int(AIM_Ints[17]) == 131074:
                    Format = "short"
                    print("     -> format " + Format)
                elif int(AIM_Ints[17]) == 65537:
                    Format = "char"
                    print("     -> format " + Format)
                elif int(AIM_Ints[17]) == 1376257:
                    Format = "bin compressed"
                    print("     -> format " + Format + " not supported! Exiting!")
                    exit(1)
                else:
                    Format = "unknown"
                    print("     -> format " + Format + "! Exiting!")
                    exit(1)
                Header = f.read(AIM_Ints[8])
                Header_Length = len(Header) + 280
                Extents = (0, AIM_Ints[24] - 1, 0, AIM_Ints[26] - 1, 0, AIM_Ints[28] - 1)

        # collect data from header if existing
        # header = re.sub('(?i) +', ' ', header)
        Header = Header.split('\n'.encode())
        Header.pop(0)
        Header.pop(0)
        Header.pop(0)
        Header.pop(0)
        Scaling = None
        Slope = None
        Intercept = None
        IPLPostScanScaling = 1
        for Line in Header:
            if Line.find(b"Orig-ISQ-Dim-p") > -1:
                origdimp = ([int(s) for s in Line.split(b" ") if s.isdigit()])

            if Line.find("Orig-ISQ-Dim-um".encode()) > -1:
                origdimum = ([int(s) for s in Line.split(b" ") if s.isdigit()])

            if Line.find("Orig-GOBJ-Dim-p".encode()) > -1:
                origdimp = ([int(s) for s in Line.split(b" ") if s.isdigit()])

            if Line.find("Orig-GOBJ-Dim-um".encode()) > -1:
                origdimum = ([int(s) for s in Line.split(b" ") if s.isdigit()])

            if Line.find("Scaled by factor".encode()) > -1:
                Scaling = float(Line.split(" ".encode())[-1])
            if Line.find("Density: intercept".encode()) > -1:
                Intercept = float(Line.split(" ".encode())[-1])
            if Line.find("Density: slope".encode()) > -1:
                Slope = float(Line.split(" ".encode())[-1])
            # if el_size scale was applied, the above still takes the original voxel size. This function works
            # only if an isotropic scaling was applied!!!!
            if Line.find("scale".encode()) > -1:
                IPLPostScanScaling = float(Line.split(" ".encode())[-1])
        # Spacing is calculated from Original Dimensions. This is wrong, when the images were coarsened and
        # the voxel size is not anymore corresponding to the original scanning resolution!

        try:
            Spacing = IPLPostScanScaling * (
                np.around(np.asarray(origdimum) / np.asarray(origdimp) / 1000, 5)
            )
        except:
            pass

        # read AIM with vtk
        Reader = vtk.vtkImageReader2()
        Reader.SetFileName(File)
        Reader.SetDataByteOrderToLittleEndian()
        Reader.SetFileDimensionality(3)
        Reader.SetDataExtent(Extents)
        Reader.SetHeaderSize(Header_Length)
        if Format == "short":
            Reader.SetDataScalarTypeToShort()
        elif Format == "char":
            Reader.SetDataScalarTypeToChar()
        Reader.SetDataSpacing(Spacing)
        Reader.Update()
        VTK_Image = Reader.GetOutput()


        # Convert VTK to numpy
        Data = VTK_Image.GetPointData().GetScalars()
        Dimension = VTK_Image.GetDimensions()
        Numpy_Image = vtk_to_numpy(Data)
        Numpy_Image = Numpy_Image.reshape(Dimension[2], Dimension[1], Dimension[0])

        # Y symmetry (thanks Michi for notified this!)
        Numpy_Image = Numpy_Image[:,::-1,:]

        # Converty numpy to ITK image
        Image = sitk.GetImageFromArray(Numpy_Image)
        Image.SetSpacing(Spacing)
        Image.SetOrigin([0.0, 0.0, 0.0])

        AdditionalData = {'Scaling':Scaling,
                        'Slope':Slope,
                        'Intercept':Intercept,
                        'Header':Header}

        Toc = time.time()
        PrintTime(Tic, Toc)

        return Image, AdditionalData

    def ISQ(self, File, BMD=False, Info=False):

        """
        This function read an ISQ file from Scanco and return an ITK image and additional data.
        
        Adapted from https://github.com/mdoube/BoneJ/blob/master/src/org/bonej/io/ISQReader.java
        
        Little endian byte order (the least significant bit occupies the lowest memory position.
        00   char    check[16];              // CTDATA-HEADER_V1
        16   int     data_type;
        20   int     nr_of_bytes;
        24   int     nr_of_blocks;
        28   int     patient_index;          //p.skip(28);
        32   int     scanner_id;				//p.skip(32);
        36   int     creation_date[2];		//P.skip(36);
        44   int     dimx_p;					//p.skip(44);
        48   int     dimy_p;
        52   int     dimz_p;
        56   int     dimx_um;				//p.skip(56);
        60   int     dimy_um;
        64   int     dimz_um;
        68   int     slice_thickness_um;		//p.skip(68);
        72   int     slice_increment_um;		//p.skip(72);
        76   int     slice_1_pos_um;
        80   int     min_data_value;
        84   int     max_data_value;
        88   int     mu_scaling;             //p.skip(88);  /* p(x,y,z)/mu_scaling = value [1/cm]
        92	int     nr_of_samples;
        96	int     nr_of_projections;
        100  int     scandist_um;
        104  int     scanner_type;
        108  int     sampletime_us;
        112  int     index_measurement;
        116  int     site;                   //coded value
        120  int     reference_line_um;
        124  int     recon_alg;              //coded value
        128  char    name[40]; 		 		//p.skip(128);
        168  int     energy;        /* V     //p.skip(168);
        172  int     intensity;     /* uA    //p.skip(172);
        ...
        508 int     data_offset;     /* in 512-byte-blocks  //p.skip(508);
        * So the first 16 bytes are a string 'CTDATA-HEADER_V1', used to identify
        * the type of data. The 'int' are all 4-byte integers.
        *
        * dimx_p is the dimension in pixels, dimx_um the dimensions in micrometer
        *
        * So dimx_p is at byte-offset 40, then dimy_p at 44, dimz_p (=number of
        * slices) at 48.
        *
        * The microCT calculates so called 'x-ray linear attenuation' values. These
        * (float) values are scaled with 'mu_scaling' (see header, e.g. 4096) to
        * get to the signed 2-byte integers values that we save in the .isq file.
        *
        * e.g. Pixel value 8192 corresponds to lin. att. coeff. of 2.0 [1/cm]
        * (8192/4096)
        *
        * Following to the headers is the data part. It is in 2-byte short integers
        * (signed) and starts from the top-left pixel of slice 1 to the left, then
        * the next line follows, until the last pixel of the last sclice in the
        * lower right.
        """

        print('\nRead ISQ file')
        Tic = time.time()

        try:
            f = open(File, 'rb')
        except IOError:
            print("\n **ERROR**: ISQReader: intput file ' % s' not found!\n\n" % File)
            print('\n E N D E D  with ERRORS \n\n')

        for Index in np.arange(0, 200, 4):
            f.seek(Index)
            print('Index %s :          %s' % (Index, struct.unpack('i', f.read(4))[0]))
            f.seek(Index)
            try:
                print('Index %s :          %s' % (Index, struct.unpack('c', f.read(4))[0]))
            except:
                print('')

        f.seek(32)
        CT_ID = struct.unpack('i', f.read(4))[0]
        print('scanner ID:                 ', CT_ID)

        if CT_ID != 6020:
            print('!!! unknown muCT -> no Slope and Intercept known !!!')

        f.seek(28)
        #    sample_nb = struct.unpack('i', f.read(4))[0]

        f.seek(108)
        Scanning_time = struct.unpack('i', f.read(4))[0] / 1000
        print('Scanning time in ms:         ', Scanning_time)

        f.seek(168)
        Energy = struct.unpack('i', f.read(4))[0] / 1000.
        print('Energy in keV:              ', Energy)

        f.seek(172)
        Current = struct.unpack('i', f.read(4))[0]
        print('Current in muA:             ', Current)

        f.seek(44)
        X_pixel = struct.unpack('i', f.read(4))[0]
        print('Nb X pixel:                 ', X_pixel)

        f.seek(48)
        Y_pixel = struct.unpack('i', f.read(4))[0]
        print('Nb Y pixel:                 ', Y_pixel)

        f.seek(52)
        Z_pixel = struct.unpack('i', f.read(4))[0]
        print('Nb Z pixel:                 ', Z_pixel)

        f.seek(56)
        Res_General_X = struct.unpack('i', f.read(4))[0]
        print('Resolution general X in mu: ', Res_General_X)

        f.seek(60)
        Res_General_Y = struct.unpack('i', f.read(4))[0]
        print('Resolution general Y in mu: ', Res_General_Y)

        f.seek(64)
        Res_General_Z = struct.unpack('i', f.read(4))[0]
        print('Resolution general Z in mu: ', Res_General_Z)

        Res_X = Res_General_X / float(X_pixel)
        print('Pixel resolution X in mu:    %.2f' % Res_X)

        Res_Y = Res_General_Y / float(Y_pixel)
        print('Pixel resolution Y in mu:    %.2f' % Res_Y)

        Res_Z = Res_General_Z / float(Z_pixel)
        print('Pixel resolution Z in mu:    %.2f' % Res_Z)

        Header_Txt = ['scanner ID:                 %s' % CT_ID,
                    'scaning time in ms:         %s' % Scanning_time,
                    'scaning time in ms:         %s' % Scanning_time,
                    'Energy in keV:              %s' % Energy,
                    'Current in muA:             %s' % Current,
                    'nb X pixel:                 %s' % X_pixel,
                    'nb Y pixel:                 %s' % Y_pixel,
                    'nb Z pixel:                 %s' % Z_pixel,
                    'resolution general X in mu: %s' % Res_General_X,
                    'resolution general Y in mu: %s' % Res_General_Y,
                    'resolution general Z in mu: %s' % Res_General_Z,
                    'pixel resolution X in mu:   %.2f' % Res_X,
                    'pixel resolution Y in mu:   %.2f' % Res_Y,
                    'pixel resolution Z in mu:   %.2f' % Res_Z]
        #    np.savetxt(inFileName.split('.')[0]+'.txt', Header_Txt)

        if Info:
            Write_File = open(File.split('.')[0] + '_info.txt', 'w')
            for Item in Header_Txt:
                Write_File.write("%s\n" % Item)
            Write_File.close()

        f.seek(44)
        Header = np.zeros(6)
        for i in range(0, 6):
            Header[i] = struct.unpack('i', f.read(4))[0]
        print(Header)

        ElementSpacing = [Header[3] / Header[0] / 1000, Header[4] / Header[1] / 1000, Header[5] / Header[2] / 1000]
        f.seek(508)

        HeaderSize = 512 * (1 + struct.unpack('i', f.read(4))[0])
        f.seek(HeaderSize)


        VoxelModel = np.fromfile(f, dtype='i2')
        # VoxelModel = np.fromfile(f, dtype=np.float)

        NDim = [int(Header[0]), int(Header[1]), int(Header[2])]
        LDim = [float(ElementSpacing[0]), float(ElementSpacing[1]), float(ElementSpacing[2])]

        AdditionalData = {'-LDim': LDim,
                        '-NDim': NDim,
                        'ElementSpacing': LDim,
                        'DimSize': NDim,
                        'HeaderSize': HeaderSize,
                        'TransformMatrix': [1, 0, 0, 0, 1, 0, 0, 0, 1],
                        'CenterOfRotation': [0.0, 0.0, 0.0],
                        'Offset': [0.0, 0.0, 0.0],
                        'AnatomicalOrientation': 'LPS',
                        'ElementType': 'int16',
                        'ElementDataFile': File}

        Toc = time.time()
        PrintTime(Tic, Toc)

        print('\nReshape data')
        Tic = time.time()

        try:
            VoxelModel = VoxelModel.reshape((NDim[2], NDim[1], NDim[0]))
            f.close()
            del f

        except:
            # if the length does not fit the dimensions (len(VoxelModel) != NDim[2] * NDim[1] * NDim[0]),
            # add an offset with seek to reshape the image -> actualise length, delta *2 = seek

            Offset = (len(VoxelModel) - (NDim[2] * NDim[1] * NDim[0]))
            f.seek(0)
            VoxelModel = np.fromfile(f, dtype='i2')

            print('len(VoxelModel) = ', len(VoxelModel))
            print('Should be ', (NDim[2] * NDim[1] * NDim[0]))
            print('Delta:', len(VoxelModel) - (NDim[2] * NDim[1] * NDim[0]))

            f.seek((len(VoxelModel) - (NDim[2] * NDim[1] * NDim[0])) * 2)
            VoxelModel = np.fromfile(f, dtype='i2')
            f.close()
            del f

            VoxelModel = VoxelModel.reshape((NDim[2], NDim[1], NDim[0]))
            # the image is flipped by the Offset --> change the order to obtain the continuous image:
            VoxelModel = np.c_[VoxelModel[:, :, -Offset:], VoxelModel[:, :, :(VoxelModel.shape[2] - Offset)]]


        Toc = time.time()
        PrintTime(Tic, Toc)

        if CT_ID == 6020 and BMD is True:
            # BE CAREFULL, THIS IS FOR BMD CONVERSION:
            print('muCT 100 of ISTB detected, IS IT CORRECT?')
            Slope = 369.154  # ! ATTENTION, dependent on voltage, Current and time!!!
            Intercept = -191.56
            try:
                VoxelModel = VoxelModel.astype('i4')
                VoxelModel *= Slope
                VoxelModel += Intercept
            except:
                print('\n********* memory not sufficient for BMD values ************\n')

        # Convert numpy array to image
        Image = sitk.GetImageFromArray(VoxelModel)
        Image.SetSpacing(LDim[::-1])
        Image.SetOrigin([0.0, 0.0, 0.0])

        return Image, AdditionalData

#%% Loading
# Load data

CW, Data, Scripts, Results = SetDirectories('FRACTIB')

FileNames = {}
# Common area
FileNames['Common'] = str(Results / '03_hFE' / '432_L_77_F' / 'CommonMask.mhd')
FileNames['Common_uCT'] = str(Results / '05_Localizations' / '432_L_77_F' / 'CommonMask.mhd')
Common_hFE = sitk.ReadImage(FileNames['Common'])
Common_uCT = sitk.ReadImage(FileNames['Common_uCT'])


# Transform parameters
FileNames['InitialTransform'] = str(Results / '05_Localizations' / '432_L_77_F' / 'InitialTransform.txt')
FileNames['Transform'] = str(Results / '05_Localizations' / '432_L_77_F' / 'TransformParameters.0.txt')


uCT_File = 'C0001901_reso_0.098_DOWNSCALED_FULLMASK.mhd'
uCT_Mask = sitk.ReadImage(str(Data / '02_uCT' / '432_L_77_F' / uCT_File))

HRpQCT_File = 'C0003103'
Cort_File = str(Data / '01_HRpQCT' / '432_L_77_F' / (HRpQCT_File + '_CORT_MASK_UNCOMP.AIM'))
Trab_File = str(Data / '01_HRpQCT' / '432_L_77_F' / (HRpQCT_File + '_TRAB_MASK_UNCOMP.AIM'))
Reader = Read()
HRpQCT_Cort_Mask, AdditionalData = Reader.AIM(Cort_File)
HRpQCT_Trab_Mask, AdditionalData = Reader.AIM(Trab_File)
HRpQCT_Mask = HRpQCT_Cort_Mask + HRpQCT_Trab_Mask


#%% Meshing
# Mesh and resample

Image = Common_uCT
Spacing = np.array(Image.GetSpacing())
CoarseFactor = int(round(1.2747 / Spacing[0]))
FEelSize = np.copy(Spacing) * CoarseFactor

ShowSlice(Image)
Resampled = Resample(Image, Factor=CoarseFactor)
ShowSlice(Resampled)
MHD(Resampled, 'Resampled')

Array = sitk.GetArrayFromImage(Resampled)
# Adjusted_uCT = Adjust_Image_Size(Array_uCT, CoarseFactor, CropZ='Crop')

Coords = np.array(np.where(Array.transpose((2,1,0))))
Points = Coords * np.array([Resampled.GetSpacing()]).T
Points = Points.T

#%% Inverse transform 1
# Transform point

FT = GetParameterMap(FileNames['Transform'])
C3 = np.array(FT['CenterOfRotationPoint'], 'float')
P3 = np.array(FT['TransformParameters'],'float')
R3 = RotationMatrix(-P3[0], -P3[1], -P3[2])
T3 = -np.array(P3[3:])

R3 = np.linalg.inv(R3)

ProcessTiming(1,'Transform points ')    
TransformedPoints = []

for iP, Point in enumerate(Points):

    TP = np.dot(R3, Point - C3) + C3 - T3
    TransformedPoints.append(TP)

    Progress = iP / len(Points) * 20
    ProgressNext(Progress)
ProcessTiming(0)
TP3 = np.array(TransformedPoints)



#%% Inverse transform 2
# Transform Point
IT = sitk.ReadTransform(FileNames['InitialTransform'])
C2 = np.array(IT.GetFixedParameters()[:-1], 'float')
P2 = IT.GetParameters()
R2 = RotationMatrix(-P2[0], -P2[1], -P2[2])
T2 = -np.array(P2[3:])

R2 = np.linalg.inv(R2)

ProcessTiming(1,'Transform points ')

TransformedPoints = []

for iP, Point in enumerate(TP3):

    TP = np.dot(R2, Point - C2) + C2 - T2
    TransformedPoints.append(TP)

    Progress = iP / len(Points) * 20
    ProgressNext(Progress)
ProcessTiming(0)
TP2 = np.array(TransformedPoints)

# Transform Image
# Transformed2 = sitk.Resample(Transformed1, IT)


#%% Transform 1
# Transform Point
I = sitk.ReadImage(FileNames['Common'])
Center = np.array(I.GetSize()) / 2 * np.array(I.GetSpacing())
C1 = Center + np.array(I.GetOrigin())
R1 = np.array([[-1, 0, 0],[0, 1, 0],[0, 0, -1]])
T1 = [0, 0, 0]

R1 = np.linalg.inv(R1)

ProcessTiming(1,'Transform points ')

TransformedPoints = []

for iP, Point in enumerate(TP2):
    
    TP = np.dot(R1, Point - C1) + C1 - T1
    TransformedPoints.append(TP)

    Progress = iP / len(Points) * 20
    ProgressNext(Progress)
ProcessTiming(0)
TP1 = np.array(TransformedPoints)

# Transform Image
# Rotation = sitk.VersorRigid3DTransform()
# M = RotationMatrix(0, sp.pi, 0)
# M = [v for v in M.flatten()]

# Rotation.SetMatrix(M)
# Center = np.array(Resampled.GetSize()) / 2 * np.array(Resampled.GetSpacing())
# CO = Center + np.array(Resampled.GetOrigin())
# Rotation.SetCenter(CO)
# Transformed1 = sitk.Resample(Resampled, Rotation)


#%% Rebuild image and write MHD

T_Coords = np.round(TP1 / np.array(Transformed2.GetSpacing())).astype('int')

# Compute padding
MinP = np.abs([min(v,0) for v in np.min(T_Coords,axis=0)])
MaxP = [1 + max(v, 0) for v in (np.max(T_Coords, axis=0) - Array.shape[::-1])]

P_Array = np.zeros(Array.shape)
P_Array = np.pad(P_Array,((MinP[2], MaxP[2]), (MinP[1], MaxP[1]), (MinP[0], MaxP[0])))
P_Array[T_Coords[:,2] + MinP[2], T_Coords[:,1] + MinP[1], T_Coords[:,0] + MinP[0]] = 1

# P_Array = P_Array[MinP[2]:-MaxP[2],
#                   MinP[1]:-MaxP[1],
#                   MinP[0]:-MaxP[0]]

T_Image = sitk.GetImageFromArray(P_Array)
T_Image.SetSpacing(Transformed2.GetSpacing())

NewOrigin = np.array(Image.GetOrigin()) + MinP
T_Image.SetOrigin(NewOrigin)


MHD(T_Image, 'Point1')
# MHD(Transformed3, 'Transformed3')


#%% Eigen vector rotation
# Eigenvector

e1 = np.array([1,0,0])
e2 = np.array([0,1,0])
e3 = np.array([0,0,1])

R1 = RotationMatrix(Gamma=sp.pi/4)
R2 = RotationMatrix(Beta=-sp.pi/4)

Re1 = np.dot(R2, np.dot(R1, e1))
Re2 = np.dot(R2, np.dot(R1, e2))
Re3 = np.dot(R2, np.dot(R1, e3))

R1_inv = np.linalg.inv(R1)
R2_inv = np.linalg.inv(R2)

Re1_inv = np.dot(R1_inv, np.dot(R2_inv, Re1))
Re2_inv = np.dot(R1_inv, np.dot(R2_inv, Re2))
Re3_inv = np.dot(R1_inv, np.dot(R2_inv, Re3))

Figure = plt.figure(figsize=(5.5, 4))
Axis = Figure.add_subplot(111, projection='3d')
Axis.quiver(0, 0, 0, e1[0], e1[1], e1[2], color=(1,0,0,0.5))
Axis.quiver(0, 0, 0, e2[0], e2[1], e2[2], color=(0,1,0,0.5))
Axis.quiver(0, 0, 0, e3[0], e3[1], e3[2], color=(0,0,1,0.5))
Axis.quiver(0, 0, 0, Re1[0], Re1[1], Re1[2], color=(1,0,0))
Axis.quiver(0, 0, 0, Re2[0], Re2[1], Re2[2], color=(0,1,0))
Axis.quiver(0, 0, 0, Re3[0], Re3[1], Re3[2], color=(0,0,1))
Axis.quiver(0, 0, 0, Re1_inv[0], Re1_inv[1], Re1_inv[2], color=(1,0,1))
Axis.quiver(0, 0, 0, Re2_inv[0], Re2_inv[1], Re2_inv[2], color=(1,1,0))
Axis.quiver(0, 0, 0, Re3_inv[0], Re3_inv[1], Re3_inv[2], color=(0,1,1))

# scaling hack
Bbox_min = np.min([e1, e2, e3])
Bbox_max = np.max([e1, e2, e3])
Axis.auto_scale_xyz([Bbox_min, Bbox_max], [Bbox_min, Bbox_max], [Bbox_min, Bbox_max])
# make the panes transparent
Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
# modify ticks
MinX, MaxX = -1, 1
MinY, MaxY = -1, 1
MinZ, MaxZ = -1, 1
Axis.set_xticks([MinX, 0, MaxX])
Axis.set_yticks([MinY, 0, MaxY])
Axis.set_zticks([MinZ, 0, MaxZ])
Axis.xaxis.set_ticklabels([MinX, 0, MaxX])
Axis.yaxis.set_ticklabels([MinY, 0, MaxY])
Axis.zaxis.set_ticklabels([MinZ, 0, MaxZ])

Axis.set_xlabel('X')
Axis.set_ylabel('Y')
Axis.set_zlabel('Z')
plt.show()

# %%
