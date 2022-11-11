#%%
#! /usr/bin/env python3

import os
import vtk
import time
import struct
import argparse
import numpy as np
from pathlib import Path
import SimpleITK as sitk
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy

#%%
Version = '01'

Description = """
    This script aims to provide utility functions
    to read and write different image types.
    Most functions are taken from different sources
    and adapted here
    
    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel
            University of Bern

    Date: November 2022
    """

#%%
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

def Get_AIM_Ints(File):

    """
    Function by Glen L. Niebur, University of Notre Dame (2010)
    reads the integer data of an AIM file to find its header length
    """

    nheaderints = 32
    File.seek(0)
    binints = File.read(nheaderints * 4)
    header_int = struct.unpack("=32i", binints)

    return header_int

def ReadAIM(File):

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
    print('\n\nRead AIM header of file: ' + File)
    with open(File, 'rb') as f:
        AIM_Ints = Get_AIM_Ints(f)
        # check AIM version
        if int(AIM_Ints[5]) == 16:
            print("     -> version 020")
            if int(AIM_Ints[10]) == 131074:
                Format = "short"
                print("     -> format " + Format)
            elif int(AIM_Ints[10]) == 65537:
                Format = "char"
                print("     -> format " + Format)
            elif int(AIM_Ints[10]) == 1376257:
                Format = "bin compressed"
                print("     -> format " + Format + " not supported! Exiting!")
                exit(1)
            else:
                Format = "unknown"
                print("     -> format " + Format + "! Exiting!")
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

def ReadISQ(File, BMD=False, Echo=False, Info=False):

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

    if Echo:
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

    if os.name == 'posix':
        File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape, order='C')
    else:
        File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape, order='F')

    File[:] = CastedImageArray[:]
    del File

    return

def WriteMHD(Image, Path, FileName, PixelType='uint'):

    print('\nWrite MHD')
    Tic = time.time()

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

    Toc = time.time()
    PrintTime(Tic, Toc)

    return

def WriteVTK(VectorField,FilePath,FileName,SubSampling=1):

    print('\nWrite VTK')
    Tic = time.time()

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

    Toc = time.time()
    PrintTime(Tic, Toc)

    return

#%%
if __name__ == '__main__':

    # Initiate the parser with a description
    Parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add required arguments
    Parser.add_argument('File', help='File to read (required)', type=str)

    # Add long and short optional arguments
    ScriptVersion = Parser.prog + ' version ' + Version
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=ScriptVersion)
    Parser.add_argument('--BMD', default=False, help='Convert gray values to BMD (bool), !!! Depends on voltage, current and time !!!', type=bool)
    Parser.add_argument('--Echo', default=False, help='Print out current operation and results (bool)', type=bool)
    Parser.add_argument('--Info', default=False, help='Write file info into text file (bool)', type=bool)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    File = Arguments.File

    try:
        BMD = Arguments.BMD
    except:
        BMD = False
    
    try:
        Echo = Arguments.Echo
    except:
        Echo = False

    try:
        Info = Arguments.Info
    except:
        Info = False

    if File.endswith('.ISQ'):
        ReadISQ(File, BMD, Echo, Info)

    elif File.endswith('.AIM'):
        ReadAIM(File)