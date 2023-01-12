#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    This script aims to provide utility functions to
    manipulate image (reading, ploting, registering,
    and writing) or signal processing. Most functions
    are taken from different sources and adapted here
    
    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel
            University of Bern

    Date: November 2022
    """

#%% Imports
# Modules import

import os
import vtk
import sys
import time
import struct
import argparse
import numpy as np
import sympy as sp
import pandas as pd
import SimpleITK as sitk
from pathlib import Path
import scipy.signal as sig
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%% Tuning
# Tune diplay settings

DWidth = 320 # display width in number of character
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', DWidth)
pd.set_option('display.width', DWidth)
np.set_printoptions(linewidth=DWidth,suppress=True,formatter={'float_kind':'{:3}'.format})

plt.rc('font', size=12) # increase slightly plot font size for readability


#%% Functions
# Define some general functions

def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results

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

def GetSlice(Image, Slice=None, Axis='Z', Slice2D=True):

    # Get image attributes
    Size = Image.GetSize()
    Spacing = np.array(Image.GetSpacing())
    Origin = np.array(Image.GetOrigin())
    Direction = np.array(Image.GetDirection()).reshape((3,3))

    if type(Slice) == list:
        Start, Stop = Slice

        if Axis == 'Z':
            Sliced = sitk.Slice(Image, (0, 0, Start), (Size[0], Size[1], Stop))
        elif Axis == 'Y':
            Sliced = sitk.Slice(Image, (0, Start, 0), (Size[0], Stop, Size[2]))
        elif Axis == 'X':
            Sliced = sitk.Slice(Image, (Start, 0, 0), (Stop, Size[1], Size[2]))

    elif type(Slice) == int:

        if Axis == 'Z':
            Sliced = sitk.Slice(Image, (0, 0, Slice), (Size[0], Size[1], Slice+1))

            if Slice2D:
                Array = sitk.GetArrayFromImage(Sliced)[0]
                Sliced2D = sitk.GetImageFromArray(Array)
                Sliced2D.SetOrigin(Origin[:2])
                Sliced2D.SetSpacing(Spacing[:2])
                Sliced2D.SetDirection(Direction[:2,:2].flatten())

        elif Axis == 'Y':
            Sliced = sitk.Slice(Image, (0, Slice, 0), (Size[0], Slice+1, Size[2]))

            if Slice2D:
                Array = sitk.GetArrayFromImage(Sliced)[1]
                Sliced2D = sitk.GetImageFromArray(Array)
                Sliced2D.SetOrigin(Origin[::2])
                Sliced2D.SetSpacing(Spacing[::2])
                Sliced2D.SetDirection(Direction[::2,::2].flatten())

        elif Axis == 'X':
            Sliced = sitk.Slice(Image, (Slice, 0, 0), (Slice+1, Size[1], Size[2]))

            if Slice2D:
                Array = sitk.GetArrayFromImage(Sliced)[2]
                Sliced2D = sitk.GetImageFromArray(Array)
                Sliced2D.SetOrigin(Origin[1:])
                Sliced2D.SetSpacing(Spacing[1:])
                Sliced2D.SetDirection(Direction[1:,1:].flatten())

    else:

        if Axis == 'Z':
            Sliced = sitk.Slice(Image, (0, 0, int(Size[2]/2)), (Size[0], Size[1], int(Size[2]/2)+1))
        
            if Slice2D:
                Array = sitk.GetArrayFromImage(Sliced)[2]
                Sliced2D = sitk.GetImageFromArray(Array)
                Sliced2D.SetOrigin(Origin[:2])
                Sliced2D.SetSpacing(Spacing[:2])
                Sliced2D.SetDirection(Direction[:2,:2].flatten())


        elif Axis == 'Y':
            Sliced = sitk.Slice(Image, (0, int(Size[1]/2), 0), (Size[0], int(Size[1]/2)+1, Size[2]))
        
            if Slice2D:
                Array = sitk.GetArrayFromImage(Sliced)[1]
                Sliced2D = sitk.GetImageFromArray(Array)
                Sliced2D.SetOrigin(Origin[::2])
                Sliced2D.SetSpacing(Spacing[::2])
                Sliced2D.SetDirection(Direction[::2,::2].flatten())


        elif Axis == 'X':
            Sliced = sitk.Slice(Image, (int(Size[0]/2), 0, 0), (int(Size[0]/2)+1, Size[1], Size[2]))

            if Slice2D:
                Array = sitk.GetArrayFromImage(Sliced)[2]
                Sliced2D = sitk.GetImageFromArray(Array)
                Sliced2D.SetOrigin(Origin[1:])
                Sliced2D.SetSpacing(Spacing[1:])
                Sliced2D.SetDirection(Direction[1:,1:].flatten())


    if Sliced2D:
        Sliced = Sliced2D
        
    return Sliced

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

def Resample(Image, Factor=None, Size=[None], Spacing=[None]):

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
    
    elif Size[0]:
        NewSize = Size
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Spacing[0]:
        NewSpacing = Spacing
        NewSize = [round(Size/Spacing) + 1 for Size,Spacing in zip(PhysicalSize, NewSpacing)]
    
    NewArray = np.zeros(NewSize[::-1],'int')
    NewImage = sitk.GetImageFromArray(NewArray)
    NewImage.SetOrigin(Origin)
    NewImage.SetDirection(Direction)
    NewImage.SetSpacing(NewSpacing)
  
    Transform = sitk.TranslationTransform(Dimension)
    
    return sitk.Resample(Image, NewImage, Transform, sitk.sitkNearestNeighbor)

def ProgressStart(Text):
    global CurrentProgress
    sys.stdout.write(Text + '|')
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


#%% Ploting functions
class Show:

    def __init__(self):
        self.FName = None
        self.ShowPlot = True
        self.IRange = [0.8, 1.2]

    def Normalize(self, Array):

        Min = np.min(Array)
        Max = np.max(Array)
        N_Array = (Array - Min) / (Max - Min) * 255

        return np.array(N_Array).astype('uint8') 

    def Slice(self, Image, Slice=None, Title=None, Axis='Z'):

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
        
        if (Title):
            Axis.set_title(Title)

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0)

        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

        return

    def Overlay(self, Fixed, Moving, Slice=None, Title=None, Axis='Z', AsBinary=False):

        FixedArray = sitk.GetArrayFromImage(Fixed)
        MovingArray = sitk.GetArrayFromImage(Moving)

        if AsBinary:
            Otsu = sitk.OtsuMultipleThresholdsImageFilter()
            Otsu.SetNumberOfThresholds(2)
            
            if len(np.unique(FixedArray)) > 2:
                Fixed_Bin = Otsu.Execute(Fixed)
                FixedArray = sitk.GetArrayFromImage(Fixed_Bin)
                FixedArray = (FixedArray == 2) * 1
            
            if len(np.unique(MovingArray)) > 2:
                Moving_Bin = Otsu.Execute(Moving)
                MovingArray = sitk.GetArrayFromImage(Moving_Bin)
                MovingArray = (MovingArray == 2) * 1

        FixedArray = self.Normalize(FixedArray)
        MovingArray = self.Normalize(MovingArray)


        if Fixed.GetDimension() == 3:
            Array = np.zeros((Fixed.GetSize()[2], Fixed.GetSize()[1], Fixed.GetSize()[0], 3), 'uint8')
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
            Array = np.zeros((Fixed.GetSize()[1], Fixed.GetSize()[0], 3), 'uint8')
            Array[:,:,0] = FixedArray
            Array[:,:,1] = MovingArray
            Array[:,:,2] = MovingArray

        Figure, Axis = plt.subplots()
        Axis.imshow(Array,interpolation=None)
        Axis.axis('Off')
        
        if (Title):
            Axis.set_title(Title)

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0)

        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

        return

    def Intensity(self, Structure, Deformations, Mask=None, Slice=None, Axis='Z', Title=None):

        Array = sitk.GetArrayFromImage(Structure)
        Values = sitk.GetArrayFromImage(Deformations)

        if Mask:
            MaskArray = sitk.GetArrayFromImage(Mask).astype('bool')

        if Structure.GetDimension() == 3:
            
            if Axis == 'Z':
                if Slice:
                    Array = Array[Slice,:,:]
                    Values = Values[Slice,:,:]
                    if Mask:
                        MaskArray = MaskArray[Slice,:,:]
                else:
                    Array = Array[Array.shape[0]//2,:,:]
                    Values = Values[Values.shape[0]//2,:,:]
                    if Mask:
                        MaskArray = MaskArray[MaskArray.shape[0]//2,:,:]

            if Axis == 'Y':
                if Slice:
                    Array = Array[:,Slice,:]
                    Values = Values[:,Slice,:]
                    if Mask:
                        MaskArray = MaskArray[:,Slice,:]
                else:
                    Array = Array[:,Array.shape[1]//2,:]
                    Values = Values[:,Values.shape[1]//2,:]
                    if Mask:
                        MaskArray = MaskArray[:,MaskArray.shape[1]//2,:]

            if Axis == 'X':
                if Slice:
                    Array = Array[:,:,Slice]
                    Values = Values[:,:,Slice]
                    if Mask:
                        MaskArray = MaskArray[:,:,Slice]
                else:
                    Array = Array[:,:,Array.shape[2]//2]
                    Values = Values[:,:,Values.shape[2]//2]
                    if Mask:
                        MaskArray = MaskArray[:,:,MaskArray.shape[2]//2]

        Structure = np.zeros((Array.shape[0], Array.shape[1], 4))
        Structure[:,:,3] = self.Normalize(Array) / 255
        
        if Mask:
            Values[~MaskArray] = np.nan
        else:
            Values[Values == 0] = np.nan

        Figure, Axis = plt.subplots(1,1)
        Plot = Axis.imshow(Values, cmap='jet', vmin=self.IRange[0], vmax=self.IRange[1], interpolation=None)
        Axis.imshow(Structure)
        Axis.axis('Off')

        # Colorbar hack
        Divider = make_axes_locatable(Axis)
        CAxis = Divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(Plot, cax=CAxis, orientation='vertical')

        if (Title):
            Axis.set_title(Title)

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0)

        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

        return

    def Signal(self, X, Y, Points=[], Normalize=False, Axes=[], Labels=[]):

        Colors = [(1,0,0), (0,0,1), (0,0,0), (0,1,0), (0,1,1), (1,0,1)]

        Figure, Axis = plt.subplots(1,1)

        for i in range(len(X)):

            if Normalize:
                Xi = self.Normalize(X[i])
                Yi = self.Normalize(Y[i])
            else:
                Xi, Yi = X[i], Y[i]
            
            if len(Labels) > 0:
                Axis.plot(Xi, Yi, color=Colors[i], label=Labels[i])
            else:
                Axis.plot(Xi, Yi, color=Colors[i], label='Signal ' + str(i+1))

            if len(Points) > 0:
                Px, Py = Xi[Points[i]], Yi[Points[i]]
                Axis.plot(Px, Py, marker='o', color=(0.7,0.7,0.7), fillstyle='none', linestyle='none')

        if len(Axes) > 0:
            Axis.set_xlabel(Axes[0])
            Axis.set_ylabel(Axes[1])

        Cols = i+1 if i < 3 else (1+i)//2
        plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.12), ncol=Cols)

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0.02)

        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

        return

Show = Show()
#%% Reading functions
class Read:

    def __init__(self):
        self.Echo = True
    
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

        if self.Echo:
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
            # only if an isotropic scaling was applied!!!
            if Line.find('downscaled'.encode()) > -1:
                pass
            elif Line.find("scale".encode()) > -1:
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

        # Y symmetry (thanks Michi for notifying this!)
        Numpy_Image = Numpy_Image[:,::-1,:]
        
        # Converty numpy to ITK image
        Image = sitk.GetImageFromArray(Numpy_Image)
        Image.SetSpacing(Spacing)
        Image.SetOrigin([0.0, 0.0, 0.0])

        AdditionalData = {'Scaling':Scaling,
                        'Slope':Slope,
                        'Intercept':Intercept,
                        'Header':Header}

        if self.Echo:
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

Read = Read()
#%% Writing functions
class Write:

    def __init__(self):
        pass

    def Raw(self, Image, OutputFileName, PixelType):

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

    def MHD(self, Image, FileName, PixelType='uint'):

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
        Offset = str(np.array(Image.GetOrigin()))[1:-1]
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

        outs.write('ElementDataFile = %s\n' % (FileName.split('/')[-1] + '.raw'))
        outs.close()

        self.Raw(Image, FileName + '.raw', PixelType)

        Toc = time.time()
        PrintTime(Tic, Toc)

        return

    def VTK(self, VectorField,FilePath,FileName,SubSampling=1):

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

Write = Write()
#%% Registration funtions
class Registration:

    def __init__(self):
        self.Echo = True

    def Register(self, FixedImage, MovingImage, Type, FixedMask=None, MovingMask=None, Path=None, Dictionary={}):

        if self.Echo:
            print('\nPerform ' + Type + ' registration')
            Tic = time.time()
        PM = sitk.GetDefaultParameterMap(Type)

        # Set standard parameters if not specified otherwise
        if 'MaximumNumberOfIterations' not in Dictionary.keys():
            PM['MaximumNumberOfIterations'] = ['2000']

        if 'FixedImagePyramidSchedule' not in Dictionary.keys():
            PM['FixedImagePyramidSchedule'] = ['50', '20', '10']

        if 'MovingImagePyramidSchedule' not in Dictionary.keys():
            PM['MovingImagePyramidSchedule'] = ['50', '20', '10']

        if 'SP_alpha' not in Dictionary.keys():
            PM['SP_alpha'] = ['0.6']

        if 'SP_A' not in Dictionary.keys():
            PM['SP_A'] = ['1000']

        # Set other defined parameters
        for Key in Dictionary.keys():
            PM[Key] = [str(Item) for Item in Dictionary[Key]]


        # Set Elastix and perform registration
        EIF = sitk.ElastixImageFilter()
        EIF.SetParameterMap(PM)
        EIF.SetFixedImage(FixedImage)
        EIF.SetMovingImage(MovingImage)

        if FixedMask:
            FixedMask = sitk.Cast(FixedMask, sitk.sitkUInt8)
            EIF.SetFixedMask(FixedMask)

        if MovingMask:
            MovingMask = sitk.Cast(MovingMask, sitk.sitkUInt8)
            EIF.SetMovingMask(MovingMask)

        if Path:
            EIF.SetOutputDirectory(Path)
            EIF.LogToConsoleOff()
            EIF.LogToFileOn()

        if os.name == 'posix':
            EIF.SetNumberOfThreads(8)

        EIF.Execute()

        # Get results
        Result_Image = EIF.GetResultImage()
        TransformParameters = EIF.GetTransformParameterMap()

        # Print elapsed time
        if self.Echo:
            Toc = time.time()
            PrintTime(Tic, Toc)

        return Result_Image, TransformParameters
        
    def ComputeInverse(self, FixedImage, TPMFileName, MovingImage=None, Path=None):

        """
        Compute inverse of rigid elastix transform. Manual 6.1.6
        """

        if self.Echo:
            print('\nCompute registration inverse transform')
            Tic = time.time()

        # Set Elastix and perform registration
        EF = sitk.ElastixImageFilter()
        EF.SetFixedImage(FixedImage)
        EF.SetMovingImage(FixedImage)
        EF.SetInitialTransformParameterFileName(TPMFileName)

        EF.SetParameter('HowToCombineTransforms','Compose')
        EF.SetParameter('MaximumNumberOfIteration','2000')
        EF.SetParameter('FixedImagePyramidSchedule', ['50', '20', '10'])
        EF.SetParameter('MovingImagePyramidSchedule', ['50', '20', '10'])
        EF.SetParameter('SP_alpha', '0.6')
        EF.SetParameter('SP_A', '1000')

        if MovingImage:
            EF.SetParameter('Size', '%i %i %i'%MovingImage.GetSize())
            EF.SetParameter('Spacing', '%f %f %f'%MovingImage.GetSpacing())
            EF.SetParameter('Origin', '%f %f %f'%MovingImage.GetOrigin())
        
        if Path:
            EF.SetOutputDirectory(Path)
            EF.LogToConsoleOff()
            EF.LogToFileOn()

        EF.Execute()
        InvertedTransform = EF.GetTransformParameterMap()[0]
        del InvertedTransform['InitialTransformParametersFileName']

        # Print elapsed time
        if self.Echo:
            Toc = time.time()
            PrintTime(Tic, Toc)

        return InvertedTransform

    def Apply(self, Image,TransformParameterMap,Path=None,Jacobian=None):

        """
        Apply transform parametermap from elastix to an image
        """

        if self.Echo:
            print('\nApply transform using Transformix')
            Tic = time.time()

        TIF = sitk.TransformixImageFilter()
        TIF.ComputeDeterminantOfSpatialJacobianOff()
        TIF.SetTransformParameterMap(TransformParameterMap)

        if Jacobian:
            TIF.ComputeDeformationFieldOn()
            TIF.ComputeSpatialJacobianOn()

        else:
            TIF.ComputeDeformationFieldOff()
            TIF.ComputeSpatialJacobianOff()


        if Path:
            TIF.SetOutputDirectory(Path)

        TIF.SetMovingImage(Image)
        TIF.Execute()

        ResultImage = TIF.GetResultImage()

        ResultImage.SetOrigin(np.array(TransformParameterMap[0]['Origin'], float))
        ResultImage.SetSpacing(np.array(TransformParameterMap[0]['Spacing'], float))

        # Print elapsed time
        if self.Echo:
            Toc = time.time()
            PrintTime(Tic, Toc)

        return ResultImage
 
    def ApplyCustom(self, Image, TransformParameterMap):

        """
        Apply costum parameter map to an image
        """

        # Get transform parameters
        D = np.array(TransformParameterMap['FixedImageDimension'], 'int')[0]
        TP = np.array(TransformParameterMap['TransformParameters'],'float')
        COR = np.array(TransformParameterMap['CenterOfRotationPoint'], 'float')

        # Apply rotation
        R = sitk.VersorRigid3DTransform()
        RM = RotationMatrix(Alpha=TP[0], Beta=TP[1], Gamma=TP[2])
        R.SetMatrix([Value for Value in RM.flatten()])
        R.SetCenter(COR)
        Image_R = sitk.Resample(Image, R)

        # Apply translation
        T = sitk.TranslationTransform(int(D), -TP[3:])
        Image_T = sitk.Resample(Image_R, T)

        return Image_T
           
    def ApplyInverse(self, Image, TransformParameterMap):

        """
        Apply inverse rigid transform from transform parameter map
        """

        # Get transform parameters
        D = np.array(TransformParameterMap['FixedImageDimension'], 'int')[0]
        TP = np.array(TransformParameterMap['TransformParameters'],'float')
        COR = np.array(TransformParameterMap['CenterOfRotationPoint'], 'float')

        # Apply inverse translation
        T = sitk.TranslationTransform(int(D), -COR - TP[3:])
        Image_T = sitk.Resample(Image, T)

        # Apply inverse rotation
        RM = RotationMatrix(Alpha=TP[0], Beta=TP[1], Gamma=TP[2])
        RM_Inv = np.linalg.inv(RM)
        R = sitk.VersorRigid3DTransform()
        R.SetMatrix([Value for Value in RM_Inv.flatten()])
        R.SetCenter([0, 0, 0])

        Image_R = sitk.Resample(Image_T, R)

        # Translate again for center of rotation
        T = sitk.TranslationTransform(int(D), COR)
        Image_T = sitk.Resample(Image_R, T)

        return Image_T

Registration = Registration()
#%% Signal treatment functions
class Signal:

    def FFT(Signal, Sampling, Show=False):

        """
        Analyze signal spectrum in frequency domain
        
        :param Signal: the signal to analyze
        :param Sampling: signal sampling interval (in /s or /m)
        :param Show: Plot the frequency spectrum
        """

        SamplingFrequency = 1 / Sampling
        NormalizedSpectrum = np.fft.fft(Signal) / len(Signal)
        Frequencies = np.fft.fftfreq(Signal.shape[-1], SamplingFrequency)

        RealHalfSpectrum = np.abs(NormalizedSpectrum.real[Frequencies >= 0])
        HalfFrequencies = Frequencies[Frequencies >= 0]

        if Show:
            Figure, Axis = plt.subplots(1,1)
            Axis.semilogx(HalfFrequencies, RealHalfSpectrum, color=(1,0,0))
            Axis.set_xlabel('Frequencies [Hz]')
            Axis.set_ylabel('Amplitude [-]')
            plt.show()

        return HalfFrequencies, RealHalfSpectrum

    def DesignFilter(Frequency, Order=2):

        """
        Design Butterworth filter according to cut-off frequency and order

        :param Frequency: cut-off frequency of the filter
        :param Order: order of the filter
        """
        
        b, a = sig.butter(Order, Frequency, 'low', analog=True)
        w, h = sig.freqs(b, a)

        Figure, Axis = plt.subplots(1,1)
        Axis.semilogx(w, 20 * np.log10(abs(h)))
        Axis.set_xlabel('Frequency [radians / second]')
        Axis.set_ylabel('Amplitude [dB]')
        Axis.grid(which='both', axis='both')
        Axis.axvline(Frequency, color='green') # cutoff frequency
        Axis.set_ylim([-50,5])
        plt.show()

        return

    def Filter(Signal, Sampling, Frequency, Order=2, Show=False):
        
        """
        Filter signal and look filtering effect

        :param Signal: signal to filter
        :param Sampling: signal sampling interval (in /s or /m)
        :param Frequency: cut-off frequency
        :param Order: filter order
        :param Show: plot results
        """

        SOS = sig.butter(Order, Frequency / Sampling, output='sos')
        FilteredSignal = sig.sosfiltfilt(SOS, Signal)

        if Show:
            Figure, Axis = plt.subplots(1,1)
            Axis.plot(Signal, color=(0,0,0))
            Axis.plot(FilteredSignal, color=(1,0,0))
            plt.show()

        return FilteredSignal

#%% Abaqus functions
class Abaqus:

    def __init__(self):
        pass

    def ReadDAT(self, File):

        try:
            with open(File) as F:
                Text = F.read()

                Values = []
                Condition = Text.find('U1') + 1
                Start = Text.find('U1')
                Steps, Increments = [], []
                Step, Increment = 1, 1
                while Condition:

                    Inc = int(Text[Start-349:Start-347])
                    if Inc < Increment:
                        Step += 1
                    Increment = Inc
                    Increments.append(Increment)
                    Steps.append(Step)

                    for i in range(6):
                        iStart = Start + 125 + 15*i
                        iStop = iStart + 14
                        Values.append(float(Text[iStart : iStop]))

                        iStart += Text[Start:].find('RF1')
                        iStop += Text[Start:].find('RF1')
                        Values.append(float(Text[iStart : iStop]))

                    Start = iStop + Text[iStop:].find('U1')
                    Condition = Text[iStop:].find('U1') + 1

                Values = np.array(Values)
                Cols = 12
                Rows = Values.size // Cols
                Values = np.reshape(Values,(Rows,Cols))

                Data = pd.DataFrame()
                ColNames = []
                for i in range(3):
                    for V in ['U', 'F']:
                        ColNames.append(V + str(i+1))
                for i in range(3):
                    for V in ['R', 'M']:
                        ColNames.append(V + str(i+1))

                for iName, Name in enumerate(ColNames):
                    Data[Name] = Values[:,iName]
                
                Data.columns = ['X', 'FX', 'Y', 'FY', 'Z', 'FZ', 'Phi', 'MX', 'Theta', 'MY', 'Psi', 'MZ']

                Data['Step'] = Steps
                Data['Increment'] = Increments

            return Data

        except FileNotFoundError:
            print('File' + File + 'does not exist')

            return

Abaqus = Abaqus()
#%%
if __name__ == '__main__':

    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)

    # Add required arguments
    Parser.add_argument('File', help='File to process', type=str)

    # Add long and short optional arguments
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-v', '--Version', help='Show script version', action='version', version=SV)
    # Read arguments from the command line
    Arguments = Parser.parse_args()
    
    # This part to finish
    # Do something