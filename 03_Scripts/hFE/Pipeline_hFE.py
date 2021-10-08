"""
This script runs ACCURATE pipeline without PSL, converted from Denis's pipeline.
Meant to be run locally from the FRACTIB root folder

Author: Mathieu Simon, ARTORG Center for Biomedical Engineering Research, SITEM Insel, University of Bern
Date: October 2021
"""

import os
import yaml
import resource
import struct
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
import sys


def ReadConfigFile(Filename):

    """ Read configuration file and store to dictionary """

    print('\n\nReading initialization file', Filename)
    with open(Filename, 'r') as File:
        Configuration = yaml.load(File, Loader=yaml.FullLoader)

    return Configuration
def Set_FileNames(Config, Sample):

    """
    Adapted from Denis's io_utils_SA.py
    Set filenames for each grayscale file
    Filenames depend on pipeline (fast/accurate)
    Always:
        - Native image for header
        - BMD or Native image cropped to ROI
    Additional for fast pipeline:
        - Periosteal mask
    Additional for accurate pipeline:
        - Trabecular mask
        - Cortical mask
        - Two-phase segmentation
    """

    # Always, not depending on phase (accurate or fast)
    Version = Config['Version']
    Folder_IDs = Config['Folder_IDs']
    Folder = Folder_IDs[Sample]
    AIMDir = os.path.join(Config['AIMDir'], Folder)

    FileName = {}
    FileName['BCs'] = Config['BCs']

    # Additional for accurate pipeline
    Postfix_CortMask = Config['Postfix_CortMask']
    Postfix_TrabMask = Config['Postfix_TrabMask']
    Postfix_BMD = Config['Postfix_BMD']
    Postfix_SEG = Config['Postfix_SEG']

    FileName['FILEMASKCORT'] = "{}{}".format(Sample, Postfix_CortMask)
    FileName['FILEMASKTRAB'] = "{}{}".format(Sample, Postfix_TrabMask)
    FileName['FILEBMD'] = "{}{}".format(Sample, Postfix_BMD)
    FileName['FILESEG'] = "{}{}".format(Sample, Postfix_SEG)

    # Always, not depending on phase (accurate or fast)
    FileName['FILEGRAY'] = FileName['FILEBMD']
    FileName['RAWname'] = os.path.join(AIMDir, FileName['FILEGRAY'])
    FileName['BMDname'] = os.path.join(AIMDir, FileName['FILEBMD'])

    FileName['CORTMASKname'] = os.path.join(AIMDir, FileName['FILEMASKCORT'])
    FileName['TRABMASKname'] = os.path.join(AIMDir, FileName['FILEMASKTRAB'])
    FileName['SEGname'] = os.path.join(AIMDir, FileName['FILESEG'])

    # General filenames
    print(FileName['BMDname'])
    New_FileName = "{}_{}.inp".format(Sample, Version)
    FileName["INPname"] = os.path.join(AIMDir, New_FileName)

    New_FileName = "{}_{}_summary.txt".format(Sample, Version)
    FileName["SUMname"] = os.path.join(AIMDir, New_FileName)

    New_FileName = "{}_{}_BPVb".format(Sample, Version)
    FileName["VER_BPVname"] = os.path.join(AIMDir, New_FileName)

    return FileName
def Print_Memory_Usage():
    Memory_Used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1e-6
    print("Memory usage: {:.2f} (GB)".format(Memory_Used))
    return
def Get_AIM_Ints(f):

    """
    Function by Glen L. Niebur, University of Notre Dame (2010)
    reads the integer data of an AIM file to find its header length
    """

    nheaderints = 32
    nheaderfloats = 8
    f.seek(0)
    binints = f.read(nheaderints * 4)
    header_int = struct.unpack("=32i", binints)

    return header_int
def AIMReader(File, Spacing):

    """
    Reads an AIM file and provides
    the corresponding vtk image with spacing,
    calibration data and header
    """

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
    # read AIM
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
    return VTK_Image, Spacing, Scaling, Slope, Intercept, Header
def Read_Image_Parameters(FileNames, Bone):

    """
    Read image parameters from AIM image header.
    Input: AIM image (Scanco Medical)
    Output: (bone dictionary)
    - Spacing
    - Scaling
    - Slope
    - Intercept
    """

    print("\n\nRead AIM files")

    VTK_Image, Spacing, Scaling, Slope, Intercept, Header = AIMReader(FileNames["RAWname"], 0)

    Bone['Spacing'] = Spacing
    Bone['Scaling'] = Scaling
    Bone['Slope'] = Slope
    Bone['Intercept'] = Intercept

    return Bone
def VTK2Numpy(VTK_Image):

    """ Turns a vtk image data into a numpy array """

    Dimension = VTK_Image.GetDimensions()
    Data = VTK_Image.GetPointData().GetScalars()
    Numpy_Image = vtk_to_numpy(Data)
    # vtk and numpy have different array conventions
    Numpy_Image = Numpy_Image.reshape(Dimension[2], Dimension[1], Dimension[0])
    Numpy_Image = Numpy_Image.transpose(2, 1, 0)
    return Numpy_Image
def Read_AIM(Name, FileNames, Bone):

    """
    Read AIM image
    Adapted from Denis's io_utils_SA.py
    --------------
    All necessary AIM files are imported and stored in Bone dict
    Input: name specifier, FileNames dict, Bone dict
    Output: Bone dict
    - numpy array containing AIM image
    """

    print("\n\nRead AIM file :" + Name)

    Spacing = Bone["Spacing"]
    # Read image as vtk
    VTK_Image = AIMReader(FileNames[Name + 'name'], Spacing)[0]
    # convert AIM files to numpy arrays
    IMG_Array = VTK2Numpy(VTK_Image)
    if Name == 'SEG':
        IMG_Array[IMG_Array == 127] = 2
        IMG_Array[IMG_Array == 126] = 1
        Bone[Name + "_array"] = IMG_Array
    else:
        Bone[Name + "_array"] = IMG_Array

    return Bone
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
def Adjust_Image(Name, Bone, Config, CropType='Crop'):

    """
    Adapted from Denis's preprocessing_SA.py
    Adjust image size to current FE element size,
    that no layers of the image are removed, size in z direction has to fit.
    If empty layers are added in this dimension,
    this will create a weak layer at the bottom of the image.
    Expansions in x and y dimension,
    will probably not affect strength, but will lower stiffness.

    """

    # get bone values
    IMG_array = Bone[Name + '_array']
    Spacing = Bone['Spacing']

    # coarsening factor = FE element size / CT voxel size
    CoarseFactor = int(round(Config['ElementSize'] / Spacing[0]))
    FEelSize = np.copy(Spacing) * CoarseFactor

    # Adjustment for BMD image and Mask
    IMG_Array_Adjusted = Adjust_Image_Size(IMG_array, CoarseFactor, CropZ=CropType)

    # set bone values
    # copy old IMG_array to IMG_array_original and store new adjusted IMG as IMG_array
    Bone[Name + '_array_original'] = IMG_array
    Bone[Name + '_array'] = IMG_Array_Adjusted
    Bone['FEelSize'] = FEelSize
    Bone['CoarseFactor'] = CoarseFactor

    return Bone
def Write_TestImages(Bone, Config, Sample):

    """
    Write images to test quality of pipeline / FE mesh
    """

    Image_Path = Config['AIMDir'] + '/' + Config['Folder_IDs'][Sample] + '/'

    matplotlib.image.imsave(Image_Path + 'bmd_image.png',
                            Bone['BMD_array'][int(Bone['BMD_array'].shape[0] / 2), :, :])
    matplotlib.image.imsave(Image_Path + 'seg_image.png',
                            Bone['SEG_array'][int(Bone['SEG_array'].shape[0] / 2), :, :])
    matplotlib.image.imsave(Image_Path + 'trabmask_image.png',
                            Bone['TRABMASK_array'][int(Bone['TRABMASK_array'].shape[0] / 2), :, :])
    matplotlib.image.imsave(Image_Path + 'cortmask_image.png',
                            Bone['CORTMASK_array'][int(Bone['CORTMASK_array'].shape[0] / 2), :, :])

    return
def Calculate_BVTV(Bone, Config, ImageType):

    """
    Adapted from Denis's preprocessing_SA.py
    Calculate BVTV and mask images
    ------------------------------
    Scaling, slope and intercept are printed out for review
    If image is already in BMD units, Scaling, Slope and
    Intercept are not applied. BVTVraw is scaled according to
    Hosseini et al. 2017
    (This scaling function could as well be defined externally).

    Parameters
    ----------
    bone    bone results dictionary
    config  configuration parameters dictionary
    IMTYPE  string defining the type of image (BMD/NATIVE)

    Returns
    -------
    bone    bone results dictionary
    """

    # get bone values
    Scaling = Bone["Scaling"]
    Slope = Bone["Slope"]
    Intercept = Bone["Intercept"]
    BMD_array = Bone["BMD_array"]

    print("\n ... prepare mask and BVTV images")
    print("     -> Scaling   = ", Scaling)
    print("     -> Slope     = ", Slope)
    print("     -> Intercept = ", Intercept)

    if ImageType.find("BMD") > -1:
        # if image is already in BMD units (e.g. Hosseinis data)
        BVTVraw = BMD_array / 1200.0
    elif ImageType.find("NATIVE") > -1:
        BMD_array = (BMD_array / Scaling) * Slope + Intercept
        BVTVraw = BMD_array / 1200.0  # if image is in native units

    # Scaling of BVTV 61um to BVTV 11.4um [Hosseini2017]
    Seg_Scaling_Slope = 0.963
    Seg_Scaling_Intercept = 0.03814

    if Config['BVTV_Scaling'] == 1:
        BVTV_Scaled = Seg_Scaling_Slope * BVTVraw + Seg_Scaling_Intercept
    else:
        BVTV_Scaled = BVTVraw

    # Set bone values
    Bone['BVTV_Scaled'] = BVTV_Scaled
    Bone['BMD_Scaled'] = BVTV_Scaled * 1200

    return Bone
def ProgressStart(text):
    global curProgress
    sys.stdout.write(text + '|')
    curProgress = 0
    sys.stdout.flush()
def ProgressNext(progress):
    global curProgress
    if progress > curProgress:
        curProgress += 1
        sys.stdout.write('=')
        sys.stdout.flush()
def ProgressEnd():
    sys.stdout.write('|\n')
    sys.stdout.flush()
def WriteAbaqusGeneral(outFileName, curVoxelModel, dimList):

    """
    Modified function of Medtool
    General Abaqus *.inp file writer. For these materials a default material will be
    applied. Supported commands:
      *USER NODE
      *USER ELEMENT
      *USER NSET, type=point, location=arbitrary
        generate NSET: ARB_NODE_S, ARB_NODE_N, ARB_NODE_E, ARB_NODE_W, ARB_NODE_T, ARB_NODE_B
      *USER NSET, type=point, location=addcorner
        generate NSET: ACOR_NODE_SWB, ACOR_NODE_SEB, ACOR_NODE_NEB, ACOR_NODE_NWB,
                       ACOR_NODE_SWT, ACOR_NODE_SET, ACOR_NODE_NET, ACOR_NODE_NWT
      *USER NSET, type=face
        generate NSET: ALL_NODE_S, ALL_NODE_N, ALL_NODE_E, ALL_NODE_W, ALL_NODE_T, ALL_NODE_B
      *USER ELSET, type=face
        generate ELSET: ALL_S, ALL_ELEM_N, ALL_ELEM_E, ALL_ELEM_W, ALL_ELEM_T, ALL_ELEM_B
      *USER PROPERTY, file=property_temp.inp, range=5:367
        generate multiple material cards, internal variables are "SetName, CardName, GrayValue"
        for the given example: GrayValues > 5 and GrayValues <= 367 are written
        This card can be used multiple times
        If range=... is not given, material cards for all GrayValues are written
    Elements are only written for the given ranges in *USER PROPERTY

    @param outFileName: name of the output file
    @param curVoxelModel: voxel model of the RVE
        - TYPE: np.array[iX, jY, kZ] = grayValue
        - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
        - int grayValue  ... value of voxel
    @param  dimList: list of voxel dimension's
        - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
        - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
    @param  templateFile: name of the template file
    @param  smoothParam: taubin voxel model smoothing parameter
        - TYPE: list[0] = iter, list[1] = lambda, list[2] = kPB,
                list[3] = nearIntf, list[4] = bcid, list[5] = shrink
        - int iter, float lambda, float kPB, int nearIntf

    @return:
      no return value
    """
    sys.stdout.write('\n\nSetup ABAQUS *.inp file from template')
    sys.stdout.write("    -> recast model from '%s' to 'i'" % curVoxelModel.dtype.char)
    sys.stdout.flush()
    curVoxelModel = castType(curVoxelModel, 'i')
    time1 = time.time()
    if dimList.all() == None:  # 12.01.01 change: if dimList == None:     to if dimList.all() == None
        print('\n **ERROR** writeAbaqusGeneral(): Voxel size not optional for this function!\n')
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    xvox = dimList[0]
    yvox = dimList[1]
    zvox = dimList[2]
    nx, ny, nz = get_Shape(curVoxelModel)
    minVox, maxVox = computeNumpyMinMax(curVoxelModel, 0)
    minVox = int(minVox + 0.5)
    maxVox = int(maxVox + 0.5)
    activeNodes = {}
    nodeSets = {}
    nodeSets['ALL_NODE_S'] = []
    nodeSets['ALL_NODE_N'] = []
    nodeSets['ALL_NODE_E'] = []
    nodeSets['ALL_NODE_W'] = []
    nodeSets['ALL_NODE_T'] = []
    nodeSets['ALL_NODE_B'] = []
    elemSets = {}
    elemSets['ALL_ELEM_S'] = []
    elemSets['ALL_ELEM_N'] = []
    elemSets['ALL_ELEM_E'] = []
    elemSets['ALL_ELEM_W'] = []
    elemSets['ALL_ELEM_T'] = []
    elemSets['ALL_ELEM_B'] = []

    tempflag = True
    OS = open('temp.inp', 'w')
    OS.write('*USER NODE\n*USER ELEMENT\n*USER PROPERTY, file=prop.inp, range=1:255\n')
    OS.close()
    OS = open('prop.inp', 'w')
    OS.write('*SOLID SECTION, ELSET=SetName, MATERIAL=CardName\n1.\n')
    OS.write('*MATERIAL,NAME=CardName\n')
    OS.write('*ELASTIC\n')
    OS.write('20000., 0.3\n')
    OS.close()
    templateFile = 'temp.inp'
    OS = open(outFileName, 'w')
    try:
        osTempFile = open(templateFile, 'r')
    except IOError:
        sys.stdout.write(
            "\n **ERROR** mic.writeAbaqusGeneral(): Abaqus Template file '%s' not found!\n\n" % templateFile)
        sys.stdout.flush()
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)

    lines = osTempFile.readlines()
    elsetNodes = {}
    ranges = {}
    thresList = []
    rangeMin = 0
    rangeMax = 255
    outFlag = False
    overlap = np.zeros(rangeMax + 1, np.int)
    for line in lines:
        line = line.replace('\n', '')
        if line.upper().find('*USER PROPERTY') == 0:
            line = line.replace(' ', '')
            args = line.split(',')
            matTemplateFilename = ''
            for arg in args:
                if arg.upper().find('RANGE') == 0:
                    dummy, rangeStr = arg.split('=')
                    rangeMin, rangeMax = userSplit(rangeStr)
                    rangeMin = int(rangeMin)
                    rangeMax = int(rangeMax)
                    if rangeMin < 1:
                        sys.stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): Minimum Range < 1!\n\n')
                        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
                        sys.stdout.flush()
                        exit(1)
                    if rangeMax > maxVox:
                        outFlag = True
                    for ii in range(rangeMax - rangeMin + 1):
                        overlap[rangeMin + ii] += 1

                if arg.upper().find('FILE') == 0:
                    dummy, matTemplateFilename = arg.split('=')

            ranges[matTemplateFilename] = (
                rangeMin, rangeMax)

    if len(ranges) == 0:
        sys.stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): *USER PROPERTY: keyword missing!\n\n')
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    if rangeMax > maxVox:
        sys.stdout.write(
            '\n **WARNING** mic.writeAbaqusGeneral(): *USER PROPERTY: Max GV Range (%i) > Max Image GV (%i)!\n\n' % (
                rangeMax, maxVox))
    if np.sum(np.greater(overlap, 1)) > 0:
        sys.stdout.write(
            '\n **ERROR** mic.writeAbaqusGeneral(): *USER PROPERTY: Ranges in property template overlap!\n\n')
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    for crange in ranges:
        for matId in range(ranges[crange][0], ranges[crange][1] + 1):
            # print('matID', matId)
            elsetNodes[repr(matId)] = []
            thresList.append(matId)

    elid = 0
    nx1 = nx + 1
    nxy1 = (ny + 1) * (nx + 1)
    sum = 0
    ProgressStart('     -> setup Element Data  : ')
    for k in range(nz):
        sum += 1
        progress = float(sum) / float(nz) * 10.0
        for j in range(ny):
            for i in range(nx):
                grayValue = curVoxelModel[k, j, i]
                if repr(grayValue) in elsetNodes:
                    # if elsetNodes.has_key(repr(grayValue)):
                    elid = elid + 1
                    elnds = [nxy1 * k + nx1 * j + (i + 1),
                             nxy1 * k + nx1 * j + (i + 2),
                             nxy1 * k + nx1 * (j + 1) + (i + 2),
                             nxy1 * k + nx1 * (j + 1) + (i + 1),
                             nxy1 * (k + 1) + nx1 * j + (i + 1),
                             nxy1 * (k + 1) + nx1 * j + (i + 2),
                             nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2),
                             nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)]
                    elsetNodes[repr(grayValue)].append((elid, elnds))
                    if k == 0:
                        elemSets['ALL_ELEM_B'].append(elid)
                    if k == nz - 1:
                        elemSets['ALL_ELEM_T'].append(elid)
                    if j == 0:
                        elemSets['ALL_ELEM_S'].append(elid)
                    if j == ny - 1:
                        elemSets['ALL_ELEM_N'].append(elid)
                    if i == 0:
                        elemSets['ALL_ELEM_W'].append(elid)
                    if i == nx - 1:
                        elemSets['ALL_ELEM_E'].append(elid)

        ProgressNext(progress)

    ProgressEnd()
    sys.stdout.write('     -> setup Node Data     :')
    for matid in thresList:
        if len(elsetNodes[repr(matid)]) > 0:
            matidStr = 'SET' + repr(matid)
            for elnds in elsetNodes[repr(matid)]:
                elid = elnds[0]
                for elnd in elnds[1]:
                    activeNodes[elnd] = 1

    noid = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                noid = noid + 1
                if noid in activeNodes:
                    # if activeNodes.has_key(noid):
                    if k == 0:
                        nodeSets['ALL_NODE_B'].append(noid)
                    if k == nz:
                        nodeSets['ALL_NODE_T'].append(noid)
                    if j == 0:
                        nodeSets['ALL_NODE_S'].append(noid)
                    if j == ny:
                        nodeSets['ALL_NODE_N'].append(noid)
                    if i == 0:
                        nodeSets['ALL_NODE_W'].append(noid)
                    if i == nx:
                        nodeSets['ALL_NODE_E'].append(noid)

    sys.stdout.write(' Done\n')
    sys.stdout.flush()
    nodeCoord = {}
    nodeCoordOrig = {}


    noid = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                noid = noid + 1
                if noid in activeNodes:
                    # if activeNodes.has_key(noid):
                    nodeCoord[noid] = (
                        float(xvox * i), float(yvox * j), float(zvox * k))

    curPathFilename, ext = getFilenameAndExtension(outFileName)
    curFilename, ext = getShortFilenameAndExtension(outFileName)
    sys.stdout.write(' ... write ABAQUS *.inp file from template\n')
    for line in lines:
        line = line.replace('\n', '')
        line = line.replace('$filename', curFilename)
        line = line.replace('$pathfilename', curPathFilename)
        if line.upper().find('*USER NODE') > -1:
            OS.write('*NODE\n')
            noid2 = 0
            noid = 0
            ProgressStart('     -> process Node IDs    : ')
            for k in range(nz + 1):
                progress = float(k + 1) / float(nz + 1) * 10.0
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        noid = noid + 1
                        if noid in activeNodes:
                            # if activeNodes.has_key(noid):
                            noid2 = noid2 + 1
                            OS.write('%12i,%13.6g,%13.6g,%13.6g\n' % (
                                noid, nodeCoord[noid][0], nodeCoord[noid][1], nodeCoord[noid][2]))

                ProgressNext(progress)

            ProgressEnd()
            sys.stdout.write('     -> write Nodes         : %10i \n' % noid2)
            sys.stdout.flush()
        elif line.upper().find('*USER ELEMENT') > -1:
            count = 0
            ProgressStart('     -> process Elements    : ')
            for matid in thresList:
                count += 1
                progress = count / float(len(thresList)) * 10.0
                if len(elsetNodes[repr(matid)]) > 0:
                    matidStr = 'SET' + repr(matid)
                    OS.write('*ELEMENT, TYPE=C3D8, ELSET=%s\n' % matidStr)
                    for elnds in elsetNodes[repr(matid)]:
                        elid = elnds[0]
                        OS.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (
                            elid, elnds[1][0], elnds[1][1], elnds[1][2], elnds[1][3], elnds[1][4], elnds[1][5],
                            elnds[1][6], elnds[1][7]))
                        for elnd in elnds[1]:
                            activeNodes[elnd] = 1

                ProgressNext(progress)

            ProgressEnd()
            sys.stdout.write('     -> write Elements      : %10i             \n' % elid)
            sys.stdout.flush()
        elif line.upper().find('*USER NSET') > -1:
            if line.upper().find('TYPE=FACE') > -1:
                sys.stdout.write('     -> write BCs Node Sets     \n')
                sys.stdout.flush()
                for nsetName in nodeSets:
                    OS.write('*NSET, NSET=%s\n' % nsetName)
                    entry = 0
                    for noid in nodeSets[nsetName]:
                        entry = entry + 1
                        if entry == 16:
                            OS.write('%s' % repr(noid))
                            entry = 0
                            OS.write('\n')
                        else:
                            OS.write('%s,' % repr(noid))

                    OS.write('\n')

            if line.upper().find('TYPE=POINT') > -1:
                if line.upper().find('LOCATION=ARBITRARY') > -1:
                    for nsetName in nodeSets:
                        if len(nodeSets[nsetName]) > 0:
                            nid = nodeSets[nsetName][0]
                            name = nsetName.replace('ALL_NODE_', 'ARB_NODE_')
                            OS.write('*NSET, NSET=%s\n' % name)
                            OS.write('%s\n' % repr(nid))

                if line.upper().find('LOCATION=ADDCORNER') > -1:
                    nid = (nx + 1) * (ny + 1) * (nz + 1)
                    OS.write('*NODE, NSET=ACOR_NODE_SWB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 1, 0.0, 0.0, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_SEB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 2, nx * xvox, 0.0, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_NEB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 3, nx * xvox, ny * yvox, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_NWB\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 4, 0.0, ny * yvox, 0.0))
                    OS.write('*NODE, NSET=ACOR_NODE_SWT\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 5, 0.0, 0.0, nz * zvox))
                    OS.write('*NODE, NSET=ACOR_NODE_SET\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 6, nx * xvox, 0.0, nz * zvox))
                    OS.write('*NODE, NSET=ACOR_NODE_NET\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 7, nx * xvox, ny * yvox, nz * zvox))
                    OS.write('*NODE, NSET=ACOR_NODE_NWT\n')
                    OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 8, 0.0, ny * yvox, nz * zvox))
        elif line.upper().find('*USER ELSET') > -1:
            if line.upper().find('TYPE=FACE') > -1:
                sys.stdout.write('     -> Write BCs Elem Sets          \n')
                sys.stdout.flush()
                for elsetName in elemSets:
                    OS.write('*ELSET, ELSET=%s\n' % elsetName)
                    entry = 0
                    for elid in elemSets[elsetName]:
                        entry = entry + 1
                        if entry == 16:
                            OS.write('%s' % repr(elid))
                            entry = 0
                            OS.write('\n')
                        else:
                            OS.write('%s,' % repr(elid))

                    OS.write('\n')

        elif line.upper().find('*USER PROPERTY') > -1:
            line = line.replace(' ', '')
            args = line.split(',')
            rangeMin = minVox
            rangeMax = maxVox
            matTemplateFilename = ''
            for arg in args:
                if arg.upper().find('RANGE') == 0:
                    dummy, rangeStr = arg.split('=')
                    rangeMin, rangeMax = userSplit(rangeStr)
                    rangeMin = int(rangeMin)
                    rangeMax = int(rangeMax)
                if arg.upper().find('FILE') == 0:
                    dummy, matTemplateFilename = arg.split('=')

            sys.stdout.write('     -> Write Property      : %s \n' % matTemplateFilename)
            try:
                osMatCard = open(matTemplateFilename, 'r')
            except IOError:
                sys.stdout.write(
                    "\n **ERROR** writeAbaqusGeneral(): Material template file '%s' not found!\n\n" % matTemplateFilename)
                sys.stdout.flush()
                sys.stdout.write('\n E N D E D  with ERRORS \n\n')
                sys.stdout.flush()
                exit(1)

            lines = osMatCard.readlines()
            for matid in thresList:
                GrayValue = matid
                if len(elsetNodes[repr(matid)]) > 0:
                    if matid >= rangeMin and matid <= rangeMax:
                        matidStr = 'SET' + repr(matid)
                        GrayValue = matid
                        for line in lines:
                            line = line.replace('\n', '')
                            if line.find('SetName') > -1:
                                line = line.replace('SetName', matidStr)
                            if line.find('CardName') > -1:
                                line = line.replace('CardName', 'MAT' + matidStr)
                            if line.find('GrayValue') > -1:
                                exprList = line.split(',')
                                count = 0
                                for expr in exprList:
                                    if expr.find('GrayValue') > -1:
                                        compValue = eval(expr)
                                        OS.write('%s' % repr(compValue))
                                    else:
                                        OS.write('%s' % expr)
                                    if count < len(exprList) - 1:
                                        OS.write(',')
                                    count += 1

                                OS.write('\n')
                            else:
                                OS.write('%s\n' % line)

            osMatCard.close()
        else:
            OS.write('%s\n' % line)

    osTempFile.close()
    if tempflag:
        os.remove('temp.inp')
        os.remove('prop.inp')
    time2 = time.time()
    sys.stdout.write('     -> Write finished in   :   %8.1f sec  \n' % (time2 - time1))
    sys.stdout.flush()
    OS.close
    return
def Generate_Mesh(Bone, Config, FileNames):

    """
    Adapted from Denis's preprocessing_SA.py -> PSL_generate_full_block_mesh_accurate
    Creates Abaqus mesh from coarsened BVTV image and writes MESH to input file (.inp)
    Elements, nodes and element sets are read from the input file and stored in bone.
    Extended to add artificial layers at top and bottom of the image, for reducing
    influences of boundary conditions on homogeneity of strain measures.
    Debugged and checked for right orientations

    Parameters
    ----------
    bone
    config
    filenames

    Returns
    -------
    bone
    """

    # get bone values
    BVTV_Scaled = Bone['BVTV_Scaled']
    CORTMASK_Array = Bone['CORTMASK_array']
    TRABMASK_Array = Bone['TRABMASK_array']
    FEelSize = Bone['FEelSize']
    Spacing = Bone['Spacing']
    CoarseFactor = Bone['CoarseFactor']

    BVTV_Masked = np.copy(BVTV_Scaled)
    MASK_Array = np.add(CORTMASK_Array, TRABMASK_Array)
    BVTV_Masked[MASK_Array == 0] = 0

    # Create array for MESH (no padding)
    MESH = np.ones(([int(dim) for dim in np.rint(np.array(BVTV_Scaled.shape) / CoarseFactor)]))
    MESH = MESH.transpose(2, 1, 0)  # weird numpy array convention (z,y,x)

    Print_Memory_Usage()
    print('Spacing = ' + str(Spacing))
    print('FEelSize = ' + str(FEelSize))
    print(FEelSize[0] / 0.082)

    # Write MESH to Abaqus input file
    print('\n\nGenerate full block mesh (Abaqus inp file)')
    Input_FilneName = FileNames['INPname']

    WriteAbaqusGeneral(Input_FilneName, MESH, FEelSize, templateFile=None, smooth=None)
    inp = mf.readAbaqus(Input_FilneName)
    # title = inp[0]
    nodes = inp[1]
    # nsets = inp[2]
    elems = inp[3]
    elsets = inp[4]
    # element set "SET1" correspondes to the full block including padding layers

    # New element sets "BONE" and "GHOST" are created after material mapping
    elems = {elem: elems[elem] for elem in elsets["SET1"]}
    nodes = {node: nodes[node] for elem in elems for node in elems[elem].get_nodes()}
    elsets = {}
    print(" - finished")
    # *****************************************************************

    # set bone values
    bone["elems"] = elems
    bone["nodes"] = nodes
    bone["elsets"] = elsets
    bone["MESH"] = MESH

    return bone
def Calculate_Iso_Fabric(bone, config):

    """
    Adapted from Denis's preprocessing_SA.py
    Compute isotropic fabric
    """

    FTYPE = config["ftype"]
    if FTYPE.find("Iso") > -1 or FTYPE.find("iso") > -1 or FTYPE.find("test") > -1:
        print("\n ... compute isotropic fabric")
        evalue, evect = utils.compute_isoFAB()
        print(" - finished")

    io_utils.print_mem_usage()

    bone["evalue"] = evalue
    bone["evect"] = evect

    return bone
def AIM2FE_SA_PSL(Config, Sample):

    """
    Convert an AIM file to a HDF5 file
    Adapted from Denis's aim2fe_SA_PSL.py
    """

    print('\n\nPerform material mapping for sample: ', Sample)
    FileNames = Set_FileNames(Config, Sample)

    print(yaml.dump(FileNames, default_flow_style=False))
    Print_Memory_Usage()

    # Create bone dictionary, storing arrays and variables
    Bone = {}

    # Read AIM images and image parameters
    Bone = Read_Image_Parameters(FileNames, Bone)
    Image_List = ['BMD', 'SEG', 'CORTMASK', 'TRABMASK']

    for Index, Item in enumerate(Image_List):
        Bone = Read_AIM(Item, FileNames, Bone)
        Bone = Adjust_Image(Item, Bone, Config, 'Crop')

    Print_Memory_Usage()

    # Write out test images to check quality
    Write_TestImages(Bone, Config, Sample)


    # Prepare material mapping
    ImageType = Config['ImageType']
    Bone = Calculate_BVTV(Bone, Config, ImageType)
    Bone = Generate_Mesh(Bone, Config, FileNames)
    Bone = Calculate_Iso_Fabric(Bone, Config)

    # 4 Material mapping
    # Compute MSL kernel list
    bmd = bone['BMDscaled']
    cort_mask = bone['CORTMASK_array']
    trab_mask = bone['TRABMASK_array']
    np.mean(bmd[trab_mask == 127]) / 1200.0

    if config['fabric_type'] == 'local':
        bone = preprocessing.compute_local_MSL(bone, config)
    io_utils.print_mem_usage()
    # store MSL_kernel_list as np array for shortcut during debugging
    # with open(config["feadir"] + '/' + config['folder_id'][sample] + '/' + 'MSL_kernel_list_cort_shortcut.npy', 'wb') as f:
    #     np.save(f, bone['MSL_kernel_list_cort'])
    # with open(config["feadir"] + '/' + config['folder_id'][sample] + '/' + 'MSL_kernel_list_trab_shortcut.npy', 'wb') as f:
    #     np.save(f, bone['MSL_kernel_list_trab'])
    # #
    # with open(config["feadir"] + '/' + config['folder_id'][sample] + '/' + 'MSL_kernel_list_cort_shortcut.npy', 'rb') as f:
    #     bone["MSL_kernel_list_cort"] = np.load(f)
    # with open(config["feadir"] + '/' + config['folder_id'][sample] + '/' + 'MSL_kernel_list_trab_shortcut.npy', 'rb') as f:
    #     bone["MSL_kernel_list_trab"] = np.load(f)

    # # Write out VTK decimate STL file (can't be used anymore, STL file is deleted in compute_MSL
    # STLname = config["feadir"] + '/' + config['folder_id'][sample] + '/' + 'STL_deci.stl'
    # writer = vtk.vtkSTLWriter()
    # STLdeci = bone['STLdeci']
    # writer.SetInputConnection(STLdeci.GetOutputPort())
    # writer.SetFileName(STLname)
    # writer.Write()

    # ---------------------------------------------------------------------------------
    # Ghost layer mode
    if config["mode_ghost_layer"] == 1:
        if config['isotropic_cortex']:
            bone = preprocessing.PSL_material_mapping_copy_layers_accurate_iso_cort(bone, config, umat_parameters,
                                                                                    filenames)
        else:
            bone = preprocessing.PSL_material_mapping_copy_layers_accurate(bone, config, umat_parameters, filenames)

    elif config["mode_ghost_layer"] == 2:
        bone = preprocessing.PSL_material_mapping_predefined_properties(bone, config, umat_parameters, filenames)
    else:
        raise TypeError("Ghost layer mode was not set correctly")
    io_utils.print_mem_usage()
    # 5 Compute and store summary and performance variables
    # ---------------------------------------------------------------------------------
    reload(preprocessing)
    reload(io_utils)
    summary_variables = preprocessing.set_summary_variables(bone)
    io_utils.log_summary(bone, config, filenames, summary_variables)
    bone = dict(list(bone.items()) + list(summary_variables.items()))

    preprocessing.plot_MSL_fabric_fast(config, bone, sample)
    io_utils.print_mem_usage()
    return bone

WorkingDirectory = os.getcwd()
ConfigFile = os.path.join(WorkingDirectory,'03_Scripts/hFE/ConfigFile.yaml')

# Read config and store to dictionary
Config = ReadConfigFile(ConfigFile)

print(yaml.dump(Config, default_flow_style=False))

# Current version of the pipeline/run
Version = Config['Version']

# Directories
AIMDir = Config['AIMDir']
FEADir = Config['FEADir']

# File names and folders
GrayScale_FileNames = Config['GrayScale_FileNames']
Folder_IDs = Config['Folder_IDs']


Sample = GrayScale_FileNames[0]
###############################################################
# for Sample in GrayScale_FileNames:

# Set paths
Folder = Folder_IDs[Sample]
InputFileName = "{}_{}.inp".format(Sample, Version)
InputFile = os.path.join(AIMDir, Folder, InputFileName)
SampleDir = os.path.join(FEADir, Folder)

# Perform material mapping
Bone = aim2fe_SA_PSL(Config, sample)

