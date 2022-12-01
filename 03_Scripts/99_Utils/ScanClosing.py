import struct
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats.distributions import t
import SimpleITK as sitk
import os
import pandas as pd

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

## Define functions
def ReadMHD(FileName):

    # Reads the image using SimpleITK
    Itkimage = sitk.ReadImage(FileName)

    # Convert the image to a  numpy array first and then shuffle the dimensions to get axis in the order z,y,x
    CT_Scan = sitk.GetArrayFromImage(Itkimage)

    # Read the origin of the ct_scan, will be used to convert the coordinates from world to voxel and vice versa.
    Origin = np.array(list(reversed(Itkimage.GetOrigin())))

    # Read the spacing along each dimension
    Spacing = np.array(list(reversed(Itkimage.GetSpacing())))

    # Read the dimension of the CT scan (in voxel)
    Size = np.array(list(reversed(Itkimage.GetSize())))

    return CT_Scan, Origin, Spacing, Size
def Closing(img, radius=20):
    """
    morphologic closing operation (fills holes)
    img: numpy array or sitk image
    radius: (optional) closing sphere radius, int or float
    return: closed sitk image
    """
    print('    ... start closing sitk with radius %s mm' % radius)
    if type(img) == np.ndarray:
        img = sitk.GetImageFromArray(img.transpose(2, 1, 0))
    else:
        pass
    if img.GetPixelIDTypeAsString() != '16-bit int':
        img = sitk.Cast(img, sitk.sitkInt16)
    vectorRadius = (int(radius), int(radius), int(radius))
    kernel = sitk.sitkBall
    closed = sitk.BinaryMorphologicalClosing(img, vectorRadius, kernel)
    return closed


# 01 Set variables
WorkingDirectory = os.getcwd()
DataPath = os.path.join(WorkingDirectory,'02_Data/03_uCT/')

## List sample directories
SamplesDirectories = [Dir for Dir in os.listdir(DataPath) if os.path.isdir(DataPath+Dir)]
SamplesDirectories.sort()


DataFrame = pd.DataFrame()
for Sample in SamplesDirectories[:-1]:

    SamplesDirectory = os.path.join(DataPath, Sample)

    SampleFiles = [File for File in os.listdir(SamplesDirectory) if File.endswith("DOWNSCALED_BMD.mhd")]
    SampleFiles.sort()
    SampleFile = SampleFiles[0]

    CT_Scan, Origin, Spacing, Size = ReadMHD(os.path.join(SamplesDirectory, SampleFile))


Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Image = Axes.imshow(CT_Scan[int(Size[0]/2),:,:], cmap='viridis')
ColorBar = plt.colorbar(Image)
plt.show()
plt.close(Figure)


a = Closing(CT_Scan, radius=12)
b = sitk.GetArrayFromImage(a)

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Image = Axes.imshow(b[:,:,int(Size[0]/2)], cmap='viridis')
ColorBar = plt.colorbar(Image)
plt.show()
plt.close(Figure)

sitk.WriteImage(a,'Test.mhd')