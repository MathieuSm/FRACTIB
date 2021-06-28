# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)



# 01 Define functions
def Load_Itk(Filename):

    # Reads the image using SimpleITK
    ItkImage = sitk.ReadImage(Filename)

    # Convert the image to a  numpy array first and then shuffle the dimensions to get axis in the order z,y,x
    CT_Scan = sitk.GetArrayFromImage(ItkImage)

    # Read the origin of the ct_scan, will be used to convert the coordinates from world to voxel and vice versa.
    Origin = np.array(list(reversed(ItkImage.GetOrigin())))

    # Read the spacing along each dimension
    Spacing = np.array(list(reversed(ItkImage.GetSpacing())))

    # Read the dimension of the CT scan (in voxel)
    Size = np.array(list(reversed(ItkImage.GetSize())))

    return CT_Scan, Origin, Spacing, Size


## Set variables
WorkingDirectory = os.getcwd()
DataPath = os.path.join(WorkingDirectory,'02_Data/')
ResultsPath = os.path.join(WorkingDirectory,'04_Results/02_MorphoAnalysis/')



## Load files
SampleList = pd.read_csv(DataPath + 'SampleList.csv')

for Index in SampleList.index:
    Sample = SampleList['Internal ID'].loc[Index]
    FileName = 'C000' + str(SampleList['MicroCT pretest file number'].loc[Index]) + '_reso_0.098_DOWNSCALED_SEG.mhd'
    File = os.path.join(DataPath,'03_uCT',Sample,FileName)
    CT_Scan, Origin, Spacing, Size = Load_Itk(File)
    MaskFile = os.path.join(DataPath,'03_uCT',Sample,FileName[:-7]+'MASK.mhd')
    CT_Mask, Mask_Origin, Mask_Spacing, Mask_Size = Load_Itk(MaskFile)

    Morpho_Scan = CT_Scan + CT_Mask

    # ## Plot mid-planes
    # Ratio = Size[2] / Size[1]
    # GS = gridspec.GridSpec(1,2, width_ratios=[1,Ratio])
    # Figure = plt.figure(figsize=(11,4.5), dpi=100)
    # Axes1 = plt.subplot(GS[0])
    # Axes2 = plt.subplot(GS[1])
    # Axes = [Axes1, Axes2]
    #
    # X, Y = int(Size[2]/2), int(Size[1]/2)
    #
    # for Plane in range(2):
    #
    #     if Plane == 0:
    #         Axes[Plane].imshow(Morpho_Scan[:, :, X], cmap='bone', clim=[0,2])
    #         Axes[Plane].set_xlabel('Direction 2 (voxel)')
    #         Axes[Plane].set_ylabel('Direction 3 (voxel)')
    #
    #     else:
    #         Axes[Plane].imshow(Morpho_Scan[:, Y, :], cmap='bone', clim=[0,2])
    #         Axes[Plane].set_yticks([])
    #         Axes[Plane].set_xlabel('Direction 1 (voxel)')
    #
    #     Axes[Plane].set_xlim([0, Size[Plane + 1]])
    #     Axes[Plane].set_ylim([0, Size[0]])
    #
    # plt.show()
    # plt.close(Figure)


    ## Write binarized image
    Image = sitk.GetImageFromArray(Morpho_Scan)
    ImagePath = ResultsPath + FileName[:-29] + 'Binarized.mhd'
    sitk.WriteImage(Image, ImagePath)

    # read input file
    MHD = open(File, "rt")
    # read file contents to string
    MHDText = MHD.read()
    # replace all occurrences of the required string
    OldName = FileName[:-4]
    MHDText = MHDText.replace(OldName, FileName[:-29] + 'Binarized')
    MHDText = MHDText.replace('MET_CHAR','MET_UCHAR')
    # Delete last line
    MHDText = MHDText[:-73]
    # Rewrite resolution
    Start = MHDText.find('ElementSpacing')
    Stop = MHDText.find('DimSize')
    MHDText = MHDText.replace(MHDText[Start+17:Stop-1], '0.098 0.098 0.098')
    # close the input file
    MHD.close()
    # open the input file in write mode
    MHD = open(ImagePath, "wt")
    # overrite the input file with the resulting data
    MHD.write(MHDText)
    # close the file
    MHD.close()



# 05 Write the analysis batch
BatchFile = open(os.path.join(WorkingDirectory, '03_Scripts', '01_ScanBatch.bash'),'w')
Binarized_ROIs = [File for File in os.listdir(ResultsPath) if File.endswith('.mhd')]

PipelineFile  = open(os.path.join(WorkingDirectory, '03_Scripts', '01_MorphometryPipeline.txt'), 'r')
PipelineText  = PipelineFile.read()

for Binarized_ROI in Binarized_ROIs:

    ROIText = PipelineText
    ROIText = ROIText.replace('$Scan',Binarized_ROI[:-4])

    BatchFile.write(ROIText+'\n\n')

BatchFile.close()
