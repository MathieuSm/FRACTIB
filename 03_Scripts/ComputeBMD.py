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
def ReadIsqFile(FileName, info=False):

    # Read ISQ from Scanco

    File = open(FileName, 'rb')
    File.seek(44)
    Header = np.zeros(6)
    for i in range(0, 6):
        Header[i] = struct.unpack('i', File.read(4))[0]

    Elsp = [round(Header[3] / Header[0] / 1000.0, 4),
            round(Header[4] / Header[1] / 1000.0, 4),
            round(Header[5] / Header[2] / 1000.0, 4)]
    File.seek(508)
    Headersize = 512 * (1 + struct.unpack('i', File.read(4))[0])
    File.seek(Headersize)
    Ndim = [
        int(Header[0]), int(Header[1]), int(Header[2])]
    NoEntries = Ndim[0] * Ndim[1] * Ndim[2]
    if info == True:
        CT_Scan = None
    else:
        CT_Scan = np.fromfile(File, dtype='i2', count=NoEntries)
    File.close()
    Ldim = [
        float(Elsp[0]), float(Elsp[1]), float(Elsp[2])]
    AdditionalData = {}
    AdditionalData['-Ldim'] = Ldim
    AdditionalData['-Ndim'] = Ndim
    AdditionalData['ElementSpacing'] = Ldim
    AdditionalData['DimSize'] = Ndim
    AdditionalData['HeaderSize'] = Headersize
    AdditionalData['TransformMatrix'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    AdditionalData['CenterOfRotation'] = [0.0, 0.0, 0.0]
    AdditionalData['Offset'] = [0.0, 0.0, 0.0]
    AdditionalData['AnatomicalOrientation'] = 'LPS'
    AdditionalData['ElementType'] = 'int16'
    AdditionalData['ElementDataFile'] = FileName
    if info == False:
        CT_Scan = CT_Scan.reshape((Ndim[2], Ndim[1], Ndim[0]))

    return CT_Scan, AdditionalData
def EdgesDetection(Slice):
    Edges = np.zeros(Slice.shape)
    for X in range(1, Slice.shape[1] - 1):
        for Y in range(1, Slice.shape[0] - 1):
            XMag = Slice[X + 1, Y] - Slice[X - 1, Y]
            YMag = Slice[X, Y + 1] - Slice[X, Y - 1]

            # Draw in black and white the magnitude
            Color = np.sqrt(XMag ** 2 + YMag ** 2)
            Edges[X, Y] = Color
    return Edges
def PlotRegressionResults(Model,Alpha=0.95):

    print(Model.summary())

    ## Plot results
    Y_Obs = Model.model.endog
    Y_Fit = Model.fittedvalues
    N = int(Model.nobs)
    C = np.matrix(Model.cov_params())
    X = np.matrix(Model.model.exog)
    X_Obs = np.sort(np.array(X[:,1]).reshape(len(X)))


    ## Compute R2 and standard error of the estimate
    E = Y_Obs - Y_Fit
    RSS = np.sum(E ** 2)
    SE = np.sqrt(RSS / Model.df_resid)
    TSS = np.sum((Model.model.endog - Model.model.endog.mean()) ** 2)
    RegSS = TSS - RSS
    R2 = RegSS / TSS
    R2adj = 1 - RSS/TSS * (N-1)/(N-X.shape[1]+1-1)

    ## Compute CI lines
    B_0 = np.sqrt(np.diag(np.abs(X * C * X.T)))
    t_Alpha = t.interval(Alpha, N - X.shape[1] - 1)
    CI_Line_u = Y_Fit + t_Alpha[0] * SE * B_0
    CI_Line_o = Y_Fit + t_Alpha[1] * SE * B_0

    t_Alpha2 = t.interval(0.9, N - X.shape[1] - 1)
    CI_Line_u2 = Y_Fit + t_Alpha2[0] * SE * B_0
    CI_Line_o2 = Y_Fit + t_Alpha2[1] * SE * B_0


    ## Plots
    DPI = 100
    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI, sharey=True, sharex=True)
    Axes.plot(X[:,1], Y_Fit, color=(1,0,0))
    Axes.fill_between(X_Obs, np.sort(CI_Line_o2), np.sort(CI_Line_u2), color=(0, 0, 0), alpha=0.1)
    Axes.plot(X_Obs, np.sort(CI_Line_u), color=(0, 0, 1), linestyle='--')
    Axes.plot(X_Obs, np.sort(CI_Line_o), color=(0, 0, 1), linestyle='--')
    Axes.annotate(r'$N$  : ' + str(N), xy=(0.7, 0.2), xycoords='axes fraction')
    Axes.annotate(r'$R^2$ : ' + format(round(R2, 5), '.5f'), xy=(0.7, 0.125), xycoords='axes fraction')
    Axes.annotate(r'$SE$ : ' + format(round(SE, 2), '.2f'), xy=(0.7, 0.05), xycoords='axes fraction')
    Axes.plot(X[:,1], Y_Obs, linestyle='none', marker='o', color=(0,0,0), fillstyle='none',label='Data')
    Axes.plot([], color=(1,0,0), label='Fit')
    Axes.fill_between([], [], color=(0, 0, 0), alpha=0.1, label='90% CI')
    Axes.plot([], color=(0, 0, 1), linestyle='--',label='95% CI')
    Axes.set_xlabel('Phantom gray values (-)')
    Axes.set_ylabel('Cylinder reference densities (mg HA/ccm)')
    plt.subplots_adjust(left=0.15, bottom=0.15)
    plt.legend(loc='upper left')
    plt.show()
    plt.close(Figure)
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
def WriteRaw(ImageArray, OutputFileName):

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
    CastedImageArray = ImageArray.astype(np.short)

    File = np.memmap(OutputFileName, dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape)
    File[:] = CastedImageArray[:]
    del File

    return
def WriteMHD(ImageArray, Spacing, Path, FileName):

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
    outs.write('ElementType = %s\n' % 'MET_SHORT')
    outs.write('ElementDataFile = %s\n' % (FileName + '.raw'))
    outs.close()
    WriteRaw(ImageArray, os.path.join(Path, FileName) + '.raw')
    return



# 01 Set variables
WorkingDirectory = os.getcwd()
DataPath = os.path.join(WorkingDirectory,'02_Data/03_uCT/')

## Load phantom scan
ISQFiles = [File for File in os.listdir(DataPath+'Phantom') if File.endswith('.ISQ')]
ISQFiles.sort()
QC_Scan, AdditionalData = ReadIsqFile(DataPath+'Phantom/'+ISQFiles[0])

## List sample directories
SamplesDirectories = [Dir for Dir in os.listdir(DataPath) if os.path.isdir(DataPath+Dir)]
SamplesDirectories.sort()


# 02 Analysis of phantom to compute BMD equation
Slice = QC_Scan[0,:,:]
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(Slice)
plt.show()
plt.close(Figure)

## Select 1st ROI
ROI_1_Indices = np.array([[1345,1485],[1130,1270]])
Margin = 100
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(Slice[ROI_1_Indices[0][0]-Margin:ROI_1_Indices[0][1]+Margin,
            ROI_1_Indices[1][0]-Margin:ROI_1_Indices[1][1]+Margin])
Axes.plot([Margin,Margin+ROI_1_Indices[1][1]-ROI_1_Indices[1][0]],
          [Margin,Margin],
          color=(1,0,0))
Axes.plot([Margin,Margin+ROI_1_Indices[1][1]-ROI_1_Indices[1][0]],
          [Margin+ROI_1_Indices[0][1]-ROI_1_Indices[0][0],Margin+ROI_1_Indices[0][1]-ROI_1_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin,Margin],
          [Margin,Margin+ROI_1_Indices[0][1]-ROI_1_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin+ROI_1_Indices[1][1]-ROI_1_Indices[1][0],Margin+ROI_1_Indices[1][1]-ROI_1_Indices[1][0]],
          [Margin,Margin+ROI_1_Indices[0][1]-ROI_1_Indices[0][0]],
          color=(1,0,0))
plt.show()
plt.close(Figure)

## Select 2nd ROI
ROI_2_Indices = np.array([[1750,1890],[1245,1385]])
Margin = 100
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(Slice[ROI_2_Indices[0][0]-Margin:ROI_2_Indices[0][1]+Margin,
            ROI_2_Indices[1][0]-Margin:ROI_2_Indices[1][1]+Margin])
Axes.plot([Margin,Margin+ROI_2_Indices[1][1]-ROI_2_Indices[1][0]],
          [Margin,Margin],
          color=(1,0,0))
Axes.plot([Margin,Margin+ROI_2_Indices[1][1]-ROI_2_Indices[1][0]],
          [Margin+ROI_2_Indices[0][1]-ROI_2_Indices[0][0],Margin+ROI_2_Indices[0][1]-ROI_2_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin,Margin],
          [Margin,Margin+ROI_2_Indices[0][1]-ROI_2_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin+ROI_2_Indices[1][1]-ROI_2_Indices[1][0],Margin+ROI_2_Indices[1][1]-ROI_2_Indices[1][0]],
          [Margin,Margin+ROI_2_Indices[0][1]-ROI_2_Indices[0][0]],
          color=(1,0,0))
plt.show()
plt.close(Figure)

## Select 3rd ROI
ROI_3_Indices = np.array([[1760,1900],[1665,1805]])
Margin = 100
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(Slice[ROI_3_Indices[0][0]-Margin:ROI_3_Indices[0][1]+Margin,
            ROI_3_Indices[1][0]-Margin:ROI_3_Indices[1][1]+Margin])
Axes.plot([Margin,Margin+ROI_3_Indices[1][1]-ROI_3_Indices[1][0]],
          [Margin,Margin],
          color=(1,0,0))
Axes.plot([Margin,Margin+ROI_3_Indices[1][1]-ROI_3_Indices[1][0]],
          [Margin+ROI_3_Indices[0][1]-ROI_3_Indices[0][0],Margin+ROI_3_Indices[0][1]-ROI_3_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin,Margin],
          [Margin,Margin+ROI_3_Indices[0][1]-ROI_3_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin+ROI_3_Indices[1][1]-ROI_3_Indices[1][0],Margin+ROI_3_Indices[1][1]-ROI_3_Indices[1][0]],
          [Margin,Margin+ROI_3_Indices[0][1]-ROI_3_Indices[0][0]],
          color=(1,0,0))
plt.show()
plt.close(Figure)

## Select 4th ROI
ROI_4_Indices = np.array([[1370,1510],[1810,1950]])
Margin = 100
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.imshow(Slice[ROI_4_Indices[0][0]-Margin:ROI_4_Indices[0][1]+Margin,
            ROI_4_Indices[1][0]-Margin:ROI_4_Indices[1][1]+Margin])
Axes.plot([Margin,Margin+ROI_4_Indices[1][1]-ROI_4_Indices[1][0]],
          [Margin,Margin],
          color=(1,0,0))
Axes.plot([Margin,Margin+ROI_4_Indices[1][1]-ROI_4_Indices[1][0]],
          [Margin+ROI_4_Indices[0][1]-ROI_4_Indices[0][0],Margin+ROI_4_Indices[0][1]-ROI_4_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin,Margin],
          [Margin,Margin+ROI_4_Indices[0][1]-ROI_4_Indices[0][0]],
          color=(1,0,0))
Axes.plot([Margin+ROI_4_Indices[1][1]-ROI_4_Indices[1][0],Margin+ROI_4_Indices[1][1]-ROI_4_Indices[1][0]],
          [Margin,Margin+ROI_4_Indices[0][1]-ROI_4_Indices[0][0]],
          color=(1,0,0))
plt.show()
plt.close(Figure)



# 03 Compute BMD equation
Total1, Total2, Total3, Total4 = 0,0,0,0
ROI_1_Dim = ROI_1_Indices[:,1]-ROI_1_Indices[:,0]
ROI_2_Dim = ROI_2_Indices[:,1]-ROI_2_Indices[:,0]
ROI_3_Dim = ROI_3_Indices[:,1]-ROI_3_Indices[:,0]
ROI_4_Dim = ROI_4_Indices[:,1]-ROI_4_Indices[:,0]
for SliceIndex in range(QC_Scan.shape[0]):

    Slice = QC_Scan[SliceIndex,:,:]

    Total1 += np.sum(Slice[ROI_1_Indices[0][0]:ROI_1_Indices[0][1],ROI_1_Indices[1][0]:ROI_1_Indices[1][1]])
    Total2 += np.sum(Slice[ROI_2_Indices[0][0]:ROI_2_Indices[0][1],ROI_2_Indices[1][0]:ROI_2_Indices[1][1]])
    Total3 += np.sum(Slice[ROI_3_Indices[0][0]:ROI_3_Indices[0][1],ROI_3_Indices[1][0]:ROI_3_Indices[1][1]])
    Total4 += np.sum(Slice[ROI_4_Indices[0][0]:ROI_4_Indices[0][1],ROI_4_Indices[1][0]:ROI_4_Indices[1][1]])
MeanROI_1 = Total1 / (QC_Scan.shape[0]*ROI_1_Dim[0]*ROI_1_Dim[1])
MeanROI_2 = Total2 / (QC_Scan.shape[0]*ROI_2_Dim[0]*ROI_2_Dim[1])
MeanROI_3 = Total3 / (QC_Scan.shape[0]*ROI_3_Dim[0]*ROI_3_Dim[1])
MeanROI_4 = Total4 / (QC_Scan.shape[0]*ROI_4_Dim[0]*ROI_4_Dim[1])
Means = np.array([MeanROI_1,MeanROI_2,MeanROI_3,MeanROI_4])

## Fit reference densities with gray values
ReferenceDensities = np.array([784.7238, 410.8419, 211.6873, 100.2155])
Data2Fit = pd.DataFrame()
Data2Fit['PhantomValues'] = Means
Data2Fit['ReferenceDensities'] = ReferenceDensities
Fit = smf.ols("ReferenceDensities ~ PhantomValues", data=Data2Fit).fit()
PlotRegressionResults(Fit)



# 04 Write BMD scan and compute values
BMDData = pd.DataFrame()
for Sample in SamplesDirectories[:-1]:

    SamplesDirectory = os.path.join(DataPath, Sample)

    SampleFiles = [File for File in os.listdir(SamplesDirectory) if File.endswith("DOWNSCALED.mhd")]
    SampleFiles.sort()
    SampleFile = SampleFiles[0]

    CT_Scan, Origin, Spacing, Size = ReadMHD(os.path.join(SamplesDirectory,SampleFile))
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # ColorBar = Axes.imshow(CT_Scan[int(Size[0]/2),:,:])
    # plt.colorbar(ColorBar)
    # plt.show()
    # plt.close(Figure)

    BMD_Scan = CT_Scan * Fit.params['PhantomValues'] + Fit.params['Intercept']
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # ColorBar = Axes.imshow(BMD_Scan[int(Size[0]/2),:,:])
    # plt.colorbar(ColorBar)
    # plt.show()
    # plt.close(Figure)

    BMD_Image = sitk.GetImageFromArray(BMD_Scan)
    Otsu_Filter = sitk.OtsuThresholdImageFilter()
    Otsu_Filter.SetInsideValue(0)
    Otsu_Filter.SetOutsideValue(1)
    Segmentation = Otsu_Filter.Execute(BMD_Image)
    Threshold = Otsu_Filter.GetThreshold()

    # Otsu_Scan = np.zeros(BMD_Scan.shape)
    # Otsu_Scan += BMD_Scan
    # Otsu_Scan[Otsu_Scan<threshold] = 0
    # Otsu_Scan[Otsu_Scan>threshold] = 1
    #
    # BMD_Otsu = Otsu_Scan * BMD_Scan
    #
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    # Image = Axes.imshow(BMD_Otsu[int(Size[0] / 2), :, :], cmap='viridis')
    # ColorBar = plt.colorbar(Image)
    # plt.show()
    # plt.close(Figure)

    # WriteMHD(BMD_Scan, Spacing, SamplesDirectory, SampleFile[:-4]+'_BMD')

    Mask, Origin, Spacing, Size = ReadMHD(os.path.join(SamplesDirectory, SampleFile[:-4] + '_MASK.mhd'))
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Axes.imshow(Mask[int(Size[0]/2),:,:])
    # plt.show()
    # plt.close(Figure)


    Volume = Mask.sum() * Spacing[0] * Spacing[1] * Spacing[2]
    MeanArea = Volume / (Size[0] * Spacing[0])


    Shift = (Origin/Spacing).astype('int')
    BMD_Masked = BMD_Scan[Shift[0]:Mask.shape[0]+Shift[0],
                          Shift[1]:Mask.shape[1]+Shift[1],
                          Shift[2]:Mask.shape[2]+Shift[2]] * Mask
    # # Plot scan
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Image = Axes.imshow(BMD_Masked[int(Size[0]/2),:,:],cmap='viridis')
    # ColorBar = plt.colorbar(Image)
    # plt.show()
    # plt.close(Figure)

    ## Plot values distribution
    # BMDValues = np.reshape(BMD_Masked,(BMD_Masked.shape[0]*BMD_Masked.shape[1]*BMD_Masked.shape[2]))
    # MaskMBDValues = pd.DataFrame({'BMD Values':np.round(BMDValues,-2)})
    # MaskMBDDistribution = MaskMBDValues['BMD Values'].value_counts()
    # MaskMBDDistribution.sort_index().plot.bar(color=(0,0,0))
    # plt.show()

    vBMD = BMD_Masked.sum()/Mask.sum()
    BMC = vBMD * (Volume / 1e3)

    SEG_Scan, Origin, Spacing, Size = ReadMHD(os.path.join(SamplesDirectory,SampleFile[:-4]+'_SEG.mhd'))
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Image = Axes.imshow(SEG_Scan[int(Size[0]/2),:,:], cmap='viridis')
    # ColorBar = plt.colorbar(Image)
    # plt.show()
    # plt.close(Figure)
    Shift = (Origin/Spacing).astype('int')
    BMD_SEG = BMD_Scan[Shift[0]:SEG_Scan.shape[0]+Shift[0],
                       Shift[1]:SEG_Scan.shape[1]+Shift[1],
                       Shift[2]:SEG_Scan.shape[2]+Shift[2]] * SEG_Scan

    # # Plot maks + segmented scan
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Image = Axes.imshow(SEG_Scan[int(Size[0]/2),:,:]+Mask[int(Size[0]/2),:,:], cmap='viridis')
    # ColorBar = plt.colorbar(Image)
    # plt.show()
    # plt.close(Figure)
    #
    #
    # # Plot scan
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Image = Axes.imshow(BMD_SEG[int(Size[0]/2),:,:], cmap='viridis')
    # ColorBar = plt.colorbar(Image)
    # plt.show()
    # plt.close(Figure)
    #
    # # Plot values distribution
    # BMDValues = np.reshape(BMD_SEG,(BMD_SEG.shape[0]*BMD_SEG.shape[1]*BMD_SEG.shape[2]))
    # SEGMBDValues = pd.DataFrame({'BMD Values':np.round(BMDValues,-2)})
    # SEGMBDDistribution = SEGMBDValues['BMD Values'].value_counts()
    # SEGMBDDistribution.sort_index().plot.bar(color=(0,0,0))
    # plt.ylim([0,1e7])
    # plt.subplots_adjust(bottom=0.2)
    # plt.show()

    Tissue_BMD = BMD_SEG.sum()/SEG_Scan.sum()
    SEG_Volume = SEG_Scan.sum() * Spacing[0] * Spacing[1] * Spacing[2]
    Tissue_BMC = Tissue_BMD * (SEG_Volume / 1e3)


    BMDData = BMDData.append({'Sample':Sample,
                              'Volume (mm3)':Volume,
                              'Mean Area (mm2)':MeanArea,
                              'BMC (HA mg)':BMC,
                              'vBMD (HA g/cm3)':vBMD,
                              'Tissue BMD (HA g/cm3)':Tissue_BMD,
                              'Tissue BMC (HA mg)':Tissue_BMC,
                              'Tissue BMD ratio (-)':BMD_Scan.max()/Threshold},ignore_index=True)
BMDData.to_csv(os.path.join(DataPath,'BMD_Values.csv'),index=False)

BMDData = pd.read_csv(os.path.join(DataPath,'BMD_Values.csv'))
BMDRatio = pd.DataFrame()
for Sample in SamplesDirectories[:-1]:

    SamplesDirectory = os.path.join(DataPath, Sample)

    SampleFiles = [File for File in os.listdir(SamplesDirectory) if File.endswith("DOWNSCALED.mhd")]
    SampleFiles.sort()
    SampleFile = SampleFiles[0]

    CT_Scan, Origin, Spacing, Size = ReadMHD(os.path.join(SamplesDirectory, SampleFile))
    BMD_Scan = CT_Scan * Fit.params['PhantomValues'] + Fit.params['Intercept']
    WriteMHD(BMD_Scan, Spacing, SamplesDirectory, SampleFile[:-4]+'_BMD')


    # # Plot maks + segmented scan
    # Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
    # Image = Axes.imshow(BMD_Scan[int(Size[0]/2),:,:], cmap='viridis')
    # ColorBar = plt.colorbar(Image)
    # plt.show()
    # plt.close(Figure)

    # Remove negative values
    SEG_Scan = np.zeros(BMD_Scan.shape)
    SEG_Scan += BMD_Scan
    SEG_Scan[SEG_Scan < 0] = 0

    SEG_Image = sitk.GetImageFromArray(SEG_Scan)
    Otsu_Filter = sitk.OtsuThresholdImageFilter()
    Otsu_Filter.SetInsideValue(0)
    Otsu_Filter.SetOutsideValue(1)
    Segmentation = Otsu_Filter.Execute(SEG_Image)
    Threshold = Otsu_Filter.GetThreshold()

    BMDRatio = BMDRatio.append({'Tissue BMD ratio (-)':BMD_Scan.max()/Threshold,
                                'Threshold':Threshold},ignore_index=True)

BMDRatio
BMDRatio['Gray value Threshold'] = (BMDRatio['Threshold'] - Fit.params['Intercept']) / Fit.params['PhantomValues']

BMDData['Tissue BMD ratio (-)'] = BMDRatio['Tissue BMD ratio (-)']
BMDData['BMD Threshold (HA g/cm3)'] = BMDRatio['Threshold']
BMDData['Gray value Threshold'] = BMDRatio['Gray value Threshold']

BMDData.to_csv(os.path.join(DataPath,'BMD_Values.csv'),index=False)


# 05 Boxplots of BMD values
C = BMDData.columns
DataLabel = C[-1]

Figure, Axes = plt.subplots(1, 1, figsize=(3.5, 4.5),dpi=100)
Axes.boxplot(BMDData[DataLabel],vert=True,
             showmeans=False,
             boxprops=dict(linestyle='-',color=(0,0,0)),
             medianprops=dict(linestyle='-',color=(1,0,0)),
             whiskerprops=dict(linestyle='--',color=(0,0,0)),
             meanprops=dict(marker='x',markeredgecolor=(0,0,1)))
Axes.plot([1],BMDData.loc[8,DataLabel],linestyle='none',marker='o',color=(0,0,1))
Axes.set_ylabel(DataLabel)
Axes.set_xticks([])
Axes.set_title('')
plt.title('')
plt.suptitle('')
plt.subplots_adjust(0.3)
plt.show()
plt.close(Figure)
