#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)

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
Data_Directory = os.path.join(WorkingDirectory,'04_Results/05_FractureLinePrediction/')

SampleList = [Dir for Dir in os.listdir(Data_Directory) if os.path.isdir(Data_Directory+Dir)]
SampleList.sort()

Index = 0
# for Index in range(len(SampleList)):

Sample = SampleList[Index]

# 02 Set Paths and Scans
SamplePath = os.path.join(Data_Directory,Sample)


## Do some tests
uCT_Mask = sitk.ReadImage((SamplePath + '/uCT_Mask_Cropped.mhd'))
J_Cropped = sitk.ReadImage(SamplePath + '/J_Cropped.mhd')
J_Registration = sitk.ReadImage(WorkingDirectory + '/04_Results/04_Registration/432_L_77_F/J.mhd')


Test = J_Cropped*uCT_Mask
TestArray = sitk.GetArrayFromImage(Test)

# Crop J registration with transverse plane
J_Array = sitk.GetArrayFromImage(J_Registration)
J_Zeros = np.zeros(uCT_Mask.GetSize()[::-1])

Shift = int(J_Cropped.GetOrigin()[2] / J_Cropped.GetSpacing()[2])
for i in range(J_Cropped.GetSize()[2]):
    for j in range(J_Array.shape[1]):
        for k in range(J_Array.shape[2]):
            J_Zeros[i,j,k] += J_Array[Shift + i, j, k]

TestArray2 = J_Zeros * sitk.GetArrayFromImage(uCT_Mask)

# Plot
Vmin = min(TestArray.min(),TestArray2.min())
Vmax = max(TestArray.max(),TestArray2.max())


import matplotlib.gridspec as gridspec

Ratio = 1

GS = gridspec.GridSpec(1,2, width_ratios=[1,Ratio])
Figure = plt.figure(figsize=(11,4.5), dpi=500)
Axes1 = plt.subplot(GS[0])
Axes2 = plt.subplot(GS[1])
Axes = [Axes1, Axes2]

Axes[0].imshow(TestArray2[:, :, int(TestArray2.shape[2]/2)], cmap='jet', vmin=Vmin, vmax=Vmax)
CMap = Axes[1].imshow(TestArray[:, :, int(TestArray.shape[2]/2)], cmap='jet', vmin=Vmin, vmax=Vmax)
Axes[0].axis('off')
Axes[1].axis('off')
Axes[0].set_title('Registration')
Axes[1].set_title('hFE')
Figure.colorbar(CMap,ax=Axes,orientation='horizontal',shrink=0.5)
plt.subplots_adjust(bottom=0.25,wspace=0.01)
plt.show()
plt.close(Figure)

Origin = np.array(uCT_Mask.GetOrigin()[::-1])
Spacing = np.array(uCT_Mask.GetSpacing()[::-1])

WriteMHD(TestArray, Spacing, Origin, SamplePath, 'J_hFE', PixelType='float')
WriteMHD(TestArray2, Spacing, Origin, SamplePath, 'J_Registration', PixelType='float')


## Look at correlations
J_hFE = TestArray.flatten()
J_R = TestArray2.flatten()
J_hFE[J_hFE == 0] = np.nan
J_R[J_R == 0] = np.nan

StackLength = TestArray.shape[2] * TestArray.shape[1]
StackColors = plt.cm.jet(np.linspace(0,1,TestArray.shape[0]))
Colors = np.tile(StackColors,StackLength).reshape((StackLength*TestArray.shape[0],4))

ColoredTibia = TestArray.copy()
ColoredTibia[ColoredTibia != 0] = 1
ColoredTibia[ColoredTibia == 0] = np.nan
for Stack in range(ColoredTibia.shape[0]):
    ScanStack = ColoredTibia[Stack,:,:]
    ScanStack[ScanStack == 1] = Stack / ColoredTibia.shape[0]
ColoredTibia = ColoredTibia / np.nanmax(ColoredTibia)




import statsmodels.formula.api as smf
Data = pd.DataFrame([J_R,J_hFE],index=['Registration','hFE']).T
Model = smf.ols("Registration ~ hFE - 1", data=Data).fit()




## Plot results
Y_Obs = Model.model.endog
Y_Fit = Model.fittedvalues
N = int(Model.nobs)
C = np.matrix(Model.cov_params())
X = np.matrix(Model.model.exog)
X_Obs = np.sort(Model.model.exog[:,0])


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
from scipy.stats.distributions import t
Alpha = 0.95
t_Alpha = t.interval(Alpha, N - X.shape[1] - 1)
CI_Line_u = Y_Fit + t_Alpha[0] * SE * B_0
CI_Line_o = Y_Fit + t_Alpha[1] * SE * B_0




import matplotlib.gridspec as gridspec

Ratio = 1

GS = gridspec.GridSpec(1,2, width_ratios=[1,Ratio])
Figure = plt.figure(figsize=(11,4.5), dpi=500)
Axes1 = plt.subplot(GS[0])
Axes2 = plt.subplot(GS[1])
Axes = [Axes1, Axes2]

Axes[0].imshow(ColoredTibia[:,:,int(ColoredTibia.shape[2]/2)],cmap='jet')
Axes[1].scatter(J_R,J_hFE,marker='o',c=Colors)
Axes[1].plot([],linestyle='none',marker='o',color=(0,0,0),label='Data')
Axes[1].fill_between(X_Obs, np.sort(CI_Line_o), np.sort(CI_Line_u), color=(0, 0, 0), alpha=0.1)
Axes[1].plot(X_Obs,np.sort(Y_Fit),linestyle='--',color=(0,0,0),label='Fit')
Axes[0].axis('off')
Axes[0].set_title('Color bar')
Axes[1].set_xlabel('Registration values')
Axes[1].set_ylabel('hFE values')
Axes[1].annotate(r'R$^2$: ' + str(round(R2,3)),(0.75,0.3),xycoords='axes fraction')
Axes[1].annotate(r'SE: ' + str(round(SE,3)),(0.748,0.25),xycoords='axes fraction')
plt.legend(loc='lower right')
plt.show()
plt.close(Figure)

