#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used to compare trabecular morphometry results
    from medtool, scanco, and pypore3d.

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: April 2023
    """

#%% Imports
# Modules import

import os
import vtk
import time
import struct
import argparse
import numpy as np
import pandas as pd
import SimpleITK as sitk
from pathlib import Path
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats.distributions import t
from vtk.util.numpy_support import vtk_to_numpy # type: ignore
from pypore3d.p3dSITKPy import py_p3dReadRaw8 as ReadRaw8
from pypore3d.p3dBlobPy import py_p3dMorphometricAnalysis as MA


#%% Functions
# Define functions

def Trabecular(Image, Mask=None):

    """
    Perform trabecular morphological analysis
    :param Image: Segmented image of trabecular bone
                    - Type: binary 0-255 sitkImage
    :param Mask : Mask for region of interest selection
                    - Type: binary 0-255 sitkImage
    :param DA   : Compute or not the degree of anisotropy
                    using mean intercept length
    :return Data: Pandas data frame with computed parameters
    """

    Data = pd.DataFrame()

    # Write temporary images to disk
    dX, dY, dZ = Image.GetSize()
    Spacing = Image.GetSpacing()[0]
    sitk.WriteImage(Image, 'TempROI.mhd')

    if Mask:
        sitk.WriteImage(Mask, 'TempMask.mhd')
    else:
        Mask = np.ones((dZ, dY, dZ), 'int') * 255
        Mask = sitk.GetImageFromArray(Mask)
        Mask.SetSpacing(Image.GetSpacing())
        sitk.WriteImage(Mask, 'TempMask.mhd')

    # Perform morphometric analysis
    ROI = ReadRaw8('TempROI.raw', dX, dY, dimz=dZ)
    ROIMask = ReadRaw8('TempMask.raw', dX, dY, dimz=dZ)
    MS = MA(ROI, ROIMask, dX, dY, dimz=dZ, resolution=Spacing)

    # Remove temporary files and write data to disk
    for F in ['TempROI','TempMask']:
        os.remove(F + '.mhd')
        os.remove(F + '.raw')

    # Store data
    Props = ['BV/TV (-)', 'Tb.Th. (mm)', 'Tb.N. (-)', 'Tb.Sp. (mm)']
    Measures = [MS.BvTv, MS.TbTh, MS.TbN, MS.TbSp]
    for Prop, Stat in zip(Props, Measures):
        Data.loc[0,Prop] = Stat

    return Data
def OLS(X, Y, Cmap=np.array(None), Labels=None, Alpha=0.95, Annotate=['N','R2','SE','Slope','Intercept'], FName=None, ShowPlot=True):

    if Labels == None:
        Labels = ['X', 'Y']
    
    # Perform linear regression
    Array = np.array([X,Y])
    if Array.shape[0] == 2:
        Array = Array.T
    Data = pd.DataFrame(Array,columns=['X','Y'])
    FitResults = smf.ols('Y ~ X', data=Data).fit()
    Slope = FitResults.params[1]

    # Build arrays and matrices
    Y_Obs = FitResults.model.endog
    Y_Fit = FitResults.fittedvalues
    N = int(FitResults.nobs)
    C = np.matrix(FitResults.normalized_cov_params)
    X = np.matrix(FitResults.model.exog)

    # Sort X values and Y accordingly
    Sort = np.argsort(np.array(X[:,1]).reshape(len(X)))
    X_Obs = np.sort(np.array(X[:,1]).reshape(len(X)))
    Y_Fit = Y_Fit[Sort]
    Y_Obs = Y_Obs[Sort]

    ## Compute R2 and standard error of the estimate
    E = Y_Obs - Y_Fit
    RSS = np.sum(E ** 2)
    SE = np.sqrt(RSS / FitResults.df_resid)
    TSS = np.sum((FitResults.model.endog - FitResults.model.endog.mean()) ** 2)
    RegSS = TSS - RSS
    R2 = RegSS / TSS
    R2adj = 1 - RSS/TSS * (N-1)/(N-X.shape[1]+1-1)

    ## Compute CI lines
    B_0 = np.sqrt(np.diag(np.abs(X * C * X.T)))
    t_Alpha = t.interval(Alpha, N - X.shape[1] - 1)
    CI_Line_u = Y_Fit + t_Alpha[0] * SE * B_0[Sort]
    CI_Line_o = Y_Fit + t_Alpha[1] * SE * B_0[Sort]

    ## Plots
    DPI = 100
    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=DPI, sharey=True, sharex=True)

    if Cmap.any():
        Colors = plt.cm.winter((Cmap-min(Cmap))/(max(Cmap)-min(Cmap)))
        Scatter = Axes.scatter(X_Obs, Y_Obs, facecolor='none', edgecolor=Colors, marker='o',)
    else:
        Axes.plot(X_Obs, Y_Obs, linestyle='none', marker='o', color=(0,0,1), fillstyle='none')

    Axes.plot(X_Obs, Y_Fit, color=(1,0,0))
    Axes.fill_between(X_Obs, CI_Line_o, CI_Line_u, color=(0, 0, 0), alpha=0.1)

    if Slope > 0:

        YPos = 0.925
        if 'N' in Annotate:
            Axes.annotate(r'$N$  : ' + str(N), xy=(0.025, YPos), xycoords='axes fraction')
            YPos -= 0.075
        if 'R2' in Annotate:
            Axes.annotate(r'$R^2$ : ' + format(round(R2, 2), '.2f'), xy=(0.025, YPos), xycoords='axes fraction')
            YPos -= 0.075
        if 'SE' in Annotate:
            Axes.annotate(r'$SE$ : ' + format(round(SE, 2), '.2f'), xy=(0.025, YPos), xycoords='axes fraction')
        
        YPos = 0.025
        if 'Intercept' in Annotate:
            Intercept = str(FitResults.params[0])
            Round = 3 - Intercept.find('.')
            Intercept = round(FitResults.params[0], Round)
            CI = FitResults.conf_int().loc['Intercept'].round(Round)
            if Round <= 0:
                Intercept = int(Intercept)
                CI = [int(v) for v in CI]
            Text = r'Intercept : ' + str(Intercept) + ' (' + str(CI[0]) + ',' + str(CI[1]) + ')'
            Axes.annotate(Text, xy=(0.425, YPos), xycoords='axes fraction')
            YPos += 0.075

        if 'Slope' in Annotate:
            Round = 3 - str(FitResults.params[1]).find('.')
            Slope = round(FitResults.params[1], Round)
            CI = FitResults.conf_int().loc['X'].round(Round)
            if Round <= 0:
                Slope = int(Slope)
                CI = [int(v) for v in CI]
            Text = r'Slope : ' + str(Slope) + ' (' + str(CI[0]) + ',' + str(CI[1]) + ')'
            Axes.annotate(Text, xy=(0.425, YPos), xycoords='axes fraction')

    elif Slope < 0:

        YPos = 0.025
        if 'N' in Annotate:
            Axes.annotate(r'$N$  : ' + str(N), xy=(0.025, YPos), xycoords='axes fraction')
            YPos += 0.075
        if 'R2' in Annotate:
            Axes.annotate(r'$R^2$ : ' + format(round(R2, 2), '.2f'), xy=(0.025, YPos), xycoords='axes fraction')
            YPos += 0.075
        if 'SE' in Annotate:
            Axes.annotate(r'$SE$ : ' + format(round(SE, 2), '.2f'), xy=(0.025, YPos), xycoords='axes fraction')
        
        YPos = 0.925
        if 'Intercept' in Annotate:
            Intercept = str(FitResults.params[0])
            Round = 3 - Intercept.find('.')
            Intercept = round(FitResults.params[0], Round)
            CI = FitResults.conf_int().loc['Intercept'].round(Round)
            if Round <= 0:
                Intercept = int(Intercept)
                CI = [int(v) for v in CI]
            Text = r'Intercept : ' + str(Intercept) + ' (' + str(CI[0]) + ',' + str(CI[1]) + ')'
            Axes.annotate(Text, xy=(0.425, YPos), xycoords='axes fraction')
            YPos -= 0.075

        if 'Slope' in Annotate:
            Round = 3 - str(FitResults.params[1]).find('.')
            Slope = round(FitResults.params[1], Round)
            CI = FitResults.conf_int().loc['X'].round(Round)
            if Round <= 0:
                Slope = int(Slope)
                CI = [int(v) for v in CI]
            Text = r'Slope : ' + str(Slope) + ' (' + str(CI[0]) + ',' + str(CI[1]) + ')'
            Axes.annotate(Text, xy=(0.425, YPos), xycoords='axes fraction')
    
    Axes.set_xlabel(Labels[0])
    Axes.set_ylabel(Labels[1])
    plt.subplots_adjust(left=0.15, bottom=0.15)

    if (FName):
        plt.savefig(FName, bbox_inches='tight', pad_inches=0.02)
    if ShowPlot:
        plt.show()
    else:
        plt.close()

    return FitResults
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
def ReadAIM(File, Echo=False):

    """
    Reads an AIM file and provides
    the corresponding itk image additional
    data (i.e. spacing, calibration data, 
    and header)

    from Denis hFE pipeline
    """

    if Echo:
        Text = 'Read AIM'
        Time.Process(1, Text)

    # read header
    with open(File, 'rb') as f:
        AIM_Ints = Get_AIM_Ints(f)
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

    if Echo:
        Time.Process(0, Text)

    return Image, AdditionalData

#%% Time functions
class Time():

    def __init__(self):
        self.Width = 15
        self.Length = 16
        self.Text = 'Process'
        self.Tic = time.time()
        pass
    
    def Set(self, Tic=None):
        
        if Tic == None:
            self.Tic = time.time()
        else:
            self.Tic = Tic

    def Print(self, Tic=None,  Toc=None):

        """
        Print elapsed time in seconds to time in HH:MM:SS format
        :param Tic: Actual time at the beginning of the process
        :param Toc: Actual time at the end of the process
        """

        if Tic == None:
            Tic = self.Tic
            
        if Toc == None:
            Toc = time.time()


        Delta = Toc - Tic

        Hours = np.floor(Delta / 60 / 60)
        Minutes = np.floor(Delta / 60) - 60 * Hours
        Seconds = Delta - 60 * Minutes - 60 * 60 * Hours

        print('\nProcess executed in %02i:%02i:%02i (HH:MM:SS)' % (Hours, Minutes, Seconds))

        return

    def Update(self, Progress, Text=''):

        Percent = int(round(Progress * 100))
        Np = self.Width * Percent // 100
        Nb = self.Width - Np

        if len(Text) == 0:
            Text = self.Text
        else:
            self.Text = Text

        Ns = self.Length - len(Text)
        if Ns >= 0:
            Text += Ns*' '
        else:
            Text = Text[:self.Length]
        
        Line = '\r' + Text + ' [' + Np*'=' + Nb*' ' + ']' + f' {Percent:.0f}%'
        print(Line, sep='', end='', flush=True)

    def Process(self, StartStop:bool, Text=''):

        if len(Text) == 0:
            Text = self.Text
        else:
            self.Text = Text

        if StartStop*1 == 1:
            print('')
            self.Tic = time.time()
            self.Update(0, Text)

        elif StartStop*1 == 0:
            self.Update(1, Text)
            self.Print()

Time = Time()
#%% Main
# Main code

def Main():

    # List data
    WD = Path.cwd()
    DD = WD / '..' / '..' / '02_Data'
    uD = DD / '02_uCT'
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))

    # Set dataframes
    Props = ['BV/TV (-)', 'Tb.Th. (mm)', 'Tb.N. (-)', 'Tb.Sp. (mm)']
    MProps = ['$BVTV_voxel', '$Tb_Th_mean', '$Tb_N_mean', '$Tb_Sp_mean']
    SProps = ['VOX-BV/TV','DT-Tb.Th','DT-Tb.N','DT-Tb.Sp']

    Columns = pd.MultiIndex.from_product([['Medtool', 'Scanco'],Props])
    PyPore = pd.DataFrame(columns=Columns, index=SampleList['Internal ID'])
    Medtool = pd.DataFrame(columns=Props)
    Scanco = pd.DataFrame(columns=Props)

    for Index, Row in SampleList.iterrows():

        Sample = Row['Internal ID']
        Time.Process(1, Sample)

        # Read and store medtool values
        uCT = 'C000' + str(Row['MicroCT pretest file number']) + '_Binarized.csv'
        MR = pd.read_csv(str(WD / 'Medtool' / uCT), sep=';')
        for iP, P in enumerate(Props):
            Medtool.loc[Sample, P] = MR.loc[0,MProps[iP]]

        # Use PyPore for morphometry and store results
        Time.Update(1/6, 'Medtool')
        MHD = sitk.ReadImage(str(WD / 'Medtool' / (uCT[:-4] + '.mhd')))
        Trab = Trabecular((MHD == 2) * 255)
        for P in Props:
            PyPore.loc[Sample,'Medtool'][P] = Trab.loc[0,P]

        # Read and store Scanco values
        uCT = uCT[:8] + '_DOWNSCALED_BONE_MORPHO.TXT'
        SR = pd.read_csv(str(WD / 'Scanco' / uCT), sep='\t')
        for iP, P in enumerate(Props):
            Scanco.loc[Sample, P] = SR.loc[0,SProps[iP]]
        Scanco.loc[Sample, 'DA (-)'] = SR.loc[0,'TRI-DA']

        # Read AIMs to compare with Scanco
        Time.Update(2/6, 'Read AIMs')
        AIM = ReadAIM(str(uD / Sample / (uCT[:19] + '.AIM')))[0]
        SEG = ReadAIM(str(uD / Sample / (uCT[:19] + '_TRAB_MASK.AIM')))[0]

        # Perform same image analysis as Scanco
        Time.Update(3/6, 'Gauss Filter')
        Filter = sitk.DiscreteGaussianImageFilter()
        Filter.SetVariance(SR.loc[0,'Sigma'])
        Filter.SetMaximumKernelWidth(int(round(SR.loc[0,'Support'] / AIM.GetSpacing()[0])))
        Filter.SetUseImageSpacing(False)
        Gauss = Filter.Execute(AIM)

        Time.Update(4/6, 'Binarize Scans')
        Filter = sitk.BinaryThresholdImageFilter()
        Filter.SetUpperThreshold(float(SR.loc[0,'Data-Threshold']))
        Filter.SetOutsideValue(0)
        Bin = (1 - Filter.Execute(Gauss)) * 255

        Filter.SetUpperThreshold(1)
        Seg = (1 - Filter.Execute(SEG)) * 255

        # Use Pypore for morphometry and store results
        Time.Update(5/6, 'Scanco')
        Trab = Trabecular(Bin, Mask=Seg)
        for P in Props:
            PyPore.loc[Sample,'Scanco'][P] = Trab.loc[0,P]

        Time.Process(0, Sample)


    # Perform linear regressions and show plot
    for P in Props:
        print(P)
        Fit = OLS(Medtool[P].astype('float'), PyPore['Medtool', P].astype('float'),
                Labels=['Medtool', 'Pypore3D'])
        Fit = OLS(Scanco[P].astype('float'), PyPore['Scanco', P].astype('float'),
                Labels=['Scanco', 'Pypore3D'])
        Fit = OLS(Scanco[P].astype('float'), Medtool[P].astype('float'),
                Labels=['Scanco', 'Medtool'])



    return

#%% Execution part
# Execution as main
if __name__ == '__main__':

    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)

    # Add long and short argument
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-V', '--Version', help='Show script version', action='version', version=SV)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main()