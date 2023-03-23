#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    This script aims to provide utility functions to
    manipulate image (reading, ploting, registering,
    and writing), signal processing, abaqus file read-
    ing, and bone analysis and tensor algebra. Most 
    functions are taken from different sources and adapted here
    
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
from numba import njit
import SimpleITK as sitk
from pathlib import Path
import scipy.signal as sig
from numba.core import types
import matplotlib.pyplot as plt
from multiprocessing import Pool
from numba.typed import Dict, List
import statsmodels.formula.api as smf
from skimage import measure, morphology
from matplotlib.colors import ListedColormap
from scipy.stats.distributions import t, norm
from vtk.util.numpy_support import vtk_to_numpy # type: ignore
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pypore3d.p3dSITKPy import py_p3dReadRaw8 as ReadRaw8
from pypore3d.p3dBlobPy import py_p3dMorphometricAnalysis as MA



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

def RotationMatrix(Phi=0.0, Theta=0.0, Psi=0.0, V=np.zeros(3), A=0):

    if (V != 0).any():
        a = np.cos(A) * np.eye(3)
        b = np.sin(A) * np.array([[0, -V[2], V[1]],[V[2], 0, -V[0]],[-V[1], V[0], 0]])
        c = (1-np.cos(A)) * np.outer(V, V)
        R = np.round(a + b + c, 15)

    else:

        # if list of angles, use numpy for speed
        try:
            len(Phi)
            Phi, Theta, Psi = np.array(Phi), np.array(Theta), np.array(Psi)
            Rx = np.array([[np.ones(len(Phi)),  np.zeros(len(Phi)), np.zeros(len(Phi))],
                           [np.zeros(len(Phi)),        np.cos(Phi),       -np.sin(Phi)],
                           [np.zeros(len(Phi)),        np.sin(Phi),        np.cos(Phi)]])

            Ry = np.array([[ np.cos(Theta),      np.zeros(len(Theta)),        np.sin(Theta)],
                           [np.zeros(len(Theta)), np.ones(len(Theta)), np.zeros(len(Theta))],
                           [-np.sin(Theta),      np.zeros(len(Theta)),        np.cos(Theta)]])

            Rz = np.array([[np.cos(Psi),              -np.sin(Psi), np.zeros(len(Psi))],
                           [np.sin(Psi),               np.cos(Psi), np.zeros(len(Psi))],
                           [np.zeros(len(Psi)), np.zeros(len(Psi)),  np.ones(len(Psi))]])

            R = np.einsum('ijl,jkl->lik',Rz, np.einsum('ijl,jkl->ikl',Ry, Rx))

        # if only float angles, use sympy for more accuracy
        except:
            Rx = sp.Matrix([[1,             0,              0],
                            [0, sp.cos(Phi), -sp.sin(Phi)],
                            [0, sp.sin(Phi),  sp.cos(Phi)]])

            Ry = sp.Matrix([[ sp.cos(Theta), 0, sp.sin(Theta)],
                            [0,             1,              0],
                            [-sp.sin(Theta), 0, sp.cos(Theta)]])

            Rz = sp.Matrix([[sp.cos(Psi), -sp.sin(Psi),     0],
                            [sp.sin(Psi),  sp.cos(Psi),     0],
                            [0,             0,              1]])

            R = Rz * Ry * Rx
    
    return np.array(R, dtype='float')

def GetAngles(R):

    # Compute Euler angles from rotation matrix
    # Assuming R = RxRyRz
    # Adapted from https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix

    if len(R.shape) == 2:
        # Special case
        if R[0,2] == 1 or R[0,2] == -1:
            E3 = 0 # Set arbitrarily
            dlta = np.arctan2(R[0,1],R[0,2])

            if R[0,2] == -1:
                E2 = np.pi/2;
                E1 = E3 + dlta

            else:
                E2 = -np.pi/2;
                E1 = -E3 + dlta

        else:
            E2 = - np.arcsin(R[0,2])
            E1 = np.arctan2(R[1,2]/np.cos(E2), R[2,2]/np.cos(E2))
            E3 = np.arctan2(R[0,1]/np.cos(E2), R[0,0]/np.cos(E2))
    
    else:
        E1, E2, E3 = np.zeros(len(R)), np.zeros(len(R)), np.zeros(len(R))
        
        M1 = R[:,0,2] == 1
        M2 = R[:,0,2] == -1
        if sum(M1 + M2) > 0:
            dlta = np.arctan2(R[:,0,1],R[:,0,2])

            if sum(M2) > 0:
                E2[M2] = np.pi/2
                E1[M2] = E3 + dlta

            else:
                E2[M1] = -np.pi/2
                E1[M1] = -E3 + dlta

        else:
            E2 = - np.arcsin(R[:,0,2])
            E1 = np.arctan2(R[:,1,2]/np.cos(E2), R[:,2,2]/np.cos(E2))
            E3 = np.arctan2(R[:,0,1]/np.cos(E2), R[:,0,0]/np.cos(E2))

    return np.array([-E1, -E2, -E3]).T

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

def Resample(Image, Factor=None, Size=[None], Spacing=[None], Order=0):

    Dimension = Image.GetDimension()
    OriginalSpacing = np.array(Image.GetSpacing())
    OriginalSize = np.array(Image.GetSize())
    PhysicalSize = OriginalSize * OriginalSpacing

    Origin = Image.GetOrigin()
    Direction = Image.GetDirection()

    if Factor:
        NewSize = [round(Size/Factor) for Size in Image.GetSize()] 
        NewSpacing = [PSize/(Size-1) for Size,PSize in zip(NewSize, PhysicalSize)]
    
    elif Size[0]:
        NewSize = Size
        NewSpacing = [PSize/Size for Size, PSize in zip(NewSize, PhysicalSize)]
    
    elif Spacing[0]:
        NewSpacing = Spacing
        NewSize = [np.floor(Size/Spacing).astype('int') + 1 for Size,Spacing in zip(PhysicalSize, NewSpacing)]

    NewArray = np.zeros(NewSize[::-1],'int')
    NewImage = sitk.GetImageFromArray(NewArray)
    NewImage.SetOrigin(Origin - OriginalSpacing/2)
    NewImage.SetDirection(Direction)
    NewImage.SetSpacing(NewSpacing)
  
    Transform = sitk.TranslationTransform(Dimension)
    Resampled = sitk.Resample(Image, NewImage, Transform, Order+1)
    
    return Resampled


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
#%% Ploting functions
class Show():

    def __init__(self):
        self.FName = None
        self.ShowPlot = True
        self.IRange = [0.8, 1.2]

    def Normalize(self, Array, uint=False):

        Min = np.min(Array)
        Max = np.max(Array)
        N_Array = (Array - Min) / (Max - Min)

        if uint:
            N_Array = np.array(N_Array * 255).astype('uint8') 

        return N_Array

    def Slice(self, Image, Slice=None, Title=None, Axis='Z'):

        try:
            Array = sitk.GetArrayFromImage(Image)
            Dimension = Image.GetDimension()
        except:
            Array = Image
            Dimension = len(Array.shape)

        if Dimension == 3:
            
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

        FixedArray = self.Normalize(FixedArray.astype(float))
        MovingArray = self.Normalize(MovingArray.astype(float))


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
        Structure[:,:,3] = self.Normalize(Array.astype(float)) / 255
        
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

    def Signal(self, X, Y=[], Points=[], Normalize=False, Axes=[], Labels=[], Legend=True):

        if len(X) > 6:
            N = len(X)
            Values = np.ones((N, 4))
            Values[:N//3*2, 0] = np.linspace(1, 0, N//3*2)
            Values[N//3*2:, 0] = 0

            Values[:N//3*2, 2] = np.linspace(0, 1, N//3*2)
            Values[-N//3:, 2] = np.linspace(1, 0, N//3+1)

            Values[:-N//3, 1] = 0
            Values[-N//3:, 1] = np.linspace(0, 1, N//3+1)
            Colors = ListedColormap(Values)(np.linspace(0,1,N))
        else:
            Colors = [(1,0,0), (0,0,1), (0,0,0), (0,1,0), (0,1,1), (1,0,1)]

        Figure, Axis = plt.subplots(1,1)

        if len(Y) == 0:
            self.Y = []
            for ix, x in enumerate(X):
                self.Y.append(x)
                X[ix] = np.arange(len(x))
            Y = self.Y
            delattr(self, 'Y')

        for i in range(len(X)):

            if Normalize:
                Yi = self.Normalize(Y[i])
            else:
                Yi = Y[i]
            Xi = X[i]
            
            if len(Labels) > 0:
                Axis.plot(Xi, Yi, color=Colors[i], label=Labels[i])
            else:
                Axis.plot(Xi, Yi, color=Colors[i], label='Signal ' + str(i+1))

            if len(Points) > 0:
                Px, Py = Xi[Points[i]], Yi[Points[i]]
                Axis.plot(Px, Py, marker='o', color=(0, 0, 0), fillstyle='none', linestyle='none')

        if len(Axes) > 0:
            Axis.set_xlabel(Axes[0])
            Axis.set_ylabel(Axes[1])

        if i < 3:
            Cols = i+1
        else:
            Cols = (1+i)//2

        if Legend:
            plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.12), ncol=Cols)

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0.02)

        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

        return

    def OLS(self, X, Y, Cmap=np.array(None), Labels=None, Alpha=0.95, Annotate=['N','R2','SE','Slope','Intercept']):

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

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0.02)
        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

        return FitResults

    def BoxPlot(self, ArraysList, Labels=['', 'Y'], SetsLabels=None, Vertical=True):

        Figure, Axis = plt.subplots(1,1)

        for i, Array in enumerate(ArraysList):
            RandPos = np.random.normal(i,0.02,len(Array))

            Axis.boxplot(Array, vert=Vertical, widths=0.35,
                        showmeans=False,meanline=True,
                        showfliers=False, positions=[i],
                        capprops=dict(color=(0,0,0)),
                        boxprops=dict(color=(0,0,0)),
                        whiskerprops=dict(color=(0,0,0),linestyle='--'),
                        medianprops=dict(color=(0,0,1)),
                        meanprops=dict(color=(0,1,0)))
            Axis.plot(RandPos - RandPos.mean() + i, Array, linestyle='none',
                      marker='o',fillstyle='none', color=(1,0,0))
        
        Axis.plot([],linestyle='none',marker='o',fillstyle='none', color=(1,0,0), label='Data')
        Axis.plot([],color=(0,0,1), label='Median')
        Axis.set_xlabel(Labels[0])
        Axis.set_ylabel(Labels[1])

        if SetsLabels:
            Axis.set_xticks(np.arange(len(SetsLabels)))
            Axis.set_xticklabels(SetsLabels, rotation=0)
        else:
            Axis.set_xticks([])
        
        plt.legend()
        plt.subplots_adjust(left=0.25, right=0.75)
        
        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0.02)
        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

    def Fabric(self, eValues, eVectors, nPoints=32, Title=None):

        # New coordinate system
        Q = np.array(eVectors)

        ## Build data for fabric plotting
        u = np.arange(0, 2 * np.pi + 2 * np.pi / nPoints, 2 * np.pi / nPoints)
        v = np.arange(0, np.pi + np.pi / nPoints, np.pi / nPoints)
        X = eValues[0] * np.outer(np.cos(u), np.sin(v))
        Y = eValues[1] * np.outer(np.sin(u), np.sin(v))
        Z = eValues[2] * np.outer(np.ones_like(u), np.cos(v))
        nNorm = np.zeros(X.shape)

        for i in range(len(X)):
            for j in range(len(X)):
                [X[i, j], Y[i, j], Z[i, j]] = np.dot([X[i, j], Y[i, j], Z[i, j]], Q)
                n = np.array([X[i, j], Y[i, j], Z[i, j]])
                nNorm[i, j] = np.linalg.norm(n)

        NormedColor = nNorm - nNorm.min()
        NormedColor = NormedColor / NormedColor.max()

        Figure = plt.figure(figsize=(5.5, 4))
        Axe = Figure.add_subplot(111, projection='3d')
        Axe.plot_surface(X, Y, Z, facecolors=plt.cm.jet(NormedColor), rstride=1, cstride=1, alpha=0.2, shade=False)
        Axe.plot_wireframe(X, Y, Z, color='k', rstride=1, cstride=1, linewidth=0.1)
        
        # scaling hack
        Bbox_min = np.min([X, Y, Z])
        Bbox_max = np.max([X, Y, Z])
        Axe.auto_scale_xyz([Bbox_min, Bbox_max], [Bbox_min, Bbox_max], [Bbox_min, Bbox_max])
        
        # make the panes transparent
        Axe.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        Axe.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        Axe.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        
        # make the grid lines transparent
        Axe.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        Axe.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        Axe.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        
        # modify ticks
        MinX, MaxX = -1, 1
        MinY, MaxY = -1, 1
        MinZ, MaxZ = -1, 1
        Axe.set_xticks([MinX, 0, MaxX])
        Axe.set_yticks([MinY, 0, MaxY])
        Axe.set_zticks([MinZ, 0, MaxZ])
        Axe.xaxis.set_ticklabels([MinX, 0, MaxX])
        Axe.yaxis.set_ticklabels([MinY, 0, MaxY])
        Axe.zaxis.set_ticklabels([MinZ, 0, MaxZ])

        Axe.set_xlabel('X')
        Axe.set_ylabel('Y')
        Axe.set_zlabel('Z')

        if (Title):
            Axe.set_title(Title)

        ColorMap = plt.cm.ScalarMappable(cmap=plt.cm.jet)
        ColorMap.set_array(nNorm)
        if not NormedColor.max() == 1:
            ColorBar = plt.colorbar(ColorMap, 
                                    ticks=[int(NormedColor.mean() - 1),
                                    int(NormedColor.mean()),
                                    int(NormedColor.mean() + 1)])
            plt.cm.ScalarMappable.set_clim(self=ColorMap,
                                           vmin=int(NormedColor.mean() - 1),
                                           vmax=int(NormedColor.mean() + 1))
        ColorBar = plt.colorbar(ColorMap)
        ColorBar.set_label('Vector norm (-)')

        plt.tight_layout()
        plt.show()
        plt.close(Figure)

    def Histogram(self, Arrays, Labels=[], Density=False, Norm=False):

        if len(Arrays) > 6:
            N = len(Arrays)
            Values = np.ones((N, 4))
            Values[:N//3*2, 0] = np.linspace(1, 0, N//3*2)
            Values[N//3*2:, 0] = 0

            Values[:N//3*2, 2] = np.linspace(0, 1, N//3*2)
            Values[-N//3:, 2] = np.linspace(1, 0, N//3+1)

            Values[:-N//3, 1] = 0
            Values[-N//3:, 1] = np.linspace(0, 1, N//3+1)
            Colors = ListedColormap(Values)(np.linspace(0,1,N))
        else:
            Colors = [(1,0,0), (0,0,1), (0,0,0), (0,1,0), (0,1,1), (1,0,1)]

        Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=300)
        for i, Array in enumerate(Arrays):

            X = pd.DataFrame(Array)
            SortedValues = np.sort(X.T.values)[0]
            N = len(X)
            X_Bar = X.mean()
            S_X = np.std(X, ddof=1)

            if i == 0:
                Bins = 20
                Div = (X.max()-X.min()) / Bins
            else:
                Bins = int(round((X.max()-X.min()) / Div))

            ## Histogram and density distribution
            Histogram, Edges = np.histogram(X, bins=Bins)
            Width = 0.9 * (Edges[1] - Edges[0])
            Center = (Edges[:-1] + Edges[1:]) / 2
            Axes.bar(Center, Histogram, align='center', width=Width, edgecolor=Colors[i], color=(1, 1, 1, 0))
            
            if Density and N < 1E3:
                KernelEstimator = np.zeros(N)
                NormalIQR = np.sum(np.abs(norm.ppf(np.array([0.25, 0.75]), 0, 1)))
                DataIQR = np.abs(X.quantile(0.75)) - np.abs(X.quantile(0.25))
                KernelHalfWidth = 0.9 * N ** (-1 / 5) * min(np.abs([S_X, DataIQR / NormalIQR]))
                for Value in SortedValues:
                    KernelEstimator += norm.pdf(SortedValues - Value, 0, KernelHalfWidth * 2)
                KernelEstimator = KernelEstimator / N
            
                Axes.plot(SortedValues, KernelEstimator, color=Colors[i], label='Kernel Density')
            
            if Norm:
                TheoreticalDistribution = norm.pdf(SortedValues, X_Bar, S_X)
                Axes.plot(SortedValues, TheoreticalDistribution, linestyle='--', color=Colors[i], label='Normal Distribution')
            
        if len(Labels) > 0:
            plt.xlabel(Labels[0])
            plt.ylabel(Labels[1])
        
        # plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.15), prop={'size': 10})
        # plt.legend(loc='upper left')

        if (self.FName):
            plt.savefig(self.FName, bbox_inches='tight', pad_inches=0.02)
        if self.ShowPlot:
            plt.show()
        else:
            plt.close()

Show = Show()
#%% Reading functions
class Read():

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
            Text = 'Read AIM'
            Time.Process(1, Text)

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
            Time.Process(0, Text)

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
        92	 int     nr_of_samples;
        96	 int     nr_of_projections;
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

        if self.Echo:
            Text = 'Read ISQ'
            Time.Process(1, Text)

        try:
            f = open(File, 'rb')
        except IOError:
            print("\n **ERROR**: ISQReader: intput file ' % s' not found!\n\n" % File)
            print('\n E N D E D  with ERRORS \n\n')

        for Index in np.arange(0, 200, 4):
            f.seek(Index)
            #print('Index %s :          %s' % (Index, struct.unpack('i', f.read(4))[0]))
            f.seek(Index)

        f.seek(32)
        CT_ID = struct.unpack('i', f.read(4))[0]

        if CT_ID != 6020:
            print('!!! unknown muCT -> no Slope and Intercept known !!!')

        f.seek(28)
        #    sample_nb = struct.unpack('i', f.read(4))[0]

        f.seek(108)
        Scanning_time = struct.unpack('i', f.read(4))[0] / 1000

        f.seek(168)
        Energy = struct.unpack('i', f.read(4))[0] / 1000.

        f.seek(172)
        Current = struct.unpack('i', f.read(4))[0]

        f.seek(44)
        X_pixel = struct.unpack('i', f.read(4))[0]

        f.seek(48)
        Y_pixel = struct.unpack('i', f.read(4))[0]

        f.seek(52)
        Z_pixel = struct.unpack('i', f.read(4))[0]

        f.seek(56)
        Res_General_X = struct.unpack('i', f.read(4))[0]
        #print('Resolution general X in mu: ', Res_General_X)

        f.seek(60)
        Res_General_Y = struct.unpack('i', f.read(4))[0]
        #print('Resolution general Y in mu: ', Res_General_Y)

        f.seek(64)
        Res_General_Z = struct.unpack('i', f.read(4))[0]
        #print('Resolution general Z in mu: ', Res_General_Z)

        Res_X = Res_General_X / float(X_pixel)
        Res_Y = Res_General_Y / float(Y_pixel)
        Res_Z = Res_General_Z / float(Z_pixel)

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
        #print(Header)

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

        #print('\nReshape data')
        #Tic = time.time()

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

            if self.Echo:
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


        if self.Echo:
            Time.Process(0,Text)
            print('\nScanner ID:                 ', CT_ID)
            print('Scanning time in ms:        ', Scanning_time)
            print('Energy in keV:              ', Energy)
            print('Current in muA:             ', Current)
            print('Nb X pixel:                 ', X_pixel)
            print('Nb Y pixel:                 ', Y_pixel)
            print('Nb Z pixel:                 ', Z_pixel)
            print('Pixel resolution X in mu:    %.2f' % Res_X)
            print('Pixel resolution Y in mu:    %.2f' % Res_Y)
            print('Pixel resolution Z in mu:    %.2f' % Res_Z)


        if CT_ID == 6020 and BMD is True:
            # BE CAREFULL, THIS IS FOR BMD CONVERSION:
            if self.Echo:
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
class Write():

    def __init__(self):
        self.Echo = True
        self.FName = 'Image'
        pass

    def Raw(self, Image, PixelType):

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

        File = np.memmap(self.FName + '.raw', dtype=CastedImageArray.dtype, mode='w+', shape=CastedImageArray.shape, order='C')
        File[:] = CastedImageArray[:]
        del File

        return

    def MHD(self, Image, PixelType='uint'):

        if self.Echo:
            Text = 'Write MHD'
            Time.Process(1, Text)

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

        outs = open(self.FName + '.mhd', 'w')
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

        if '\\' in self.FName:
            Fname = self.FName.split('\\')[-1]
        elif '/' in self.FName:
            Fname = self.FName.split('/')[-1]
        else:
            Fname = self.FName
        outs.write('ElementDataFile = %s\n' % (Fname + '.raw'))

        if PixelType == 'norm':
            Array = sitk.GetArrayFromImage(Image)
            outs.write('\n# Min Max = %i %i\n' % (Array.min(), Array.max()))

        outs.close()

        self.Raw(Image, PixelType)

        if self.Echo:
            Time.Process(0, Text)

        return

    def VectorFieldVTK(self, VectorField, SubSampling=1):

        if self.Echo:
            Text = 'Write VTK'
            Time.Process(1, Text)

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
            z = np.arange(0, Size[0], SubSampling)
            y = np.arange(0, Size[1], SubSampling)
            x = np.arange(0, Size[2], SubSampling)
            NumberOfElements = len(x) * len(y) * len(z)

            # Build vector arrays
            u = VectorField[:, :, :, 0]
            v = VectorField[:, :, :, 1]
            w = VectorField[:, :, :, 2]

        File = open(self.FName + '.vtk','w')

        # ASCII file header
        File.write('# vtk DataFile Version 4.2\n')
        File.write('VTK from Python\n')
        File.write('ASCII\n\n')
        File.write('DATASET RECTILINEAR_GRID\n')
        File.write('DIMENSIONS ' + str(len(x)) + ' ' + str(len(y)) + ' ' + str(len(z)) + '\n\n')
        File.close()

        # Append ascii x,y,z
        MaxLineWidth = Size[2]*Size[1]*Size[0]
        File = open(self.FName + '.vtk','a')
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
        File = open(self.FName + '.vtk','a')
        File.write('\nPOINT_DATA ' + str(NumberOfElements) + '\n\n')
        File.write('VECTORS Deformation float\n')
        File.close()

        # Append ascii u,v,w and build deformation magnitude array
        Magnitude = np.zeros(VectorField.shape[:-1])
        File = open(self.FName + '.vtk','a')

        if Dimension == 2:
            for j in range(0,Size[1],SubSampling):
                if self.Echo:
                    Time.Update(j / Size[1], 'Write Vectors')
                for i in range(0,Size[2],SubSampling):
                    Magnitude[j,i] = np.sqrt(u[j,i]**2 + v[j,i]**2 + w[j,i]**2)
                    File.write(str(u[j,i]) + ' ' + str(v[j,i]) + ' ' + str(w[j,i]) + ' ')

        elif Dimension == 3:
            for k in range(0,Size[0],SubSampling):
                if self.Echo:
                    Time.Update(k / Size[0], 'Write Vectors')
                for j in range(0,Size[1],SubSampling):
                    for i in range(0,Size[2],SubSampling):
                        Magnitude[k,j,i] = np.sqrt(u[k,j,i] ** 2 + v[k,j,i] ** 2 + w[k,j,i] ** 2)
                        File.write(str(u[k,j,i]) + ' ' + str(v[k,j,i]) + ' ' + str(w[k,j,i]) + ' ')

        File.close()

        # Append another subheader
        File = open(self.FName + '.vtk', 'a')
        File.write('\n\nSCALARS Magnitude float\n')
        File.write('LOOKUP_TABLE default\n')

        if Dimension == 2:
            for j in range(0, Size[1], SubSampling):
                if self.Echo:
                    Time.Update(j / Size[1], 'Write Magnitudes')
                for i in range(0, Size[2], SubSampling):
                    File.write(str(Magnitude[j,i]) + ' ')

        elif Dimension == 3:
            for k in range(0, Size[0], SubSampling):
                if self.Echo:
                    Time.Update(k / Size[0], 'Write Magnitudes')
                for j in range(0, Size[1], SubSampling):
                    for i in range(0, Size[2], SubSampling):
                        File.write(str(Magnitude[k,j,i]) + ' ')

        File.close()

        if self.Echo:
            Time.Process(0,Text)

        return

    def FabricVTK(self, eValues, eVectors, nPoints=32, Scale=1, Origin=(0,0,0)):

        if self.Echo:
            Text = 'Write Fabric VTK'
            Time.Process(1, Text)

        ## New coordinate system
        Q = np.array(eVectors)

        ## Build data for fabric plotting
        u = np.arange(0, 2 * np.pi + 2 * np.pi / nPoints, 2 * np.pi / nPoints)
        v = np.arange(0, np.pi + np.pi / nPoints, np.pi / nPoints)
        X = eValues[0] * np.outer(np.cos(u), np.sin(v))
        Y = eValues[1] * np.outer(np.sin(u), np.sin(v))
        Z = eValues[2] * np.outer(np.ones_like(u), np.cos(v))
        nNorm = np.zeros(X.shape)

        for i in range(len(X)):
            for j in range(len(X)):
                [X[i, j], Y[i, j], Z[i, j]] = np.dot([X[i, j], Y[i, j], Z[i, j]], Q)
                n = np.array([X[i, j], Y[i, j], Z[i, j]])
                nNorm[i, j] = np.linalg.norm(n)

        # Scale the arrays
        X = Scale/2 * X
        Y = Scale/2 * Y
        Z = Scale/2 * Z

        # Translate origin to the center of ROI
        X = Origin[2] + X + Scale/2
        Y = Origin[1] + Y + Scale/2
        Z = Origin[0] + Z + Scale/2

        # Write VTK file
        VTKFile = open(self.FName + '.vtk', 'w')

        # Write header
        VTKFile.write('# vtk DataFile Version 4.2\n')
        VTKFile.write('VTK from Python\n')
        VTKFile.write('ASCII\n')
        VTKFile.write('DATASET UNSTRUCTURED_GRID\n')

        # Write points coordinates
        Points = int(X.shape[0] * X.shape[1])
        VTKFile.write('\nPOINTS ' + str(Points) + ' floats\n')
        for i in range(Points):
            if self.Echo:
                Time.Update(i / len(Points), 'Write Points')
            VTKFile.write(str(X.reshape(Points)[i].round(3)))
            VTKFile.write(' ')
            VTKFile.write(str(Y.reshape(Points)[i].round(3)))
            VTKFile.write(' ')
            VTKFile.write(str(Z.reshape(Points)[i].round(3)))
            VTKFile.write('\n')

        # Write cells connectivity
        Cells = int(nPoints**2)
        ListSize = int(Cells*5)
        VTKFile.write('\nCELLS ' + str(Cells) + ' ' + str(ListSize) + '\n')

        ## Add connectivity of each cell
        Connectivity = np.array([0, 1])
        Connectivity = np.append(Connectivity,[nPoints+2,nPoints+1])

        for i in range(Cells):
            if self.Echo:
                Time.Update(i / len(Cells), 'Write Cells')

            VTKFile.write('4')

            if i > 0 and np.mod(i,nPoints) == 0:
                Connectivity = Connectivity + 1

            for j in Connectivity:
                VTKFile.write(' ' + str(j))
            VTKFile.write('\n')

            ## Update connectivity
            Connectivity = Connectivity+1

        # Write cell types
        VTKFile.write('\nCELL_TYPES ' + str(Cells) + '\n')
        for i in range(Cells):
            VTKFile.write('9\n')

        # Write MIL values
        VTKFile.write('\nPOINT_DATA ' + str(Points) + '\n')
        VTKFile.write('SCALARS MIL float\n')
        VTKFile.write('LOOKUP_TABLE default\n')

        for i in range(nPoints+1):
            if self.Echo:
                Time.Update(i / len(nPoints+1), 'Write Values')
            for j in range(nPoints+1):
                VTKFile.write(str(nNorm.reshape(Points)[j + i * (nPoints+1)].round(3)))
                VTKFile.write(' ')
            VTKFile.write('\n')
        VTKFile.close()

        if self.Echo:
            Time.Process(0, Text)

        return

Write = Write()
#%% Registration funtions
class Registration():

    def __init__(self):
        self.Echo = True

    def Register(self, FixedImage, MovingImage, Type, FixedMask=None, MovingMask=None, Path=None, Dictionary={}):

        if self.Echo:
            Text = Type + ' reg.'
            Time.Process(1, Text)
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

        EIF.Execute()

        # Get results
        Result_Image = EIF.GetResultImage()
        TransformParameters = EIF.GetTransformParameterMap()

        # Print elapsed time
        if self.Echo:
            Time.Process(0, Text)

        return Result_Image, TransformParameters
        
    def ComputeInverse(self, FixedImage, TPMFileName, MovingImage=None, Path=None):

        """
        Compute inverse of rigid elastix transform. Manual 6.1.6
        """

        if self.Echo:
            Text = 'Inverse reg.'
            Time.Process(1, Text)

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
            Time.Process(0, Text)

        return InvertedTransform

    def Apply(self, Image,TransformParameterMap,Path=None,Jacobian=None):

        """
        Apply transform parametermap from elastix to an image
        """

        if self.Echo:
            Text = 'Apply Transform'
            Time.Process(1, Text)

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
            Time.Process(0, Text)

        return ResultImage
 
    def ApplyCustom(self, Image, TransformParameterMap):

        """
        Apply costum parameter map to an image
        """

        if self.Echo:
            Text = 'Custom transform'
            Time.Process(1, Text)

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

        if self.Echo:
            Time.Process(0, Text)

        return Image_T
           
    def ApplyInverse(self, Image, TransformParameterMap):

        """
        Apply inverse rigid transform from transform parameter map
        """

        if self.Echo:
            Text = 'Inverse transform'
            Time.Process(1, Text)

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

        if self.Echo:
            Time.Process(0, Text)

        return Image_T

Registration = Registration()
#%% Signal treatment functions
class Signal():

    def __init__(self):
        pass

    def FFT(self, Signal, Sampling, Show=False):

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

    def DesignFilter(self, Frequency, Order=2):

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

    def Filter(self, Signal, Sampling, Frequency, Order=2, Show=False):
        
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

    def MaxSlope(self, X, Y=[], WindowWidth=1, StepSize=1):

        if len(Y) == 0:
            Y = X.copy()
            X = np.arange(len(X))

        Slope = NumbaMaxSlope(np.array(X), np.array(Y), int(WindowWidth), StepSize)

        return Slope

@njit
def NumbaMaxSlope(Xdata, Ydata, WindowWidth, StepSize):

    Slopes = []
    Iterations = round((len(Xdata) - WindowWidth) / StepSize - 1)
    
    for i in range(Iterations):

        Start = i * StepSize
        Stop = Start + WindowWidth
        XPoints = Xdata[Start:Stop] 
        YPoints = Ydata[Start:Stop] 
        N = len(XPoints)

        X = np.ones((N,2))
        X[:,1] = XPoints
        Y = np.reshape(YPoints, (N,1))

        X1 = np.linalg.inv(np.dot(X.T, X))
        Intercept, Coefficient = np.dot(X1, np.dot(X.T, Y))
        Slopes.append(Coefficient)

    Slope = max(Slopes)[0]

    return Slope

Signal = Signal()
#%% Abaqus functions
class Abaqus():

    def __init__(self):

        self.NLGEOM = 'YES'
        self.MaxINC = 1000

        self.StepText = """
**
*STEP,AMPLITUDE=RAMP,UNSYMM=YES,INC={MAXINC},NLGEOM={NL}
***********************************************************
**               INCLUDE
***********************************************************
*INCLUDE, input={BCsFile}
***********************************************************
*OUTPUT,FIELD
*ELEMENT OUTPUT, POSITION=CENTROIDAL
SDV2,
SDV15,
SDV16,
SDV17,
SDV18,
SDV22,
SDV23,
SDV24,
SDV25,
SDV26,
SDV27,
SDV28,
SDV29,
SDV30,
SDV31,
S,
LE,
COORD,
*NODE OUTPUT
U,
RF,
CF,
**
*OUTPUT, HISTORY
*NODE OUTPUT, NSET=REF_NODE
U,
RF,
*NODE PRINT, NSET=REF_NODE, FREQUENCY=1, SUMMARY=NO
U,
RF
CF
*END STEP"""

        return

    def ReadDAT(self, File):

        """
        Read .dat file from abaqus and extract reference point data
        """

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

                ColNames = []
                for i in range(3):
                    for V in ['U', 'F']:
                        ColNames.append(V + str(i+1))
                for i in range(3):
                    for V in ['R', 'M']:
                        ColNames.append(V + str(i+1))

                Data = pd.DataFrame()
                for iName, Name in enumerate(ColNames):
                    Data[Name] = np.concatenate([[0], Values[:,iName]])
                
                Data.columns = ['X', 'FX', 'Y', 'FY', 'Z', 'FZ', 'Phi', 'MX', 'Theta', 'MY', 'Psi', 'MZ']

                Data['Step'] = np.concatenate([[0],Steps])
                Data['Increment'] = np.concatenate([[0],Increments])

            return Data

        except FileNotFoundError:
            print('File' + File + 'does not exist')

            return

    def WriteRefNodeBCs(self, FileName, DOFs, Values, BCType='DISPLACEMENT', Parameters=[0.05, 1, 5e-05, 0.05]):

        """
        Write boundary conditions for a reference node called REF_NODE
        :param FileName: Name of the BCs file
                         Type: str
        :param DOFs: Degrees of freedom to constrain
                     Type: List of ints
        :param Values: Constrain value for each DOF listed
                       Type: List of floats
        :param BCType: Type of boundary condition
                       Type: str
                       Either 'DISPLACEMENT' (Dirichlet)
                       or 'FORCE' (Neumann)
        :param Parameters: Step parameters
                           P1: Start step size
                           P2: Time for displacement
                           P3: Min step size
                           P4: Max step size
                           Type: List of floats

        """

        with open(FileName, 'w') as F:
            F.write('*STATIC\n')

            F.write(str(Parameters[0])+ ', ')
            F.write(str(Parameters[1])+ ', ')
            F.write(str(Parameters[2])+ ', ')
            F.write(str(Parameters[3])+ '\n')

            F.write('*BOUNDARY, TYPE=' + BCType + '\n')

            for DOF, Value in zip(DOFs, Values):
                Line = str(DOF) + ', ' + str(DOF) + ', ' + str(Value)
                F.write('REF_NODE, ' + Line + '\n')

        return

    def AddStep(self, File, Step, DOFs, Values, BCType='DISPLACEMENT', Parameters=[0.05, 1, 5e-05, 0.05]):

        self.WriteRefNodeBCs(Step, DOFs, Values, BCType, Parameters)

        Context = {'MAXINC':self.MaxINC,
                   'NL':self.NLGEOM,
                   'BCsFile':Step}
        Text = self.StepText.format(**Context)

        with open(File, 'a') as F:
            F.write(Text)

        return

    def RemoveSteps(self, File, NSteps=1):

        with open(File, 'r') as F:
            Text = F.read()

        if Text.find('*STEP') == -1:
            print('\nNo step in file!')

        else:
            Start = 0
            Index = Text.find('*STEP', Start)
            Indices = []
            while bool(Index+1):
                Indices.append(Index)
                Start = Index + len('*STEP')
                Index = Text.find('*STEP', Start)

            if NSteps > len(Indices):
                Index = Indices[0]
            else:
                Index = np.array(Indices)[-NSteps]
                
            with open(File, 'w') as F:
                F.write(Text[:Index - len('*STEP') + 1])

        return

    def WriteUMAT(self, UMATPath, Name='Elastic'):

        if Name == 'Elastic':
            Text = """C================================================================================================
C
C                          UMAT FOR ISOTROPIC ELASTIC CONSTITUTIVE LAW
C
C     FOUND ON https://simplifiedfem.wordpress.com/about/tutorial-write-a-simple-umat-in-abaqus/
C
C                  ADAPTED BY MATHIEU SIMON - ISTB - UNIVERSITY OF BERN - 2021
C
C================================================================================================
C
C
    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
    1 RPL,DDSDDT,DRPLDE,DRPLDT,
    2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
    3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
    4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
C
    INCLUDE 'ABA_PARAM.INC'
C
C
    CHARACTER*80 CMNAME
    DIMENSION DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
    1 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
    2 PROPS(NPROPS),COORDS(3),DROT(3,3), JSTEP(4)
C     
C     
    INTEGER NSTATV, KSTEP, KINC
    DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
    DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)
C
C
C     ELASTIC USER SUBROUTINE
    PARAMETER (ONE=1.0D0, TWO=2.0D0)
C
C
C    Get engineering variables
        E=PROPS(1)
        ANU=PROPS(2)
C
C
C    Compute lamé parameters
        ALAMBDA=E/(ONE+ANU)/(ONE-TWO*ANU)
        BLAMBDA=(ONE-ANU)
        CLAMBDA=(ONE-TWO*ANU)
C
C
C    Initialize material stiffness matrix
            DO I=1,NTENS
            DO J=1,NTENS
            DDSDDE(I,J)=0.0D0
            ENDDO
            ENDDO
C
C
C    Update stiffness matrx
                DDSDDE(1,1)=(ALAMBDA*BLAMBDA)
                DDSDDE(2,2)=(ALAMBDA*BLAMBDA)
                DDSDDE(3,3)=(ALAMBDA*BLAMBDA)
                DDSDDE(4,4)=(ALAMBDA*CLAMBDA)
                DDSDDE(5,5)=(ALAMBDA*CLAMBDA)
                DDSDDE(6,6)=(ALAMBDA*CLAMBDA)
                DDSDDE(1,2)=(ALAMBDA*ANU)
                DDSDDE(1,3)=(ALAMBDA*ANU)
                DDSDDE(2,3)=(ALAMBDA*ANU)
                DDSDDE(2,1)=(ALAMBDA*ANU)
                DDSDDE(3,1)=(ALAMBDA*ANU)
                DDSDDE(3,2)=(ALAMBDA*ANU)
C
C
C     Update stress tensor
        DO I=1,NTENS
            DO J=1,NTENS
            STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
            ENDDO
        ENDDO
C
C        
C     Deformation gradient
    STATEV(1) = DFGRD1(1,1)
    STATEV(2) = DFGRD1(1,2)
    STATEV(3) = DFGRD1(1,3)
    STATEV(4) = DFGRD1(2,1)
    STATEV(5) = DFGRD1(2,2)
    STATEV(6) = DFGRD1(2,3)
    STATEV(7) = DFGRD1(3,1)
    STATEV(8) = DFGRD1(3,2)
    STATEV(9) = DFGRD1(3,3)
C
C
        RETURN
        END"""

        with open(str(Path(UMATPath, (Name + '.f'))), 'w') as File:
            File.write(Text)

    def uFE(self, BinImage, FileName, UMAT='Elastic', Cap=False):

        Text = 'Create uFE file'
        Time.Process(1,Text)

        # Define directory and parameters
        Compression = 0.001 # Relative compression
        E = 1E4             # Young's modulus
        Nu = 0.3            # Poisson's ratio

        # Define input file name and create it
        if os.name == 'nt':
            Last = FileName[::-1].find('\\')
            NProc = 4
        else:
            Last = FileName[::-1].find('/')
            NProc = 12
        UMATPath = FileName[:-Last]
        self.WriteUMAT(UMATPath, UMAT)

        # Compute nodes coordinates
        BinArray = sitk.GetArrayFromImage(BinImage)
        Spacing = np.round(BinImage.GetSpacing(),3)

        if len(np.unique(BinArray)) > 2:
            print('Image is not binary')
            return
        else:
            BinArray = BinArray.astype(int)
        
        # Add bottom and top cap
        if Cap:
            CapSize = 5
            BinArray = np.pad(BinArray, ((CapSize,CapSize), (0,0), (0,0)), constant_values=(125))
            BinArray = BinArray.astype('uint8')

        # Perform elements and nodes mapping
        Time.Update(1/6,'Perform mapping')
        Nodes, Coords, Elements, ElementsNodes = NumbaMapping(BinArray)
        NodesNeeded = np.unique(ElementsNodes[BinArray.astype(bool)])

        ElementsNodes = ElementsNodes[BinArray.astype(bool)]
        Elements = Elements[BinArray.astype(bool)]
        Materials = BinArray[BinArray.astype(bool)]
        Coords = Coords[np.isin(Nodes,NodesNeeded)]
        Nodes = Nodes[np.isin(Nodes,NodesNeeded)]

        # Sort nodes according to coordinates
        Indices = np.lexsort((-Coords[:,1],-Coords[:,2],-Coords[:,0]))
        Coords = np.round(Coords[Indices] * Spacing,3)
        Nodes = Nodes[Indices]

        # Sort elements according to their number
        Indices = np.argsort(Elements)
        Elements = Elements[Indices]
        ElementsNodes = ElementsNodes[Indices]
        Materials = Materials[Indices]

        # Identify top and bottom nodes
        Time.Update(2/6,'Find top-bottom')
        Margin = 0.0
        BottomNodes = [Node[0] + 1 for Node in np.argwhere(Coords[:,2] <= Coords[:,2].min()+Margin)]
        TopNodes = [Node[0] + 1 for Node in np.argwhere(Coords[:,2] >= Coords[:,2].max()-Margin)]

        # Identify top and bottom elements
        BottomElements = np.isin(ElementsNodes,Nodes[np.array(BottomNodes)-1])
        BottomElements = [Element[0] + 1 for Element in np.argwhere(BottomElements.any(axis=1))]
        TopElements = np.isin(ElementsNodes,Nodes[np.array(TopNodes)-1])
        TopElements = [Element[0] + 1 for Element in np.argwhere(TopElements.any(axis=1))]

        # Compute boundary conditions
        Displacement = round((Coords[:,2].max() - Coords[:,2].min()) * Compression,3)

        # Generate nodes text
        Time.Update(3/6,'Write nodes text')
        NodesStr = (np.arange(len(Nodes)) + 1).astype('<U32')
        CoordsStr = Coords.astype('<U32')
        Sep = np.array([', ', '\n']).astype('<U32')

        Ns = np.array_split(NodesStr,NProc)
        Cs = np.array_split(CoordsStr, NProc)
        Processes = [[N, C, Sep] for N, C in zip(Ns, Cs)]
        NText = []
        with Pool(processes=NProc) as P:
            for TaskResult in P.map(NodesText, Processes):
                NText.append(TaskResult)
        NText = ''.join(NText)

        # Generate element text
        Time.Update(4/6,'Write elem. text')
        EN = (np.arange(len(Elements)) + 1).astype('<U32')
        NS = np.argsort(Nodes)
        SM = np.searchsorted(Nodes[NS], ElementsNodes)
        SEN = (NS[SM] + 1).astype('<U32')

        ENs = np.array_split(EN,NProc)
        SENs = np.array_split(SEN, NProc)
        Processes = [[E, S, Sep] for E, S in zip(ENs, SENs)]
        EText = []
        with Pool(processes=NProc) as P:
            for TaskResult in P.map(ElementsText, Processes):
                EText.append(TaskResult)
        EText = ''.join(EText)

        # Write file
        Time.Update(5/6,'Write input file')
        with open(FileName,'w') as File:

            # Write heading
            File.write('*Heading\n')
            File.write('** Job name: ' + FileName[-Last:] + ' Model name: Bone\n')
            File.write('** Generated by: Python\n')
            File.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
            File.write('**\n')

            # Write parts
            File.write('** PARTS\n')
            File.write('*Part, name=SAMPLE\n')

            # Write nodes
            File.write('**\n')
            File.write('*Node\n')
            File.write(NText)

            # Write elements
            File.write('**\n')
            File.write('*Element, type=C3D8\n')
            File.write(EText)

            # Write node set
            File.write('**\n')
            File.write('*Nset, nset=NODESET, generate\n')
            File.write('1,  ' + str(len(Coords)) + ',  1\n')

            # Write bone elements set
            File.write('**\n')
            File.write('*Elset, elset=BONE\n')
            Line = ''
            N = 0
            for Index, Material in enumerate(Materials):
                if Material == 255:
                    Index += 1
                    N += 1
                    if np.mod((N),16) == 0:
                        File.write(Line[:-2] + '\n')
                        Line = str(Index) + ', '
                    else:
                        Line += str(Index) + ', '
            File.write(Line[:-2] + '\n')

            # Write caps element set
            if Cap:
                File.write('**\n')
                File.write('*Elset, elset=CAPS\n')
                Line = ''
                N = 0
                for Index, Material in enumerate(Materials):
                    if Material == 125:
                        Index += 1
                        N += 1
                        if np.mod((N),16) == 0:
                            File.write(Line[:-2] + '\n')
                            Line = str(Index) + ', '
                        else:
                            Line += str(Index) + ', '
                File.write(Line[:-2] + '\n')

            # Write section
            File.write('**\n')
            File.write('*Solid Section, elset=BONE, material=' + UMAT + '\n')
            if Cap:
                File.write('*Solid Section, elset=CAPS, material=Steel\n')
            File.write('*End Part\n')

            # Write assembly
            File.write('**\n')
            File.write('** ASSEMBLY\n')
            File.write('*Assembly, name=Assembly\n')

            ## Write instance
            File.write('**\n')
            File.write('*Instance, name=SAMPLE, part=SAMPLE\n')
            File.write('*End Instance\n')

            # Write bottom nodes set
            File.write('**\n')
            File.write('*Nset, nset=BOTTOMNODES, instance=SAMPLE\n')
            Line = ''
            for iBN, BN in enumerate(BottomNodes):
                if np.mod((iBN+1),16) == 0:
                    File.write(Line[:-2] + '\n')
                    Line = str(BN) + ', '
                else:
                    Line += str(BN) + ', '
            File.write(Line[:-2] + '\n')

            # Write bottom elements set
            File.write('**\n')
            File.write('*Elset, elset=BOTTOMELEMENTS, instance=SAMPLE\n')
            Line = ''
            for iBE, BE in enumerate(BottomElements):
                if np.mod((iBE+1),16) == 0:
                    File.write(Line[:-2] + '\n')
                    Line = str(BE) + ', '
                else:
                    Line += str(BE) + ', '
            File.write(Line[:-2] + '\n')

            # Write top nodes set
            File.write('**\n')
            File.write('*Nset, nset=TOPNODES, instance=SAMPLE\n')
            Line = ''
            for iTN, TN in enumerate(TopNodes):
                if np.mod((iTN+1),16) == 0:
                    File.write(Line[:-2] + '\n')
                    Line = str(TN) + ', '
                else:
                    Line += str(TN) + ', '
            File.write(Line[:-2] + '\n')

            # Write top elements set
            File.write('**\n')
            File.write('*Elset, elset=TOPELEMENTS, instance=SAMPLE\n')
            Line = ''
            for iTE, TE in enumerate(TopElements):
                if np.mod((iTE+1),16) == 0:
                    File.write(Line[:-2] + '\n')
                    Line = str(TE) + ', '
                else:
                    Line += str(TE) + ', '
            File.write(Line[:-2] + '\n')

            File.write('**\n')
            File.write('*End Assembly\n')

            # Write materials
            File.write('**\n')
            File.write('** MATERIALS\n')
            File.write('*Material, name=' + UMAT + '\n')
            File.write('*User Material, constants=2\n')
            File.write(str(E) + ', ' + str(Nu) + '\n')

            # Write dependent variables
            File.write('*DEPVAR\n')
            File.write('9\n')
            File.write('1, F11,\n')
            File.write('2, F12,\n')
            File.write('3, F13,\n')
            File.write('4, F21,\n')
            File.write('5, F22,\n')
            File.write('6, F23,\n')
            File.write('7, F31,\n')
            File.write('8, F32,\n')
            File.write('9, F33\n')
            File.write('**\n')

            if Cap:
                File.write('*Material, name=Steel\n')
                File.write('*ELASTIC\n')
                File.write('100000.0, 0.45\n')
                File.write('**\n')

            # Write step
            File.write('** STEP: Compression\n')
            File.write('*Step, name=Compression, nlgeom=YES\n')     # nlgeom=YES for F different from identity
            File.write('Static compression\n')
            File.write('*Static\n')
            File.write('1., 1., 1e-05, 1.\n')

            ## Write boundary conditions
            File.write('**\n')
            File.write('** BOUNDARY CONDITIONS\n')
            File.write('*Boundary\n')
            File.write('BOTTOMNODES, 1, 3, 0\n')
            File.write('*Boundary\n')
            File.write('TOPNODES, 1, 2, 0\n')
            File.write('*Boundary\n')
            File.write('TOPNODES, 3, 3, -' + str(Displacement) + '\n')

            ## Write output requests
            File.write('**\n')
            File.write('** OUTPUT REQUESTS\n')
            File.write('*Restart, write, frequency=0\n')
            File.write('** FIELD OUTPUT: F-Output-1\n')
            File.write('**\n')
            File.write('*Output, field\n')
            File.write('*Element Output\n')
            File.write('S,\n')
            File.write('SDV1,\n')
            File.write('SDV2,\n')
            File.write('SDV3,\n')
            File.write('SDV4,\n')
            File.write('SDV5,\n')
            File.write('SDV6,\n')
            File.write('SDV7,\n')
            File.write('SDV8,\n')
            File.write('SDV9,\n')
            File.write('*Node Output\n')
            File.write('U,\n')
            File.write('RF\n')
            File.write('**El Print\n')
            File.write('**DG\n')
            File.write('*Output, history, frequency=0\n')
            File.write('*End Step\n')

        Time.Process(0,Text)

        return

    def RunSim(self, WorkingDir, InputFile, UMAT, nCPUs=12):

        Text = 'Run ' + InputFile
        Time.Process(1, Text)

        # Set working directory
        os.chdir(WorkingDir)

        # Absolutely necessary to start abaqus job
        try:
            os.environ.pop('PYTHONIOENCODING')
        except KeyError:
            pass

        # Write script
        Script = '''
    #!/bin/bash
    abaqus interactive job={Job} inp={InputFile} user={UMAT} cpus={nCPUs} ask_delete=OFF
    '''
        
        Context = {'Job':InputFile[:-4],
                'InputFile':InputFile,
                'UMAT':UMAT,
                'nCPUs':nCPUs}

        FileName = 'Run.sh'
        with open(FileName, 'w') as File:
            File.write(Script.format(**Context))

        # Run simulation with script
        try:
            sh.bash(str(Path(WorkingDir, FileName))) # type: ignore
            Completed = True
        except:
            Completed = False
            print('Analysis not completed')

        if Completed:
            # Remove unnecessary files
            os.remove(InputFile[:-4] + '.com')
            os.remove(InputFile[:-4] + '.msg')
            os.remove(InputFile[:-4] + '.dat')
            os.remove(InputFile[:-4] + '.prt')
            os.remove(InputFile[:-4] + '.sim')
            os.remove(InputFile[:-4] + '.sta')

        # Print time
        Time.Process(0)
        return

    def ReadODB(self, WorkingDir, OdbFile, Variables=['U','RF'], Step=0, Frame=-1, NodeSet=False):

        # Change working directory
        os.chdir(WorkingDir)

        # Write odb reader
        with open('ReadOdb.py','w') as File:

            # Write heading and initial part
            Text = r"""# ReadODB.py
# A script to read deformation gradient from odb file from ABAQUS.

import csv
import sys
import numpy as np
from odbAccess import *

print \'\nOpen odb file\'

Odb = openOdb(\'{File}\')

# Create variable that refers to the last frame of the step.
Steps = Odb.steps.keys()
Frame = Odb.steps[Steps[{Step}]].frames[{Frame}]

# Create variable refering to model instance
Instances = Odb.rootAssembly.instances.keys()
Instance = Odb.rootAssembly.instances[Instances[0]]

"""
            File.write(Text.format(**{'File':OdbFile + '.odb',
                                      'Step':Step,
                                      'Frame':Frame}))
                    
            # Select fields outputs
            if 'U' in Variables:
                File.write('Displacements = Frame.fieldOutputs[\'U\'].values\n')
                
            if 'RF' in Variables:
                File.write('Forces = Frame.fieldOutputs[\'RF\'].values\n')
                
            if 'F' in Variables:
                File.write("""F11 = Frame.fieldOutputs['SDV_F11']
F12 = Frame.fieldOutputs['SDV_F12']
F13 = Frame.fieldOutputs['SDV_F13']
F21 = Frame.fieldOutputs['SDV_F21']
F22 = Frame.fieldOutputs['SDV_F22']
F23 = Frame.fieldOutputs['SDV_F23']
F31 = Frame.fieldOutputs['SDV_F31']
F32 = Frame.fieldOutputs['SDV_F32']
F33 = Frame.fieldOutputs['SDV_F33']

F = [F11,F12,F13,F21,F22,F23,F31,F32,F33]
F_Names = ['F11','F12','F13','F21','F22','F23','F31','F32','F33']
""")

            # Select nodes sets
            if NodeSet:
                Line = """\nNodes = odb.rootAssembly.nodeSets['{NodeSet}'].nodes[0]\n"""
                File.write(Line.format(**{'NodeSet':NodeSet}))
            else:
                File.write('Nodes = Instance.nodes\n')

            # Store nodes values
            if 'U' in Variables or 'RF' in Variables:
                File.write("""NodesData = []
for Node in Nodes:
    Cx = float(Node.coordinates[0])
    Cy = float(Node.coordinates[1])
    Cz = float(Node.coordinates[2])
""")
                Line = '    NodesData.append([Cx, Cy, Cz'

                if 'U' in Variables:
                    File.write("""    Ux = float(Displacements[Node.label].data[0])
    Uy = float(Displacements[Node.label].data[1])
    Uz = float(Displacements[Node.label].data[2])
""")
                    Line += ', Ux, Uy, Uz'

                if 'RF' in Variables:
                    File.write("""    Fx = float(Forces[Node.label].data[0])
    Fy = float(Forces[Node.label].data[1])
    Fz = float(Forces[Node.label].data[2])
""")
                    Line += ', Fx, Fy, Fz'
                
                File.write(Line + '])\n')
                File.write('NodesData = np.array(NodesData)\n')

            # For deformation gradient analysis (to clean/adapt)
            if 'F' in Variables:
                File.write(r"""NodeLabels = np.array([])
for iNode, Node in enumerate(Nodes):
    NodeLabels = np.append(NodeLabels,Node.label)

# Initialize loop
N = len(Instance.elements)
ElementsDG = np.zeros((12,N))

for ElementNumber, Element in enumerate(Instance.elements[:N]):
    sys.stdout.write("\r" + "Read element " + str(ElementNumber+1) + "/" + str(N))
    sys.stdout.flush()

    # Compute element central position
    XYZ = np.zeros((len(Element.connectivity),3))
    NodeNumber = 0
    for NodeLabel in Element.connectivity:
        Node = np.where(NodeLabels==NodeLabel)[0][0]
        for Axis in range(3):
            XYZ[NodeNumber,Axis] = Nodes[Node].coordinates[Axis]
        NodeNumber += 1

    # Get element mean deformation gradient
    F_IntegrationPoints = np.zeros((1,9))

    for F_Component in range(9):
        F_Value = F[F_Component].getSubset(region=Element).values[0].data
        F_IntegrationPoints[0,F_Component] = F_Value

    # Add data to arrays and increment
    ElementsDG[:3,ElementNumber] = np.mean(XYZ,axis=0)
    ElementsDG[3:,ElementNumber] = F_IntegrationPoints
""")

            # Write into csv file
            File.write('\nprint  \'' + r'\n' + 'Save to csv\'')

            if 'U' in Variables or 'RF' in Variables:
                File.write("""
with open('{File}_Nodes.csv', 'w') as File:
    Writer = csv.writer(File)
    for Row in NodesData:
        Writer.writerow(Row)  
""".format(**{'File':OdbFile}))

            if 'F' in Variables:
                File.write("""
with open('{File}_{Step}_{Frame}_DG.csv', 'w') as File:
    Writer = csv.writer(File)
    for Row in ElementsDG.T:
        Writer.writerow(Row)            
""".format(**{'File':OdbFile, 'Step':Step, 'Frame':Frame}))

            # Close odb file
            File.write('Odb.close()')

        # Run odb reader
        os.system('abaqus python ReadOdb.py')

        # Collect results
        Data = []
        if 'U' in Variables or 'RF' in Variables:
            ND = pd.read_csv(str(Path(WorkingDir, OdbFile + '_Nodes.csv')), header=None)
            if 'RF' not in Variables:
                ND.columns = ['Cx','Cy','Cz','Ux','Uy','Uz']
            if 'U' not in Variables:
                ND.columns = ['Cx','Cy','Cz','Fx','Fy','Fz']
            else:
                ND.columns = ['Cx','Cy','Cz','Ux','Uy','Uz','Fx','Fy','Fz']
            Data.append(ND)

        if 'F' in Variables:
            DG = pd.read_csv(OdbFile + '_' + str(Step) + '_' + str(Frame) + '_DG.csv', header=None)
            Data.append(DG)

        return Data

    def ODB2VTK(self, WorkingDir, File):

        """
        Read odb file and write it into VTK to visualize with Paraview

        Author: Mathieu Simon
                ARTORG Center for Biomedical Engineering Research
                SITEM Insel
                University of Bern

        Adapted from: Qingbin Liu, Jiang Li, Jie Liu, 2017
                    ParaView visualization of Abaqus output on the mechanical deformation of complex microstructures
                    Computers and Geosciences, 99: 135-144
        
        March 2023
        """

        odbname = File                # Odb file name
        mesh_name = "Hexahedron"      # Element type
        mesh_corner = 8               # Nodes per element
        input_frame = range(0, 1 + 1) # Frames
        input_step = ['0']            # Step
        input_instance = ['0']        # Instance

        Dict = {'odbname':odbname,
                'mesh_name':mesh_name,
                'mesh_corner':mesh_corner,
                'input_frame':input_frame,
                'input_step':input_step,
                'input_instance':input_instance}

        Text = r"""from odbAccess import *
from textRepr import *
from string import *
from time import *

# Define frame, step and instance
input_frame = {input_frame}
input_step = {input_step}
input_instance = {input_instance}

# display the reading result of odb2vtk file
print('\nBasic Information:')
print('Model:' + '{odbname}' + '; Mesh type:' + '{mesh_name}' + '; Number of blocks:' + str(1))
print('Convert frames: ' + str(input_frame[0]) + ' to ' + str(input_frame[-1]))
print('Step & Instance : ' + str(input_step) + ', ' + str(input_instance))

# open an ODB ( Abaqus output database )
starttime = time()
odb = openOdb('{odbname}', readOnly=True)
print('\nODB opened')

# access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
rootassembly = odb.rootAssembly
instance = rootassembly.instances
# access attribute information
step = odb.steps
# get instance & step information : Quantity and all names
allinstancestr = str(instance)
autoins = allinstancestr.split("'")
inslen = len(autoins) / 4
instance_N = range(0, inslen)
allstepstr = str(step)
autostep = allstepstr.split("'")
steplen = len(autostep) / 4
step_N = range(0, steplen)

for i in input_step:
    if (steplen < int(i)):
        print('\nInput step exceeds the range of steps')
        os._exit(0)
for i in input_instance:
    if (inslen < int(i)):
        print('\nInput instance exceeds the range of instances')
        os._exit(0)

# step cycle
for step_i in input_step:
    n = int(step_i) * 4 + 1
    stepname = autostep[n]
    print('\nStep: ' + stepname)
    # instance cycle
    for ins_i in input_instance:
        n = int(ins_i) * 4 + 1
        instancename = autoins[n]
        print('\nInstance: ' + instancename)

        # access nodes & elements
        node = instance[instancename].nodes
        element = instance[instancename].elements
        n_nodes = len(node)
        n_elements = len(element)
        # access attribute(fieldOutputs) information
        frame = step[stepname].frames

        # compute the number of element of each block
        p_elements = n_elements + 1
        lp_elements = n_elements # last block

        # match nodes' label and its order in sequence (for empty nodes in tetra mesh)
        MLN = node[n_nodes - 1].label
        TOTAL = []
        # read node in sequence, and get the largest label of node(non-empty)
        # MLN is the max label of nodeset
        for i in node:
            TOTAL.append(i.label)
            if (i.label > MLN):
                MLN = i.label
        # match (the key)
        L = []
        n = 0
        for i in range(MLN):
            L.append(0)
        for i in TOTAL:
            L[i - 1] = n
            n += 1

        # frame cycle
        for i_frame in input_frame:

            # Detect whether the input frame is out of range
            try:
                TRY = odb.steps[stepname].frames[int(i_frame)]
            except:
                print('\nInput frame exceeds the range of frames')
                os._exit(0)
                break

            # Access a frame
            N_Frame = odb.steps[stepname].frames[int(i_frame)]
            print('\nFrame:' + str(i_frame))

            # create array for store result data temporarily
            # Vector-U,A,V,RF
            L0 = []
            # Tensors-S
            L1 = []
            for i in range(MLN):
                L0.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                L1.append([0, 0])

            print('\nReading U, RF ...')
            time1 = time()
            # Access Spatial displacement
            displacement = N_Frame.fieldOutputs['U']
            fieldValues = displacement.values
            for valueX in fieldValues:
                i = valueX.nodeLabel
                L0[i - 1][0] = valueX.data[0]
                L0[i - 1][1] = valueX.data[1]
                L0[i - 1][2] = valueX.data[2]

            # Access Reaction force
            Reaction_force = N_Frame.fieldOutputs['RF']
            fieldValues = Reaction_force.values
            for valueX in fieldValues:
                i = valueX.nodeLabel
                L0[i - 1][9] = valueX.data[0]
                L0[i - 1][10] = valueX.data[1]
                L0[i - 1][11] = valueX.data[2]
            print('Time elapsed: %.3f s' % (time() - time1))
            print('\nReading Stress ...')
            time1 = time()
            # access Stress components
            Stress = N_Frame.fieldOutputs['S']
            node_Stress = Stress.getSubset(position=ELEMENT_NODAL)
            fieldValues = node_Stress.values
            for valueX in fieldValues:
                L1[valueX.nodeLabel - 1][0] += 1
                L1[valueX.nodeLabel - 1][1] += valueX.mises
            # can first ave
            print('Time elapsed: %.3f s' % (time() - time1))

            '''============================================================'''

            print('\nPartitionning model and writing vtk files ...')
            time1 = time()
            print('Frame:' + str(i_frame) + '; Block:' + str(0))
            
            # Reorganization
            # Control&Storage
            # estimate whether the node has already existed
            stg_p = []
            
            # store the reorganized node for element
            stg_e = []
            
            # store the reorganized node for node
            stg_n = []
            for i in range(MLN):
                stg_p.append(-1)
            nodecount = 0
            
            # reorganize the node and element (reconstruct the mesh)
            for i in range(n_elements):
                for j in range({mesh_corner}):
                    k = element[i].connectivity[j] - 1
                    if (stg_p[k] < 0):
                        stg_p[k] = nodecount
                        stg_n.append(L[k])
                        stg_e.append(nodecount)
                        nodecount += 1
                    else:
                        stg_e.append(stg_p[k])
                        
            # compute point quantity
            n_reop = len(stg_n)
            reop_N = range(0, len(stg_n))

            # create and open a VTK(.vtu) files
            outfile = open('{odbname}'[:-4] + '_' + stepname + '_' + instancename + 'f%03d' % int(i_frame) + '.vtu', 'w')

            # <VTKFile>, including the type of mesh, version, and byte_order
            outfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
            
            # <UnstructuredGrid>
            outfile.write('<UnstructuredGrid>' + '\n')
            
            # <Piece>, including the number of points and cells
            outfile.write('<Piece NumberOfPoints="' + str(n_reop) + '"' + ' ' + 'NumberOfCells="' + str(
                lp_elements) + '">' + '\n')

            print('Writing Nodes ...')
            # <Points> Write nodes into vtk files
            displacement = N_Frame.fieldOutputs['U']
            fieldValues = displacement.values
            outfile.write('<Points>' + '\n')
            outfile.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">' + '\n')
            for i in reop_N:
                nt = stg_n[i]
                k = node[stg_n[i]].label - 1
                X, Y, Z = node[nt].coordinates[0] + L0[k][0], node[nt].coordinates[1] + L0[k][1], \
                            node[nt].coordinates[2] + L0[k][2]
                outfile.write(' ' + '%11.8e' % X + '  ' + '%11.8e' % Y + '  ' + '%11.8e' % Z + '\n')
            outfile.write('</DataArray>' + '\n')
            outfile.write('</Points>' + '\n')
            # </Points>

            print('Writing Results data ...')
            # <PointData> Write results data into vtk files
            outfile.write("<" + "PointData" + " " + "Vectors=" + '"' + "Spatial_displacement,Reaction_force" + '"' \
                            + " " + "Scalars=" + '"' + "Stress_Mises" + '"' + ">" + '\n')

            # Spatial displacement, <DataArray>
            outfile.write(
                "<" + "DataArray" + " " + "type=" + '"' + "Float32" + '"' + " " + "Name=" + '"' + "Spatial_displacement" + '"' + " " + "NumberOfComponents=" + '"' + "3" + '"' + " " + "format=" + '"' + "ascii" + '"' + ">" + '\n')
            for i in reop_N:
                k = node[stg_n[i]].label - 1
                X, Y, Z = L0[k][0], L0[k][1], L0[k][2]
                outfile.write('%11.8e' % X + ' ' + '%11.8e' % Y + ' ' + '%11.8e' % Z + '\n')
            outfile.write("</DataArray>" + '\n')
            # </DataArray>

            # Reaction force
            outfile.write(
                "<" + "DataArray" + " " + "type=" + '"' + "Float32" + '"' + " " + "Name=" + '"' + "Reaction_force" + '"' + " " + "NumberOfComponents=" + '"' + "3" + '"' + " " + "format=" + '"' + "ascii" + '"' + ">" + '\n')
            for i in reop_N:
                k = node[stg_n[i]].label - 1
                X, Y, Z = L0[k][9], L0[k][10], L0[k][11]
                outfile.write('%11.8e' % X + ' ' + '%11.8e' % Y + ' ' + '%11.8e' % Z + '\n')
            outfile.write("</DataArray>" + '\n')
            # </DataArray>

            # Stress Mises, <DataArray>
            outfile.write(
                "<" + "DataArray" + " " + "type=" + '"' + "Float32" + '"' + " " + "Name=" + '"' + "Stress_Mises" + '"' + " " + "format=" + '"' + "ascii" + '"' + ">" + '\n')
            for i in reop_N:
                k = node[stg_n[i]].label - 1
                X = L1[k][1] / L1[k][0]
                outfile.write('%11.8e' % X + '\n')
            outfile.write('</DataArray>' + '\n')
            # </DataArray>

            outfile.write("</PointData>" + '\n')
            # </PointData>

            print('Writing Cells ...')
            # <Cells> Write cells into vtk files
            outfile.write('<Cells>' + '\n')
            # Connectivity
            outfile.write('<DataArray type="Int32" Name="connectivity" format="ascii">' + '\n')
            for i in range(len(stg_e) / 8):
                outfile.write(str(stg_e[i * 8]) + ' ' + str(stg_e[i * 8 + 1]) + ' ' + str(
                    stg_e[i * 8 + 2]) + ' ' + str(stg_e[i * 8 + 3]) + ' ' + str(
                    stg_e[i * 8 + 4]) + ' ' + str(stg_e[i * 8 + 5]) + ' ' + str(
                    stg_e[i * 8 + 6]) + ' ' + str(stg_e[i * 8 + 7]) + '\n')
            outfile.write('</DataArray>' + '\n')
            
            # Offsets
            outfile.write('<DataArray type="Int32" Name="offsets" format="ascii">' + '\n')
            for i in range(len(stg_e) / {mesh_corner}):
                outfile.write(str(i * {mesh_corner} + {mesh_corner}) + '\n')
            outfile.write('</DataArray>' + '\n')
            
            # Type
            outfile.write('<DataArray type="UInt8" Name="types" format="ascii">' + '\n')
            for i in range(len(stg_e) / {mesh_corner}):
                outfile.write(str(12) + '\n')
            outfile.write('</DataArray>' + '\n')
            outfile.write('</Cells>' + '\n')
            # </Cells>

            # </Piece>
            outfile.write('</Piece>' + '\n')
            # </UnstructuredGrid>
            outfile.write('</UnstructuredGrid>' + '\n')
            # </VTKFile>
            outfile.write('</VTKFile>' + '\n')

            outfile.close()
            print('Time elapsed: %.3f s' % (time() - time1))

            '''====================================================================='''
            print('\nCreating .pvtu file for frame ' + str(i_frame) + ' ...')

        odb.close()

print('Total time elapsed: %.3f s' % (time() - starttime))
    """

        # Write and execute ODB to VTK
        os.chdir(WorkingDir)

        with open('Odb2Vtk.py','w') as File:
            File.write(Text.format(**Dict))

        os.system('abaqus python Odb2Vtk.py')

        return

Abaqus = Abaqus()

@njit
def NumbaMapping(Array):

    X, Y, Z = Array.T.shape

    # Generate nodes map
    Index = 0
    Nodes = np.zeros((Z+1,Y+1,X+1),'int')
    Coords = np.zeros((Z+1,Y+1,X+1,3),'int')
    for Xn in range(X + 1):
        for Yn in range(Y + 1):
            for Zn in range(Z + 1):
                Index += 1
                Nodes[Zn,Yn,Xn] = Index
                Coords[Zn,Yn,Xn] = [Xn, Yn, Zn]
    Nodes = Nodes[::-1,::-1]
    Coords = Coords[:,::-1]

    # Generate elements map
    Index = 0
    Elements = np.zeros((Z, Y, X),'int')
    ElementsNodes = np.zeros((Z, Y, X, 8), 'int')
    for Xn in range(X):
            for Yn in range(Y):
                for Zn in range(Z):
                    Index += 1
                    Elements[Zn, Yn, Xn] = Index
                    ElementsNodes[Zn, Yn, Xn, 0] = Nodes[Zn+1, Yn, Xn]
                    ElementsNodes[Zn, Yn, Xn, 1] = Nodes[Zn+1, Yn+1, Xn]
                    ElementsNodes[Zn, Yn, Xn, 2] = Nodes[Zn, Yn+1, Xn]
                    ElementsNodes[Zn, Yn, Xn, 3] = Nodes[Zn, Yn, Xn]
                    ElementsNodes[Zn, Yn, Xn, 4] = Nodes[Zn+1, Yn, Xn+1]
                    ElementsNodes[Zn, Yn, Xn, 5] = Nodes[Zn+1, Yn+1, Xn+1]
                    ElementsNodes[Zn, Yn, Xn, 6] = Nodes[Zn, Yn+1, Xn+1]
                    ElementsNodes[Zn, Yn, Xn, 7] = Nodes[Zn, Yn, Xn+1]
    Elements = Elements[::-1, ::-1]

    return Nodes, Coords, Elements, ElementsNodes

#@njit
def NodesText(NCS):
    Nodes, Coords, Sep = NCS
    NText = [''.join([N, Sep[0],
                    C[0], Sep[0],
                    C[1], Sep[0],
                    C[2], Sep[1]]) for N, C in zip(Nodes, Coords)]
    NText = ''.join(NText)
    return NText

#@njit
def ElementsText(ESS):
    ElementsNumber, SortedElementsNodes, Sep = ESS
    EText = [''.join([N, Sep[0],
                    E[0], Sep[0],
                    E[1], Sep[0],
                    E[2], Sep[0],
                    E[3], Sep[0],
                    E[4], Sep[0],
                    E[5], Sep[0],
                    E[6], Sep[0],
                    E[7], Sep[1]]) for N, E in zip(ElementsNumber, SortedElementsNodes)]
    EText = ''.join(EText)

    return EText

#%% Morphometry functions
class Morphometry():

    def __init__(self):
        self.Echo = True
        pass

    def SplitTriangle(self, Tri):

        """ 
        Used in SphereTriangles for MIL computation
        Splits one triange into four triangles. 
        """

        P1 = Tri[0]
        P2 = Tri[1]
        P3 = Tri[2]
        P1x = P1[0]
        P1y = P1[1]
        P1z = P1[2]
        P2x = P2[0]
        P2y = P2[1]
        P2z = P2[2]
        P3x = P3[0]
        P3y = P3[1]
        P3z = P3[2]
        P4 = ((P1x + P2x) / 2, (P1y + P2y) / 2, (P1z + P2z) / 2)
        P5 = ((P3x + P2x) / 2, (P3y + P2y) / 2, (P3z + P2z) / 2)
        P6 = ((P1x + P3x) / 2, (P1y + P3y) / 2, (P1z + P3z) / 2)
        nTri1 = (P1, P4, P6)
        nTri2 = (P4, P2, P5)
        nTri3 = (P4, P5, P6)
        nTri4 = (P5, P3, P6)
        return nTri1, nTri2, nTri3, nTri4
    
    def CorrectValues(self, X, Y, Z, Precision=1e-06):

        """
        Used in Project2UnitSphere for MIL computation
        Ensure that current direction do not go through corner or edge  
        i.e. has an angle of 45 deg.
        """

        C1 = abs(int(X / Precision)) == abs(int(Y / Precision))
        C2 = abs(int(X / Precision)) == abs(int(Z / Precision))

        if C1 and C2:
            X += Precision
            Z += 2.0 * Precision
        elif abs(int(X / Precision)) == abs(int(Y / Precision)):
            X += Precision
        elif abs(int(X / Precision)) == abs(int(Z / Precision)):
            X += Precision
        elif abs(int(Z / Precision)) == abs(int(Y / Precision)):
            Z += Precision

        return (X, Y, Z)

    def Project2UnitSphere(self, PointRS):

        """
        Used in SphereTriangles for MIL computation
        Projects an equally sided triangle patch to a unit sphere
        """

        S45 = np.sin(np.pi / 4.0)
        XYZ = [(0.0, 0.0, 0.0),
               (1.0, 0.0, 0.0),
               (0.0, 1.0, 0.0),
               (0.0, 0.0, 1.0),
               (0.5, 0.0, 0.0),
               (S45, S45, 0.0),
               (0.0, 0.5, 0.0),
               (S45, 0.0, S45),
               (0.0, S45, S45),
               (0.0, 0.0, 0.5)]

        R = PointRS[0]
        S = PointRS[1]
        T = PointRS[2]

        N5 = 4.0 * R * (1.0 - R - S - T)
        N6 = 4.0 * R * S
        N7 = 4.0 * S * (1.0 - R - S - T)
        N8 = 4.0 * R * T
        N9 = 4.0 * S * T
        N10 = 4.0 * T * (1.0 - R - S - T)

        N1 = 1.0 - R - S - T - 0.5 * N5 - 0.5 * N7 - 0.5 * N10
        N2 = R - 0.5 * N5 - 0.5 * N6 - 0.5 * N8
        N3 = S - 0.5 * N6 - 0.5 * N7 - 0.5 * N9
        N4 = T - 0.5 * N8 - 0.5 * N9 - 0.5 * N10

        aN = [N1, N2, N3, N4, N5, N6, N7, N8, N9, N10]
        X = 0.0
        Y = 0.0
        Z = 0.0
        for Node in range(10):
            X += XYZ[Node][0] * aN[Node]
            Y += XYZ[Node][1] * aN[Node]
            Z += XYZ[Node][2] * aN[Node]

        X, Y, Z = self.CorrectValues(X, Y, Z)

        Factor = 1.0 / np.sqrt(X * X + Y * Y + Z * Z)

        return (Factor * X, Factor * Y, Factor * Z)

    def SphereTriangles(self, nDirs):

        """ 
         Used in MIL computation to setup a mesh for a unit sphere. 
         
         :param nDirs: Parameter for number of triangles on unit sphere 
                       (No of triangles = 8*4^power). 
                       - TYPE: int          
                      
         :return: Triangles: List of triangles 
                  - TYPE: list[ (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ]
                  - float xi, yi, zi ... x,y,z coordinates of triangle corner point i                          
        """

        Triangles = []
        Triangles.append(((1.0, 0.0, 0.0),
                          (0.0, 1.0, 0.0),
                          (0.0, 0.0, 1.0)))
        for cDir in range(int(nDirs)):
            NewTriangles = []
            for Triangle in Triangles:
                nTri1, nTri2, nTri3, nTri4 = self.SplitTriangle(Triangle)
                NewTriangles.append(nTri1)
                NewTriangles.append(nTri2)
                NewTriangles.append(nTri3)
                NewTriangles.append(nTri4)

            Triangles = NewTriangles

        NewTriangles = []
        for Triangle in Triangles:
            nP1 = self.Project2UnitSphere(Triangle[0])
            nP2 = self.Project2UnitSphere(Triangle[1])
            nP3 = self.Project2UnitSphere(Triangle[2])
            nTr = (nP1, nP2, nP3)
            NewTriangles.append(nTr)

        Triangles = NewTriangles
        NewTriangles2 = []
        for Triangle in Triangles:
            NewTriangles2.append(Triangle)
            T1 = (Triangle[0][0], -Triangle[0][1], Triangle[0][2])
            T2 = (Triangle[1][0], -Triangle[1][1], Triangle[1][2])
            T3 = (Triangle[2][0], -Triangle[2][1], Triangle[2][2])
            NewTriangles2.append((T1, T2, T3))

        Triangles = NewTriangles2
        NewTriangles3 = []
        for Triangle in Triangles:
            NewTriangles3.append(Triangle)
            T1 = (-Triangle[0][0], Triangle[0][1], Triangle[0][2])
            T2 = (-Triangle[1][0], Triangle[1][1], Triangle[1][2])
            T3 = (-Triangle[2][0], Triangle[2][1], Triangle[2][2])
            NewTriangles3.append((T1, T2, T3))

        Triangles = NewTriangles3
        NewTriangles4 = []
        for Triangle in Triangles:
            NewTriangles4.append(Triangle)
            T1 = (Triangle[0][0], Triangle[0][1], -Triangle[0][2])
            T2 = (Triangle[1][0], Triangle[1][1], -Triangle[1][2])
            T3 = (Triangle[2][0], Triangle[2][1], -Triangle[2][2])
            NewTriangles4.append((T1, T2, T3))

        Triangles = NewTriangles4
        return Triangles

    def AreaAndCOG(self, Triangle):

        """ 
        Used in NormalAndArea for MIL computation
        Computes area and center of gravity of a triangle
        The length of the normal is "1". 
        """
        P1 = np.array(Triangle[0])
        P2 = np.array(Triangle[1])
        P3 = np.array(Triangle[2])

        P21 = P2 - P1
        P31 = P3 - P1

        A = 0.5 * np.linalg.norm(np.cross(P21, P31))

        X = (P1[0] + P2[0] + P3[0]) / 3.0
        Y = (P1[1] + P2[1] + P3[1]) / 3.0
        Z = (P1[2] + P2[2] + P3[2]) / 3.0

        Factor = 1.0 / np.sqrt(X * X + Y * Y + Z * Z)

        return (A, (Factor * X, Factor * Y, Factor * Z))

    def NormalAndArea(self, Power):

        """
        Used in OriginalDistribution for MIL computation
        Computes the normals at COG and area (weight) of 
        a triangulated unit sphere.        
        
        :param Power: Parameter for number of triangles on unit sphere 
                      (No of triangles = 8*4^power). 
                      - TYPE: int
                     
        :return: Normals: normals from COG with unit length 
                 - TYPE: list[ (nix, niy, niz) ]  
                 - float nix, niy, niz ... components of normal vectors 
                                  
                 Area_n: area of triangles which build the surface of sphere    
                 - TYPE: dict[ (nix, niy, niz) ] = value    
                 - float nix, niy, niz ... components of normal vectors 
                 - float value         ... Area for that direction                         
        """
        Triangles = self.SphereTriangles(Power)
        Normals = []
        Area_n = {}
        ASum = 0.0
        for Triangle in Triangles:
            A, COG = self.AreaAndCOG(Triangle)
            Normals.append(COG)
            Area_n[COG] = A
            ASum += A

        k = 4.0 * np.pi / ASum
        for n in Area_n:
            Area_n[n] = Area_n[n] * k

        return Normals, Area_n

    def OriginalDistribution(self, Array, Step, Power):

        """
        Used in step 2 of MIL computation
        Function computes MIL/SLD/SVD distributions for direction vectors "n" 
        using a voxel ray field going trought the RVE. Normals n = tuple(nix,niy,niz) 
        are the directions from the midpoint of a unit sphere to the COG of triangles 
        which build the surface of the sphere. Very similar to 
        self.computeOrigDistribution_STAR(). 
        A segement voxel model with isotropic resolution is needed.   
        
        :param Array: Segmented voxel model
                      - TYPE: numpy.array[kZ, jY, iX] = grayValue 
                      - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.       
                      - int grayValue  ... gray value of voxel, 0..255                  
        :param Step: Step in the considered voxel 
                      - TYPE: int > 0 
        :param Power: Power for the number of star directions = 8*4^power
                      - TYPE: int > 1
        :param Valid: Smallest valid intersection length.
                      - TYPE: float
        :param Echo: Flag for printing the  " ... Threshold Data" on stdout
                      - TYPE: bool
                                     
        @return: MIL: Mean Intercept Length 
                      - TYPE: dict[ (nix, niy, niz) ] = value 
                      - float nix, niy, niz ... components of normal vectors 
                      - float value         ... MIL for that direction          
                 SLD: Star Length Distribution 
                      - TYPE: same as for MIL 
                 SVD: Star Volume Distribution 
                      - TYPE: same as for MIL       
                 Area: Weight (Area of triange on unit sphere) for each direction 
                      - TYPE: same as for MIL          
        """

        if self.Echo:
            Text = 'Compute MIL'
            Time.Process(1, Text)

        MIL = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64,)
        SVD = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64,)
        SLD = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64,)

        SumL = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64,)
        SumL2 = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64,)
        SumL4 = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64,)

        Corners = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:],)
        Corners['swb'] = np.asarray([0.0, 0.0, 0.0], dtype='f8')
        Corners['seb'] = np.asarray([float(self.nX), 0.0, 0.0], dtype='f8')
        Corners['neb'] = np.asarray([float(self.nX), float(self.nY), 0.0], dtype='f8')
        Corners['nwb'] = np.asarray([0.0, float(self.nY), 0.0], dtype='f8')
        Corners['swt'] = np.asarray([0.0, 0.0, float(self.nZ)], dtype='f8')
        Corners['set'] = np.asarray([float(self.nX), 0.0, float(self.nZ)], dtype='f8')
        Corners['net'] = np.asarray([float(self.nX), float(self.nY), float(self.nZ)], dtype='f8')
        Corners['nwt'] = np.asarray([0.0, float(self.nY), float(self.nZ)], dtype='f8')


        sDict = Dict.empty(key_type=types.unicode_type, value_type=types.UniTuple(types.float64, 3),)
        sDict['r-dir'] = (1.0, 0.0, 0.0)
        sDict['s-dir'] = (0.0, 0.0, 1.0)
        eDict = Dict.empty(key_type=types.unicode_type, value_type=types.UniTuple(types.float64, 3),)
        eDict['r-dir'] = (0.0, 1.0, 0.0)
        eDict['s-dir'] = (0.0, 0.0, 1.0)
        nDict = Dict.empty(key_type=types.unicode_type, value_type=types.UniTuple(types.float64, 3),)
        nDict['r-dir'] = (1.0, 0.0, 0.0)
        nDict['s-dir'] = (0.0, 0.0, 1.0)
        wDict = Dict.empty(key_type=types.unicode_type, value_type=types.UniTuple(types.float64, 3),)
        wDict['r-dir'] = (0.0, 1.0, 0.0)
        wDict['s-dir'] = (0.0, 0.0, 1.0)
        bDict = Dict.empty(key_type=types.unicode_type, value_type=types.UniTuple(types.float64, 3),)
        bDict['r-dir'] = (1.0, 0.0, 0.0)
        bDict['s-dir'] = (0.0, 1.0, 0.0)
        tDict = Dict.empty(key_type=types.unicode_type, value_type=types.UniTuple(types.float64, 3),)
        tDict['r-dir'] = (1.0, 0.0, 0.0)
        tDict['s-dir'] = (0.0, 1.0, 0.0)


        InnerDict = types.DictType(types.unicode_type, types.UniTuple(types.float64, 3))
        ModelPlanes = Dict.empty(key_type=types.unicode_type, value_type=InnerDict,)
        ModelPlanes['s'] = sDict
        ModelPlanes['e'] = eDict
        ModelPlanes['n'] = nDict
        ModelPlanes['w'] = wDict
        ModelPlanes['b'] = bDict
        ModelPlanes['t'] = tDict
        
        BaseModel = Dict.empty(key_type=types.unicode_type, value_type=types.unicode_type,)
        BaseModel['s'] = 'swb'
        BaseModel['e'] = 'seb'
        BaseModel['n'] = 'nwb'
        BaseModel['w'] = 'swb'
        BaseModel['b'] = 'swb'
        BaseModel['t'] = 'swt'
        
        ViewerAt = {}
        ViewerAt['swb'] = (1.0, 1.0, 1.0)
        ViewerAt['seb'] = (-1.0, 1.0, 1.0)
        ViewerAt['neb'] = (-1.0, -1.0, 1.0)
        ViewerAt['nwb'] = (1.0, -1.0, 1.0)

        ViewerTo = {}
        ViewerTo['swb'] = 'net'
        ViewerTo['seb'] = 'nwt'
        ViewerTo['neb'] = 'swt'
        ViewerTo['nwb'] = 'set'

        Normals, Area = self.NormalAndArea(Power)
        Dict1 = {}
        Direction = ''
        for n in Normals:
            nX = n[0]
            nY = n[1]
            nZ = n[2]
            if nX >= 0.0 and nY >= 0.0 and nZ >= 0.0:
                VoxX = 1
                VoxY = 1
                VoxZ = 1
                StepX = 1
                StepY = 1
                StepZ = 1
                VoxelRay = []
                VoxelRay.append((VoxX, VoxY, VoxZ))
                if abs(nX) > abs(nY):
                    if abs(nZ) > abs(nX):
                        Direction = 'Z'
                    else:
                        Direction = 'X'
                elif abs(nZ) > abs(nY):
                    Direction = 'Z'
                else:
                    Direction = 'Y'
                PreVoxX = 1
                PreVoxY = 1
                PreVoxZ = 1
                C1 = abs(VoxX) <= self.nX
                C2 = abs(VoxY) <= self.nY
                C3 = abs(VoxZ) <= self.nZ
                while C1 and C2 and C3:
                    TMaxX = VoxX / nX
                    TMaxY = VoxY / nY
                    TMaxZ = VoxZ / nZ
                    if abs(TMaxX) < abs(TMaxY):
                        if abs(TMaxX) < abs(TMaxZ):
                            VoxX = VoxX + StepX
                        else:
                            VoxZ = VoxZ + StepZ
                    elif abs(TMaxY) < abs(TMaxZ):
                        VoxY = VoxY + StepY
                    else:
                        VoxZ = VoxZ + StepZ

                    Cc1 = abs(VoxX) <= self.nX
                    Cc2 = abs(VoxY) <= self.nY
                    Cc3 = abs(VoxZ) <= self.nZ
                    if Cc1 and Cc2 and Cc3:
                        if Direction == 'X':
                            if VoxX > PreVoxX:
                                VoxelRay.append((VoxX, VoxY, VoxZ))
                        if Direction == 'Y':
                            if VoxY > PreVoxY:
                                VoxelRay.append((VoxX, VoxY, VoxZ))
                        if Direction == 'Z':
                            if VoxZ > PreVoxZ:
                                VoxelRay.append((VoxX, VoxY, VoxZ))
                    PreVoxX = VoxX
                    PreVoxY = VoxY
                    PreVoxZ = VoxZ

                    C1 = abs(VoxX) <= self.nX
                    C2 = abs(VoxY) <= self.nY
                    C3 = abs(VoxZ) <= self.nZ

                Dict1[n] = VoxelRay

        i = 0
        Sum = len(ViewerAt) * 4.0 ** float(Power)
        for v in ViewerAt:
            Dict2 = Dict.empty(key_type=types.UniTuple(types.float64, 3), value_type=types.float64[:,:],)
            CornVoxX = int(Corners[v][0])
            CornVoxY = int(Corners[v][1])
            CornVoxZ = int(Corners[v][2])
            if CornVoxX == 0:
                CornVoxX = 1
            if CornVoxY == 0:
                CornVoxY = 1
            if CornVoxZ == 0:
                CornVoxZ = 1
            StepX = int(ViewerAt[v][0])
            StepY = int(ViewerAt[v][1])
            StepZ = int(ViewerAt[v][2])
            for n in Dict1:
                VoxelRay = Dict1[n]
                NewVoxelRay = []
                for Voxel in VoxelRay:
                    VoxelX = CornVoxX + StepX * Voxel[0] - StepX
                    VoxelY = CornVoxY + StepY * Voxel[1] - StepY
                    VoxelZ = CornVoxZ + StepZ * Voxel[2] - StepZ
                    NewVoxelRay.append((VoxelX, VoxelY, VoxelZ))

                D = n[0] * ViewerAt[v][0], n[1] * ViewerAt[v][1], n[2] * ViewerAt[v][2]
                Dict2[tuple(D)] = np.asarray(NewVoxelRay, dtype='f8')

            EntryPlanes = List([v[0], v[1], v[2]])

            if self.Echo:
                Time.Update(i/Sum, 'Setup Data')

            CVx = List([CornVoxX, CornVoxY, CornVoxZ])
            Ns = List([self.nX, self.nY, self.nZ])
            Vars = List([MIL, SVD, SLD])
            Sums = List([SumL, SumL2, SumL4])
            i, MIL, SVD, SLD = NumbaSetupMILData(i, v, Step, Array,
                                                 CVx, Ns, Vars, Sums,
                                                 Dict2, Corners, EntryPlanes,
                                                 ModelPlanes, BaseModel)

        if self.Echo:
            Time.Process(0, Text)
        
        return MIL, SVD, SLD, Area

    def FabricTensor(self, OrgMIL):

        """ 
        Used in ApproximalDistribution for MIL computation
        Compute the fabric tensors using an ellipsoidal fit
        
         :param OrgMIL: Original distribution 
                - TYPE: dict[ (nix, niy, niz) ] = value 
                - float nix, niy, niz ... components of normal vector 
                - float value         ... value for that direction             
                      
         :return: M: fabric tensor from ellipsoidal fit 
                  - TYPE: float numpy.array[3,3]            
        """

        nDir = len(OrgMIL)
        nHat = np.array(np.zeros((nDir, 6), float))
        An = np.array(np.zeros(nDir, float))
        H = np.array(np.zeros((3, 3), float))
        d = 0
        for n in OrgMIL:
            nHat[d, 0] = n[0] * n[0]
            nHat[d, 1] = n[1] * n[1]
            nHat[d, 2] = n[2] * n[2]
            nHat[d, 3] = np.sqrt(2.0) * n[1] * n[2]
            nHat[d, 4] = np.sqrt(2.0) * n[2] * n[0]
            nHat[d, 5] = np.sqrt(2.0) * n[0] * n[1]
            MILn = np.array(OrgMIL[n], dtype=float)
            An[d] = 1.0 / MILn * (1.0 / MILn)
            d += 1

        N1 = np.dot(np.transpose(nHat), nHat)
        N2 = np.dot(np.transpose(nHat), An)
        VM = np.dot(np.linalg.inv(N1), N2)
        H[(0, 0)] = VM[0]
        H[(1, 1)] = VM[1]
        H[(2, 2)] = VM[2]
        H[(1, 2)] = VM[3] / np.sqrt(2.0)
        H[(2, 0)] = VM[4] / np.sqrt(2.0)
        H[(0, 1)] = VM[5] / np.sqrt(2.0)
        H[(2, 1)] = VM[3] / np.sqrt(2.0)
        H[(0, 2)] = VM[4] / np.sqrt(2.0)
        H[(1, 0)] = VM[5] / np.sqrt(2.0)

        return H

    def EigenValuesAndVectors(self, orgMIL):
        
        """
        Used in step 4 of MIL computation
        computes the eigenvalueS and eigenvectors by fitting an ellipsoid 
        
        :param  orgMIL: Original distribution 
                - TYPE: dict[ (nix, niy, niz) ] = value 
                - float nix, niy, niz ... components of normal vector 
                - float value         ... value for that direction             
            
        :return: evalue: Eigenvalues of fabric tensor 
                 - TYPE: numpy.array[evalID] = eval
                 - int evalID ... eigenvalues ID. 0,1, or 2
                 - float eval ... current eigenvalue 
                 evector : Eigenvectors of fabric tensor
                 - TYPE: numpy.array[evectID, evalID] = evect
                 - int evectID ... eigenvector ID. 0,1, or 2                         
                 - int evalID  ... eigenvalue ID. 0,1, or 2
                 - float evect ... component of eigenvectors, e.g. evector[0,2] = ev1_z
        """

        M = self.FabricTensor(orgMIL)
        eValue, eVector = np.linalg.eig(M)
        eValue[0] = 1.0 / np.sqrt(eValue[0])
        eValue[1] = 1.0 / np.sqrt(eValue[1])
        eValue[2] = 1.0 / np.sqrt(eValue[2])


        Norm = (eValue[0] + eValue[1] + eValue[2]) / 3.0

        eValue[0] = eValue[0] / Norm
        eValue[1] = eValue[1] / Norm
        eValue[2] = eValue[2] / Norm

        return eValue, eVector

    def MIL(self, Image, Power2=4, Step=5, Power=2):

        """
        Compute Mean Intercept Length of image
        Based on mia.py from medtool from D. H. Pahr
        """

        if hasattr(Image, 'GetSize'):
            self.nX, self.nY, self.nZ = Image.GetSize()
            Array = sitk.GetArrayFromImage(Image)
        elif hasattr(Image, 'shape'):
            self.nZ, self.nY, self.nX = Image.shape
            Array = Image
        else:
            print('Image must be either numpy array or sitk image')

        # Step 1: Compute original distribution
        OrgMIL, OrgSVD, OrgSLD, Area = self.OriginalDistribution(Array, Step, Power)

        # Step 2: Compute eigen values and eigen vectors
        eValMIL, eVectMIL = self.EigenValuesAndVectors(OrgMIL)

        return eValMIL, eVectMIL

    def Trabecular(self, Image, Mask=None, DA=True):

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

        # Compute MIL for degree of anisotropy assessment
        if DA:
            Masked = sitk.Mask(Image, Mask)
            eVal, eVect = self.MIL(Masked)
            Data.loc[0,'DA'] = max(eVal) / min(eVal)

        return Data

    def Cortical(self, Image):

        """
        Compute morphology standard parameters for cortical bone
        :param Image: Segmented image of trabecular bone
                      - Type: binary sitkImage
        :return Data: Pandas data frame with computed parameters
        """

        Data = pd.DataFrame()
        Spacing = Image.GetSpacing()[0]
        Size = Image.GetSize()

        T = []
        I = []
        D = []
        Circle = measure.CircleModel()
        for S in range(Size[2]):
            Slice = sitk.Slice(Image, (0, 0, S), (Size[0], Size[1], S+1))
            Pad = sitk.ConstantPad(Slice, (1, 1, 0), (1, 1, 0))
            Array = sitk.GetArrayFromImage(Pad)

            if Array.sum() > 0:
                # Cortical thickness
                Skeleton, Distance = morphology.medial_axis(Array[0], return_distance=True)
                T.append(2 * np.mean(Distance[Skeleton]) * Spacing)

                # Inertia
                Properties = measure.regionprops(Array)[0]
                I.append(Properties.inertia_tensor[0,0] * Properties.area * Spacing**2)

                # Fitted diameter
                Circle.estimate(Properties.coords[:,1:])
                D.append(Circle.params[2] * Spacing * 2)

        
        Props = ['C.Th (mm)', 'I (mm4)', 'D (mm)']
        Measures = [np.median(T), np.median(I), np.median(D)]
        for Prop, Stat in zip(Props, Measures):
            Data.loc[0,Prop] = Stat
        
        return Data

    def SegmentBone(self, Image, Sigma=0.02, nThresholds=2, Mask=True, CloseSize=1):

        """
        Perform segmentation of bone form gray value image
        Step 1: gaussian filter to smooth values
        Step 2: multiple Otsu's algorithm for segmentation
        Step 3: Crop image at bone limits to reduce memory
        Step 4: Pad to avoid border contacts
        Optional
        Step 5: Close contour for further hole filling
        Step 6: Mask generation by filling holes in consecutive z slices

        :param Image: Gray value image
                      - Type: sitkImage
        :param Sigma: Filter width
                      - Type: float
        :param nThresholds: Number of Otsu's threshold
                      - Type: int
        :param Mask: Generate mask or not
                      - Type: bool
        :param CloseSize: Radius used to close the contour
                      - Type: int

        :return GrayCrop: Cropped gray value image
                BinCrop: Cropped binary segmented image
                Mask: Generated mask
        """

        # Filter image to reduce noise
        Gauss = sitk.SmoothingRecursiveGaussianImageFilter()
        Gauss.SetSigma(Sigma)
        Smooth  = Gauss.Execute(Image)

        # Segment image by thresholding using Otsu's method
        Otsu = sitk.OtsuMultipleThresholdsImageFilter()
        Otsu.SetNumberOfThresholds(nThresholds)
        Seg = Otsu.Execute(Smooth)

        # Binarize image to keep bone only
        Binarize = sitk.BinaryThresholdImageFilter()
        Binarize.SetUpperThreshold(1)
        Binarize.SetOutsideValue(255)
        Binarize.SetInsideValue(0)
        Bin = Binarize.Execute(Seg)

        # Crop image to bone
        Array = sitk.GetArrayFromImage(Bin)
        Z, Y, X = np.where(Array > 0)
        X1, X2 = int(X.min()), int(X.max())
        Y1, Y2 = int(Y.min()), int(Y.max())
        Z1, Z2 = int(Z.min()), int(Z.max())
        BinCrop = sitk.Slice(Bin, (X1, Y1, Z1), (X2, Y2, Z2))
        GrayCrop = sitk.Slice(Image, (X1, Y1, Z1), (X2, Y2, Z2))

        # Pad images to avoid contact with border
        BinCrop = sitk.ConstantPad(BinCrop, (1,1,1), (1,1,1))
        GrayCrop = sitk.ConstantPad(GrayCrop, (1,1,1), (1,1,1))

        if Mask:
            # Close contour
            Close = sitk.BinaryMorphologicalClosingImageFilter()
            Close.SetForegroundValue(255)
            Close.SetKernelRadius(CloseSize)
            Closed = Close.Execute(BinCrop)

            # Generate mask slice by slice
            Size = BinCrop.GetSize()
            Mask = BinCrop
            for Start in range(Size[2]):
                Slice = sitk.Slice(Closed, (0, 0, Start), (Size[0], Size[1], Start+1))

                # Pad slice in z direction to "close" holes
                Pad = sitk.ConstantPad(Slice, (0,0,1), (0,0,1), 255)

                # Fill holes
                Fill = sitk.BinaryFillholeImageFilter()
                Fill.SetForegroundValue(255)
                Filled = Fill.Execute(Pad)

                # Get center slice
                Slice = sitk.Slice(Filled, (0, 0, 1), (Size[0], Size[1], 2))

                # Paste slice into original image
                Mask = sitk.Paste(Mask, Slice, Slice.GetSize(), destinationIndex=(0,0,Start+1))

            return GrayCrop, BinCrop, Mask
        
        else:
            return GrayCrop, BinCrop

Morphometry = Morphometry()

@njit
def NumbaSetupMILData(i, v, Step, Array,
                      CVx, Ns, Vars, Sums,
                      Dict2, Corners, EntryPlanes,
                      ModelPlanes, BaseModel):
    
    # Unpack variables
    CornVoxX, CornVoxY, CornVoxZ = CVx
    nX, nY, nZ = Ns
    MIL, SVD, SLD = Vars
    SumL, SumL2, SumL4 = Sums
    
    for n in Dict2:
        i += 1
        nL = 0
        nNotValid0 = 0
        nValid0 = 0
        nNotValid1 = 0
        nValid1 = 0
        nNotValid3 = 0
        nValid3 = 0
        SumL[n] = 0.0
        SumL2[n] = 0.0
        SumL4[n] = 0.0
        NewVoxelRay = Dict2[n]
        nn = np.array([n[0], n[1], n[2]])
        nb = np.array((0.0, 0.0, 1.0))
        ng = np.cross(nn, nb)
        ns = np.cross(ng, nn)
        nr = np.cross(ns, nn)
        ng = ng / np.linalg.norm(ng)
        ns = ns / np.linalg.norm(ns)
        nr = nr / np.linalg.norm(nr)
        rmax = 0.0
        rmin = 0.0
        smax = 0.0
        smin = 0.0
        r1c = Corners[v]
        for c in Corners:
            r0c = Corners[c]
            b = r0c - r1c
            a11 = nr[0]
            a12 = ns[0]
            a13 = -nn[0]
            a21 = nr[1]
            a22 = ns[1]
            a23 = -nn[1]
            a31 = nr[2]
            a32 = ns[2]
            a33 = -nn[2]
            DET = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13)
            x = [0.0, 0.0, 0.0]
            x[0] = 1.0 / DET * ((a33 * a22 - a32 * a23) * b[0] - (a33 * a12 - a32 * a13) * b[1] + (a23 * a12 - a22 * a13) * b[2])
            x[1] = 1.0 / DET * (-(a33 * a21 - a31 * a23) * b[0] + (a33 * a11 - a31 * a13) * b[1] - (a23 * a11 - a21 * a13) * b[2])
            if (x[0] > rmax):
                rmax = x[0]
            if x[0] < rmin:
                rmin = x[0]
            if x[1] > smax:
                smax = x[1]
            if x[1] < smin:
                smin = x[1]

        for curR in range(int(rmin), int(rmax + 1), Step):
            for curS in range(int(smin), int(smax + 1), Step):
                for Plane in EntryPlanes:
                    CutPlane = ModelPlanes[Plane]
                    r1 = Corners[BaseModel[Plane]]
                    r0 = curR * nr + curS * ns + r1c
                    at = nn
                    br = np.array(CutPlane['r-dir'])
                    cs = np.array(CutPlane['s-dir'])
                    b = r0 - r1
                    a11 = br[0]
                    a12 = cs[0]
                    a13 = -at[0]
                    a21 = br[1]
                    a22 = cs[1]
                    a23 = -at[1]
                    a31 = br[2]
                    a32 = cs[2]
                    a33 = -at[2]
                    DET = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13)
                    x = [0.0, 0.0, 0.0]
                    x[0] = 1.0 / DET * ((a33 * a22 - a32 * a23) * b[0] - (a33 * a12 - a32 * a13) * b[1] + (a23 * a12 - a22 * a13) * b[2])
                    x[1] = 1.0 / DET * (-(a33 * a21 - a31 * a23) * b[0] + (a33 * a11 - a31 * a13) * b[1] - (a23 * a11 - a21 * a13) * b[2])
                    ipt = x[0] * br + x[1] * cs + r1
                    C1 = ipt[0] >= 0.0
                    C2 = ipt[1] >= 0.0
                    C3 = ipt[2] >= 0.0
                    C4 = ipt[0] <= nX
                    C5 = ipt[1] <= nY
                    C6 = ipt[2] <= nZ
                    if C1 and C2 and C3 and C4 and C5 and C6:
                        EntryVoxX = int(ipt[0] + 0.5)
                        EntryVoxY = int(ipt[1] + 0.5)
                        EntryVoxZ = int(ipt[2] + 0.5)
                        if EntryVoxX == 0:
                            EntryVoxX = 1
                        if EntryVoxY == 0:
                            EntryVoxY = 1
                        if EntryVoxZ == 0:
                            EntryVoxZ = 1
                        StartBone = (1, 1, 1)
                        EndBone = (1, 1, 1)
                        PrevVox = (1, 1, 1)
                        StartFlag = False
                        for StartRayVox in NewVoxelRay:
                            VoxX = StartRayVox[0] - (CornVoxX - EntryVoxX)
                            VoxY = StartRayVox[1] - (CornVoxY - EntryVoxY)
                            VoxZ = StartRayVox[2] - (CornVoxZ - EntryVoxZ)
                            Xv = int(VoxX - 1)
                            Yv = int(VoxY - 1)
                            Zv = int(VoxZ - 1)
                            Cc1 = VoxX < 1
                            Cc2 = VoxY < 1
                            Cc3 = VoxZ < 1
                            Cc4 = VoxX > nX
                            Cc5 = VoxY > nY
                            Cc6 = VoxZ > nZ
                            if Cc1 or Cc2 or Cc3 or Cc4 or Cc5 or Cc6:
                                if StartFlag == True:
                                    if VoxX > nX or VoxY > nY or VoxZ > nZ:
                                        StartFlag = False
                                        EndBone = PrevVox[0],PrevVox[1], PrevVox[2]
                                        lx = StartBone[0] - EndBone[0]
                                        ly = StartBone[1] - EndBone[1]
                                        lz = StartBone[2] - EndBone[2]
                                        L2 = lx * lx + ly * ly + lz * lz
                                        if L2 > 0.0:
                                            nL += 1
                                            SumL[n] += L2 ** 0.5
                                            SumL2[n] += L2
                                            SumL4[n] += L2 * L2
                            elif Array[Zv, Yv, Xv] == 0:
                                if StartFlag == True:
                                    StartFlag = False
                                    EndBone = PrevVox[0], PrevVox[1], PrevVox[2]
                                    lx = StartBone[0] - EndBone[0]
                                    ly = StartBone[1] - EndBone[1]
                                    lz = StartBone[2] - EndBone[2]
                                    L2 = lx * lx + ly * ly + lz * lz
                                    Text2 = L2
                                    if L2 > 0.0:
                                        nL += 1
                                        SumL[n] += L2 ** 0.5
                                        SumL2[n] += L2
                                        SumL4[n] += L2 * L2
                            elif StartFlag == False:
                                StartBone = (VoxX, VoxY, VoxZ)
                                StartFlag = True
                            PrevVox = (VoxX, VoxY, VoxZ)

                        break

        n2 = (-n[0], -n[1], -n[2])
        MIL[n] = SumL[n] / float(nL)
        MIL[n2] = SumL[n] / float(nL)
        SLD[n] = SumL2[n] / SumL[n]
        SLD[n2] = SumL2[n] / SumL[n]
        SVD[n] = np.pi / 3.0 * SumL4[n] / SumL[n]
        SVD[n2] = np.pi / 3.0 * SumL4[n] / SumL[n]

    return i, MIL, SVD, SLD

#%% Tensor algebra function
class Tensor():

    def __init__(self):
        pass

    def CheckPosition(self, A, B):
        AShape = A.shape
        BShape = B.shape
        if AShape[len(AShape) - 1] < BShape[0]:
            print('\nInconsistent Shape  A=', AShape, ' B=', BShape)
        return

    def CheckShape(self, A):
        Shape = A.shape
        for Index in range(len(Shape)):
            if Shape[Index] != 3:
                print('\nOrder of Tensor is not correct: ', Shape, '\n')

    def Length(self, a):
        c = np.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
        return c

    def UnitVector(self, a):
        l = self.Length(a)
        c = a / l
        return c

    def UnitMatrix(self, n):
        I = np.zeros((n, n), float)
        for row in range(n):
            for col in range(n):
                if row == col:
                    I[col, row] = 1.0

        return I

    def CrossProduct(self, a, b):
        c1 = a[1] * b[2] - a[2] * b[1]
        c2 = a[2] * b[0] - a[0] * b[2]
        c3 = a[0] * b[1] - a[1] * b[0]
        c = np.array([c1, c2, c3])
        return c

    def DyadicProduct(self, A, B):

        self.CheckShape(A)
        self.CheckShape(B)
        self.CheckPosition(A, B)

        type = 10 * len(A.shape) + len(B.shape)
        C = np.array([])
        if type == 11:
            C = np.zeros((3, 3), float)
            for i in range(3):
                for j in range(3):
                    C[i, j] = A[i] * B[j]

        elif type == 21:
            C = np.zeros((3, 3, 3), float)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        C[i, j, k] = A[i, j] * B[k]

        elif type == 22:
            C = np.zeros((3, 3, 3, 3), float)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            C[i, j, k, l] = A[i, j] * B[k, l]

        else:
           print('Dyadic product not supported')

        return C

    def DoubleContraction(self, A, B):

        self.CheckShape(A)
        self.CheckShape(B)
        self.CheckPosition(A, B)

        type = 10 * len(A.shape) + len(B.shape)
        C = np.array([])
        if type == 22:
            C = np.zeros((1, ), float)
            for i in range(3):
                for j in range(3):
                    C[(0, )] = C[(0, )] + A[i, j] * B[i, j]

        elif type == 42:
            C = np.zeros((3, 3), float)
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):
                            C[i, j] = C[i, j] + A[i, j, m, n] * B[m, n]

        elif type == 44:
            C = np.zeros((1, ), float)
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):
                            C[(0, )] = C[(0, )] + A[i, j, m, n] * B[i, j, m, n]

        else:
            print('Double contraction not supported')

        if C.shape[0] == 1:
            return C[0]
            
        else:
            return C

    # Rotate fabric (use this trick in tensor rotation)
    # R = np.linalg.inv(Fabric[1])
    # M = np.zeros((3,3))
    # for i, m in enumerate(Fabric[0]):
    #     M += m * np.outer(Fabric[1][:,i], Fabric[1][:,i])
    # RM = np.dot(np.dot(R,M),Fabric[1])
    # np.linalg.eigvals(RM) # check that eigen values didn't change

Tensor = Tensor()
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
