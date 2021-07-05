#!/usr/bin/env python3

# 00 Initialization
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import time
from matplotlib.widgets import Slider

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


# 01 Set variables
WorkingDirectory = os.getcwd()
DataDirectory = os.path.join(WorkingDirectory,'02_Data/05_hFE/01_AIMs/')

SampleList = os.listdir(DataDirectory)
SampleList.sort()


for Index in range(len(SampleList)):
    
    Sample = SampleList[Index]

    SamplePath = os.path.join(DataDirectory,Sample)
    Scans = [Scan for Scan in os.listdir(SamplePath) if Scan.endswith('.mhd')]
    Scans.sort()

    # 03 Load HR-pQCT files
    FixedImage = sitk.ReadImage(SamplePath + '/' + Scans[0][:8] + '_UNCOMP.mhd')
    Offset = FixedImage.GetOrigin()
    Direction = FixedImage.GetDirection()


    # 04 Load uCT file
    Registration_Directory = os.path.join(WorkingDirectory, '04_Results/05_FractureLinePrediction', Sample)
    MovingImage = sitk.ReadImage(Registration_Directory + '/' + 'uCT_Registered.mhd')
    Spacing = MovingImage.GetSpacing()


    # 05 Resample image
    Resample = sitk.ResampleImageFilter()
    Resample.SetInterpolator = sitk.sitkLinear
    Resample.SetOutputDirection(Direction)
    Resample.SetOutputOrigin(Offset)
    Resample.SetOutputSpacing(Spacing)

    Orig_Size = np.array(FixedImage.GetSize(), dtype=np.int)
    Orig_Spacing = FixedImage.GetSpacing()
    New_Size = Orig_Size * (np.array(Orig_Spacing) / np.array(Spacing))
    New_Size = np.ceil(New_Size).astype(np.int)  # Image dimensions are in integers
    New_Size = [int(s) for s in New_Size]
    Resample.SetSize(New_Size)

    ResampledImage = Resample.Execute(FixedImage)


    # 06 Add pixels to get same image size
    FixedImage = sitk.GetArrayFromImage(ResampledImage)
    MovingImage = sitk.GetArrayFromImage(MovingImage)
    
    ZerosArray = np.zeros(MovingImage.shape)
    Diff_Z, Diff_Y, Diff_X = np.array(MovingImage.shape) - np.array(FixedImage.shape)
    for Z in range(min(FixedImage.shape[0],MovingImage.shape[0])):
        for Y in range(min(FixedImage.shape[1],MovingImage.shape[1])):
            for X in range(min(FixedImage.shape[2],MovingImage.shape[2])):
                ZerosArray[Z + int(Diff_Z/2), Y + int(Diff_Y/2), X + int(Diff_X/2)] = FixedImage[Z,Y,X]
    FixedImage = ZerosArray


    # 06 Compute similar mid position
    Mid_Position = int(round(FixedImage.shape[2]/2))


    # # 08 Plot images
    # Figure, Axes = plt.subplots(1, 1, figsize=(11, 9), dpi=100)
    #
    # Fixed_Show = Axes.imshow(FixedImage[:, :, Mid_Position], cmap='bone', alpha=1)
    # Moving_Show = Axes.imshow(MovingImage[:, :, Mid_Position], cmap='jet', alpha=0.5)
    #
    # SliderAxis = plt.axes([0.25, 0.15, 0.65, 0.03])
    # Mid_Position_Slider = Slider(SliderAxis, 'Y Position', 0, FixedImage.shape[1], valinit=Mid_Position)
    #
    # def update(val):
    #     Position = Mid_Position_Slider.val
    #     Fixed_Show.set_data(FixedImage[:, :, int(Position)])
    #     Moving_Show.set_data(MovingImage[:, :, int(Position)])
    #     Figure.canvas.draw_idle()
    #
    # Mid_Position_Slider.on_changed(update)
    #
    # Axes.set_title(Sample)
    # plt.show()
    # plt.close(Figure)

    # 08 Plot images
    print('Plot sample ' + Sample)
    Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
    Fixed_Show = Axes.imshow(FixedImage[:, :, Mid_Position], cmap='bone', alpha=1)
    Moving_Show = Axes.imshow(MovingImage[:, :, Mid_Position], cmap='jet', alpha=0.5)
    Axes.set_title(Sample)
    plt.savefig(Registration_Directory + '/RegistrationResults.png')
    # plt.show()
    plt.close(Figure)