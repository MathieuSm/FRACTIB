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
DataDirectory = os.path.join(WorkingDirectory,'04_Results/06_FractureLinePrediction/')

SampleList = os.listdir(DataDirectory)
SampleList.sort()


for Index in range(len(SampleList)):
    
    Sample = SampleList[Index]

    SamplePath = os.path.join(DataDirectory,Sample)

    # 03 Load files
    uCT_Scan = sitk.ReadImage(SamplePath + '/uCT.mhd')
    HRpQCT_Scan = sitk.ReadImage(SamplePath + '/HR-pQCT_Registered.mhd')


    # 06 Add pixels to get same image size
    uCT_Scan = sitk.GetArrayFromImage(uCT_Scan)
    HRpQCT_Scan = sitk.GetArrayFromImage(HRpQCT_Scan)


    # 06 Compute similar mid position
    Mid_Position = int(round(uCT_Scan.shape[2]/2))


    # 08 Plot images        
    Figure, Axes = plt.subplots(1, 1, figsize=(11, 9), dpi=100)

    uCT_Show = Axes.imshow(uCT_Scan[:, :, Mid_Position], cmap='bone', alpha=1)
    HRpQCT_Show = Axes.imshow(HRpQCT_Scan[:, :, Mid_Position], cmap='jet', alpha=0.5)

    SliderAxis = plt.axes([0.25, 0.15, 0.65, 0.03])
    Mid_Position_Slider = Slider(SliderAxis, 'Y Position', 0, uCT_Scan.shape[1], valinit=Mid_Position)
    
    def update(val):
        Position = Mid_Position_Slider.val
        uCT_Show.set_data(uCT_Scan[:, :, int(Position)])
        HRpQCT_Show.set_data(HRpQCT_Scan[:, :, int(Position)])
        Figure.canvas.draw_idle()

    Mid_Position_Slider.on_changed(update)

    Axes.set_title(Sample)
    plt.show()
    plt.close(Figure)

