# 00 Initialization
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

desired_width = 500
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', desired_width)
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width,suppress=True,formatter={'float_kind':'{:3}'.format})
plt.rc('font', size=12)


# 01 Set variables
WorkingDirectory = os.getcwd()
DataDirectory = os.path.join(WorkingDirectory,'04_Results/04_Registration/')

## Load file
Data = pd.read_csv(DataDirectory+'RegistrationStudy.csv')


# 02 Build plot data
DeformationsLevel = 'High'
XAxis, YAxis = 'N iterations', 'A'
Filter1 = Data['Deformations']==DeformationsLevel
Filter2 = Data['Pyramid Schedule']=='[50, 20, 10]'
Slice = Data[Filter1&Filter2]
XData = Slice[XAxis]
YData = Slice[YAxis]
Color = Slice['DSC'].values

## Build table
Table = Slice.pivot(XAxis, YAxis, 'DSC')
XPositions = Table.index.values
YPositions = Table.columns.values
Contourlevels = Table.transpose().values

## Levels for contour
SortedSlice = Slice.sort_values('DSC')['DSC'].values
Vmax = SortedSlice[-1]
Vmin = 0.7

## Plot
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(XData,YData,marker='o',fillstyle='none',linestyle='none',color=(0,0,0))
Axes.contour(XPositions,YPositions,Contourlevels,vmax=Vmax, vmin=Vmin, cmap='viridis')
CostContour = Axes.contourf(XPositions,YPositions,Contourlevels,vmax=Vmax,vmin=Vmin, alpha=0.85,cmap='viridis')
ColorBar = Figure.colorbar(CostContour,ax=Axes)
ColorBar.set_label('DSC (-)')
Axes.set_xlim([XPositions[0]*0.9,XPositions[-1]*1.1])
Axes.set_ylim([YPositions[0]*0.9,YPositions[-1]*1.1])
Axes.set_xscale('log')
Axes.set_yscale('log')
Axes.set_xlabel(XAxis + ' (-)')
Axes.set_ylabel(YAxis + ' (-)')
# Axes.set_title(DeformationsLevel + ' Deformations')
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.show()
plt.close(Figure)



DeformationsLevel = 'High'
XAxis, YAxis = 'N iterations', 'Non-Rigid'
Filter1 = Data['Deformations']==DeformationsLevel
Filter2 = Data['Pyramid Schedule']=='[50, 20, 10]'
Slice = Data[Filter1&Filter2]
XData = Slice[XAxis]
YData = Slice[YAxis]
Color = Slice['DSC'].values


## Plot
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
CMap = Axes.scatter(XData, YData, marker='o', c=Color, cmap='binary', linewidth=1, edgecolor=(0,0,0), s=50)
ColorBar = Figure.colorbar(CMap,ax=Axes)
ColorBar.set_label('DSC (-)')
Axes.set_xscale('log')
Axes.set_xlabel(XAxis + ' (-)')
Axes.set_ylabel(YAxis + ' registration time (s)')
Axes.set_ylim([250,700])
# Axes.set_title(DeformationsLevel + ' Deformations')
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.show()
plt.close(Figure)



DeformationsLevel = 'Low'
XAxis, YAxis = 'Pyramid Schedule', 'Non-Rigid'
Filter = Data['Deformations']==DeformationsLevel
Slice = Data[Filter]
XData = Slice[XAxis]
YData = Slice[YAxis]
Color = Slice['DSC'].values


## Plot
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
CMap = Axes.scatter(XData, YData, marker='o', c=Color, cmap='binary', linewidth=1, edgecolor=(0,0,0), s=50)
ColorBar = Figure.colorbar(CMap,ax=Axes)
ColorBar.set_label('DSC (-)')
Axes.set_xlabel(XAxis + ' (-)')
Axes.set_ylabel(YAxis + ' registration time (s)')
Axes.set_ylim([250,700])
# Axes.set_title(DeformationsLevel + ' Deformations')
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.show()
plt.close(Figure)







DeformationsLevel = 'Low'
XAxis, YAxis = 'A', 'DSC'
Filter = Data['Deformations']==DeformationsLevel
Slice = Data[Filter]
Color = 'N iterations'
Color_Categories = Slice[Color].unique()
XData = Slice[[XAxis,Color]]
YData = Slice[[YAxis,Color]]

Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(XData[XData[Color]==Color_Categories[0]][XAxis], YData[YData[Color]==Color_Categories[0]][YAxis],
          marker='o', fillstyle='none', color=(1,0,0), linestyle='none', label=Color_Categories[0])
Axes.plot(XData[XData[Color]==Color_Categories[1]][XAxis], YData[YData[Color]==Color_Categories[1]][YAxis],
          marker='o', fillstyle='none', color=(0,1,0), linestyle='none', label=Color_Categories[1])
Axes.plot(XData[XData[Color]==Color_Categories[2]][XAxis], YData[YData[Color]==Color_Categories[2]][YAxis],
          marker='o', fillstyle='none', color=(0,0,1), linestyle='none', label=Color_Categories[2])
Axes.set_xlabel(XAxis + ' (-)')
Axes.set_ylabel(YAxis + ' (-)')
# Axes.set_xscale('log')
Axes.set_ylim([0.5,0.8])
# Axes.set_title(DeformationsLevel + ' Deformations')
plt.legend(loc='lower left')
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.show()
plt.close(Figure)

X = np.log10(XData['A'].sort_values())
p = np.polyfit(X,YData.loc[X.index,'DSC'].values,deg=2)
Fit = p[2] + p[1] * X + p[0] * X ** 2


Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5),dpi=100)
Axes.plot(XData[XData[Color]==Color_Categories[0]][XAxis], YData[YData[Color]==Color_Categories[0]][YAxis],
          marker='o', fillstyle='none', color=(1,0,0), linestyle='none', label='  ' + str(Color_Categories[0]) + ' iterations')
Axes.plot(XData[XData[Color]==Color_Categories[1]][XAxis], YData[YData[Color]==Color_Categories[1]][YAxis],
          marker='o', fillstyle='none', color=(0,1,0), linestyle='none', label=str(Color_Categories[1]) + ' iterations')
Axes.plot(XData[XData[Color]==Color_Categories[2]][XAxis], YData[YData[Color]==Color_Categories[2]][YAxis],
          marker='o', fillstyle='none', color=(0,0,1), linestyle='none', label=str(Color_Categories[2]) + ' iterations')
Axes.plot(10**X,Fit,color=(0,0,0),linestyle='--',label='$p_1 + p_{2}x + p_{3}x^2$ fit')
Axes.set_xlabel(XAxis + ' (-)')
Axes.set_ylabel(YAxis + ' (-)')
# Axes.set_xscale('log')
Axes.set_ylim([0.4,0.8])
# Axes.set_title(DeformationsLevel + ' Deformations')
plt.legend(loc='lower right')
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.show()
plt.close(Figure)


X = np.arange(0,3,0.001)
Fit = p[2] + p[1] * X + p[0] * X ** 2

10**X[np.where(Fit==np.max(Fit))]
