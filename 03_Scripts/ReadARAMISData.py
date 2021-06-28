import os
import sys
import pandas as pd
import matplotlib.pyplot as plt


## List csv files
ARAPath = '/home/mathieu/Documents/MscThesis/02_Data/FractureZoneAssessment/04_Experiment/2_ARAMIS/'
SamplesFiles = [File for File in os.listdir(ARAPath) if File.endswith('csv')]
SamplesFiles.sort()

## Read a defined sample
SampleNumber = 1
Time = 'Time UTC'
ARASampleData = pd.read_csv(os.path.join(ARAPath,SamplesFiles[SampleNumber-1]),sep=';',header=1,parse_dates=[Time])


## Compare the 3 Coordinate systems
CoordinateSystems = ['TopPlate_Csys','AnatomicalCsys','Global coordinate system']
Distances = ['L [mm]','LX [mm]','LY [mm]','LZ [mm]']
Angles = ['Phi(X) [°]','Theta(Y) [°]','Psi(Z) [°]']

'TopPlate_Csys→AnatomicalCsys.LZ [mm]'


Time = ARASampleData['Time UTC'] - ARASampleData['Time UTC'].loc[0]

RelativeMovement = CoordinateSystems[1] + '→' + CoordinateSystems[2]
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[0]],color=(0,0,0),label=Distances[0][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[1]],color=(0,0,1),label=Distances[1][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[2]],color=(0,1,0),label=Distances[2][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[3]],color=(1,0,0),label=Distances[3][:-5])
Axes.set_xlabel('Time [s]')
Axes.set_ylabel('Displacement [mm]')
Axes.set_title(RelativeMovement)
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)

RelativeMovement = CoordinateSystems[0] + '→' + CoordinateSystems[2]
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[0]],color=(0,0,0),label=Distances[0][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[1]],color=(0,0,1),label=Distances[1][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[2]],color=(0,1,0),label=Distances[2][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[3]],color=(1,0,0),label=Distances[3][:-5])
Axes.set_xlabel('Time [s]')
Axes.set_ylabel('Displacement [mm]')
Axes.set_title(RelativeMovement)
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)

RelativeMovement = CoordinateSystems[0] + '→' + CoordinateSystems[1]
Figure, Axes = plt.subplots(1, 1, figsize=(5.5, 4.5), dpi=100)
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[0]],color=(0,0,0),label=Distances[0][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[1]],color=(0,0,1),label=Distances[1][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[2]],color=(0,1,0),label=Distances[2][:-5])
Axes.plot(Time/1e9,ARASampleData[RelativeMovement + '.' + Distances[3]],color=(1,0,0),label=Distances[3][:-5])
Axes.set_xlabel('Time [s]')
Axes.set_ylabel('Displacement [mm]')
Axes.set_title(RelativeMovement)
plt.legend(loc='upper left')
plt.show()
plt.close(Figure)