#%%
#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
from Utils import *
from ReadDAT import Main as ReadDAT
from multiprocess import Process
import multiprocessing

plt.rc('font', size=12)

def GetAngles(R):

    # Compute Euler angles from rotation matrix
    # Assuming R = RxRyRz
    # Adapted from https://stackoverflow.com/questions/15022630/how-to-calculate-the-angle-from-rotation-matrix

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

    return [E1, E2, E3]

CW, Data, Scripts, Results = SetDirectories('FRACTIB')

SampleList = pd.read_csv(str(Data / 'SampleList.csv'))

Index = 0
Sample = SampleList.loc[Index, 'Internal ID']
uCT_ID = SampleList.loc[Index, 'MicroCT pretest file number']

#%% Load data

ExpData = pd.read_csv(str(Results / '02_Experiment' / Sample / 'MatchedSignals.csv'))
SimData5 = ReadDAT('C000' + str(uCT_ID) + '_5.dat')
SimData2 = ReadDAT('C000' + str(uCT_ID) + '_2.dat')
SimData3 = ReadDAT('C000' + str(uCT_ID) + '_3.dat')
SimData1 = ReadDAT('C000' + str(uCT_ID) + '_1.dat')
SimDataD = ReadDAT('C000' + str(uCT_ID) + '_D.dat')
SimData = ReadDAT('C000' + str(uCT_ID) + '.dat')

#%% Force displacement

X, Y = ['Z', 'FZ']

Figure, Axis = plt.subplots(1,1)
Axis.plot(ExpData[X], -ExpData[Y], color=(0,0,0), label='Experiment')
# Axis.plot(SimData5[X], SimData5[Y], color=(1,0,0), label='Original UMAT')
# Axis.plot(SimData2[X], SimData2[Y], color=(1,0,1), label='Translations locked')
# Axis.plot(SimData3[X], SimData3[Y], color=(0,0,1), label='Rotations locked')
# Axis.plot(SimData1[X], SimData1[Y], color=(0,1,1), label='All locked')
Axis.plot(SimDataD[X], SimDataD[Y], color=(0,0,1), label='All Free')
Axis.plot(SimData[X], SimData[Y], color=(1,0,0), label='Motion Imposed')
Axis.set_xlabel('Displacement (mm)')
Axis.set_ylabel('Force (N)')
plt.legend()
plt.show()

#%% Peak detection

Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=0.02)
Start = Peaks[6]
Stop = np.argmin(np.abs(ExpData['Z'] - SimDataD['Z'].max()))

Figure, Axis = plt.subplots(1,1)
Axis.plot(ExpData['T'], ExpData['Z'] / ExpData['Z'].max(), color=(1,0,0), label='Displacement')
Axis.plot(ExpData['T'], ExpData['FZ'] / ExpData['FZ'].min(), color=(0,0,1), label='Force')
Axis.fill_between([ExpData['T'][Start],ExpData['T'][Stop]],[1,1], color=(0,1,0,0.2), label='Simulation')
Axis.set_xlabel('Time (s)')
Axis.set_ylabel('Normalized signal (-)')
plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.12))
plt.subplots_adjust(top=0.86)
plt.show()

# %%

Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=False, sharey=True)
Axis[0].plot(ExpData['X'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,0), label='Experiment')
# Axis[0].plot(SimData1['X'], -SimData1['Z'], color=(0,1,1), label='All locked')
# Axis[0].plot(SimData2['X'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
# Axis[0].plot(SimData3['X'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
# Axis[0].plot(SimData5['X'], -SimData5['Z'], color=(0,0,1), label='All Free')
Axis[0].plot(SimDataD['X'], -SimDataD['Z'], color=(1,0,0), label='All Free')
Axis[0].plot(SimData['X'], -SimData['Z'], color=(0,0,1), label='Simulation')
Axis[1].plot(ExpData['Y'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,0))
# Axis[1].plot(SimData1['Y'], -SimData1['Z'], color=(0,1,1))
# Axis[1].plot(SimData2['Y'], -SimData2['Z'], color=(0,0,1))
# Axis[1].plot(SimData3['Y'], -SimData3['Z'], color=(1,0,1))
# Axis[1].plot(SimData5['Y'], -SimData5['Z'], color=(0,0,1))
Axis[1].plot(SimDataD['Y'], -SimDataD['Z'], color=(1,0,0))
Axis[1].plot(SimData['Y'], -SimData['Z'], color=(0,0,1))
Axis[0].set_xlabel('X (mm)')
Axis[0].set_ylabel('Z (mm)')
Axis[1].set_xlabel('Y (mm)')
Figure.legend(loc='upper center', ncol=3)
plt.show()

#%% Rotate coordinate system to align them
from scipy.optimize import minimize

InterpX = np.interp(ExpData['Z'][Start:Stop], SimDataD['Z'], SimDataD['X'])
InterpY = np.interp(ExpData['Z'][Start:Stop], SimDataD['Z'], SimDataD['Y'])
Ref = np.array([InterpX, InterpY, ExpData['Z'][Start:Stop]]).T

Sys = np.array(ExpData[['X','Y','Z']][Start:Stop])


def RotateSystem(Angle, Ref, Sys):
    M = RotationMatrix(Gamma=Angle[0])
    rSys = np.dot(M, Sys.T).T
    Cost = np.linalg.norm(Ref - rSys, axis=1).sum()
    return Cost

Result = minimize(RotateSystem, [3*np.pi/2], args=(Ref, Sys), bounds=([0, 2*np.pi],))
Angle = Result.x[0]

RotatedSystem = np.dot(RotationMatrix(Gamma=Angle), Sys.T).T

Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=False, sharey=True)
Axis[0].plot(RotatedSystem[:,0], -RotatedSystem[:,2], color=(0,0,0), label='Experiment')
# Axis[0].plot(SimData1['X'], -SimData1['Z'], color=(0,1,1), label='All locked')
# Axis[0].plot(SimData2['X'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
# Axis[0].plot(SimData3['X'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
# Axis[0].plot(SimData5['X'], -SimData5['Z'], color=(1,0,0), label='All Free')
Axis[0].plot(SimDataD['X'], -SimDataD['Z'], color=(0,0,1), label='All Free')
Axis[0].plot(SimData['X'], -SimData['Z'], color=(1,0,0), label='Imposed Motion')
Axis[1].plot(RotatedSystem[:,1], -RotatedSystem[:,2], color=(0,0,0))
# Axis[1].plot(SimData1['Y'], -SimData1['Z'], color=(0,1,1))
# Axis[1].plot(SimData2['Y'], -SimData2['Z'], color=(0,0,1))
# Axis[1].plot(SimData3['Y'], -SimData3['Z'], color=(1,0,1))
# Axis[1].plot(SimData5['Y'], -SimData5['Z'], color=(1,0,0))
Axis[1].plot(SimDataD['Y'], -SimDataD['Z'], color=(0,0,1))
Axis[1].plot(SimData['Y'], -SimData['Z'], color=(1,0,0))
Axis[0].set_xlabel('X (mm)')
Axis[0].set_ylabel('Z (mm)')
Axis[1].set_xlabel('Y (mm)')
Figure.legend(loc='upper center', ncol=3)
plt.show()

# %%

AnglesData = ExpData.loc[Start:Stop-1, ['Phi', 'Theta', 'Psi']]/180 * np.pi
# CPUs = multiprocessing.cpu_count()//2
# Splits = np.array_split(AnglesData, CPUs)

# Processes = [Process(target=RotateAngles, args=(Data, Angle)) for Data in Splits]

# # Start the processes
# for P in Processes:
#     P.start()

# # Wait for all processes to finish
# for P in Processes:
#     P.join()


def RotateAngles(Data, Angle):

    R_Data = pd.DataFrame(columns=Data.columns, index=Data.index)

    for Index in Data.index:
        Angles = Data.loc[Index]
        R = RotationMatrix(Alpha=-Angles[0], Beta=-Angles[1], Gamma=-Angles[2])
        R_R = np.dot(RotationMatrix(Gamma=Angle), R)
        Angles = GetAngles(R_R)
        R_Data.loc[Index] = Angles

    return R_Data

R_Data = RotateAngles(AnglesData, Angle)
R_Phi, R_Theta, R_Psi = R_Data.values.T


Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=False, sharey=True)
Axis[0].plot(-R_Phi, -ExpData['Z'][Start:Stop], color=(0,0,0), label='Experiment')
# Axis[0].plot(SimData1['Phi'], -SimData1['Z'], color=(0,1,1), label='All locked')
# Axis[0].plot(SimData2['Phi'], -SimData2['Z'], color=(0,0,1), label='Translations locked')
# Axis[0].plot(SimData3['Phi'], -SimData3['Z'], color=(1,0,1), label='Rotations locked')
# Axis[0].plot(SimData5['Phi'], -SimData5['Z'], color=(1,0,0), label='All Free')
Axis[0].plot(SimDataD['Phi'], -SimDataD['Z'], color=(0,0,1), label='All Free')
Axis[0].plot(SimData['Phi'], -SimData['Z'], color=(1,0,0), label='Imposed Motion')
Axis[1].plot(-R_Theta, -ExpData['Z'][Start:Stop], color=(0,0,0))
# Axis[1].plot(SimData1['Theta'], -SimData1['Z'], color=(0,1,1))
# Axis[1].plot(SimData2['Theta'], -SimData2['Z'], color=(0,0,1))
# Axis[1].plot(SimData3['Theta'], -SimData3['Z'], color=(1,0,1))
# Axis[1].plot(SimData5['Theta'], -SimData5['Z'], color=(1,0,0))
Axis[1].plot(SimDataD['Theta'], -SimDataD['Z'], color=(0,0,1))
Axis[1].plot(SimData['Theta'], -SimData['Z'], color=(1,0,0))
Axis[0].set_xlabel('Psi')
Axis[0].set_ylabel('Z')
Axis[1].set_xlabel('Theta')
Figure.legend(loc='upper center', ncol=3)
plt.show()


#%% Extract points for simulation

Peaks, Properties = sig.find_peaks(-ExpData['FZ'][Start:Stop], prominence=0.02)
Break = Peaks[1]

print('Displacements at break')
print(RotatedSystem[Break])

print('Rotations at break')
print([-R_Phi[Break], -R_Theta[Break]])


# End point
for P in [1.0, 1.5, 2.0, 2.5, 2.78]:
    End = np.argmin(np.abs(ExpData['Z'][Start:Stop] - P))

    print('Displacements at %.2f' % (P))
    print(RotatedSystem[End])

    print('Rotations at %.2f' % (P))
    print([-R_Phi[End], -R_Theta[End]])

#%% Verification

TestA = ExpData.loc[Stop,'Phi'] / 180 * np.pi
TestR = RotationMatrix(Alpha=TestA)
TestR_R = np.dot(RotationMatrix(Gamma=Angle), TestR)
TestAngles = GetAngles(TestR_R)
TestRA = RotationMatrix(Alpha=-TestAngles[0], Beta=TestAngles[1], Gamma=-TestAngles[2])


Figure = plt.figure(figsize=(5.5, 4))
Axis = Figure.add_subplot(111, projection='3d')
Axis.quiver(0,0,0,1,0,0,color=(0,0,0, 0.5))
Axis.quiver(0,0,0,0,1,0,color=(0,0,0, 0.5))
Axis.quiver(0,0,0,0,0,1,color=(0,0,0, 0.5))
# Axis.quiver(0,0,0,TestR[0,0], TestR[0,1], TestR[0,2],color=(1,0,0))
# Axis.quiver(0,0,0,TestR[1,0], TestR[1,1], TestR[1,2],color=(0,0,1))
# Axis.quiver(0,0,0,TestR[2,0], TestR[2,1], TestR[2,2],color=(0,1,0))
Axis.quiver(0,0,0,TestR_R[0,0], TestR_R[0,1], TestR_R[0,2],color=(1,0,0,0.5))
Axis.quiver(0,0,0,TestR_R[1,0], TestR_R[1,1], TestR_R[1,2],color=(0,0,1,0.5))
Axis.quiver(0,0,0,TestR_R[2,0], TestR_R[2,1], TestR_R[2,2],color=(0,1,0,0.5))
Axis.quiver(0,0,0,TestRA[0,0], TestRA[0,1], TestRA[0,2],color=(1,0,0))
Axis.quiver(0,0,0,TestRA[1,0], TestRA[1,1], TestRA[1,2],color=(0,0,1))
Axis.quiver(0,0,0,TestRA[2,0], TestRA[2,1], TestRA[2,2],color=(0,1,0))
# scaling hack
Bbox_min = np.min([-1, -1, -1])
Bbox_max = np.max([1, 1, 1])
Axis.auto_scale_xyz([Bbox_min, Bbox_max], [Bbox_min, Bbox_max], [Bbox_min, Bbox_max])
# make the panes transparent
Axis.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
Axis.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
Axis.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
Axis.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
Axis.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
Axis.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
# modify ticks
MinX, MaxX = -1, 1
MinY, MaxY = -1, 1
MinZ, MaxZ = -1, 1
Axis.set_xticks([MinX, 0, MaxX])
Axis.set_yticks([MinY, 0, MaxY])
Axis.set_zticks([MinZ, 0, MaxZ])
Axis.xaxis.set_ticklabels([MinX, 0, MaxX])
Axis.yaxis.set_ticklabels([MinY, 0, MaxY])
Axis.zaxis.set_ticklabels([MinZ, 0, MaxZ])

Axis.set_xlabel('X')
Axis.set_ylabel('Y')
Axis.set_zlabel('Z')

plt.show()

# %%
