#%% #!/usr/bin/python python3
# Initialization
import numpy as np
import sympy as sp
import SimpleITK as sitk
import matplotlib.pyplot as plt

def RotationMatrix(Alpha=0, Beta=0, Gamma=0):

    Rx = sp.Matrix([[1,             0,              0],
                    [0, sp.cos(Alpha), -sp.sin(Alpha)],
                    [0, sp.sin(Alpha),  sp.cos(Alpha)]])

    Ry = sp.Matrix([[ sp.cos(Beta), 0, sp.sin(Beta)],
                    [0,             1,              0],
                    [-sp.sin(Beta), 0, sp.cos(Beta)]])

    Rz = sp.Matrix([[sp.cos(Gamma), -sp.sin(Gamma), 0],
                    [sp.sin(Gamma),  sp.cos(Gamma), 0],
                    [0,             0,              1]])

    R = Rz * Ry * Rx

    return np.array(R, dtype='float')

#%% Array
# Create artificial image

Spacing = (0.5, 0.5, 0.5)
Center = (3.5, 3.5, 3.5)

Xl, Yl, Zl = 15, 15, 15

Array = np.zeros((Xl, Yl, Zl),'float')
Array[5:10,5:10,2:7] = 1.0

Figure, Axis = plt.subplots(1,1)
Axis.imshow(Array[7,:,:], cmap='binary_r')
Axis.set_xticks(np.arange(Xl), np.arange(Xl)*Spacing[0])
Axis.set_yticks(np.arange(Yl), np.arange(Yl)*Spacing[1])
plt.show()

Image = sitk.GetImageFromArray(Array)
Image.SetSpacing(Spacing)

#%% Transform
# Define artificial transform

T = sitk.Euler3DTransform()
T.SetParameters((0, 0, 0, -5, 0, 0))
T.SetFixedParameters((Center[0], Center[1], Center[2], 0))
T_Image = sitk.Resample(Image, Image, T, sitk.sitkLinear, 0.0)

#%% Plot
# Plot results

T_Array = sitk.GetArrayFromImage(T_Image).round(0)

Figure, Axis = plt.subplots(1,1)
Axis.imshow(T_Array[7,:,:], cmap='binary_r')
Axis.set_xticks(np.arange(Xl), np.arange(Xl)*Spacing[0])
Axis.set_yticks(np.arange(Yl), np.arange(Yl)*Spacing[1])
plt.show()

#%% Point
# Perform point transformation
Coords = np.array(np.where(Array.transpose(2,1,0)))
Points = Coords * np.array([Spacing]).T
Points = Points.T

P = T.GetParameters()
R = RotationMatrix(-P[0], -P[1], -P[2])
C = np.array(T.GetFixedParameters()[:-1], 'float')
t = -np.array(P[3:])

T_Points = []
for iP, Point in enumerate(Points):
    T_Point = np.dot(R, Point + t - C) + C
    T_Points.append(T_Point)
T_Points = np.array(T_Points)
T_Coords = np.round(T_Points / np.array(Spacing)).astype('int')

# Compute padding
MinP = np.abs([min(v,0) for v in np.min(T_Coords,axis=0)])
MaxP = [1 + max(v, 0) for v in (np.max(T_Coords, axis=0) - Array.shape)]

P_Array = np.zeros(Array.shape)
P_Array = np.pad(P_Array,((MinP[2], MaxP[2]), (MinP[1], MaxP[1]), (MinP[0], MaxP[0])))
P_Array[T_Coords[:,2] + MinP[2], T_Coords[:,1] + MinP[1], T_Coords[:,0] + MinP[0]] = 1

P_Array = P_Array[MinP[2]:-MaxP[2],
                  MinP[1]:-MaxP[1],
                  MinP[0]:-MaxP[0]]

Figure, Axis = plt.subplots(1,1)
Axis.imshow(P_Array[7,:,:], cmap='binary_r')
Axis.set_xticks(np.arange(Xl), np.arange(Xl)*Spacing[0])
Axis.set_yticks(np.arange(Yl), np.arange(Yl)*Spacing[1])
plt.show()

# %%
