#%%
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import mpld3

mpld3.enable_notebook()
np.set_printoptions(linewidth=500,suppress=True,formatter={'float_kind':'{:3}'.format})

def RotationMatrix(Phi=0, Theta=0, Psi=0, V=np.zeros(3), A=0):

    if sum((Phi!=0, Theta!=0, Psi!=0)) > 0:
        Rx = sp.Matrix([[1,             0,              0],
                        [0, sp.cos(Phi), -sp.sin(Phi)],
                        [0, sp.sin(Phi),  sp.cos(Phi)]])

        Ry = sp.Matrix([[ sp.cos(Theta), 0, sp.sin(Theta)],
                        [0,             1,              0],
                        [-sp.sin(Theta), 0, sp.cos(Theta)]])

        Rz = sp.Matrix([[sp.cos(Psi), -sp.sin(Psi), 0],
                        [sp.sin(Psi),  sp.cos(Psi), 0],
                        [0,             0,              1]])

        R = Rz * Ry * Rx

    elif (V != 0).any():
        a = np.cos(A) * np.eye(3)
        b = np.sin(A) * np.array([[0, -V[2], V[1]],[V[2], 0, -V[0]],[-V[1], V[0], 0]])
        c = (1-np.cos(Angle)) * np.outer(V, V)
        R = np.round(a + b + c, 15)

    return np.array(R, dtype='float')

def Show3DSys(Sys, Vector=np.zeros(3)):

    Figure = plt.figure(figsize=(5.5, 4))
    Axis = Figure.add_subplot(111, projection='3d')
    Axis.quiver(-1,0,0,2,0,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(0,-1,0,0,2,0,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(0,0,-1,0,0,2,color=(0,0,0, 0.5), linewidth=0.5, arrow_length_ratio=0)
    Axis.quiver(0,0,0,Sys[0,0], Sys[0,1], Sys[0,2],color=(1,0,0))
    Axis.quiver(0,0,0,Sys[1,0], Sys[1,1], Sys[1,2],color=(0,0,1))
    Axis.quiver(0,0,0,Sys[2,0], Sys[2,1], Sys[2,2],color=(0,1,0))

    if (Vector != 0).any():
        Axis.quiver(0, 0, 0, Vector[0], Vector[1], Vector[2], color=(0,0,0))

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

    return

def GetVectorAndAngle(R):
    """
    https://en.wikipedia.org/wiki/Rotation_matrix
    """

    V = np.array([R[2,1] - R[1,2], R[0,2] - R[2,0],R[1,0] - R[0,1]])

    Angle = np.arccos((np.trace(R)- 1)/2)

    return V / np.linalg.norm(V), Angle

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


#%% Settings

# Coordinate systems
Sys = np.eye(3)                            # Original system
Rsys = RotationMatrix(Psi=sp.pi/4)         # Rotation from original to new
Sys2 = np.dot(Rsys, Sys)                   # New system rotated

# Partial rotations
Rx = RotationMatrix(Phi=sp.pi/4)
Ry = RotationMatrix(Theta=-sp.pi/4)

# Combined rotations
R =RotationMatrix(Phi=sp.pi/4, Theta=-sp.pi/4)

# Eigen vector of rotation
Vector, Angle = GetVectorAndAngle(R)

# Check correct angle sign
RCheck = RotationMatrix(V=Vector, A=Angle)
print('Check angle sign')
print(np.round(R - RCheck,12))

#%% Partial rotations
%matplotlib widget

Show3DSys(np.dot(Rx, Sys))
Show3DSys(np.dot(Ry, np.dot(Rx, Sys)))

#%% Combined rotations
%matplotlib widget

Show3DSys(np.dot(R, Sys), Vector)

#%% Rotated coordinate system
%matplotlib widget

Show3DSys(Sys2)

#%% Match final configurations
%matplotlib widget

Show3DSys(np.dot(R, np.dot(Rsys.T, Sys2)))

#%% Match final rotation


#%% Matrices analysis

print('Rotation matrix in original configuration')
print(R)
# or np.dot(R, np.dot(Rsys.T, Sys2))
print('Angles in original configuration')
Angles = GetAngles(R)
print('Phi: %i' % (Angles[0]/np.pi*180))
print('Theta: %i' % (Angles[1]/np.pi*180))
print('Psi: %i' % (Angles[2]/np.pi*180))

print('Projected angles of the eigen vector in original configuration')
x = np.dot(Sys[0,:][:-1], Vector[:-1])

print('Phi: %i' % (np.arccos(a / (np.linalg.norm()))      /np.pi * 180))
print('Theta: %i' % (np.arccos(np.dot(Sys[1,:], Vector))/np.pi * 180))
print('Psi: %i' % (np.arccos(np.dot(Sys[2,:], Vector))/np.pi * 180))

print('Angles of the eigen vector')
print('Phi: %i' % (np.arccos(np.dot(Sys2[0,:], Vector))/np.pi * 180))
print('Theta: %i' % (np.arccos(np.dot(Sys2[1,:], Vector))/np.pi * 180))
print('Psi: %i' % (np.arccos(np.dot(Sys2[2,:], Vector))/np.pi * 180))

# %%
