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


def WriteVTK(VectorField,FilePath,FileName):

    # Determine 2D of 3D vector field
    Dimension = VectorField.shape[-1]
    Size = VectorField.shape[:-1]

    if Dimension == 2:

        Size = np.append(1,Size)
        NumberOfElements = Size[0] * Size[1] * Size[2]

        # Build vector arrays
        u = VectorField[:, :, 0]
        v = VectorField[:, :, 1]
        w = np.zeros(NumberOfElements).reshape(Size[1:])

    elif Dimension == 3:

        NumberOfElements = Size[0] * Size[1] * Size[2]

        # Build vector arrays
        u = VectorField[:, :, 0]
        v = VectorField[:, :, 1]
        w = VectorField[:, :, 2]


    # Build Grid point
    z, y, x = np.arange(0, Size[0]), np.arange(0, Size[1]), np.arange(0, Size[2])

    File = open(os.path.join(FilePath,FileName + '.vtk'),'w')

    # ASCII file header
    File.write('# vtk DataFile Version 4.2\n')
    File.write('VTK from Python\n')
    File.write('ASCII\n\n')
    File.write('DATASET RECTILINEAR_GRID\n')
    File.write('DIMENSIONS ' + str(Size[2]) + ' ' + str(Size[1]) + ' ' + str(Size[0]) + '\n\n')
    File.close()

    # Append ascii x,y,z
    File = open(os.path.join(FilePath,FileName + '.vtk'),'a')
    File.write('X_COORDINATES ' + str(Size[2]) + ' int\n')
    File.write(np.array2string(x.astype('int'),
                               max_line_width=NumberOfElements,
                               threshold=NumberOfElements)[1:-1] + '\n')
    File.write('\nY_COORDINATES ' + str(Size[1]) + ' int\n')
    File.write(np.array2string(y.astype('int'),
                               max_line_width=NumberOfElements,
                               threshold=NumberOfElements)[1:-1] + '\n')
    File.write('\nZ_COORDINATES ' + str(Size[0]) + ' int\n')
    File.write(np.array2string(z.astype('int'),
                               max_line_width=NumberOfElements,
                               threshold=NumberOfElements)[1:-1] + '\n\n')
    File.close()

    # Append another subheader
    File = open(os.path.join(FilePath,FileName + '.vtk'),'a')
    File.write('\nPOINT_DATA ' + str(NumberOfElements) + '\n\n')
    File.write('VECTORS Deformation float\n')
    File.close()

    # Append ascii u,v,w and build deformation magnitude array
    Magnitude = np.zeros(VectorField.shape[:-1])
    File = open(os.path.join(FilePath,FileName + '.vtk'),'a')

    if Dimension == 2:
        for j in range(Size[1]):
            for i in range(Size[2]):
                Magnitude[j,i] = np.sqrt(u[j,i]**2 + v[j,i]**2 + w[j,i]**2)
                File.write(str(u[j,i]) + ' ' + str(v[j,i]) + ' ' + str(w[j,i]) + ' ')

    elif Dimension == 3:
        for k in range(Size[0]):
            for j in range(Size[1]):
                for i in range(Size[2]):
                    Magnitude[k,j,i] = np.sqrt(u[k,j,i] ** 2 + v[k,j,i] ** 2 + w[k,j,i] ** 2)
                    File.write(str(u[k,j,i]) + ' ' + str(v[k,j,i]) + ' ' + str(w[k,j,i]) + ' ')

    File.close()

    # Append another subheader
    File = open(os.path.join(FilePath, FileName + '.vtk'), 'a')
    File.write('\n\nSCALARS Magnitude float\n')
    File.write('LOOKUP_TABLE default\n')

    if Dimension == 2:
        for j in range(Size[1]):
            for i in range(Size[2]):
                File.write(str(Magnitude[j,i]) + ' ')

    elif Dimension == 3:
        for k in range(Size[0]):
            for j in range(Size[1]):
                for i in range(Size[2]):
                    File.write(str(Magnitude[k,j,i]) + ' ')

    File.close()

    return

# 01 Set variables
WorkingDirectory = os.getcwd()
DataDirectory = os.path.join(WorkingDirectory,'02_Data/')

## Load Files
SampleList = pd.read_csv(DataDirectory+'SampleList.csv')
Sample = SampleList.loc[0,'Internal ID']

SampleDirectory = os.path.join(DataDirectory,'03_uCT',Sample+'/')
Files = [File for File in os.listdir(SampleDirectory) if File.endswith('DOWNSCALED.mhd')]
Files.sort()

FixedImage = sitk.ReadImage(SampleDirectory+Files[0])
MovingImage = sitk.ReadImage(SampleDirectory+Files[1])

## Load Test images
DataDirectory = os.path.join(WorkingDirectory,'06_Problems/01_Registration/3D_Tests/Sample443/')
FixedImage = sitk.ReadImage(DataDirectory+'FixedImage.mhd')
MovingImage = sitk.ReadImage(DataDirectory+'MovingImage.mhd')


FixedImage = sitk.GetArrayFromImage(FixedImage)
MovingImage = sitk.GetArrayFromImage(MovingImage)

ResultImage = sitk.ReadImage(DataDirectory+'ResultImage.mhd')
ResultImage = sitk.GetArrayFromImage(ResultImage)

AbsMax = 1e4

Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5),dpi=100)
Axes[0].imshow(FixedImage[:,:,225] - MovingImage[:,:470,225], cmap='jet', vmin=-AbsMax, vmax=AbsMax)
ColorBar = Axes[1].imshow(ResultImage[:,:,225]-FixedImage[:,:,225], cmap='jet', vmin=-AbsMax, vmax=AbsMax)
for i in range(2):
    Axes[i].set_xlim([0,FixedImage.shape[1]])
    Axes[i].set_ylim([0,FixedImage.shape[0]])
plt.subplots_adjust(bottom=0.3)
ColorBarAxis = Figure.add_axes([0.25, 0.1, 0.5, 0.05])
Figure.colorbar(ColorBar, cax=ColorBarAxis, orientation='horizontal')
plt.show()
plt.close(Figure)

a = ResultImage[:,:,225]-FixedImage[:,:,225]



FixedImage = FixedImage[:,100:570,150:600]
MovingImage = MovingImage[:,120:600,140:670]

Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5),dpi=100)
Axes[0].imshow(FixedImage,cmap='binary')
Axes[1].imshow(MovingImage,cmap='binary')
for i in range(2):
    Axes[i].set_xlim([0,FixedImage.shape[1]])
    Axes[i].set_ylim([0,FixedImage.shape[0]])
plt.show()
plt.close(Figure)


PreTestScan = FixedImage
PostTestScan = MovingImage

SEG_PreTestScan = np.zeros(PreTestScan.shape[:-1])
for i in range(PreTestScan.shape[0]):
    for j in range(PreTestScan.shape[1]):
        R, G, B, A = PreTestScan[i,j]
        if G < 250 or B < 250:
            SEG_PreTestScan[i,j] = 1

SEG_PostTestScan = np.zeros(PostTestScan.shape[:-1])
for i in range(PostTestScan.shape[0]):
    for j in range(PostTestScan.shape[1]):
        R, G, B, A = PostTestScan[i,j]
        if G < 250 or B < 250:
            SEG_PostTestScan[i,j] = 1

FixedImage = SEG_PreTestScan
MovingImage = SEG_PostTestScan

Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5),dpi=100)
Axes[0].imshow(FixedImage,cmap='binary')
Axes[1].imshow(MovingImage,cmap='binary')
for i in range(2):
    Axes[i].set_xlim([0,FixedImage.shape[1]])
    Axes[i].set_ylim([0,FixedImage.shape[0]])
# plt.axis('off')
plt.show()
plt.close(Figure)

Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5),dpi=100)
Axes[0].imshow(SEG_PreTestScan+SEG_PostTestScan,cmap='jet')
Axes[1].imshow(sitk.GetArrayFromImage(ResultImage),cmap='binary')
plt.show()
plt.close(Figure)

Threshold = 4000
SEG_PreTestScan = np.zeros(PreTestScan.shape)
SEG_PreTestScan += PreTestScan
SEG_PreTestScan[SEG_PreTestScan<Threshold] = 0
SEG_PreTestScan[SEG_PreTestScan>=Threshold] = 1

SEG_PostTestScan = np.zeros(PostTestScan.shape)
SEG_PostTestScan += PostTestScan
SEG_PostTestScan[SEG_PostTestScan<Threshold] = 0
SEG_PostTestScan[SEG_PostTestScan>=Threshold] = 1
SEG_PostTestScan = SEG_PostTestScan[:,60:-65,60:-65]

FixedImage = sitk.GetImageFromArray(SEG_PreTestScan)
MovingImage = sitk.GetImageFromArray(SEG_PostTestScan)

Figure, Axes = plt.subplots(1, 2, figsize=(11, 4.5),dpi=100)
Axes[0].imshow(SEG_PreTestScan[50,:,:],cmap='binary')
Axes[1].imshow(SEG_PostTestScan[50,:,:],cmap='binary')
plt.show()
plt.close(Figure)