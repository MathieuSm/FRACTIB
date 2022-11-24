import numpy as np
import sympy as sp
import SimpleITK as sitk

def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results
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



CW, Data, Scripts, Results = SetDirectories('FRACTIB')

FileNames = {}
# Common area
FileNames['Common'] = str(Results / '03_hFE' / '432_L_77_F' / 'CommonMask.mhd')
FileNames['Common_uCT'] = str(Results / '05_Localizations' / '432_L_77_F' / 'CommonMask.mhd')

# Transform parameters
FileNames['InitialTransform'] = str(Results / '05_Localizations' / '432_L_77_F' / 'InitialTransform.txt')
FileNames['Transform'] = str(Results / '05_Localizations' / '432_L_77_F' / 'TransformParameters.0.txt')


Array = sitk.GetArrayFromImage(HRpQCT_Mask)
Coords = np.array(np.where(Array))
print(Array[Coords[0,0], Coords[1,0], Coords[2,0]])
Points = Coords.T * np.array(HRpQCT_Mask.GetSpacing()[::-1])

Point = Points[0]

# Transform Center of gravity from uCT to HRpQCT space
I = sitk.ReadImage(FileNames['Common'])
Center = np.array(I.GetSize()) / 2 * np.array(I.GetSpacing())
C1 = Center + np.array(I.GetOrigin())
R1 = np.array([[-1, 0, 0],[0, -1, 0],[0, 0, -1]])
T1 = [0, 0, 0]

IT = sitk.ReadTransform(FileNames['InitialTransform'])
C2 = np.array(IT.GetFixedParameters()[:-1], 'float')
P2 = IT.GetParameters()
R2 = RotationMatrix(P2[0], P2[1], P2[2])
T2 = P2[3:]

FT = GetParameterMap(FileNames['Transform'])
C3 = np.array(FT['CenterOfRotationPoint'], 'float')
P3 = np.array(FT['TransformParameters'],'float')
R3 = np.linalg.inv(RotationMatrix(P3[0], P3[1], P3[2]))
T3 = P3[3:]

# First transform
TP1 = np.dot(R1, Point - C1) + C1 + T1
R_TP1 = Rotation.TransformPoint(Point)
print(TP1)
print(R_TP1)

# Second transform
TP2 = np.dot(R2, TP1 - C2) + C2 + T2
R_TP2 = T.TransformPoint(R_TP1)
print(TP2)
print(R_TP2)

# Third transform
TP3 = np.dot(R3, TP2 - C3) + C3 + T3

TP

