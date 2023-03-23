#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used for the assessment of localization prediction accuracy.
    This is done by performing ordinary least square regression on the
    values (det(F) and |F~|) obtained with the registration vs the simulation

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: January 2023
    """

#%% Imports
# Modules import

import yaml
import argparse
from Utils import *
Show.ShowPlot = False
Read.Echo = False
Write.Echo = False
from skimage import morphology

#%% Functions
# Define functions
def ReadConfigFile(Filename):

    """ Read configuration file and store to dictionary """

    print('\n\nReading configuration file', Filename)
    with open(Filename, 'r') as File:
        Configuration = yaml.load(File, Loader=yaml.FullLoader)

    return Configuration
def AdjustImageSize(Image, CoarseFactor, CropZ='Crop'):

    """
    Adapted from Denis's utils_SA.py
    Images are adjusted according to CropType:
    0 = CropType.expand     (Expand image by copying layers)
    1 = CropType.crop       (Crop image)
    2 = CropType.variable   (Either crop or expand, depending on what includes less layers)
    """

    # Get array
    Array = sitk.GetArrayFromImage(Image)
    Array = Array.transpose(2, 1, 0)

    # Measure image shape
    IMDimX = np.shape(Array)[0]
    IMDimY = np.shape(Array)[1]
    IMDimZ = np.shape(Array)[2]

    AddDimX = CoarseFactor - (IMDimX % CoarseFactor)
    AddDimY = CoarseFactor - (IMDimY % CoarseFactor)

    # adjust in x and y direction
    Shape_Diff = [AddDimX, AddDimY]
    IMG_XY_Adjusted = np.lib.pad(Array,
                                 ((0, Shape_Diff[0]), (0, Shape_Diff[1]), (0, 0)),
                                 'constant', constant_values=(0),)

    if CropZ == 'Crop':
        Image_Adjusted = IMG_XY_Adjusted

    if CropZ == 'Expand':
        AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
        Shape_Diff = [AddDimX, AddDimY, AddDimZ]
        Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                    ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                    'edge')

    if CropZ == 'Variable':
        Limit = CoarseFactor / 2.0
        if IMDimZ % CoarseFactor > Limit:
            AddDimZ = CoarseFactor - (IMDimZ % CoarseFactor)
            Shape_Diff = [AddDimX, AddDimY, AddDimZ]
            Image_Adjusted = np.lib.pad(IMG_XY_Adjusted,
                                        ((0, 0), (0, 0), (0, Shape_Diff[2])),
                                        'edge')
        if IMDimZ % CoarseFactor < Limit:
            Image_Adjusted = IMG_XY_Adjusted

    Image_Adjusted = sitk.GetImageFromArray(Image_Adjusted.transpose(2, 1, 0))
    Image_Adjusted.SetSpacing(Image.GetSpacing())
    Image_Adjusted.SetOrigin(Image.GetOrigin())
    Image_Adjusted.SetDirection (Image.GetDirection())

    return Image_Adjusted
def ReadInputFile(File):

    """
    Writen in order to read the BV/TV values stored in the main
    .inp file generated using 'hFE_1Preprocessing' and returning
    it as an image where each voxel contain the computed BVTV value
    """

    # Store lines offsets, start and stop
    N = 0
    with open(File) as F:
        Offset = 0
        Offsets = [Offset]
        for Line in F:
            Offset += len(Line)
            Offsets.append(Offset)
            if '*ELEMENT,' in Line:
                N += 1
                if N == 2:
                    Start = Offset

            if 'TOPNODES' in Line:
                Stop = Offset
                break
    
    # Read and store elements data
    Elements = []
    Step = 29
    iStart = Offsets.index(Start)
    iStop = Offsets.index(Stop)
    with open(File) as F:

        for i in range(iStart, iStop, Step):
            F.seek(Offsets[i + 1])
            P = F.readline().split()

            F.seek(Offsets[i + 9])
            V = F.readline().split()

            EData = [P[3], P[6], P[9]]
            for Value in V[:4]:
                EData.append(Value[:-1])
            Elements.append(EData)
    
    # Compute BVTV
    Columns = ['X', 'Y', 'Z', 'RC', 'RT', 'PC', 'PT']
    Elements = pd.DataFrame(Elements, columns=Columns, dtype=float)
    Elements['BVTV'] = Elements['RC'] * Elements['PC'] + Elements['RT'] * Elements['PT']

    # Create array
    X = Elements['X'].unique()
    X.sort()
    Y = Elements['Y'].unique()
    Y.sort()
    Z = Elements['Z'].unique()
    Z.sort()
    Mesh = np.zeros((len(Z), len(Y), len(X)))

    for Index, Element in Elements.iterrows():
        iX = np.argmin(abs(X - Element['X']))
        iY = np.argmin(abs(Y - Element['Y']))
        iZ = np.argmin(abs(Z - Element['Z']))
        Mesh[iZ, iY, iX] = Element['BVTV']

    # Create image
    Spacing = np.array([X[1]-X[0],Y[1]-Y[0],Z[1]-Z[0]])
    Origin = np.array([X.min(), Y.min(), Z.min()])

    BVTV = sitk.GetImageFromArray(Mesh)
    BVTV.SetSpacing(Spacing)
    BVTV.SetOrigin(Origin)

    return BVTV

#%% Classes
# Define classes

class Arguments():

    def __init__(self):
        pass

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    # Set directories
    WD, DD, SD, RD = SetDirectories('FRACTIB')
    RegDir = RD / '04_Registration'
    hFEDir  = RD / '03_hFE'
    ResultsDir = RD / '05_Comparison'

    # Read sample list and configuration file
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))
    ConfigFile = str(SD / '3_hFE' / 'ConfigFile.yaml')
    Config = ReadConfigFile(ConfigFile)

    Data = pd.DataFrame(index=SampleList['Internal ID'].values, columns=['SC','ID','Dice 1','Dice 2'])
    for Index, Sample in enumerate(SampleList['Internal ID']):

        Text = 'Sample ' + str(Index+1) + '/' + str(len(SampleList))
        Time.Process(1,Text)

        # Generate images with same scales for hFE and registration
        uCTDir = DD / '02_uCT' / Sample
        Files = [File for File in os.listdir(uCTDir) if File.endswith('DOWNSCALED.AIM')]
        Files.sort()

        Image = Read.AIM(str(uCTDir / Files[0]))[0]
        Spacing = Image.GetSpacing()
        CF = int(round(Config['ElementSize'] / Spacing[0]))
        PreI = AdjustImageSize(Image, CF)
      
        # Read BVTV values from .inp file
        BVTV = ReadInputFile(str(hFEDir / Sample / 'Simulation.inp'))
        Spacing = BVTV.GetSpacing()
        Origin = BVTV.GetOrigin()
        Write.FName = 'BVTV'
        Write.MHD(BVTV, PixelType='float')
        BVTV = sitk.GetArrayFromImage(BVTV)

        # Read BSpline registration results (for local Dice computation)
        Time.Update(1/6,'Compute Dice')
        Reg = sitk.ReadImage(str(RegDir / Sample / 'NonRigid.mhd'))
        Org = sitk.ReadImage(str(RegDir / Sample / 'Rigid.mhd'))
        Otsu = sitk.OtsuMultipleThresholdsImageFilter()
        Otsu.SetNumberOfThresholds(2)
        BinReg = Otsu.Execute(Reg)
        BinOrg = Otsu.Execute(Org)
        BinIm = Otsu.Execute(PreI)

        Binary = sitk.BinaryThresholdImageFilter()
        Binary.SetUpperThreshold(1.5)
        Binary.SetInsideValue(0)
        Binary.SetOutsideValue(1)
        BinReg = Binary.Execute(BinReg)
        BinOrg = Binary.Execute(BinOrg)
        BinIm = Binary.Execute(BinIm)

        Measure = sitk.LabelOverlapMeasuresImageFilter()
        Measure.Execute(BinIm, BinOrg)
        Data.loc[Sample, 'Dice 1'] = Measure.GetDiceCoefficient()
        Measure.Execute(BinIm, BinReg)
        Data.loc[Sample, 'Dice 2'] = Measure.GetDiceCoefficient()

        Dice = BinReg * BinIm
        Dice = sitk.GetArrayFromImage(Dice)

        # Pad or crop to adjust XY size
        Time.Update(2/6,'Convolve Dice')
        Delta = np.ceil(np.array(BVTV.shape) * CF - (np.array(Dice.shape))).astype(int)
        bDelta = np.floor(Delta / 2).astype(int)
        aDelta = np.ceil(Delta / 2).astype(int)
        
        bDelta = np.abs(bDelta)
        if Delta[0] < 0:
            Dice = Dice[bDelta[0]:aDelta[0],:,:]
        elif Delta[0] > 0:
            Dice = np.pad(Dice, ((bDelta[0], aDelta[0]), (0, 0), (0, 0)), 'reflect')
        if Delta[1] < 0:
            Dice = Dice[:,bDelta[1]:aDelta[1],:]
        elif Delta[1] > 0:
            Dice = np.pad(Dice, ((0, 0), (bDelta[1], aDelta[1]), (0, 0)), 'reflect')
        if Delta[2] < 0:
            Dice = Dice[:,:,bDelta[2]:aDelta[2]]
        elif Delta[2] > 0:
            Dice = np.pad(Dice, ((0, 0), (0, 0), (bDelta[2], aDelta[2])), 'reflect')
        
        Shape = (BVTV.shape[0], BVTV.shape[1], BVTV.shape[2], CF, CF, CF)
        Convolve = np.ones(Shape) / CF**3

        for i in range(BVTV.shape[0]):
            for j in range(BVTV.shape[1]):
                for k in range(BVTV.shape[2]):
                    Convolve[i,j,k] = Convolve[i,j,k] * Dice[i*CF:(i+1)*CF,j*CF:(j+1)*CF,k*CF:(k+1)*CF]

        Convolve = np.sum(Convolve, axis=(-1,-2,-3))
        ImageConv = sitk.GetImageFromArray(Convolve)
        ImageConv.SetSpacing(Spacing)
        ImageConv.SetOrigin(Origin)
        Write.MHD(ImageConv, str(RegDir / Sample / 'Dice'))
        # Test = Show.OLS(BVTV[Mask], Convolve[Mask], Cmap=Y)

        # Load decompositions
        Time.Update(3/6, 'Load gradients')
        SCs, IDs = [], []
        for Dir in [RegDir, hFEDir]:
            SC = sitk.ReadImage(str(Dir / Sample / 'J.mhd'))
            SC = sitk.GetArrayFromImage(SC)
            SCs.append(SC)
            
            ID = sitk.ReadImage(str(Dir / Sample / 'F_Tilde.mhd'))
            ID = sitk.GetArrayFromImage(ID)
            IDs.append(ID)

        # Plot correlations (use BVTV as mask)
        Time.Update(4/6, 'DetF analysis')
        Mask = BVTV > 0
        X = SCs[0][Mask]
        Y = SCs[1][Mask]
        Show.FName = str(ResultsDir / (Sample + '_SC'))
        SC_R = Show.OLS(X, Y)
        X = np.cbrt(X)
        Y = np.cbrt(Y)
        SC_R = np.mean(np.abs(X - Y) / X)

        Time.Update(5/6, 'F~ analysis')
        X = IDs[0][Mask]
        Y = IDs[1][Mask]
        # For sample 446
        if 0 in X:
            Y = Y[X > 0]
            X = X[X > 0]
        Show.FName = str(ResultsDir / (Sample + '_ID'))
        ID_R = Show.OLS(X, Y)
        ID_R = np.mean(np.abs(X - Y) / X)
        Show.FName = None

        # Add results to data frame
        Data.loc[Sample, 'SC'] = SC_R
        Data.loc[Sample, 'ID'] = ID_R

        Index = SampleList[SampleList == Sample].index[0]
        Time.Process(0,Text)

    # Boxplot and save data
    Labels = [r'det($\mathbf{F})$', r'$||\widetilde{\mathbf{F}}||$']
    Show.FName = str(RD / 'AccuracySummary')
    Show.BoxPlot([Data['SC'], Data['ID']], Labels=['', 'Relative difference (-)'], SetsLabels=[Labels[0], Labels[1]])
    Data.to_csv(str(RD / 'Correlations.csv'))

    return

#%% Execution part
# Execution as main
if __name__ == '__main__':

    # Initiate the parser with a description
    FC = argparse.RawDescriptionHelpFormatter
    Parser = argparse.ArgumentParser(description=Description, formatter_class=FC)

    # Add long and short argument
    SV = Parser.prog + ' version ' + Version
    Parser.add_argument('-V', '--Version', help='Show script version', action='version', version=SV)
    Parser.add_argument('--Folder', help='Root folder of the project', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)