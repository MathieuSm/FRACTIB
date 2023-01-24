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
    Y = Elements['Y'].unique()
    Z = Elements['Z'].unique()
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
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    # Set directories
    WD, DD, SD, RD = SetDirectories(Arguments.Folder)
    RegDir = RD / '04_Registration'
    hFEDir  = RD / '03_hFE'
    ResultsDir = RD / '05_Comparison'

    # Read sample list and configuration file
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))
    ConfigFile = str(SD / '3_hFE' / 'ConfigFile.yaml')
    Config = ReadConfigFile(ConfigFile)

    Data = pd.DataFrame(index=SampleList['Internal ID'].values, columns=['SC','ID'])
    for Index, Sample in enumerate(SampleList['Internal ID']):

        Text = 'Compute accuracy for sample ' + str(Index+1) + '/' + str(len(SampleList))
        ProcessTiming(1,Text)

        # Generate images with same scales for hFE and registration
        uCTDir = DD / '02_uCT' / Sample
        Files = [File for File in os.listdir(uCTDir) if File.endswith('DOWNSCALED.AIM')]
        Files.sort()

        Image = Read.AIM(str(uCTDir / Files[0]))[0]
        Seg = Read.AIM(str(uCTDir / (Files[0][:-4] + '_SEG.AIM')))[0]
        Spacing = Image.GetSpacing()
        CoarseFactor = int(round(Config['ElementSize'] / Spacing[0]))

        PreI = AdjustImageSize(Image, CoarseFactor)
        Seg = AdjustImageSize(Seg, CoarseFactor)
        Seg = sitk.GetArrayFromImage(Seg)
        Seg[Seg == 127] = 1
        Seg[Seg == 126] = 1

        # Read BVTV values from .inp file
        FileName = SampleList.loc[Index, 'MicroCT pretest file number']
        FileName = 'C000' + str(FileName) + '_DOWNSCALED_00_FZ_MAX.inp'
        BVTV = ReadInputFile(str(hFEDir / Sample / FileName))
        ProgressNext(1)

        # Load decompositions
        SCs, IDs = [], []
        for Dir in [RegDir, hFEDir]:
            SC = sitk.ReadImage(str(Dir / Sample / 'J.mhd'))
            SC = sitk.GetArrayFromImage(SC)
            SCs.append(SC)
            
            ID = sitk.ReadImage(str(Dir / Sample / 'F_Tilde.mhd'))
            ID = sitk.GetArrayFromImage(ID)
            IDs.append(ID)
        ProgressNext(2)

        # Plot correlations (use hFE as mask)
        Mask = SCs[1] > 0
        X = SCs[0][Mask].flatten()
        Y = SCs[1][Mask].flatten()
        Show.FName = str(ResultsDir / (Sample + 'SC'))
        SC_R = Show.OLS(X, Y)
        ProgressNext(3)

        Mask = IDs[1] > 0
        X = IDs[0][Mask].flatten()
        Y = IDs[1][Mask].flatten()
        Show.FName = str(ResultsDir / (Sample + 'ID'))
        ID_R = Show.OLS(X, Y)
        Show.FName = None
        ProgressNext(4)

        # Add results to data frame
        Data.loc[Sample, 'SC'] = SC_R.rsquared
        Data.loc[Sample, 'ID'] = ID_R.rsquared

        Index = SampleList.index[SampleList == Sample][0]
        ProcessTiming(0)

    # Boxplot and save data
    Labels = [r'$|\mathbf{F}|$', r'$|\widetilde{\mathbf{F}}|$']
    Show.FName = str(RD / 'AccuracySummary')
    Show.BoxPlot([Data['SC'],Data['ID']], Labels=['', r'R$^2$'], SetsLabels=Labels)
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