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
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))['Internal ID']
    ConfigFile = str(SD / '3_hFE' / 'ConfigFile.yaml')
    Config = ReadConfigFile(ConfigFile)

    Data = pd.DataFrame(index=SampleList.values, columns=['SC','ID'])
    for Index, Sample in enumerate(SampleList):

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
        Seg[Seg == 127] = 2
        Seg[Seg == 126] = 1

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