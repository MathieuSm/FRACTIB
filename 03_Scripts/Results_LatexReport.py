#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    This script uses pylatex to generate a report with the images
    resulting from the hFE and registration analyses

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import yaml
import argparse
from pylatex import Document, Section, Figure, SubFigure, NoEscape, NewPage
from pylatex.package import Package
from Utils import *
Show.ShowPlot = False

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
    WD, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    RegDir = Results / '04_Registration'
    hFEDir  = Results / '03_hFE'
    Report = Results / 'Report'

    # Start document
    SampleList = pd.read_csv(str(Data / 'SampleList.csv'))['Internal ID']
    Doc = Document(default_filepath=str(Report))

    for Sample in SampleList:

        # Generate images with same scales for hFE and registration
        uCTDir = Data / '02_uCT' / Sample
        Files = [File for File in os.listdir(uCTDir) if File.endswith('DOWNSCALED.AIM')]
        Files.sort()

        Image = Read.AIM(str(uCTDir / Files[0]))[0]
        Spacing = Image.GetSpacing()

        ConfigFile = str(Scripts / '3_hFE' / 'ConfigFile.yaml')
        Config = ReadConfigFile(ConfigFile)
        CoarseFactor = int(round(Config['ElementSize'] / Spacing[0]))

        PreI = AdjustImageSize(Image, CoarseFactor)

        # Load decompositions
        SCs, IDs = [], []
        for Dir in [RegDir, hFEDir]:
            SC = sitk.ReadImage(str(Dir / Sample / 'J.mhd'))
            SC = Resample(SC, Size=PreI.GetSize())
            SCs.append(SC)
            
            ID = sitk.ReadImage(str(Dir / Sample / 'F_Tilde.mhd'))
            ID = Resample(ID, Size=PreI.GetSize())
            IDs.append(ID)

        # Use hFE as mask
        I1 = sitk.GetArrayFromImage(SCs[0])
        I2 = Resample(SCs[1], Size=SCs[0].GetSize())
        I2 = sitk.GetArrayFromImage(I2)
        I1[I2 == 0] = 0
        I1 = sitk.GetImageFromArray(I1)
        I1.SetSpacing(SCs[0].GetSpacing())
        SCs[0] = I1

        I1 = sitk.GetArrayFromImage(IDs[0])
        I2 = Resample(IDs[1], Size=IDs[0].GetSize())
        I2 = sitk.GetArrayFromImage(I2)
        I1[I2 == 0] = 0
        I1 = sitk.GetImageFromArray(I1)
        I1.SetSpacing(IDs[0].GetSpacing())
        IDs[0] = I1

        I1 = sitk.GetArrayFromImage(PreI)
        I2 = Resample(SCs[1], Size=PreI.GetSize())
        I2 = sitk.GetArrayFromImage(I2)
        I1[I2 == 0] = 0
        I1 =sitk.GetImageFromArray(I1)
        I1.SetSpacing(PreI.GetSpacing())
        PreI = I1


        # Compute values ranges
        SCRange = []
        for SC in SCs:
            I = sitk.GetArrayFromImage(SC)
            S = I[:,:,I.shape[2] // 2]
            SCRange.append([S[S > 0].min(), S[S > 0].max()])
        SCLow = min([SCRange[0][0], SCRange[1][0]])
        SCHigh = max([SCRange[0][1], SCRange[1][1]])

        IDRange = []
        for ID in IDs:
            I = sitk.GetArrayFromImage(ID)
            S = I[:,:,I.shape[2] // 2]
            IDRange.append([S[S > 0].min(), S[S > 0].max()])
        IDLow = min([IDRange[0][0], IDRange[1][0]])
        IDHigh = max([IDRange[0][1], IDRange[1][1]])

        # Plot
        Show.IRange = [SCLow, SCHigh]
        for iDir, Dir in enumerate([RegDir, hFEDir]):
            Show.FName = str(Dir / Sample / 'DetF')
            Show.Intensity(PreI, SCs[iDir], Axis='X')

        Show.IRange = [IDLow, IDHigh]
        for iDir, Dir in enumerate([RegDir, hFEDir]):
            Show.FName = str(Dir / Sample / 'Ftilde')
            Show.Intensity(PreI, IDs[iDir], Axis='X')
        Show.FName = None

        # Generate report
        Image1 = str(RegDir / Sample / 'RigidRegistration')
        Image2 = str(RegDir / Sample / 'BSplineRegistration')
        Image3 = str(RegDir / Sample / 'DetF')
        Image4 = str(RegDir / Sample / 'Ftilde')
        Image5 = str(hFEDir / Sample / 'DetF')
        Image6 = str(hFEDir / Sample / 'FTilde')
        Images = [Image1, Image2, Image3, Image4, Image5, Image6]

        SubCaptions = ['Rigid',
                       'B-spline',
                       NoEscape(r'$|\mathbf{F}|$'),
                       NoEscape(r'$\tilde{\mathbf{F}}$'),
                       NoEscape(r'$|\mathbf{F}|$'),
                       NoEscape(r'$\tilde{\mathbf{F}}$')]

        Captions = ['Registration results', 'Registration analysis', 'hFE analysis']

        # Create section and add pictures
        with Doc.create(Section(Sample, numbering=False)):

            for j in range(3):
                with Doc.create(Figure(position='h!')) as Fig:
                    
                    for i in range(2):
                        SubFig = SubFigure(position='b', width=NoEscape(r'0.5\linewidth'))
                        with Doc.create(SubFig) as SF:
                            SF.add_image(Images[2*j+i],
                                        width=NoEscape(r'0.9\linewidth'))
                            SF.add_caption(SubCaptions[2*j+i])

                    Fig.add_caption(NoEscape(Captions[j]))

        Doc.append(NewPage())

    Doc.packages.append(Package('subcaption', options='aboveskip=0pt'))
            
    Doc.generate_pdf(clean_tex=False)

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
    Parser.add_argument('-F', '--Folder', help='Root folder name', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)