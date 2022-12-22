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

import argparse
import pandas as pd
from pylatex import Document, Section, Figure, SubFigure, NoEscape, NewPage
from pylatex.package import Package
from Utils import *


#%% Classes
# Define classes

class Arguments():

    def __init__(self):
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(File):

    # Set directories
    WD, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    ResultsDir = Results / '04_Registration'
    Report = Results / 'Report'

    # Start document
    SampleList = pd.read_csv(str(Data / 'SampleList.csv'))['Internal ID']
    Doc = Document(default_filepath=str(Report))

    for Sample in SampleList[:2]:

        Image1 = str(ResultsDir / Sample / 'RigidRegistration')
        Image2 = str(ResultsDir / Sample / 'BSplineRegistration')
        Image3 = str(ResultsDir / Sample / 'DetF')
        Image4 = str(ResultsDir / Sample / 'Ftilde')
        Image5 = str(ResultsDir / Sample / 'RigidRegistration')
        Image6 = str(ResultsDir / Sample / 'BSplineRegistration')
        Images = [Image1, Image2, Image3, Image4, Image5, Image6]

        SubCaptions = ['Rigid',
                    'B-spline',
                    r'$\lvert \mathbf{F} \rvert$',
                    r'$\tilde{\mathbf{F}$',
                    r'$\lvert \mathbf{F} \rvert$',
                    r'$\tilde{\mathbf{F}$']

        Captions = ['Registration results', 'Registration analysis', 'hFE analysis']

        # Create section and add pictures
        with Doc.create(Section(Sample, numbering=False)):

            for j in range(3):
                with Doc.create(Figure(position='h!')) as Fig:
                    
                    for i in range(2):
                        SubFig = SubFigure(position='b', width=NoEscape(r'0.5\linewidth'))
                        with Doc.create(SubFig) as SF:
                            SF.add_image(Images[2*j+i],
                                        width=NoEscape(r'\linewidth'))
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