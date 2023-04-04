#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script used to perform morphometric analysis using PyPore3D

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: April 2023
    """

#%% Imports
# Modules import

import argparse
import numpy as np
import pandas as pd
import SimpleITK as sitk
from Utils import SetDirectories, Read, Morphometry


#%% Functions
# Define functions

def AFunction(Argument):

    return

#%% Main
# Main code

def Main():

    # Set directories and read sample list
    WD, DD, SD, RD = SetDirectories('FRACTIB')
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))

    # Compute calibration curve (numbers from uCT standard evaluation)
    Phantom = DD / '02_uCT' / 'Phantom' / 'C0001905.ISQ'
    Phantom = Read.ISQ(str(Phantom))[0]

    ROIs = np.array([[[1130,1270],[1345,1485]],
                     [[1245,1385],[1750,1890]],
                     [[1665,1805],[1760,1900]],
                     [[1810,1950],[1370,1510]]])
    
    RefValues = np.array([784.7238, 410.8419, 211.6873, 100.2155])

    GrayValues = []
    for Coords in ROIs:
        ROI = Phantom[Coords[0][0]:Coords[0][1], Coords[1][0]:Coords[1][1]]
        Values = sitk.GetArrayFromImage(ROI)
        GrayValues.append(Values.mean())
    GrayValues = np.array(GrayValues)

    Ones = np.ones((len(GrayValues),1))
    X = np.concatenate([Ones, np.matrix(GrayValues).T], axis=1)
    Y = np.matrix(RefValues).T
    Intercept, Slope = np.array(np.linalg.inv(X.T * X) * X.T * Y)

    # Create data frame to store metrics
    Metrics = ['BMD (mgHA/cm3)',
               'TMD (mgHA/cm3)',
               'BMC (mgHA)',
               'C. Th (mm)',
               'BV/TV (-)',
               'Tb.Th. (mm)',
               'Tb.N. (-)',
               'Tb.Sp. (mm)',
               'DA (-)']
    MorphoData = pd.DataFrame(index=SampleList['Internal ID'], columns=Metrics)

    for iS, Sample in enumerate(SampleList['Internal ID']):

        # Read AIM gray file and compute BMD values
        FilePath = DD / '02_uCT' / Sample
        FileNumber = SampleList.loc[iS, 'MicroCT pretest file number']
        FileName = 'C000' + str(FileNumber) + '_DOWNSCALED.AIM'
        Gray = Read.AIM(str(FilePath / FileName))[0]
        BMD = sitk.GetArrayFromImage(Gray) * Slope[0] + Intercept[0]

        # Read masks to compute mean BMD
        FileName = FileName[:-4] + '_CORT_MASK.AIM'
        Cort = Read.AIM(str(FilePath / FileName))[0]
        FileName = FileName[:-14] + '_TRAB_MASK.AIM'
        Trab = Read.AIM(str(FilePath / FileName))[0]
        Mask = sitk.GetArrayFromImage(Cort+Trab)
        MorphoData.loc[Sample,'BMD (mgHA/cm3)'] = np.mean(BMD[Mask > 0])
        
        # Compute BMC
        S = np.round(Gray.GetSpacing(),6)
        Volume = np.sum(Mask > 0) * S[0] * S[1] * S[2] / 1E3
        MorphoData.loc[Sample, 'BMC mgHA'] = MorphoData.loc[Sample,'BMD (mgHA/cm3)'] * Volume

        # Read segmented image to compute mean TMD
        FileName = FileName[:-14] + '_SEG.AIM'
        Seg = Read.AIM(str(FilePath / FileName))[0]
        Seg = sitk.GetArrayFromImage(Seg)
        MorphoData.loc[Sample,'TMD (mgHA/cm3)'] = np.mean(BMD[Seg > 0])

        # Cortical morphometry usgin cort mask
        Bin = (sitk.GetArrayFromImage(Cort) == 127) * 255
        Bin = sitk.GetImageFromArray(Bin)
        Bin.SetSpacing(S)
        Data = Morphometry.Cortical(Bin)
        MorphoData.loc[Sample,'C.Th (mm)'] = Data['C.Th (mm)']

        # Trabecular morphometry using trab mask and trab seg
        Bin = (Seg == 1) * 255
        Bin = sitk.GetImageFromArray(Bin)
        Bin.SetSpacing(S)

        EMask = sitk.GetArrayFromImage(Trab) * 255
        EMask = sitk.GetImageFromArray(EMask)
        EMask.SetSpacing(S)

        Data = Morphometry.Trabecular(Bin, EMask, DA=True)
        for C in Data.columns:
            MorphoData.loc[Sample,C] = Data[C]

    
    MorphoData.to_csv(str(RD / 'Morphometry.csv'))

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
    Parser.add_argument('File', help='File to process (required)', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main()