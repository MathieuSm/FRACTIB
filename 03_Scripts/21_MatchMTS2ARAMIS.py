#%% Script initialization
# Initialization

Version = '01'

Description = """
    Match data recorded by MTS with records of the ARAMIS device
    using peaks alignment

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import json
import argparse

from Utils import *


#%% Main
# Main code

def Main(Arguments):

    print('\nMatch ARAMIS and MTS signals')
    TIC = time.time()

    # Set directories and load data
    WD, Data, Scripts, Results = SetDirectories(Arguments.Folder)
    DataDir = Data / '03_Experiment'
    ResultsDir = Results / '02_Experiment' / Arguments.Sample
    os.makedirs(ResultsDir, exist_ok=True)

    print('\nLoad data')
    Tic = time.time()

    MTSFile = str(DataDir / '1_MTS' / 'AO_data_MTS.json')
    with open(MTSFile) as File:
        Data = json.load(File)
    MTSData = pd.DataFrame(Data[Arguments.Sample])
    MTSData.columns = ['T', 'D', 'F'] # rename time displacement and force

    ARAMISFile = str(DataDir / '2_ARAMIS' / (Arguments.Sample + '.csv'))
    ARAMISData = pd.read_csv(ARAMISFile, index_col=0, sep=';', header=1, parse_dates=['Time UTC'])

    Toc = time.time()
    PrintTime(Tic, Toc)


    # Preprocessing of ARAMIS data
    print('\nSignals preprocessing')
    Tic = time.time()

    Variable = 'TopPlate_Csys→AnatomicalCsys.' 
    Z = Variable + 'LZ [mm]'
    ARAMISData['D'] = ARAMISData[Z][0] - ARAMISData[Z]
    ARAMISData['T'] = (ARAMISData['Time UTC'] - ARAMISData.loc[0,'Time UTC']).dt.total_seconds()

    # Interpolate signal for regular sampling frequency and match MTS frequency
    print('\tARAMIS signal interpolation')
    Sampling = 1 / (MTSData.loc[1,'T'] - MTSData.loc[0,'T'])
    NPoints = int(round(ARAMISData['T'].max() * Sampling))
    Time = np.linspace(0, (NPoints-1) / Sampling, NPoints)
    Interpolated = np.interp(Time, ARAMISData['T'], ARAMISData['D'])

    # Signals filtering
    print('\tFiltering')
    FilteredMTS = pd.DataFrame()
    CutOff = 2.5 # Cut-off frequency in Hz
    FilteredMTS['D'] = Signal.Filter(MTSData['D'], Sampling, CutOff)
    FilteredMTS['F'] = Signal.Filter(MTSData['F'], Sampling, CutOff)

    FilteredARAMIS = pd.DataFrame()
    CutOff = 1 # Cut-off frequency in Hz
    FilteredARAMIS['D'] = Signal.Filter(Interpolated, Sampling, CutOff)

    Toc = time.time()
    PrintTime(Tic, Toc)


    # Peaks detection (use prominence for stability)
    print('\nSignals alignment')
    Tic = time.time()

    ARAMISPeaks, Properties = sig.find_peaks(-FilteredARAMIS['D'], prominence=0.02)
    MTSPeaks, Properties = sig.find_peaks(-FilteredMTS['D'], prominence=0.02)

    # Signals alignment of protocol 1
    Protocol1_Shift = int(round(np.mean(ARAMISPeaks[:5] - MTSPeaks[:5])))
    Time_Shifted = Time[Protocol1_Shift:] - Time[Protocol1_Shift]

    # Artificial junction in MTS data
    Junction = sig.find_peaks(MTSData['T'])[0][0]

    # Artificial data points to compensate shift
    LinkD = [FilteredMTS.loc[Junction, 'D'], FilteredMTS.loc[Junction+1, 'D']]
    LinkF = [FilteredMTS.loc[Junction, 'F'], FilteredMTS.loc[Junction+1, 'F']]

    # Signals alignment for protocol 2
    Protocol2_Shift = int(round(np.mean(ARAMISPeaks[6:] - MTSPeaks[6:])))
    Protocol2_Shift -= Protocol1_Shift

    NewForces = np.linspace(LinkF[0], LinkF[1], Protocol2_Shift+1)
    Forces = [FilteredMTS['F'][:Junction], NewForces, FilteredMTS['F'][Junction+1:]]
    Forces = np.concatenate(Forces)

    NewDiplacements = np.linspace(LinkD[0], LinkD[1], Protocol2_Shift+1)
    Displacements = [FilteredMTS['D'][:Junction], NewDiplacements, FilteredMTS['D'][Junction+1:]]
    Displacements = np.concatenate(Displacements)

    NewTimes = np.linspace(0, Protocol2_Shift / Sampling, Protocol2_Shift+1)
    NewTimes += MTSData.loc[Junction, 'T']
    Times = [MTSData['T'][:Junction], NewTimes, NewTimes[-1] + MTSData['T'][Junction+1:]]
    Times = np.concatenate(Times)

    Toc = time.time()
    PrintTime(Tic, Toc)

    # Interpolate and filter relevant signals and store in dataframe
    Variables = ['LX [mm]', 'LY [mm]', 'LZ [mm]', 'Phi(X) [°]', 'Theta(Y) [°]', 'Psi(Z) [°]']

    Matched = pd.DataFrame({'T':Time_Shifted})
    for iV, V in enumerate(Variables):
        Data2Interpolate = ARAMISData[Variable + V]
        Interpolated = np.interp(Time, ARAMISData['T'], Data2Interpolate)
        NewData = Signal.Filter(Interpolated, Sampling, CutOff)

        if iV < 3:
            Matched[V[1]] = NewData[Protocol1_Shift:]
        else:
            Matched[V[:-7]] = NewData[Protocol1_Shift:]

    # Add force data and write csv
    MissingPoints = len(Matched) - len(Forces)
    Forces = np.concatenate([Forces, np.ones(MissingPoints)* np.nan])
    Matched['F'] = Forces
    Matched.to_csv(str(ResultsDir / 'MatchedSignals.csv'), index=False)

    Figure, Axis = plt.subplots(1,1)
    Axis.plot(Matched['T'], Matched['Z'] / Matched['Z'].max(), color=(1,0,0), label='Displacement')
    Axis.plot(Matched['T'], Matched['F'] / Matched['F'].min(), color=(0,0,1), label='Force')
    Axis.set_xlabel('Time (s)')
    Axis.set_ylabel('Normalized signals (-)')
    plt.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.12))
    plt.subplots_adjust(top=0.9)
    plt.savefig(str(ResultsDir / 'Signals.png'))
    plt.close((Figure))

    print('\nSignals matched!')
    TOC = time.time()
    PrintTime(TIC, TOC)

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
    Parser.add_argument('Sample', help='Sample to process (required)', type=str)

    # Add defaults arguments
    Parser.add_argument('-F', '--Folder', help='Root folder name', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)