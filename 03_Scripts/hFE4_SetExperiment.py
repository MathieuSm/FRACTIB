#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Compute the best coordinate system rotation in order to
    rewrite the Abaqus input file to impose the 6 DOFs of the
    top node according to the experiment.

    Version Control:
        01 - Original script
        02 - Rewrite steps from start instead of adding them at
             the end

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: December 2022
    """

#%% Imports
# Modules import

import argparse
from Utils import *
from scipy.optimize import minimize

#%% Functions
# Define functions

def ComputeCost(Angle, Paths):
    R = RotationMatrix(Psi=Angle[0])
    RSys = np.dot(R, Paths[0].T).T
    Delta = np.abs(RSys - Paths[1])
    Cost = Delta.sum()
    return Cost


#%% For testing purpose
class Arguments:

    def __init__(self):
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    Time.Process(1, 'Set simulation')

    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)

    Data = pd.read_csv(str(DD / 'SampleList.csv'))

    # Step 1: compute best rotation for each sample
    MinData = pd.DataFrame(columns=['Angle','Cost'])
    for Index in Data.index:
        Sample = Data.loc[Index, 'Internal ID']
        SamplePath = RD / '03_hFE' / Sample
        
        # Read data
        FileName = str(SamplePath / (Sample + '.dat'))
        Simulation = Abaqus.ReadDAT(FileName)

        FileName = str(RD / '02_Experiment' / Sample / 'MatchedSignals.csv')
        Experiment = pd.read_csv(FileName)

        # Truncate experiment to monotonic loading part
        Peaks, Properties = sig.find_peaks(Experiment['FZ'], prominence=100)
        Start = Peaks[4]
        Stop = Experiment['FZ'].idxmin()
        Experiment = Experiment[Start:Stop].reset_index(drop=True)

        # Interpolate data to match hFE increments
        InterpData = pd.DataFrame()
        Coords = np.concatenate([np.array([0]),Experiment['Z']])
        for C in Experiment.columns:
            Values = np.concatenate([np.array([0]),Experiment[C]])
            Interpolated = np.interp(Simulation['Z'], Coords, Values)
            InterpData[C] = Interpolated

        # Compute rotation to minimize difference
        Paths = [InterpData[['X','Y','Z']].values,
                 Simulation[['X','Y','Z']].values]

        # Compute best initial guess
        C = 1E3
        for A in np.arange(0, 2*np.pi, np.pi/5):
            Ct = ComputeCost([A], Paths)
            if Ct < C:
                C = Ct
                Guess = A

        
        Optim = minimize(ComputeCost, [Guess], args=(Paths), bounds=([0, 2*np.pi],))
        
        # Store results
        MinData.loc[Index, 'Cost'] = Optim.fun
        MinData.loc[Index, 'Angle'] = Optim.x[0] / np.pi * 180


        # Plot results
        R = RotationMatrix(Psi=Optim.x[0])
        RSys = np.dot(R, Paths[0].T).T
        
        Figure, Axis = plt.subplots(1,2, sharex=True, sharey=True, figsize=(11,4.5))
        Axis[0].plot(Simulation['X'], -Simulation['Z'],color=(0,0,0))
        Axis[0].plot(InterpData['X'], -InterpData['Z'],color=(1,0,0))
        Axis[0].plot(RSys[:,0], -RSys[:,2],color=(0,0,1))
        Axis[0].set_xlabel('X (mm)')
        Axis[0].set_ylabel('Z (mm)')
        Axis[1].plot(Simulation['Y'], -Simulation['Z'],color=(0,0,0),label='hFE')
        Axis[1].plot(InterpData['Y'], -InterpData['Z'],color=(1,0,0), label='Experiment')
        Axis[1].plot(RSys[:,1], -RSys[:,2],color=(0,0,1), label='Rotated')
        Axis[1].set_xlabel('Y (mm)')
        Figure.legend(loc='upper center', ncol=3)
        plt.savefig(str(RD / '05_Comparison' / (Sample + '_XY')))
        plt.close()

    # Step 2: compute general best rotation angle



    for Index in Data.index:

        Sample = Data.loc[Index, 'Internal ID']
        SamplePath = RD / '03_hFE' / Sample
        
        # Read max displacement reached
        FileName = str(SamplePath / (Sample + '.dat'))
        Simulation = Abaqus.ReadDAT(FileName)
        MaxDisp = Simulation.loc[Simulation['Z'].idxmax(),'Z']

        # Replace in original loading BCs file
        FileName = str(SamplePath / 'boundary_conditions_FZ_MAX.inp')
        with open(FileName) as File:

            Text = File.read()
            Line = Text.split('\n')[7]
            
        NewLine = Line[:16] + str(round(MaxDisp - 0.005,2))
        NewText = Text.replace(Line, NewLine)

        with open(FileName, 'w') as File:
            File.write(NewText)

        # Write new BCs file for unloading
        FileName = str(SamplePath / 'boundary_conditions_FZ_MIN.inp')
            
        NewLine = Line[:16] + str(-round(MaxDisp - 0.005,2))
        NewText = Text.replace(Line, NewLine)

        with open(FileName, 'w') as File:
            File.write(NewText)

        # Add unloading step to simulation input file
        uCTFile = 'C000' + str(Data.loc[Index, 'MicroCT pretest file number']) + '_DOWNSCALED_00_FZ_MAX.inp'
        SimFile = str(SamplePath / uCTFile)

        with open(SimFile) as File:
            Text = File.read()
            Start = Text.find('*STEP')
            NewText = Text[Start:]
        
        NewText = '**\n' + NewText
        NewText = NewText.replace('boundary_conditions_FZ_MAX.inp',FileName[-30:])
        with open(SimFile,'a') as File:
            File.write(NewText)


        ProgressNext(Index / len(Data) * 10)

    ProcessTiming(0)

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
    Parser.add_argument('-F', '--Folder', help='Root folder of the project', default='FRACTIB', type=str)

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)