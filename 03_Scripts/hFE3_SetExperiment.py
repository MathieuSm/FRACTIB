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
from matplotlib import cm
from scipy.optimize import minimize

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

    # Step 1: compute best rotation angle
    Samples = Data['Internal ID']
    Angles = np.arange(0, 360, 10)
    Indices = pd.MultiIndex.from_product([Samples,['dX','dY']])
    Costs = pd.DataFrame(index=Indices, columns=Angles)
    for Sample in Samples:
        SamplePath = RD / '03_hFE' / Sample
        
        # Read data
        FileName = str(SamplePath / 'MaxLoad.dat')
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

        # Compute cost after rotation
        for A in Angles:
            R = RotationMatrix(Psi=A/180*np.pi)
            RSys = np.dot(R, InterpData[['X','Y','Z']].values.T).T
            Costs.loc[(Sample,'dX'), A] = np.array(np.abs(RSys[:,0] - Simulation['X']))
            Costs.loc[(Sample,'dY'), A] = np.array(np.abs(RSys[:,1] - Simulation['Y']))

        # Plot costs
        Show.Signal([Costs.loc[Sample,A]['dX'] for A in Angles])

    # Sum costs for each sample
    Sums = pd.DataFrame(index=Indices, columns=Angles)
    for I in Sums.index:
        for C in Sums.columns:
            Sums.loc[I,C] = sum(Costs.loc[I,C])
    Sums = Sums.groupby(level=0).sum()

    BestAngle = Sums.sum(axis=0).idxmin()
    # BestAngle = 90

    # Plot results
    Figure, Axis = plt.subplots(2,2, figsize=(12,10))
    Axis[0,0].imshow(Sums.T, cmap='viridis_r', aspect='auto')
    Axis[0,0].set_ylim([-0.5, len(Costs.columns)-0.5])
    Axis[0,0].axis('off')
    C = Sums.sum(axis=0).values
    CNorm = (C-C.min()) / (C.max()-C.min())
    Axis[0,1].barh(np.arange(len(Angles)), C, color=cm.viridis_r(CNorm), height=1)
    Axis[0,1].set_yticks(np.arange(len(Angles))[::2],Costs.columns[::2])
    Axis[0,1].set_xlim([C.min()*0.95, C.max()*1.025])
    Axis[0,1].set_ylim([-0.5, len(C)-0.5])
    Axis[0,1].set_xlabel('Cumulative Cost (mm)')
    Axis[0,1].set_ylabel('Angle (°)')
    Axis[0,1].yaxis.tick_right()
    Axis[0,1].yaxis.set_label_position('right')
    C = Sums.sum(axis=1).values
    Axis[1,0].bar(Samples, C, color=cm.viridis_r(C/C.max()), width=1)
    Axis[1,0].set_xlim([-0.5, len(C)-0.5])
    Axis[1,0].set_xticks(Samples, Samples, rotation=90)
    Axis[1,0].set_xlabel('Sample (-)')
    Axis[1,0].set_ylabel('Cumulative Cost (mm)')
    Axis[1,0].yaxis.tick_right()
    Axis[1,0].yaxis.set_label_position('right')
    plt.delaxes(Axis[1,1])
    # plt.savefig(str(RD / '05_Comparison' / 'Costs'))
    # plt.close()
    plt.show()


    for Index in Data.index:

        Sample = Data.loc[Index, 'Internal ID']
        SamplePath = RD / '03_hFE' / Sample

        Text = 'Sample #' + str(Sample[:3])
        Time.Update(Index / len(Data), Text)
        
        # Read data
        FileName = str(SamplePath / 'MaxLoad.dat')
        Simulation = Abaqus.ReadDAT(FileName)

        FileName = str(RD / '02_Experiment' / Sample / 'MatchedSignals.csv')
        Experiment = pd.read_csv(FileName)

        # Truncate experiment
        Peaks, Properties = sig.find_peaks(Experiment['FZ'], prominence=100)
        Start = Peaks[4]
        Stop = Peaks[5]
        Experiment = Experiment[Start:Stop].reset_index(drop=True)
        MaxForce = Experiment['FZ'].idxmin()

        # Rotate coordinate system
        R = RotationMatrix(Psi=BestAngle/180*np.pi)
        XYZ = np.dot(R, Experiment[['X','Y','Z']].values.T).T
        PTP = Experiment[['Phi', 'Theta', 'Psi']].values / 180 * np.pi
        Phi, Theta, Psi = PTP[:,0], PTP[:,1], PTP[:,2]
        Rs = RotationMatrix(Phi, Theta, Psi)
        rPTP = np.einsum('ij,ljk->lik',R,Rs)
        PTP = GetAngles(rPTP)
        PTP[:,2] -= PTP[0,2]

        # Plot rotation results
        Figure, Axis = plt.subplots(1,2, sharex=True, sharey=True, figsize=(11,4.5))
        Axis[0].plot(Simulation['X'], -Simulation['Z'],color=(0,0,0))
        Axis[0].plot(Experiment['X'][:MaxForce], -Experiment['Z'][:MaxForce],color=(1,0,0))
        Axis[0].plot(XYZ[:MaxForce,0], -XYZ[:MaxForce,2],color=(0,0,1))
        Axis[0].set_xlabel('X (mm)')
        Axis[0].set_ylabel('Z (mm)')
        Axis[1].plot(Simulation['Y'], -Simulation['Z'],color=(0,0,0),label='hFE')
        Axis[1].plot(Experiment['Y'][:MaxForce], -Experiment['Z'][:MaxForce],color=(1,0,0), label='Experiment')
        Axis[1].plot(XYZ[:MaxForce,1], -XYZ[:MaxForce,2],color=(0,0,1), label='Rotated')
        Axis[1].set_xlabel('Y (mm)')
        Figure.legend(loc='upper center', ncol=3)
        plt.savefig(str(RD / '05_Comparison' / (Sample + '_XY')))
        plt.close()

        # Modify abaqus input file
        uCTFile = 'C000' + str(Data.loc[Index, 'MicroCT pretest file number']) + '_DOWNSCALED_00_FZ_MAX.inp'
        FileName = str(SamplePath / uCTFile)
        FileName = str(SamplePath / 'Simulation.inp')
        Abaqus.RemoveSteps(FileName,1)

        NSteps = 50
        StepSize = len(XYZ) / NSteps
        DOFs = np.arange(6)+1
        Values = np.zeros(6)
        os.makedirs(str(SamplePath / 'Steps' ), exist_ok=True)
        for i in range(NSteps):
            StepName = str(SamplePath / 'Steps' / ('Step%02i.inp'%(i+1)))
            NewValues = np.concatenate([XYZ[int(StepSize * (i+1)-1)],
                                        PTP[int(StepSize * (i+1)-1)]])
            Delta = NewValues - Values
            Abaqus.AddStep(FileName,StepName,DOFs,np.round(Delta,3))
            Values = NewValues


    Time.Process(0)

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