#%% #!/usr/bin/env python3
# Initialization

Version = '01'

Description = """
    Script aimed to analyse the relationship of structural variables
    (force and stiffness namely) between experiment and simulation

    Version Control:
        01 - Original script

    Author: Mathieu Simon
            ARTORG Center for Biomedical Engineering Research
            SITEM Insel, University of Bern

    Date: January 2023
    """

#%% Imports
# Modules import

import argparse
from Utils import *
Show.ShowPlot = False
from scipy.optimize import minimize

#%% Functions
# Define functions

def RotateSystem(Angle, Ref, Sys):
    M = RotationMatrix(Gamma=Angle[0])
    rSys = np.dot(M, Sys.T).T
    Cost = np.linalg.norm(Ref - rSys, axis=0).sum()
    return Cost
def ComputeStiffness(Force, Displacement, RelativeRange=1/3, StepSize=None):

    Length = int(round(len(Force) * RelativeRange))
    if Length < 3:
        Length = 3

    if StepSize == None:
        StepSize = Length // 100
    
    if StepSize == 0:
        StepSize = 1

    Slopes = []
    for i in range(0, len(Force)-Length-StepSize+1, StepSize):
        Y = Force[i:i+Length]
        X = Displacement[i:i+Length]
        Show.FName = None
        Fit = Show.OLS(X, Y)
        Slopes.append(Fit.params['X'])

    Stiffness = max(Slopes)

    return Stiffness

#%% Classes
# Define classes

class Arguments():

    def __init__(self):
        self.Folder = 'FRACTIB'

Arguments = Arguments()

#%% Main
# Main code

def Main(Arguments):

    # Set directories and read sample list
    CWD, DD, SD, RD = SetDirectories(Arguments.Folder)
    SampleList = pd.read_csv(str(DD / 'SampleList.csv'))
    ResultsDir = RD / '05_Comparison'

    # Set assessed variables
    Variables = ['Force','Displacement at Fmax','Stiffness']
    DataSets = ['hFE', 'Experiment']
    Columns = pd.MultiIndex.from_product([Variables, DataSets])
    Data = pd.DataFrame(index=SampleList.index, columns=Columns)
    Data['Cost'] = np.nan
    Data['Angle'] = np.nan

    # Compute general best rotation angle
    Angles = pd.DataFrame(columns=np.arange(0, 360, 10))
    Time.Process(1, 'Compute rotations')
    for Index in SampleList.index:

        Sample = SampleList.loc[Index, 'Internal ID']
        FEADir = RD / '03_hFE' / Sample
        ExpDir = RD / '02_Experiment' / Sample

        # Read files
        FEAData = Abaqus.ReadDAT(str(FEADir / (Sample + '.dat')))
        ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))

        # Rotate system to align coordinates
        Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
        Start = Peaks[4]
        Stop = np.argmin(np.abs(ExpData['Z'] - FEAData['Z'].max()))
        RefData = ExpData['Z'][Start:Stop]
        FEStop = FEAData['Z'].idxmax()
        Interp = FEAData['Z'][:FEStop]
        InterpX = np.interp(RefData, Interp, FEAData['X'][:FEStop])
        InterpY = np.interp(RefData, Interp, FEAData['Y'][:FEStop])

        Ref = np.array([-InterpX, InterpY, RefData]).T
        Sys = np.array(ExpData[['X','Y','Z']][Start:Stop])

        for Angle in range(36):
            Angle = Angle * np.pi/18
            RotatedSystem = np.dot(RotationMatrix(Gamma=Angle), Sys.T).T
            Cost = np.abs(Ref - RotatedSystem).sum()
            Angle = int(round(Angle / np.pi * 180))
            Angles.loc[Index, Angle] = Cost

            # print(Angle)

            # Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=False, sharey=True)
            # Axis[0].plot(ExpData['X'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,1), label='Original')
            # Axis[0].plot(RotatedSystem[:,0], -RotatedSystem[:,2], color=(0,0,0), label='Experiment')
            # Axis[0].plot(FEAData['X'][:FEStop], -FEAData['Z'][:FEStop], color=(1,0,0), label='hFE')
            # Axis[1].plot(ExpData['Y'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,1))
            # Axis[1].plot(RotatedSystem[:,1], -RotatedSystem[:,2], color=(0,0,0))
            # Axis[1].plot(FEAData['Y'][:FEStop], -FEAData['Z'][:FEStop], color=(1,0,0))
            # Axis[0].set_xlabel('X (mm)')
            # Axis[0].set_ylabel('Z (mm)')
            # Axis[1].set_xlabel('Y (mm)')
            # Figure.legend(loc='upper center', ncol=3)
            # plt.show()

        # Angles.loc[Index] -= Angles.loc[Index].max()
        # Angles.loc[Index] = Angles.loc[Index].abs() / Angles.loc[Index].abs().max()
        Time.Update(Index / len(SampleList))
    
    Angle = Angles.T.sum(axis=1).idxmin()/180 * np.pi
    # Angle = 230/180 * np.pi
    Time.Process(0)

    # Compute variables from experiment and simulation    
    for Index in SampleList.index:

        Text = '\nGet structural results for sample ' + str(Index+1) + '/' + str(len(SampleList))
        ProcessTiming(1, Text)

        Sample = SampleList.loc[Index, 'Internal ID']
        FEADir = RD / '03_hFE' / Sample
        ExpDir = RD / '02_Experiment' / Sample

        # Read files
        FEAData = Abaqus.ReadDAT(str(FEADir / (Sample + '.dat')))
        ExpData = pd.read_csv(str(ExpDir / 'MatchedSignals.csv'))

        # Rotate system to align coordinates
        Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
        Start = Peaks[4]
        Stop = np.argmin(np.abs(ExpData['Z'] - FEAData['Z'].max()))
        RefData = ExpData['Z'][Start:Stop]
        FEStop = FEAData['Z'].idxmax()
        Interp = FEAData['Z'][:FEStop]
        InterpX = np.interp(RefData, Interp, FEAData['X'][:FEStop])
        InterpY = np.interp(RefData, Interp, FEAData['Y'][:FEStop])

        Ref = np.array([-InterpX, InterpY, RefData]).T
        Sys = np.array(ExpData[['X','Y','Z']][Start:Stop])

        Angle = 10 * np.argmin(Angles.loc[Index]) / 180 * np.pi
        Result = minimize(RotateSystem, [Angle], args=(Ref, Sys), bounds=([0, 2*np.pi],))
        Angle = Result.x[0]

        RotatedSystem = np.dot(RotationMatrix(Gamma=Angle), Sys.T).T
        Data.loc[Index, 'Cost'] = np.linalg.norm(Ref - RotatedSystem, axis=0).sum()
        Data.loc[Index, 'Angle'] = Angle / np.pi * 180

        UFD = -RotatedSystem[ExpData['FZ'].idxmin() - Start]
        X1 = np.min([ExpData['X'][Start:Stop].min(), RotatedSystem[:,0].min(), FEAData['X'][:FEStop].min()])
        X2 = np.max([ExpData['X'][Start:Stop].max(), RotatedSystem[:,0].max(), FEAData['X'][:FEStop].max()])
        Y1 = np.min([ExpData['Y'][Start:Stop].min(), RotatedSystem[:,1].min(), FEAData['Y'][:FEStop].min()])
        Y2 = np.max([ExpData['Y'][Start:Stop].max(), RotatedSystem[:,1].max(), FEAData['Y'][:FEStop].max()])
        Figure, Axis = plt.subplots(1,2, figsize=(11,4.5), sharex=False, sharey=True)
        Axis[0].plot([X1, X2],[UFD[2],UFD[2]],color=(0.7,0.7,0.7), linestyle='--', label='Ultimate force')
        Axis[0].plot(ExpData['X'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,1), label='Original')
        Axis[0].plot(RotatedSystem[:,0], -RotatedSystem[:,2], color=(0,0,0), label='Experiment')
        Axis[0].plot(FEAData['X'][:FEStop], -FEAData['Z'][:FEStop], color=(1,0,0), label='hFE')
        Axis[1].plot([Y1, Y2],[UFD[2],UFD[2]],color=(0.7,0.7,0.7), linestyle='--')
        Axis[1].plot(ExpData['Y'][Start:Stop], -ExpData['Z'][Start:Stop], color=(0,0,1))
        Axis[1].plot(RotatedSystem[:,1], -RotatedSystem[:,2], color=(0,0,0))
        Axis[1].plot(FEAData['Y'][:FEStop], -FEAData['Z'][:FEStop], color=(1,0,0))
        Axis[0].set_xlabel('X (mm)')
        Axis[0].set_ylabel('Z (mm)')
        Axis[1].set_xlabel('Y (mm)')
        Figure.legend(loc='upper center', ncol=4)
        plt.savefig(str(ResultsDir / (Sample + '_XY')))
        plt.close(Figure)

        # Plot force-displacement and min force point
        Force = FEAData['FZ'][1:]
        MinForceIdx = Force[Force > 0].idxmin()

        Xs = [FEAData['Z'][:MinForceIdx+1], ExpData['Z']]
        Ys = [FEAData['FZ'][:MinForceIdx+1], -ExpData['FZ']]
        Axes = ['Displacement (mm)', 'Force(N)']
        Labels = ['hFE', 'Experiment']

        Show.FName = str(ResultsDir / (Sample + '_FD'))
        Show.Signal(Xs, Ys, Axes=Axes, Labels=Labels, Points=[MinForceIdx,[]])
        Show.FName = None
        ProgressNext(1)

        # Max force
        Data.loc[Index, Variables[0]] = [max(FEAData['FZ']), max(ExpData['FZ'].abs())]

        # Displacement at max force
        FEAIndex = FEAData['FZ'].idxmax()
        hFE_F = FEAData.loc[FEAIndex,'Z']
        ExpIndex = ExpData['FZ'].idxmin()
        Exp_F = ExpData.loc[ExpIndex,'Z']
        Data.loc[Index, Variables[1]] = [hFE_F, Exp_F]
        ProgressNext(2)

        # Stiffness
        F = FEAData['FZ'].values[:FEAIndex]
        D = FEAData['Z'].values[:FEAIndex]
        FEAStiff = ComputeStiffness(F, D)
        ProgressNext(3)
        
        Peaks, Properties = sig.find_peaks(ExpData['FZ'], prominence=100)
        F = ExpData['FZ'].values[Peaks[4]:ExpIndex]
        D = ExpData['Z'].values[Peaks[4]:ExpIndex]
        ExpStiff = ComputeStiffness(-F, D)
        Data.loc[Index, Variables[2]] = [FEAStiff, ExpStiff]
        ProcessTiming(0)

    # Assess relationships (linear regression)
    X = Data[Variables[0]]['Experiment'].values.astype(float) / 1E3
    Y = Data[Variables[0]]['hFE'].values.astype(float) / 1E3
    Labels = ['hFE (kN)', 'Experiment (kN)']
    Show.FName = str(RD / 'UltimateForce')
    Results = Show.OLS(X, Y, Labels=Labels)

    X = Data[Variables[1]]['Experiment'].values.astype(float)
    Y = Data[Variables[1]]['hFE'].values.astype(float)
    Labels = ['hFE (mm)', 'Experiment (mm)']
    Show.FName = str(RD / 'DispUltForce')
    Results = Show.OLS(X, Y, Labels=Labels)

    X = Data[Variables[2]]['Experiment'].values.astype(float) / 1E3
    Y = Data[Variables[2]]['hFE'].values.astype(float) / 1E3
    Labels = ['hFE (kN/mm)', 'Experiment (kN/mm)']
    Show.FName = str(RD / 'Stiffness')
    Results = Show.OLS(X, Y, Labels=Labels)

    # Show rotation angles to fis systems (must be constant)
    Show.FName = str(RD / 'Cost')
    Show.BoxPlot([Data['Cost']], Labels=['','Cost (-)'])
    Show.FName = None

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
    Parser.add_argument('--Folder', help='Root folder of the project (required)', type=str, default='FRACTIB')

    # Read arguments from the command line
    Arguments = Parser.parse_args()

    Main(Arguments)