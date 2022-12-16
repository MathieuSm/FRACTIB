#%%
#! /usr/bin python
import pandas as pd
from pathlib import Path

def SetDirectories(Name):

    CWD = str(Path.cwd())
    Start = CWD.find(Name)
    WD = Path(CWD[:Start], Name)
    Data = WD / '02_Data'
    Scripts = WD / '03_Scripts'
    Results = WD / '04_Results'

    return WD, Data, Scripts, Results

CWD, DD, SD, RD = SetDirectories('FRACTIB')

Samples = pd.read_csv(str(DD / 'SampleList.csv'))

# Location of the script to execute
FileLocation = '@DISK2:[MICROCT.DATA.MATHIEU]'

# Different scripts that can be executed
FileNames = ['IPL_DOWNSCALE.COM', 'IPL_STD_EVAL.COM']

# Text to add at the end of the line (parameters) according to the script
Additional = [' 3\n', '_DOWNSCALED 2 5\n']


iFile = 1


File = str(CWD / ('SUBMIT_' + FileNames[iFile]))
print('\n\nWrite submit file')
with open(File, 'w') as F:
    for Index in Samples.index:
        Line = FileLocation + FileNames[iFile]
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT ID'])
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT pretest folder'])
        Line += ' ' + 'C%07d'%(Samples.loc[Index,'MicroCT pretest file number'])
        F.write(Line + Additional[iFile])

        Line = FileLocation + FileNames[iFile]
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT ID'])
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT posttest folder'])
        Line += ' ' + 'C%07d'%(Samples.loc[Index,'MicroCT posttest file number'])

        F.write(Line + Additional[iFile])


F.close()
print('Done!')

# %%
