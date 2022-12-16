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

FileLocation = '@DISK2:[MICROCT.DATA.MATHIEU]'
FileName = 'IPL_DOWNSCALE.COM'


File = str(CWD / ('SUBMIT_' + FileName))
print('\n\nWrite submit file')
with open(File, 'w') as F:
    for Index in Samples.index:
        Line = FileLocation + FileName
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT ID'])
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT pretest folder'])
        Line += ' ' + 'C%07d'%(Samples.loc[Index,'MicroCT pretest file number'])
        F.write(Line + ' 3\n')

        Line = FileLocation + FileName
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT ID'])
        Line += ' ' + '%08d'%(Samples.loc[Index,'MicroCT posttest folder'])
        Line += ' ' + 'C%07d'%(Samples.loc[Index,'MicroCT posttest file number'])
        F.write(Line + ' 3\n')

F.close()
print('Done!')

# %%
