#%%
#! /usr/bin python
import os
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

Samples = pd.read_csv(str(DD / 'SampleList.csv'))['Internal ID']
Scripts = [File for File in os.listdir(str(SD)) if File.endswith('.py')]

SList = [Scripts[1]]
Arguments = [{'Type':'BSpline', 'Jac':True}]

File = str(CWD / 'Submit.py')
print('\n\nWrite submit file')
with open(File, 'w') as F:

    F.write('#! /usr/bin python\n\n')
    
    for Script in SList:
        F.write('from ' + Script[:-3] + ' import Main as ' + Script[:-3] + '\n')

    F.write('\nSamples = [')
    for Sample in Samples[:-1]:
        F.write(Sample + ',\n')
    F.write(Samples[-1] + ']\n')
    
    F.write('\nclass Arguments():\n')
    F.write('\n\tdef __init__(self):\n')
    F.write('\t\tpass\n')
    F.write('\nArguments = Arguments()')

    F.write('for Sample in Samples:\n')
    for Sample in Samples:

        F.write('\n\tArguments.Sample = ' + Sample + '\n')

        for iScript, Script in enumerate(SList):

            if Arguments[iScript]:
                
                F.write()

        Line = 'python ' + str(SD / Script) + ' ' + Sample



        F.write(Line + '\n')

F.close()
print('Done!')

# %%
