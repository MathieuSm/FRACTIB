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
for i, S in enumerate(Scripts):
    print(str(i) + ' ' + S)


#%%
SList = [Scripts[4]]
Arguments = {
             'hFERunSimulation':{'UMAT':'UMAT.f','nCPUs':12},
             'Registration':{'Show':False,'Type':'BSpline','Jac':True},
            }


#%%
File = str(SD / 'Submit.py')
print('\n\nWrite submit file:\n')
with open(File, 'w') as F:

    F.write('#! /usr/bin python\n\n')
    
    for Script in SList:
        F.write('from ' + Script[:-3] + ' import Main as ' + Script[:-3] + '\n')

    F.write('\nSamples = [\'' + Samples[0] + '\',\n')
    for Sample in Samples[1:-1]:
        F.write('           \'' + Sample + '\',\n')
    F.write('           \'' + Samples.values[-1] + '\']\n')
    
    F.write('\nclass Arguments():\n')
    F.write('\n\tdef __init__(self):\n')
    F.write('\t\tself.Folder = \'FRACTIB\'\n')
    F.write('\nArguments = Arguments()\n')

    F.write('\nfor Sample in Samples:\n')
    F.write('\n\tArguments.Sample = Sample\n')

    for Script in SList:
        
        if Script[:-3] in Arguments.keys():
            SArguments = Arguments[Script[:-3]]
        else:
            SArguments = False

        if SArguments:
            
            for iKey, Key in enumerate(SArguments):
                Item = SArguments[Key]
                Line = '\tArguments.' + Key + ' = '
                if type(Item) == str:
                    Line += '\'' + Item + '\'\n'
                else:
                    Line += str(Item) + '\n'
                F.write(Line)

        F.write('\n\t' + Script[:-3] + '(Arguments)\n\n')

with open(File) as F:
    print(F.read())

# %%
