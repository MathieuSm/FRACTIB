#%%
#! /usr/bin python
import pandas as pd
from Utils import *

CWD, DD, SD, RD = SetDirectories('FRACTIB')

Samples = pd.read_csv(str(DD / 'SampleList.csv'))['Internal ID']
Scripts = [File for File in os.listdir(str(SD)) if File.endswith('.py')]
Arguments = [None]

File = str(CWD / 'Submit.py')
print('\n\nWrite submit file')
with open(File, 'w') as F:
    F.write('#!/usr/bin python\n')
    F.write('import os\n\n')
    for Sample in [Samples[0]]:
        for iScript, Script in enumerate([Scripts[1]]):

            Line = 'os.system(\'' + str(SD / Script) + ' ' + Sample

            if Arguments[iScript]:
                Line += ' ' + Arguments[iScript]

            F.write(Line + '\')\n')

F.close()
print('Done!')

# %%
