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

# lftp -u microct,***** 130.92.125.44 -e "cd DISK2:[MICROCT.DATA.00000269.00002159]; get C0001901.AIM; cd DISK2:[MICROCT.DATA.00000289.00002187]; get C0001929.AIM; cd DISK2:[MICROCT.DATA.00000277.00002169]; get C0001911.AIM; cd DISK2:[MICROCT.DATA.00000279.00002173]; get C0001915.AIM; cd DISK2:[MICROCT.DATA.00000276.00002168]; get C0001910.AIM; cd DISK2:[MICROCT.DATA.00000281.00002176]; get C0001918.AIM; cd DISK2:[MICROCT.DATA.00000288.00002185]; get C0001927.AIM; cd DISK2:[MICROCT.DATA.00000290.00002188]; get C0001930.AIM; cd DISK2:[MICROCT.DATA.00000275.00002167]; get C0001909.AIM; cd DISK2:[MICROCT.DATA.00000280.00002174]; get C0001916.AIM; cd DISK2:[MICROCT.DATA.00000267.00002157]; get C0001899.AIM; cd DISK2:[MICROCT.DATA.00000271.00002161]; get C0001903.AIM; cd DISK2:[MICROCT.DATA.00000273.00002164]; get C0001906.AIM; cd DISK2:[MICROCT.DATA.00000272.00002162]; get C0001904.AIM; cd DISK2:[MICROCT.DATA.00000268.00002158]; get C0001900.AIM; cd DISK2:[MICROCT.DATA.00000278.00002170]; get C0001912.AIM; cd DISK2:[MICROCT.DATA.00000274.00002165]; get C0001907.AIM; cd DISK2:[MICROCT.DATA.00000292.00002191]; get C0001933.AIM; cd DISK2:[MICROCT.DATA.00000282.00002177]; get C0001919.AIM; cd DISK2:[MICROCT.DATA.00000291.00002189]; get C0001931.AIM; cd DISK2:[MICROCT.DATA.00000284.00002180]; get C0001922.AIM; cd DISK2:[MICROCT.DATA.00000285.00002181]; get C0001923.AIM; cd DISK2:[MICROCT.DATA.00000286.00002183]; get C0001925.AIM; cd DISK2:[MICROCT.DATA.00000283.00002178]; get C0001920.AIM; cd DISK2:[MICROCT.DATA.00000287.00002184]; get C0001926.AIM; bye"

# %%

Command = 'lftp -u'
User = 'microct'
Password = 'made1in2ch'
Address = '130.92.125.44'

MainDir = 'DISK2:[MICROCT.DATA]'
ServerDir = '/home/ms20s284/FRACTIB/02_Data/02_uCT/'

File = str(SD / 'Scanco' / 'CopyFiles.sh')
print('\n\nWrite submit file')
with open(File, 'w') as F:

    F.write('#!/bin/sh\n\n')

    for Index in Samples.index:
        
        F.write('cd ' + ServerDir + Samples.loc[Index,'Internal ID'] + ' \n')
        
        FullCommand = Command + ' ' + User + ',' + Password + ' ' + Address + ' -e '
        
        OpenVMS = '\"cd ' + MainDir[:-1]
        OpenVMS += '.' + '%08d'%(Samples.loc[Index,'MicroCT ID'])
        OpenVMS += '.' + '%08d'%(Samples.loc[Index,'MicroCT pretest folder'])
        OpenVMS += ']; get'
        OpenVMS += ' ' + 'C%07d'%(Samples.loc[Index,'MicroCT pretest file number'])
        OpenVMS += '.AIM; '

        OpenVMS += 'cd ' + MainDir[:-1]
        OpenVMS += '.' + '%08d'%(Samples.loc[Index,'MicroCT ID'])
        OpenVMS += '.' + '%08d'%(Samples.loc[Index,'MicroCT posttest folder'])
        OpenVMS += ']; get'
        OpenVMS += ' ' + 'C%07d'%(Samples.loc[Index,'MicroCT posttest file number'])
        OpenVMS += '.AIM; bye\"\n'

        F.write(FullCommand + OpenVMS)

        F.write('cd ..\n')


F.close()
print('Done!')

# %%
