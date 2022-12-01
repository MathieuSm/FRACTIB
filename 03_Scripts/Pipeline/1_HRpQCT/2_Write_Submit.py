import os

# Set path
CurrentDirectory = os.getcwd()
DataFolder = os.path.join(CurrentDirectory,'02_Data/02_HRpQCT/')
ScriptFolder = os.path.join(CurrentDirectory,'03_Scripts/')

Command = '@DISK2:[MICROCT.DATA.MATHIEU.FRACTIB]10_IPL_COMMON_IMAGE_SIZE.COM'

SampleFolders = os.listdir(DataFolder)

SubmitScript = open(ScriptFolder + '12_SUBMIT_IPL_COMMON_IMAGE_SIZE.COM','w')
for SampleFolder in SampleFolders:

    Folder = os.path.join(DataFolder,SampleFolder)
    Files = os.listdir(Folder)
    SampleNumber = Files[0][:8]

    FullCommand = Command + ' ' + SampleFolder + ' ' + SampleNumber
    SubmitScript.write(FullCommand + '\n')
SubmitScript.close()
