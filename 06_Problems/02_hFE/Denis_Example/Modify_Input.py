import os

# Define directory and parameters
WorkingDirectory = os.getcwd()
FilePath = os.path.join(WorkingDirectory,'06_Problems/02_hFE/Denis_Example/')

# Read file
FileName = 'C0003103_V_11_iso_cortex_FZ_MAX.inp'
File = open(FilePath + FileName,'r')
Text = File.read()

# Find first DEPVAR
DEPVAR_Start = Text.find('*DEPVAR') + len('*DEPVAR') + 1

Start = DEPVAR_Start
Stop = Start + 2

# Change DEPVAR to add the 9 components of F
NewText = Text[:Start]
DEPVAR = int(Text[Start:Stop])
NewText += str(DEPVAR + 9)

# Add original and new DEPVAR
Start = Stop
Stop = Text[Start:].find('OF\n') + len('OF\n')
NewText += Text[Start: Start + Stop]
NewDEPVAR = DEPVAR
k = 1
for i in range(1,4):
    for j in range(1,4):
        Text2Add = str(NewDEPVAR+k) + ', F' + str(i) + str(j)
        if i < 3 or j < 3:
            Text2Add += ',\n'
        else:
            Text2Add += '\n'
        k += 1
        NewText += Text2Add

# Add element info and modified text (loop for every element)
ModifiedText = NewText[DEPVAR_Start:]
Start = DEPVAR_Start + Text[DEPVAR_Start:].find('*')
Test = Text[Start:].find('*DEPVAR')

while Test > 0:

    Stop = Start + Text[Start:].find('*DEPVAR') + len('*DEPVAR') + 1
    NewText += Text[Start:Stop]
    NewText += ModifiedText

    Start += Stop - Start + 98
    Test = Text[Start:].find('*DEPVAR')

# Add note sets and element sets
STATEV_Start = Text[Start:].find('SDV22,\n') + Start + len('SDV22,\n')
NewText += Text[Start:STATEV_Start]

# Add SDV output request
for k in range(9):
    NewText += 'SDV' + str(23 + k) + ',\n'
NewText += Text[STATEV_Start:]

# Replace username
NewText = NewText.replace('schenk','ms20s284')

# Write result
NewFile = open(FilePath + FileName[:-4] + '_New.inp','w')
NewFile.write(NewText)
NewFile.close()
File.close()