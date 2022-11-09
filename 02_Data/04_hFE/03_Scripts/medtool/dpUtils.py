
# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /home/ben/Software/medtool41/bin/dpUtils.py
# Compiled at: 2017-02-22 09:45:58
"""
########################################################################
File   : dpUtils.py
Author : D.H.Pahr
Date   : 15.12.2015

Copyright (C) Dr. Pahr Ingenieurs e.U. All rights reserved.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage: This is a library only, ready for Python 3

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requirements: None

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Description: This is a utility library and implemented following functions:

 - initialize & creates data dictionary from command line arguments
 - creates a doc string from a medtool gui strings 
 - throwing formated error and warning strings 
 - handles progress output on a shell 
 - splits strings by ":" or ";" 
  
########################################################################
"""
from sys import stdout
import os.path
import re

def getDataDict(guiInit):
    """
    Function converts gui string into a dictonary
    """
    data = {}
    dataDict = {}
    splitLineList = guiInit.split('\n')
    subwindow = 'base'
    pos = 0
    for splitLine in splitLineList:
        splitList = splitLine.split("'")
        if len(splitList) > 0:
            if len(splitList[0]) == 0:
                pass
            elif splitList[0].find('*modulName') == 0:
                pass
            elif splitList[0].find('*subWindow') == 0:
                if splitList[0].find('*subWindowStart') == 0:
                    subwindow = splitList[1]
            elif len(splitList) == 3:
                split1 = splitList[0].split()
                split2 = splitList[1]
                split3 = splitList[2].split()
                data = {}
                pos += 1
                description = split2
                paramText = ''
                paramTextEnd = ''
                if split2.find('|') > -1:
                    description, paramText = split2.split('|')
                    description = description.rstrip()
                    paramText = paramText.lstrip()
                    paramText = paramText + ' ('
                    paramTextEnd = ')'
                data['description'] = description
                data['entry'] = split1[0]
                data['default'] = split3[0]
                data['optional'] = split3[1]
                data['info'] = split3[2]
                data['value'] = None
                data['position'] = pos
                data['type'] = paramText + 'str' + paramTextEnd
                if data['entry'] == '*entry' or data['entry'] == '*splitEntry' or data['entry'] == '*colorEntry':
                    data['lvalue'] = convertValues(data['default'], split1[1], valueType=data['info'])
                    if data['default'].find(':') > -1 or data['default'].find(';') > -1:
                        data['type'] = paramText + 'list[' + data['info'] + ']' + paramTextEnd
                    else:
                        data['type'] = paramText + data['info'] + paramTextEnd
                    if data['type'] == 'strfil':
                        data['type'] = paramText + 'str' + paramTextEnd
                else:
                    data['lvalue'] = data['default']
                data['subwindow'] = subwindow
                dataDict[split1[1]] = data
            else:
                stdout.write('\n **ERROR** gui string not correct!\n\n')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)

    return dataDict


def replaceDocString(filename):
    try:
        fobj = open(filename, 'r')
    except IOError:
        print 'cannot open', filename
    else:
        lines = fobj.readlines()
        fobj.close()

    print type(lines)
    if lines[0].find('"""') == 0:
        pline = lines.pop(0)
        newlines = lines[:]
        for line in lines:
            if line.find('"""') == 0:
                break
            else:
                print 'POP', newlines.pop(0)

        print '======='
        print newlines


def readDocTxt(filename):
    """
    Read extra doc data from external file. 
    """
    text = {}
    fobj = open(filename, 'r')
    readArg = None
    for line in fobj:
        if line.find(':') == 0:
            readArg = line.strip().replace(':', '')
            text[readArg] = []
            continue
        if readArg:
            text[readArg].append(line)

    fobj.close()
    return text


def getDocString(guiInit, scriptName, scriptAuthor):
    """
    Create an automated docstring, add additional info from an external file.
    """
    fileName = scriptName.replace('Src.py', 'Doc.rst')
    addText = {}
    if os.path.isfile(fileName):
        addText = readDocTxt(fileName)
    dataDict = getDataDict(guiInit)
    dataDict['-gui'] = {}
    dataDict['-gui']['default'] = ''
    dataDict['-gui']['description'] = 'Print gui string on stdout'
    dataDict['-gui']['type'] = ''
    dataDict['-gui']['optional'] = 'yes'
    dataDict['-gui']['position'] = len(dataDict)
    dataDict['-doc'] = {}
    dataDict['-doc']['default'] = ''
    dataDict['-doc']['description'] = 'Print default doc string (__doc__) on stdout.'
    dataDict['-doc']['type'] = ''
    dataDict['-doc']['optional'] = 'yes'
    dataDict['-doc']['position'] = len(dataDict)
    dataDict['-help'] = {}
    dataDict['-help']['default'] = ''
    dataDict['-help']['description'] = 'Print help string (extended __doc__) on stdout.'
    dataDict['-help']['type'] = ''
    dataDict['-help']['optional'] = 'yes'
    dataDict['-help']['position'] = len(dataDict)
    posList = []
    for i in range(len(dataDict)):
        posList.append('')

    maxLenKey = 0
    maxLenDef = 0
    for key in dataDict:
        posList[dataDict[key]['position'] - 1] = key
        if len(key) > maxLenKey:
            maxLenKey = len(key)
        if len(dataDict[key]['default']) > maxLenDef:
            maxLenDef = len(dataDict[key]['default'])

    maxSpaces = ''
    for i in range(maxLenKey + 1):
        maxSpaces += ' '

    optional = [
     'no', 'yes']
    docstrUsage = ''
    docstrDescr = ''
    for opt in optional:
        bracketL = ' '
        bracketR = ' '
        if opt == 'no':
            docstrDescr += '\nParameters\n'
            docstrDescr += '~~~~~~~~~~\n\n'
        else:
            docstrDescr += '\nOptional Parameters\n'
            docstrDescr += '~~~~~~~~~~~~~~~~~~~\n\n'
            bracketL = '['
            bracketR = ']'
        for key in posList:
            if dataDict[key]['optional'] == opt:
                spaces1 = ''
                spaces2 = ''
                for i in range(maxLenKey - len(key)):
                    spaces1 += ' '

                for i in range(maxLenDef - len(dataDict[key]['default'])):
                    spaces2 += ' '

                docstrUsage += '    ' + bracketL + key + spaces1 + ' ' + dataDict[key]['default'] + spaces2 + bracketR + '\n'
                docstrDescr += '' + key + spaces1 + '  : ' + dataDict[key]['type'] + '\n\n'
                docstrDescr += '   ' + maxSpaces + dataDict[key]['description'] + '.\n'
                if key in addText:
                    for line in addText[key]:
                        docstrDescr += '    ' + maxSpaces + line

                    docstrDescr += '\n\n'
                else:
                    docstrDescr += '\n\n'

    path, scriptNameShort = os.path.split(scriptName.replace('Src', ''))
    docStr = ''
    docStr += '"""\n'
    title = 'Script ' + scriptNameShort.replace('.py', '')
    docStr += title + '\n'
    docStr += createUnderline(title, '-')
    docStr += '\n\n'
    if 'DESCRIPTION' in addText:
        for line in addText['DESCRIPTION']:
            docStr += line

        docStr += '\n'
    else:
        docStr += 'This scripts is ...\n\n'
    docStr += '\nUsage\n'
    docStr += '~~~~~\n\n'
    docStr += 'Module: ::\n\n'
    docStr += '  import ' + scriptNameShort.replace('.py', '') + '\n\n'
    docStr += 'Command line: ::\n\n'
    docStr += '  python ' + scriptNameShort + ' ... \n\n'
    docStr += docstrUsage
    docStr += '\n'
    docStr += '\n'
    docStr += docstrDescr
    docStr += '\nInfo\n'
    docStr += '~~~~\n\n'
    docStr += '- File:   ' + scriptNameShort + '\n'
    docStr += '- Author: ' + scriptAuthor + '\n\n'
    docStr += '"""'
    return docStr


def createUnderline(string, char):
    underline = ''
    for i in string:
        underline += char

    return underline


def initializeArguments(argv, dataDict, guiInit, doc, doc2='!!!Docstring was not defined!!!'):
    """
    Functions initializes arguments
    """
    argList = argv
    argc = len(argList)
    i = 0
    while i < argc:
        if argList[i] in dataDict:
            arg = argList[i]
            i += 1
            if dataDict[arg]['entry'] == '*entry' or dataDict[arg]['entry'] == '*splitEntry' or dataDict[arg]['entry'] == '*colorEntry':
                dataDict[arg]['lvalue'] = convertValues(argList[i], arg, valueType=dataDict[arg]['info'])
            dataDict[arg]['value'] = argList[i]
        elif argList[i][:4] == '-gui':
            stdout.write('%s' % guiInit)
            stdout.flush()
            exit(0)
        elif argList[i][:5] == '-help':
            stdout.write('%s\n' % doc)
            stdout.flush()
            exit(0)
        elif argList[i][:4] == '-doc':
            stdout.write('%s\n' % doc2)
            stdout.flush()
            exit(0)
        i += 1

    for arg in dataDict:
        if dataDict[arg]['optional'].lower() == 'no' and dataDict[arg]['value'] == None:
            stdout.write(doc)
            stdout.flush()
            stdout.write("\n **ERROR** Option '%s' is required!\n\n" % arg)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

    return


def userSplit(oldString):
    if oldString.find(':') > -1 and oldString.find(';') > -1:
        throwError("Option value '%s' shows a not allowed mixture of  ':' and ';' delimiters!" % oldString)
    newString = oldString.replace(':', ';')
    findList = re.findall('[A-Z];\\\\', newString)
    for val in findList:
        newString = newString.replace(val[0] + ';', val[0] + ':')

    findList = re.findall('[A-Z];/', newString)
    for val in findList:
        newString = newString.replace(val[0] + ';', val[0] + ':')

    return newString.split(';')


def convertValues(value, option, valueType=None):
    """
    Function splits value, and converts to format if given
    """
    lvalue = None
    if value != None:
        valueList = userSplit(value)
        if len(valueList) == 1:
            if valueType.lower() == 'float':
                try:
                    lvalue = float(valueList[0])
                except:
                    throwError("Option '%s' can not convert value '%s' to float" % (option, valueList[0]))

            elif valueType == 'int':
                try:
                    lvalue = int(valueList[0])
                except:
                    throwError("Option '%s' can not convert value '%s' to int" % (option, valueList[0]))

            elif valueType == 'strfil':
                from os import path
                if path.isfile(valueList[0]):
                    try:
                        file = open(valueList[0])
                        filecontents = file.read().rstrip('\n')
                        file.close()
                        lvalue = convertValues(filecontents, option, '?')
                    except IOError:
                        throwError("Option '%s' file '%s' exists but cannot be read" % (option, val))
                    except:
                        throwError('Unkonwn error??')

                else:
                    lvalue = convertValues(valueList[0], option, '?')
            else:
                lvalue = valueList[0]
        else:
            lvalue = []
            for val in valueList:
                if valueType.lower() == 'float':
                    try:
                        val1 = float(val)
                    except:
                        throwError("Option '%s' can not convert value '%s' to float" % (option, val))

                    lvalue.append(val1)
                elif valueType == 'int':
                    try:
                        val1 = int(val)
                    except:
                        throwError("Option '%s' can not convert value '%s' to float" % (option, val))

                    lvalue.append(val1)
                else:
                    lvalue.append(val)

    return lvalue


def initValues(value, valName='unknown', minVal=0, valType=''):
    """
    Function initilizes argument values i.e. check the number of values
    and returns a list of the right format.
    """
    if value != None:
        sValList = []
        sValStr = userSplit(value)
        if minVal != 0 and len(sValStr) != minVal:
            throwError("dpUtils.initValues - Option '%s' needs '%i' values!" % (valName, minVal))
        for sVal1 in sValStr:
            if valType == 'float':
                sValList.append(float(sVal1))
            elif valType == 'int':
                sValList.append(int(sVal1))
            else:
                sValList.append(sVal1)

        return sValList
    else:
        return
        return


def throwWarning(string):
    stdout.write('\n **WARNING** %s\n\n' % string)
    stdout.flush()


def throwError(string):
    stdout.write('\n **ERROR** %s\n\n' % string)
    stdout.write('\n E N D E D  with ERRORS \n\n')
    stdout.flush()
    exit(1)


def progressStart(text):
    global curProgress
    stdout.write(text + '|')
    curProgress = 0
    stdout.flush()


def progressNext(progress):
    global curProgress
    if progress > curProgress:
        curProgress += 1
        stdout.write('=')
        stdout.flush()


def progressIncrease():
    stdout.write('=')
    stdout.flush()


def progressEnd():
    stdout.write('|\n')
    stdout.flush()


if __name__ == '__main__':
    print
    print ' THIS IS A PYTHON MODULE ... ONLY'
    print
    print ' Not callable from the commandline ...'
    print
# okay decompiling dpUtils.pyc