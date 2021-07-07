# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/misc/dpStiffness.py
# Compiled at: 2015-08-20 15:14:50
"""
#########################################################################
File   : dpStiffness.py
Author : D.H.Pahr
Date   : 1.12.2005
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage: python dpStiffness.py   
                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requirements:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Description: 
  This program is a stiffness class.  
  
#########################################################################
"""
from string import split, replace, atof
from sys import argv, exit, stdout
from time import *
import numpy
import numpy.linalg as linalg
import dpTensor
import fec
import dpUtils

class stiffness():

    def __init__(self):
        self.stiffnessMatrix = numpy.zeros((6, 6), numpy.float)
        self.coeffOrder = numpy.zeros((6, 1))
        self.symmFactor = 0
        self.orthoDirections = None
        return

    def setStiffness(self, stiffness):
        self.stiffnessMatrix = stiffness.getStiffnessMatrix()
        self.coeffOrder = stiffness.getCoeffOrder()
        self.symmFactor = stiffness.getSymmFactor()

    def setStiffnessParts(self, matrix, symmFkt, coeffOrder):
        self.stiffnessMatrix = matrix
        self.coeffOrder = coeffOrder
        self.symmFactor = symmFkt

    def readFromFile(self, infile, startLineNo=1):
        inStream = open(infile, 'r')
        lines = inStream.readlines()
        rflag = -1
        coeffOrder = numpy.zeros((6, 1))
        symFact = 0
        lineNo = 0
        for line in lines:
            lineNo += 1
            if lineNo >= startLineNo:
                if line.find('**') == 0:
                    continue
                if line.find('**') > 0:
                    line = line[0:line.find('**') - 1]
                line = replace(line, '\n', '')
                Pos = line.upper().find('*STIFFNESS')
                if Pos == 0:
                    rflag = 0
                if rflag == 1:
                    compList = split(line)
                    self.coeffOrder = numpy.array([int(compList[0]), int(compList[1]),
                     int(compList[2]), int(compList[3]),
                     int(compList[4]), int(compList[5])])
                if rflag == 2:
                    self.symmFactor = float(line)
                if rflag == 3:
                    self.readStiffLine(line, 1)
                if rflag == 4:
                    self.readStiffLine(line, 2)
                if rflag == 5:
                    self.readStiffLine(line, 3)
                if rflag == 6:
                    self.readStiffLine(line, 4)
                if rflag == 7:
                    self.readStiffLine(line, 5)
                if rflag == 8:
                    self.readStiffLine(line, 6)
                if rflag > -1:
                    rflag = rflag + 1
                if rflag > 8:
                    break

        if rflag == -1:
            print '\n ** ERROR ** :  No Stiffness Matrix found in file: ', infile
            exit(1)
        stdout.write(' ... read data from file : %s\n' % infile)
        stdout.flush()

    def setOrthoStiffness(self, E1, E2, E3, G23, G13, G12, nu23, nu13, nu12, stype='eng'):
        self.coeffOrder = numpy.array([11, 22, 33, 23, 13, 12])
        S = numpy.zeros((6, 6), numpy.float)
        S[(0, 0)] = 1.0 / E1
        S[(1, 1)] = 1.0 / E2
        S[(2, 2)] = 1.0 / E3
        if stype == 'eng':
            self.symmFactor = 2.0
            S[(3, 3)] = 1.0 / G23
            S[(4, 4)] = 1.0 / G13
            S[(5, 5)] = 1.0 / G12
        elif stype == 'voight':
            self.symmFactor = numpy.sqrt(2.0)
            S[(3, 3)] = 1.0 / 2.0 / G23
            S[(4, 4)] = 1.0 / 2.0 / G13
            S[(5, 5)] = 1.0 / 2.0 / G12
        else:
            stdout.write('\n **ERROR** : dpStiffness.setOrthoStiffness: Type "%s" not supported!' % stype)
            stdout.flush()
        S[(0, 1)] = -nu12 / E1
        S[(1, 0)] = -nu12 / E1
        S[(0, 2)] = -nu13 / E1
        S[(2, 0)] = -nu13 / E1
        S[(1, 2)] = -nu23 / E2
        S[(2, 1)] = -nu23 / E2
        self.stiffnessMatrix = linalg.inv(S)

    def setIsoStiffness(self, E, nu, stype='eng'):
        self.coeffOrder = numpy.array([11, 22, 33, 23, 13, 12])
        S = numpy.zeros((6, 6), numpy.float)
        G = E / 2.0 / (1.0 + nu)
        S[(0, 0)] = 1.0 / E
        S[(1, 1)] = 1.0 / E
        S[(2, 2)] = 1.0 / E
        if stype == 'eng':
            self.symmFactor = 2.0
            S[(3, 3)] = 1.0 / G
            S[(4, 4)] = 1.0 / G
            S[(5, 5)] = 1.0 / G
        elif stype == 'voight':
            self.symmFactor = numpy.sqrt(2.0)
            S[(3, 3)] = 1.0 / 2.0 / G
            S[(4, 4)] = 1.0 / 2.0 / G
            S[(5, 5)] = 1.0 / 2.0 / G
        else:
            stdout.write('\n **ERROR** : dpStiffness.setIsoStiffness: Type "%s" not supported!' % stype)
            stdout.flush()
        S[(0, 1)] = -nu / E
        S[(1, 0)] = -nu / E
        S[(0, 2)] = -nu / E
        S[(2, 0)] = -nu / E
        S[(1, 2)] = -nu / E
        S[(2, 1)] = -nu / E
        self.stiffnessMatrix = linalg.inv(S)

    def setOrthoStiffnessModel(self, E0, nu0, G0, k, l, rho, m1, m2, m3, stype='eng'):
        self.coeffOrder = numpy.array([11, 22, 33, 23, 13, 12])
        S = numpy.zeros((6, 6), numpy.float)
        S[(0, 0)] = 1.0 / E0 / m1 ** (2 * l) / rho ** k
        S[(1, 1)] = 1.0 / E0 / m2 ** (2 * l) / rho ** k
        S[(2, 2)] = 1.0 / E0 / m3 ** (2 * l) / rho ** k
        if stype == 'eng':
            self.symmFactor = 2.0
            S[(3, 3)] = 1.0 / G0 / (m2 * m3) ** l / rho ** k
            S[(4, 4)] = 1.0 / G0 / (m1 * m3) ** l / rho ** k
            S[(5, 5)] = 1.0 / G0 / (m1 * m2) ** l / rho ** k
        elif stype == 'voight':
            self.symmFactor = numpy.sqrt(2.0)
            S[(3, 3)] = 1.0 / 2.0 / G0 / (m2 * m3) ** l / rho ** k
            S[(4, 4)] = 1.0 / 2.0 / G0 / (m1 * m3) ** l / rho ** k
            S[(5, 5)] = 1.0 / 2.0 / G0 / (m1 * m2) ** l / rho ** k
        else:
            stdout.write('\n **ERROR** : dpStiffness.setOrthoStiffnessModel: Type "%s" not supported!' % stype)
            stdout.flush()
        S[(0, 1)] = -nu0 / E0 / (m1 * m2) ** l / rho ** k
        S[(1, 0)] = -nu0 / E0 / (m1 * m2) ** l / rho ** k
        S[(0, 2)] = -nu0 / E0 / (m1 * m3) ** l / rho ** k
        S[(2, 0)] = -nu0 / E0 / (m1 * m3) ** l / rho ** k
        S[(1, 2)] = -nu0 / E0 / (m2 * m3) ** l / rho ** k
        S[(2, 1)] = -nu0 / E0 / (m2 * m3) ** l / rho ** k
        self.stiffnessMatrix = linalg.inv(S)

    def setIsoStiffnessModel(self, E0, nu0, k, rho, stype='eng'):
        self.coeffOrder = numpy.array([11, 22, 33, 23, 13, 12])
        S = numpy.zeros((6, 6), numpy.float)
        S[(0, 0)] = 1.0 / E0 / rho ** k
        S[(1, 1)] = 1.0 / E0 / rho ** k
        S[(2, 2)] = 1.0 / E0 / rho ** k
        G0 = E0 / 2.0 / (1.0 + nu0)
        if stype == 'eng':
            self.symmFactor = 2.0
            S[(3, 3)] = 1.0 / G0 / rho ** k
            S[(4, 4)] = 1.0 / G0 / rho ** k
            S[(5, 5)] = 1.0 / G0 / rho ** k
        elif stype == 'voight':
            self.symmFactor = numpy.sqrt(2.0)
            S[(3, 3)] = 1.0 / 2.0 / G0 / rho ** k
            S[(4, 4)] = 1.0 / 2.0 / G0 / rho ** k
            S[(5, 5)] = 1.0 / 2.0 / G0 / rho ** k
        else:
            stdout.write('\n **ERROR** : dpStiffness.setIsoStiffnessModel: Type "%s" not supported!' % stype)
            stdout.flush()
        S[(0, 1)] = -nu0 / E0 / rho ** k
        S[(1, 0)] = -nu0 / E0 / rho ** k
        S[(0, 2)] = -nu0 / E0 / rho ** k
        S[(2, 0)] = -nu0 / E0 / rho ** k
        S[(1, 2)] = -nu0 / E0 / rho ** k
        S[(2, 1)] = -nu0 / E0 / rho ** k
        self.stiffnessMatrix = linalg.inv(S)

    def getStiffnessMatrix(self):
        return self.stiffnessMatrix

    def getComplianceMatrix(self):
        return linalg.inv(self.stiffnessMatrix)

    def getCoeffOrder(self):
        return self.coeffOrder

    def getSymmFactor(self):
        return self.symmFactor

    def getOrthoStiffness(self):
        sm = numpy.zeros((6, 6), numpy.float)
        for i in range(6):
            for j in range(6):
                sm[i, j] = 0.5 * (self.stiffnessMatrix[i, j] + self.stiffnessMatrix[j, i])

        sm[(0, 3)] = 0.0
        sm[(0, 4)] = 0.0
        sm[(0, 5)] = 0.0
        sm[(1, 3)] = 0.0
        sm[(1, 4)] = 0.0
        sm[(1, 5)] = 0.0
        sm[(2, 3)] = 0.0
        sm[(2, 4)] = 0.0
        sm[(2, 5)] = 0.0
        sm[(3, 0)] = 0.0
        sm[(3, 1)] = 0.0
        sm[(3, 2)] = 0.0
        sm[(4, 0)] = 0.0
        sm[(4, 1)] = 0.0
        sm[(4, 2)] = 0.0
        sm[(5, 0)] = 0.0
        sm[(5, 1)] = 0.0
        sm[(5, 2)] = 0.0
        sm[(3, 4)] = 0.0
        sm[(3, 5)] = 0.0
        sm[(4, 5)] = 0.0
        sm[(4, 3)] = 0.0
        sm[(5, 3)] = 0.0
        sm[(5, 4)] = 0.0
        return sm

    def getOrthoCompliance(self):
        CompMatrix = self.getInverseStiffnessMatrix()
        sm = numpy.zeros((6, 6), numpy.float)
        for i in range(6):
            for j in range(6):
                sm[i, j] = 0.5 * (CompMatrix[i, j] + CompMatrix[i, j])

        sm[(0, 3)] = 0.0
        sm[(0, 4)] = 0.0
        sm[(0, 5)] = 0.0
        sm[(1, 3)] = 0.0
        sm[(1, 4)] = 0.0
        sm[(1, 5)] = 0.0
        sm[(2, 3)] = 0.0
        sm[(2, 4)] = 0.0
        sm[(2, 5)] = 0.0
        sm[(3, 0)] = 0.0
        sm[(3, 1)] = 0.0
        sm[(3, 2)] = 0.0
        sm[(4, 0)] = 0.0
        sm[(4, 1)] = 0.0
        sm[(4, 2)] = 0.0
        sm[(5, 0)] = 0.0
        sm[(5, 1)] = 0.0
        sm[(5, 2)] = 0.0
        sm[(3, 4)] = 0.0
        sm[(3, 5)] = 0.0
        sm[(4, 5)] = 0.0
        sm[(4, 3)] = 0.0
        sm[(5, 3)] = 0.0
        sm[(5, 4)] = 0.0
        return sm

    def getEngConstants(self):
        CompMatrix = self.getInverseStiffnessMatrix()
        engConst = {}
        engConst['E_' + repr(self.coeffOrder[0])[0]] = 1.0 / CompMatrix[(0, 0)]
        engConst['E_' + repr(self.coeffOrder[1])[0]] = 1.0 / CompMatrix[(1, 1)]
        engConst['E_' + repr(self.coeffOrder[2])[0]] = 1.0 / CompMatrix[(2, 2)]
        complCorr = 2.0 * self.symmFactor / 4.0
        engConst['G_' + repr(self.coeffOrder[3])] = 1.0 / (CompMatrix[(3, 3)] * complCorr)
        engConst['G_' + repr(self.coeffOrder[4])] = 1.0 / (CompMatrix[(4, 4)] * complCorr)
        engConst['G_' + repr(self.coeffOrder[5])] = 1.0 / (CompMatrix[(5, 5)] * complCorr)
        engConst['nu_' + repr(self.coeffOrder[3])] = -engConst['E_' + repr(self.coeffOrder[1])[0]] * CompMatrix[(1,
                                                                                                                 2)]
        engConst['nu_' + repr(self.coeffOrder[4])] = -engConst['E_' + repr(self.coeffOrder[0])[0]] * CompMatrix[(0,
                                                                                                                 2)]
        engConst['nu_' + repr(self.coeffOrder[5])] = -engConst['E_' + repr(self.coeffOrder[0])[0]] * CompMatrix[(0,
                                                                                                                 1)]
        if 'nu_13' not in engConst:
            engConst['nu_13'] = engConst['nu_31'] / engConst['E_3'] * engConst['E_1']
        if 'nu_32' not in engConst:
            engConst['nu_32'] = engConst['nu_23'] / engConst['E_2'] * engConst['E_3']
        if 'nu_21' not in engConst:
            engConst['nu_21'] = engConst['nu_12'] / engConst['E_1'] * engConst['E_2']
        if 'nu_31' not in engConst:
            engConst['nu_31'] = engConst['nu_13'] / engConst['E_1'] * engConst['E_3']
        if 'nu_23' not in engConst:
            engConst['nu_23'] = engConst['nu_32'] / engConst['E_3'] * engConst['E_2']
        if 'nu_12' not in engConst:
            engConst['nu_12'] = engConst['nu_21'] / engConst['E_2'] * engConst['E_1']
        KMODI = 0.0
        for i in range(3):
            for j in range(3):
                KMODI += CompMatrix[i, j]

        engConst['K'] = 1.0 / KMODI
        return engConst

    def getInverseStiffnessMatrix(self):
        return linalg.inv(self.stiffnessMatrix)

    def getOrthoStiffnessCoefficient(self):
        stiffList = []
        stiffList.append(self.stiffnessMatrix[(0, 0)])
        stiffList.append(self.stiffnessMatrix[(1, 1)])
        stiffList.append(self.stiffnessMatrix[(2, 2)])
        stiffList.append(self.stiffnessMatrix[(0, 1)])
        stiffList.append(self.stiffnessMatrix[(0, 2)])
        stiffList.append(self.stiffnessMatrix[(1, 2)])
        stiffList.append(self.stiffnessMatrix[(1, 0)])
        stiffList.append(self.stiffnessMatrix[(2, 0)])
        stiffList.append(self.stiffnessMatrix[(2, 1)])
        stiffList.append(self.stiffnessMatrix[(3, 3)])
        stiffList.append(self.stiffnessMatrix[(4, 4)])
        stiffList.append(self.stiffnessMatrix[(5, 5)])
        return stiffList

    def getOrthoComplianceCoefficient(self):
        CompMatrix = self.getInverseStiffnessMatrix()
        stiffList = []
        stiffList.append(CompMatrix[(0, 0)])
        stiffList.append(CompMatrix[(1, 1)])
        stiffList.append(CompMatrix[(2, 2)])
        stiffList.append(CompMatrix[(0, 1)])
        stiffList.append(CompMatrix[(0, 2)])
        stiffList.append(CompMatrix[(1, 2)])
        stiffList.append(CompMatrix[(1, 0)])
        stiffList.append(CompMatrix[(2, 0)])
        stiffList.append(CompMatrix[(2, 1)])
        stiffList.append(CompMatrix[(3, 3)])
        stiffList.append(CompMatrix[(4, 4)])
        stiffList.append(CompMatrix[(5, 5)])
        return stiffList

    def disturb(self, prec):
        self.stiffnessMatrix[(0, 0)] = self.stiffnessMatrix[(0, 0)] + prec
        self.stiffnessMatrix[(1, 1)] = self.stiffnessMatrix[(1, 1)] - prec
        self.stiffnessMatrix[(4, 4)] = self.stiffnessMatrix[(4, 4)] + prec
        self.stiffnessMatrix[(5, 5)] = self.stiffnessMatrix[(5, 5)] - prec

    def changeSymmFactor(self, symmFactor):
        oldSymmFactor = self.symmFactor
        self.symmFactor = symmFactor
        corr = oldSymmFactor
        corr2 = symmFactor
        T2A = self.getInverseStiffnessMatrix()
        T4A = dpTensor.Matrix6_to_Tensor4(T2A, self.coeffOrder, corr)
        T2A = dpTensor.Tensor4_to_Matrix6(T4A, self.coeffOrder, corr2)
        self.stiffnessMatrix = linalg.inv(T2A)

    def changeCoeffOrder(self, coeffOrder):
        oldCoeffOrder = self.coeffOrder
        self.coeffOrder = coeffOrder
        corr = self.symmFactor
        T2A = self.getInverseStiffnessMatrix()
        T4A = dpTensor.Matrix6_to_Tensor4(T2A, oldCoeffOrder, corr)
        T2A = dpTensor.Tensor4_to_Matrix6(T4A, coeffOrder, corr)
        self.stiffnessMatrix = linalg.inv(T2A)

    def symmetrize(self):
        snew = numpy.zeros((6, 6), numpy.float)
        for i in range(6):
            for j in range(6):
                snew[i, j] = (self.stiffnessMatrix[i, j] + self.stiffnessMatrix[j, i]) / 2.0

        for i in range(6):
            for j in range(6):
                self.stiffnessMatrix[i, j] = snew[i, j]

        del snew

    def computeEngConstDIST_slow(self, triaArray, basis=None):
        stdout.write(' ... compute distribution slow\n')
        stdout.flush()
        time1 = clock()
        if basis == None:
            e1a = numpy.array([1.0, 0.0, 0.0])
            e2a = numpy.array([0.0, 1.0, 0.0])
            e3a = numpy.array([0.0, 0.0, 1.0])
            Ba = numpy.array([e1a, e2a, e3a])
        else:
            Ba = basis
        mStiff_I = self.getComplianceMatrix()
        comp = self.getCoeffOrder()
        symFkt = self.getSymmFactor()
        Cold = dpTensor.Matrix6_to_Tensor4(mStiff_I, comp, symFkt)
        engConstDIST = {}
        checkdouble = {}
        i = 0
        sum = float(len(triaArray))
        dpUtils.progressStart('     -> Progress            : ')
        stdout.flush()
        for tria in triaArray:
            i += 1
            progress = float(i) / float(sum) * 10.0
            dpUtils.progressNext(progress)
            stdout.flush()
            for point in tria:
                if not checkdouble.has_key(point):
                    checkdouble[point] = True
                    e1n = numpy.array(point)
                    C1111 = 0.0
                    for r in range(3):
                        for s in range(3):
                            for t in range(3):
                                for u in range(3):
                                    C1111 = C1111 + dpTensor.SimpleContraction(Ba[r], e1n) * dpTensor.SimpleContraction(Ba[s], e1n) * dpTensor.SimpleContraction(Ba[t], e1n) * dpTensor.SimpleContraction(Ba[u], e1n) * Cold[r, s, t, u]

                    engConstDIST[point] = 1.0 / C1111

        dpUtils.progressEnd()
        stdout.flush()
        time2 = clock()
        stdout.write('     -> Compute in          :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        return engConstDIST

    def computeEngConstDIST(self, triaArray, basis=None, typ='E'):
        stdout.write(' ... compute distribution\n')
        stdout.flush()
        time1 = clock()
        mStiff_I = self.getComplianceMatrix()
        comp = self.getCoeffOrder()
        symFkt = self.getSymmFactor()
        C = dpTensor.Matrix6_to_Tensor4(mStiff_I, comp, symFkt)
        if basis != None:
            Ba = basis
            e1 = numpy.array([1.0, 0.0, 0.0])
            e2 = numpy.array([0.0, 1.0, 0.0])
            e3 = numpy.array([0.0, 0.0, 1.0])
            Bn = numpy.array([e1, e2, e3])
            C = dpTensor.TensorTransformBas(C, Ba, Bn)
        engConstDIST = {}
        checkdouble = {}
        i = 0
        sum = float(len(triaArray))
        I = dpTensor.UnitMatrix(3)
        dpUtils.progressStart('     -> Progress            : ')
        stdout.flush()
        oldprogress = 0
        for tria in triaArray:
            i += 1
            progress = int(float(i) / float(sum) * 10.0)
            if progress > oldprogress:
                dpUtils.progressNext(progress)
                stdout.flush()
            oldprogress = progress
            for point in tria:
                if not checkdouble.has_key(point):
                    checkdouble[point] = True
                    e1n = numpy.array(point)
                    N = dpTensor.DyadicProduct(e1n, e1n)
                    if typ == 'E':
                        CN = dpTensor.DoubleContraction(C, N)
                        NCN = dpTensor.DoubleContraction(N, CN)
                        engConstDIST[point] = 1.0 / NCN
                    elif typ == 'K':
                        CN = dpTensor.DoubleContraction(C, N)
                        ICN = dpTensor.DoubleContraction(I, CN)
                        if ICN < 1e-06:
                            stdout.write('\n **Warning**: dpStiffness.computeEngConstDIST() Bulk Modulus nearly infinite!')
                            stdout.flush()
                            ICN = 1e-06
                        engConstDIST[point] = 1.0 / ICN

        dpUtils.progressEnd()
        time2 = clock()
        stdout.write('     -> Compute in          :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        stdout.flush()
        return engConstDIST

    def writeMatrix(self, OS, name='STIFFNESS', value=None):
        if value == None:
            value = self.stiffnessMatrix
        OS.write('\n*%s Matrix\n' % name)
        for comp in self.coeffOrder:
            OS.write('%2i ' % comp)

        OS.write('     ** stress/strain component order\n')
        OS.write('%12.6g           ** symmetrization factor (2..engineering)' % self.symmFactor)
        OS.write('\n******************************************************************************')
        for row in range(6):
            OS.write('\n')
            for col in range(6):
                OS.write('%12.6g ' % value[col, row])

        OS.write('\n\n')
        return

    def writeAbaqusMaterial(self, OS, name, matType, density=None):
        OS.write('*MATERIAL, NAME=%s\n' % name)
        OS.write('*ELASTIC, TYPE=%s\n' % matType)
        corr = self.symmFactor
        corr2 = 2.0
        CompMatrix = self.getInverseStiffnessMatrix()
        T4A = dpTensor.Matrix6_to_Tensor4(CompMatrix, self.coeffOrder, corr)
        CompMatrix = dpTensor.Tensor4_to_Matrix6(T4A, self.coeffOrder, corr2)
        if matType == 'ENGINEERING CONSTANTS':
            engConst = {}
            engConst['E' + repr(self.coeffOrder[0])] = 1.0 / CompMatrix[(0, 0)]
            engConst['E' + repr(self.coeffOrder[1])] = 1.0 / CompMatrix[(1, 1)]
            engConst['E' + repr(self.coeffOrder[2])] = 1.0 / CompMatrix[(2, 2)]
            engConst['G' + repr(self.coeffOrder[3])] = 1.0 / CompMatrix[(3, 3)]
            engConst['G' + repr(self.coeffOrder[4])] = 1.0 / CompMatrix[(4, 4)]
            engConst['G' + repr(self.coeffOrder[5])] = 1.0 / CompMatrix[(5, 5)]
            engConst['nu' + repr(self.coeffOrder[3])] = -engConst['E22'] * CompMatrix[(1,
                                                                                       2)]
            engConst['nu' + repr(self.coeffOrder[4])] = -engConst['E11'] * CompMatrix[(0,
                                                                                       2)]
            engConst['nu' + repr(self.coeffOrder[5])] = -engConst['E11'] * CompMatrix[(0,
                                                                                       1)]
            OS.write('%g, %g, %g, %g, %g, %g, %g, %g\n' % (
             engConst['E11'], engConst['E22'], engConst['E33'],
             engConst['nu12'], engConst['nu13'], engConst['nu23'],
             engConst['G12'], engConst['G13']))
            OS.write('%g,\n' % engConst['G23'])
        else:
            stdout.write('\n **ERROR** : dpStiffness.writeAbaqusMaterial: Type "%s" not supported!' % matType)
            stdout.flush()
            exit(1)
        if density != None:
            OS.write('*density\n')
            OS.write('%g\n' % float(density))
        return

    def writeValue(self, OS, text, value):
        OS.write('%s = %s\n' % (text, repr(value)))

    def writeAbaqusOrientation(self, OS, name, a, b):
        OS.write('*ORIENTATION, NAME=%s\n' % name)
        OS.write('%g, %g, %g, %g, %g, %g\n' % (a[0], a[1], a[2], b[0], b[1], b[2]))
        OS.write('1, 0.\n')

    def writeAbaqusSolidSection(self, OS, elset, material, orientation):
        if orientation == None:
            OS.write('*solid section, elset=%s,  material=%s\n' % (elset, material))
        else:
            OS.write('*solid section, elset=%s,  material=%s,  orientation=%s \n' % (elset, material, orientation))
        return

    def writeEngConst(self, OS):
        CompMatrix = self.getInverseStiffnessMatrix()
        OS.write('\n*ENGINEERING Constants                       ')
        OS.write('\n******************************************************************************\n')
        complCorr = 2.0 * self.symmFactor / 4.0
        corr = self.symmFactor
        corr2 = 2.0
        CompMatrix = self.getInverseStiffnessMatrix()
        T4A = dpTensor.Matrix6_to_Tensor4(CompMatrix, self.coeffOrder, corr)
        CompMatrix = dpTensor.Tensor4_to_Matrix6(T4A, self.coeffOrder, corr2)
        engConst = {}
        engConst['E' + repr(self.coeffOrder[0])] = 1.0 / CompMatrix[(0, 0)]
        engConst['E' + repr(self.coeffOrder[1])] = 1.0 / CompMatrix[(1, 1)]
        engConst['E' + repr(self.coeffOrder[2])] = 1.0 / CompMatrix[(2, 2)]
        engConst['G' + repr(self.coeffOrder[3])] = 1.0 / CompMatrix[(3, 3)]
        engConst['G' + repr(self.coeffOrder[4])] = 1.0 / CompMatrix[(4, 4)]
        engConst['G' + repr(self.coeffOrder[5])] = 1.0 / CompMatrix[(5, 5)]
        engConst['nu' + repr(self.coeffOrder[3])] = -engConst['E22'] * CompMatrix[(1,
                                                                                   2)]
        engConst['nu' + repr(self.coeffOrder[4])] = -engConst['E11'] * CompMatrix[(0,
                                                                                   2)]
        engConst['nu' + repr(self.coeffOrder[5])] = -engConst['E11'] * CompMatrix[(0,
                                                                                   1)]
        OS.write('%5s = %12.5g       %5s = %12.5g       %5s = %12.5g\n' % (
         'E_1', engConst['E11'], 'G_23', engConst['G23'], 'nu_23', engConst['nu23']))
        OS.write('%5s = %12.5g       %5s = %12.5g       %5s = %12.5g\n' % (
         'E_2', engConst['E22'], 'G_13', engConst['G13'], 'nu_13', engConst['nu13']))
        OS.write('%5s = %12.5g       %5s = %12.5g       %5s = %12.5g\n' % (
         'E_3', engConst['E33'], 'G_12', engConst['G12'], 'nu_12', engConst['nu12']))
        OS.write('\n')

    def writeParaviewEngConstDist(self, name, approxDIST, triaArray):
        fecModel = fec.fec()
        mul = 1.0
        approxNodes, approxElems, results = myMiaModel.setupSphereFem(triaArray, approxDIST)
        sResults = [results]
        fecModel.write(name + '_EngConst' + '.geo', 'Stiffness Tensor', approxNodes, None, approxElems, None, NscaResults=sResults)
        return

    def rotateToMaxMidMin(self):
        CompMatrix = self.getInverseStiffnessMatrix()
        EA = 1.0 / CompMatrix[(0, 0)]
        EB = 1.0 / CompMatrix[(1, 1)]
        EC = 1.0 / CompMatrix[(2, 2)]
        e1 = numpy.array([0.0, 0.0, 0.0])
        e2 = numpy.array([0.0, 0.0, 0.0])
        e3 = numpy.array([0.0, 0.0, 0.0])
        if EA > EB and EA > EC:
            e1[0] = 1.0
            if EB > EC:
                e2[1] = 1.0
                e3[2] = 1.0
            else:
                e2[2] = 1.0
                e3[1] = 1.0
        if EB > EA and EB > EC:
            e2[0] = 1.0
            if EA > EC:
                e1[1] = 1.0
                e3[2] = 1.0
            else:
                e1[2] = 1.0
                e3[1] = 1.0
        if EC > EA and EC > EB:
            e3[0] = 1.0
            if EA > EB:
                e1[1] = 1.0
                e2[2] = 1.0
            else:
                e1[2] = 1.0
                e2[1] = 1.0
        Ba = numpy.array([e1, e2, e3])
        Bn = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        comp = self.getCoeffOrder()
        symFkt = self.getSymmFactor()
        C = dpTensor.Matrix6_to_Tensor4(self.stiffnessMatrix, comp, symFkt)
        CR = dpTensor.TensorTransformBas(C, Ba, Bn)
        self.stiffnessMatrix = dpTensor.Tensor4_to_Matrix6(CR, comp, symFkt)

    def rotateFromEvec(self, Evec):
        Ba = numpy.array([Evec[0], Evec[1], Evec[2]])
        Bn = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        comp = self.getCoeffOrder()
        symFkt = self.getSymmFactor()
        C2 = self.getInverseStiffnessMatrix()
        C4 = dpTensor.Matrix6_to_Tensor4(C2, comp, symFkt)
        C4R = dpTensor.TensorTransformBas(C4, Ba, Bn)
        C2R = dpTensor.Tensor4_to_Matrix6(C4R, comp, symFkt)
        self.stiffnessMatrix = linalg.inv(C2R)

    def rotateAngle(self, order, angle):
        RO = numpy.array([int(order[0]), int(order[1]), int(order[2])])
        RA = numpy.array([float(angle[0]), float(angle[1]), float(angle[2])])
        comp = self.getCoeffOrder()
        symFkt = self.getSymmFactor()
        C2 = self.getInverseStiffnessMatrix()
        C4 = dpTensor.Matrix6_to_Tensor4(C2, comp, symFkt)
        C4R = dpTensor.TensorTransformRot(C4, RO, RA)
        C2R = dpTensor.Tensor4_to_Matrix6(C4R, comp, symFkt)
        self.stiffnessMatrix = linalg.inv(C2R)

    def computeOrthoError(self, matrix='C'):
        comp = self.getCoeffOrder()
        symFkt = self.getSymmFactor()
        T2A = numpy.zeros((6, 6), numpy.float)
        T2O = numpy.zeros((6, 6), numpy.float)
        if matrix == 'S' or matrix == 's':
            corr = 2.0 / self.getSymmFactor()
            T2A = self.stiffnessMatrix
            T2O = self.getOrthoStiffness()
            T4A = dpTensor.Matrix6_to_Tensor4(T2A, comp, corr)
            T4O = dpTensor.Matrix6_to_Tensor4(T2O, comp, corr)
        elif matrix == 'C' or matrix == 'c':
            corr = self.getSymmFactor()
            T2A = self.getInverseStiffnessMatrix()
            T2O = self.getOrthoCompliance()
            T4A = dpTensor.Matrix6_to_Tensor4(T2A, comp, corr)
            T4O = dpTensor.Matrix6_to_Tensor4(T2O, comp, corr)
        diffT4 = T4A - T4O
        zaehler = dpTensor.DoubleContraction(diffT4, diffT4)
        nenner = dpTensor.DoubleContraction(T4A, T4A)
        return numpy.sqrt(zaehler / nenner)

    def computeLogOrthoError(self, matrix='C', orthoMatrix=None):
        prec = 1e-05
        if not (self.getSymmFactor() > 2.0 - prec and self.getSymmFactor() < 2.0 + prec):
            if not (self.getSymmFactor() > numpy.sqrt(2.0) - prec and self.getSymmFactor() < numpy.sqrt(2.0) + prec):
                stdout.write('\n **ERROR** : dpStiffness.computeLogOrthoError: Symmetrization factor not set!')
                stdout.flush()
                exit(1)
        comp = self.getCoeffOrder()
        T2A = numpy.zeros((6, 6), numpy.float)
        T2O = numpy.zeros((6, 6), numpy.float)
        if matrix == 'S' or matrix == 's':
            corr = 2.0 / self.getSymmFactor()
            corr2 = 2.0 / numpy.sqrt(2.0)
            T2A = self.stiffnessMatrix
            if orthoMatrix != None:
                T2O = orthoMatrix.getOrthoStiffness()
            else:
                T2O = self.getOrthoStiffness()
            T4A = dpTensor.Matrix6_to_Tensor4(T2A, comp, corr)
            T2A = dpTensor.Tensor4_to_Matrix6(T4A, comp, corr2)
            T4O = dpTensor.Matrix6_to_Tensor4(T2O, comp, corr)
            T2O = dpTensor.Tensor4_to_Matrix6(T4O, comp, corr2)
        elif matrix == 'C' or matrix == 'c':
            corr = self.getSymmFactor()
            corr2 = numpy.sqrt(2.0)
            T2A = self.getInverseStiffnessMatrix()
            if orthoMatrix != None:
                T2O = linalg.inv(orthoMatrix.getOrthoStiffness())
            else:
                T2O = linalg.inv(self.getOrthoStiffness())
            T4A = dpTensor.Matrix6_to_Tensor4(T2A, comp, corr)
            T2A = dpTensor.Tensor4_to_Matrix6(T4A, comp, corr2)
            T4O = dpTensor.Matrix6_to_Tensor4(T2O, comp, corr)
            T2O = dpTensor.Tensor4_to_Matrix6(T4O, comp, corr2)
        evalue, evect = linalg.eig(T2A)
        EVEC = evect
        EVAL = numpy.zeros((6, 6), numpy.float)
        for i in range(6):
            EVAL[i, i] = numpy.log(evalue[i])

        T2A = numpy.dot(numpy.dot(EVEC, EVAL), numpy.transpose(EVEC))
        evalue, evect = linalg.eig(T2O)
        EVEC = evect
        EVAL = numpy.zeros((6, 6), numpy.float)
        for i in range(6):
            EVAL[i, i] = numpy.log(evalue[i])

        T2O = numpy.dot(numpy.dot(EVEC, EVAL), numpy.transpose(EVEC))
        diffT2 = T2A - T2O
        zaehler = 0.0
        nenner = 0.0
        for i in range(6):
            for j in range(6):
                zaehler = zaehler + diffT2[i, j] * diffT2[i, j]
                nenner = nenner + T2A[i, j] * T2A[i, j]

        return numpy.sqrt(zaehler / nenner)

    def computeClosetIsoLogEuclidean(self):
        prec = 1e-05
        if not (self.getSymmFactor() > 2.0 - prec and self.getSymmFactor() < 2.0 + prec):
            if not (self.getSymmFactor() > numpy.sqrt(2.0) - prec and self.getSymmFactor() < numpy.sqrt(2.0) + prec):
                stdout.write('\n **ERROR** : dpStiffness.computeClosetIsoLogEuclidean: Symmetrization factor not set!')
                stdout.flush()
                exit(1)
        comp = self.getCoeffOrder()
        T2A = numpy.zeros((6, 6), numpy.float)
        corr = 2.0 / self.getSymmFactor()
        corr2 = 2.0 / numpy.sqrt(2.0)
        T2A = self.stiffnessMatrix
        T4A = dpTensor.Matrix6_to_Tensor4(T2A, comp, corr)
        T2A = dpTensor.Tensor4_to_Matrix6(T4A, comp, corr2)
        evalue, EVEC = linalg.eig(T2A)
        EVAL = numpy.zeros((6, 6), numpy.float)
        for i in range(6):
            EVAL[i, i] = numpy.log(evalue[i])

        logC = numpy.zeros((6, 6), numpy.float)
        logC = numpy.dot(numpy.dot(EVEC, EVAL), numpy.transpose(EVEC))
        aux = 1.0 / numpy.sqrt(3.0)
        uvec = numpy.array([[aux], [aux], [aux], [0.0], [0.0], [0.0]])
        Jhat = numpy.dot(uvec, numpy.transpose(uvec))
        Ihat = numpy.eye(6)
        Khat = Ihat - Jhat
        K = 1.0 / 3.0 * numpy.exp(numpy.trace(numpy.dot(Jhat, logC)))
        G = 1.0 / 2.0 * numpy.exp(1.0 / 5.0 * numpy.trace(numpy.dot(Khat, logC)))
        E = 9.0 * K * G / (3.0 * K + G)
        nu = (3.0 * K - 2.0 * G) / (6.0 * K + 2.0 * G)
        return (
         K, G, E, nu)

    def readStiffLine(self, line, no):
        Str = split(line)
        for i in range(6):
            self.stiffnessMatrix[no - 1, i] = float(Str[i])


if __name__ == '__main__':
    inName = None
    addName = None
    argList = argv
    argc = len(argList)
    i = 0
    while i < argc:
        if argList[i][:3] == '-in':
            i += 1
            inName = argList[i]
        elif argList[i][:4] == '-add':
            i += 1
            addName = argList[i]
        i += 1

    if not inName:
        print __doc__
        print ' **ERROR** input file name is not provided'
        exit(1)
    myStiffness = stiffness()
    myStiffness.readFromFile(inName)
    myStiffness.writeMatrix(stdout, 'STIFFNESS')
# okay decompiling dpStiffness.pyc
