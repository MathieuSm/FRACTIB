# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/misc/dpTensor.py
# Compiled at: 2017-02-10 11:21:53
"""
#########################################################################
File   : dpPTensor.py
Author : D.H.Pahr
Date   : 17.12.2004
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage: python dpTensor.py   

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requirements:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Description: 
  This program is a tensor algebra library.  
  
#########################################################################
"""
from string import split, replace, atof
from sys import argv, exit
import numpy
import numpy.linalg as linalg
import unittest

def CrossProduct(a, b):
    c1 = a[1] * b[2] - a[2] * b[1]
    c2 = a[2] * b[0] - a[0] * b[2]
    c3 = a[0] * b[1] - a[1] * b[0]
    c = numpy.array([c1, c2, c3])
    return c


def Length(a):
    c = numpy.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    return c


def UnitVector(a):
    l = Length(a)
    c = a / l
    return c


def checkInverse(MI, M, prec=1e-06):
    I = numpy.dot(MI, M)
    dim = I.shape
    ok = True
    for i in range(dim[0]):
        for j in range(dim[1]):
            if i == j:
                if abs(I[i, i] - 1.0) > prec:
                    ok = False
            elif abs(I[i, j]) > prec:
                ok = False

    return ok


def DyadicProduct(A, B):
    CheckShape(A)
    CheckShape(B)
    CheckPosition(A, B)
    type = 10 * len(A.shape) + len(B.shape)
    C = numpy.array([])
    if type == 11:
        C = numpy.zeros((3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                C[i, j] = A[i] * B[j]

    elif type == 21:
        C = numpy.zeros((3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    C[i, j, k] = A[i, j] * B[k]

    elif type == 22:
        C = numpy.zeros((3, 3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        C[i, j, k, l] = A[i, j] * B[k, l]

    else:
        SupportError('DyadicProduct')
    return C


def SymmetricProduct(A, B):
    CheckShape(A)
    CheckShape(B)
    CheckPosition(A, B)
    type = 10 * len(A.shape) + len(B.shape)
    C = numpy.array([])
    if type == 22:
        C = numpy.zeros((3, 3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        C[i, j, k, l] = 0.5 * (A[i, k] * B[j, l] + A[i, l] * B[j, k])

    else:
        SupportError('SymmetricProduct')
    return C


def SimpleContraction(A, B):
    CheckShape(A)
    CheckShape(B)
    CheckPosition(A, B)
    type = 10 * len(A.shape) + len(B.shape)
    C = numpy.array([])
    if type == 11:
        C = numpy.zeros((1, ), numpy.float)
        for m in range(3):
            C[(0, )] = C[(0, )] + A[m] * B[m]

    elif type == 21:
        C = numpy.zeros((3, ), numpy.float)
        for i in range(3):
            for m in range(3):
                C[i] = C[i] + A[i, m] * B[m]

    elif type == 22:
        C = numpy.zeros((3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    C[i, j] = C[i, j] + A[i, m] * B[m, j]

    else:
        SupportError('SimpleContraction')
    if C.shape[0] == 1:
        return C[0]
    else:
        return C


def DoubleContraction(A, B):
    CheckShape(A)
    CheckShape(B)
    CheckPosition(A, B)
    type = 10 * len(A.shape) + len(B.shape)
    C = numpy.array([])
    if type == 22:
        C = numpy.zeros((1, ), numpy.float)
        for m in range(3):
            for n in range(3):
                C[(0, )] = C[(0, )] + A[m, n] * B[m, n]

    elif type == 42:
        C = numpy.zeros((3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    for n in range(3):
                        C[i, j] = C[i, j] + A[i, j, m, n] * B[m, n]

    elif type == 44:
        C = numpy.zeros((1, ), numpy.float)
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    for n in range(3):
                        C[(0, )] = C[(0, )] + A[i, j, m, n] * B[i, j, m, n]

    else:
        SupportError('DoubleConstraction')
    if C.shape[0] == 1:
        return C[0]
    else:
        return C


def TensorTransformBas(A, BO, BN):
    CheckShape(A)
    CheckShape(BO)
    CheckShape(BN)
    type = len(A.shape)
    C = numpy.array([])
    if type == 1:
        C = numpy.zeros((3, ), numpy.float)
        for i in range(3):
            for r in range(3):
                C[i] = C[i] + numpy.dot(BO[r], BN[i]) * A[r]

    elif type == 2:
        C = numpy.zeros((3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for r in range(3):
                    for s in range(3):
                        C[i, j] = C[i, j] + numpy.dot(BO[r], BN[i]) * numpy.dot(BO[s], BN[j]) * A[r, s]

    elif type == 3:
        C = numpy.zeros((3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for r in range(3):
                        for s in range(3):
                            for t in range(3):
                                C[i, j, k] = C[i, j, k] + numpy.dot(BO[r], BN[i]) * numpy.dot(BO[s], BN[j]) * numpy.dot(BO[t], BN[k]) * A[r, s, t]

    elif type == 4:
        C = numpy.zeros((3, 3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for r in range(3):
                            for s in range(3):
                                for t in range(3):
                                    for u in range(3):
                                        C[i, j, k, l] = C[i, j, k, l] + numpy.dot(BO[r], BN[i]) * numpy.dot(BO[s], BN[j]) * numpy.dot(BO[t], BN[k]) * numpy.dot(BO[u], BN[l]) * A[r, s, t, u]

    else:
        SupportError('TensorTransformBas')
    return C


def TensorTransformRot(A, RO, RA):
    CheckShape(A)
    CheckShape(RO)
    CheckShape(RA)
    type = len(A.shape)
    C = numpy.array([])
    pi = 3.141592654
    a1 = RA[0] * pi / 180.0
    a2 = RA[1] * pi / 180.0
    a3 = RA[2] * pi / 180.0
    sa1 = numpy.sin(a1)
    sa2 = numpy.sin(a2)
    sa3 = numpy.sin(a3)
    ca1 = numpy.cos(a1)
    ca2 = numpy.cos(a2)
    ca3 = numpy.cos(a3)
    R = numpy.zeros((3, 3, 3), numpy.float)
    R[0] = numpy.array([[1.0, 0.0, 0.0],
     [
      0.0, ca1, sa1],
     [
      0.0, -sa1, ca1]])
    R[1] = numpy.array([[ca2, 0.0, -sa2],
     [
      0.0, 1.0, 0.0],
     [
      sa2, 0.0, ca2]])
    R[2] = numpy.array([[ca3, sa3, 0.0],
     [
      -sa3, ca3, 0.0],
     [
      0.0, 0.0, 1.0]])
    Rot = numpy.dot(R[RO[2] - 1], numpy.dot(R[RO[1] - 1], R[RO[0] - 1]))
    if type == 1:
        C = numpy.zeros((3, ), numpy.float)
        for i in range(3):
            for r in range(3):
                C[i] = C[i] + Rot[i, r] * A[r]

    elif type == 2:
        C = numpy.zeros((3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for r in range(3):
                    for s in range(3):
                        C[i, j] = C[i, j] + Rot[i, r] * Rot[j, s] * A[r, s]

    elif type == 3:
        C = numpy.zeros((3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for r in range(3):
                        for s in range(3):
                            for t in range(3):
                                C[i, j, k] = C[i, j, k] + Rot[i, r] * Rot[j, s] * Rot[k, t] * A[r, s, t]

    elif type == 4:
        C = numpy.zeros((3, 3, 3, 3), numpy.float)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for r in range(3):
                            for s in range(3):
                                for t in range(3):
                                    for u in range(3):
                                        C[i, j, k, l] = C[i, j, k, l] + Rot[i, r] * Rot[j, s] * Rot[k, t] * Rot[l, u] * A[r, s, t, u]

    else:
        SupportError('TensorTransformRot')
    return C


def Matrix6_to_Tensor4(A, CO, CR):
    if A.shape[0] != 6 and A.shape[1] != 6:
        print '\n** ERROR: Matrix6_to_Tensor4 Shape of A=',
        print A.shape, 'Required (6x6)!'
        exit(1)
    if CO.shape[0] != 6 and CO.shape[1] != 1:
        print '\n** ERROR: Matrix6_to_Tensor4 Shape of Comp=',
        print A.shape, 'Required (6x1)!'
        exit(1)
    CO2 = numpy.zeros((6, 2), numpy.int)
    for i in range(6):
        CO2[i, 0] = int(CO[i] / 10)
        CO2[i, 1] = int(CO[i] % 10)

    C = numpy.array([])
    C = numpy.zeros((3, 3, 3, 3), numpy.float)
    for m in range(6):
        for n in range(6):
            ccr = 1.0
            if m > 2 and n < 3:
                ccr = CR
            if m > 2 and n > 2:
                ccr = CR * CR
            if m < 3 and n > 2:
                ccr = CR
            i = CO2[m, 0] - 1
            j = CO2[m, 1] - 1
            k = CO2[n, 0] - 1
            l = CO2[n, 1] - 1
            C[i, j, k, l] = A[m, n] / ccr
            C[k, l, i, j] = C[i, j, k, l]
            C[j, i, k, l] = C[i, j, k, l]
            C[i, j, l, k] = C[i, j, k, l]
            C[j, i, l, k] = C[i, j, k, l]

    return C


def Tensor4_to_Matrix9(A, CO):
    CheckShape(A)
    if CO.shape[0] != 9 and CO.shape[1] != 1:
        print '\n** ERROR: Tensor4_to_Matrix9 Shape of Comp=',
        print A.shape, 'Required (9x1)!'
        exit(1)
    type = len(A.shape)
    C = numpy.array([])
    CO2 = numpy.zeros((9, 2), numpy.int)
    for i in range(9):
        CO2[i, 0] = int(CO[i] / 10)
        CO2[i, 1] = int(CO[i] % 10)

    if type == 4:
        C = numpy.zeros((9, 9), numpy.float)
        for m in range(9):
            for n in range(9):
                i = CO2[m, 0] - 1
                j = CO2[m, 1] - 1
                k = CO2[n, 0] - 1
                l = CO2[n, 1] - 1
                C[m, n] = A[i, j, k, l]

    else:
        SupportError('Tensor4_to_Matrix9')
    return C


def Tensor4_to_Matrix6(A, CO, CR):
    CheckShape(A)
    if CO.shape[0] != 6 and CO.shape[1] != 1:
        print '\n** ERROR: Tensor4_to_Matrix6 Shape of Comp=',
        print A.shape, 'Required (6x1)!'
        exit(1)
    type = len(A.shape)
    C = numpy.array([])
    CO2 = numpy.zeros((6, 2), numpy.int)
    for i in range(6):
        CO2[i, 0] = int(CO[i] / 10)
        CO2[i, 1] = int(CO[i] % 10)

    if type == 4:
        C = numpy.zeros((6, 6), numpy.float)
        for m in range(6):
            for n in range(6):
                ccr = 1.0
                if m > 2 and n < 3:
                    ccr = CR
                if m > 2 and n > 2:
                    ccr = CR * CR
                if m < 3 and n > 2:
                    ccr = CR
                i = CO2[m, 0] - 1
                j = CO2[m, 1] - 1
                k = CO2[n, 0] - 1
                l = CO2[n, 1] - 1
                C[m, n] = A[i, j, k, l] * ccr

    else:
        SupportError('Tensor4_to_Matrix6')
    return C


def MatrixNorm(A, Type):
    CheckSquare(A)
    Norm = 0.0
    if Type == 'Euklid':
        charMatrix = numpy.dot(A, numpy.transpose(A))
        evalu = linalg.eigvals(charMatrix)
        Norm = numpy.sqrt(abs(max(evalu)))
    elif Type == 'Frobenius':
        shapeMatrix = A.shape
        Norm = 0.0
        for m in range(shapeMatrix[0]):
            for n in range(shapeMatrix[1]):
                Norm = Norm + A[m, n] * A[m, n]

        Norm = numpy.sqrt(Norm)
    else:
        print '\n** ERROR: Norm Type "', Type, '" not known!'
        exit(1)
    return Norm


def UnitMatrix(n):
    I = numpy.zeros((n, n), numpy.float)
    for row in range(n):
        for col in range(n):
            if row == col:
                I[col, row] = 1.0

    return I


def CheckPosition(A, B):
    ash = A.shape
    bsh = B.shape
    if ash[len(ash) - 1] < bsh[0]:
        print '\n** ERROR: Inconsistent Shape  A=', ash, ' B=', bsh
        exit(1)


def CheckShape(A):
    ash = A.shape
    for index in range(len(ash)):
        if ash[index] != 3:
            print '\n** ERROR: Order of Tensor', ash, 'is not correct\n'
            exit(1)


def CheckSquare(A):
    ash = A.shape
    if ash[0] != ash[1]:
        print '\n** ERROR: Matrix is not a square Matrix\n'
        exit(1)


def SupportError(Operation):
    print '\n** ERROR:', Operation, 'Option not supported!\n'
    exit(1)


class TensorTest(unittest.TestCase):
    R1 = numpy.array([11, 12, 13])
    S1 = numpy.array([14, 15, 16])
    T1 = numpy.array([2, 3, 4])
    U1 = numpy.array([5, 6, 7])
    S2 = numpy.array([[7, 8, 9],
     [
      1, 2, 3],
     [
      4, 5, 6]])
    T2 = numpy.array([[1, 2, 3],
     [
      4, 5, 6],
     [
      7, 8, 9]])
    T3 = numpy.array([[[11, 21, 31], [41, 51, 61], [71, 81, 91]],
     [
      [
       12, 22, 32], [42, 52, 62], [72, 82, 92]],
     [
      [
       13, 23, 33], [43, 53, 63], [73, 83, 93]]])
    T4 = numpy.array([
     [
      [
       [
        11, 21, 31], [41, 51, 61], [71, 81, 91]],
      [
       [
        12, 22, 32], [42, 52, 62], [72, 82, 92]],
      [
       [
        13, 23, 33], [43, 53, 63], [73, 83, 93]]],
     [
      [
       [
        116, 216, 316], [416, 516, 616], [716, 816, 916]],
      [
       [
        126, 226, 326], [426, 526, 626], [726, 826, 926]],
      [
       [
        136, 236, 336], [436, 536, 636], [736, 836, 936]]],
     [
      [
       [
        119, 219, 319], [419, 519, 619], [719, 819, 919]],
      [
       [
        129, 229, 329], [429, 529, 629], [729, 829, 929]],
      [
       [
        139, 239, 339], [439, 539, 639], [739, 839, 939]]]])
    e1n = numpy.array([1, 0, 0])
    e2n = numpy.array([0, 1, 0])
    e3n = numpy.array([0, 0, 1])
    BN = numpy.array([e1n, e2n, e3n])
    pi = numpy.pi
    phi = 40.0 * pi / 180.0
    det = 30.0 * pi / 180.0
    sphi = numpy.sin(phi)
    cphi = numpy.cos(phi)
    sdet = numpy.sin(det)
    cdet = numpy.cos(det)
    e1o = numpy.array([cdet * cphi, cdet * sphi, sdet])
    e2o = numpy.array([-sphi, cphi, 0])
    e3o = numpy.array([-sdet * cphi, -sdet * sphi, cdet])
    BO = numpy.array([e1o, e2o, e3o])
    RO1 = numpy.array([3, 2, 1])
    RO2 = numpy.array([1, 2, 3])
    RA1 = numpy.array([0.0, 30.0, -40.0])
    RA2 = numpy.array([0.0, -30.0, 40.0])
    checkDyad1 = numpy.array([[154.0, 165.0, 176.0],
     [
      168.0, 180.0, 192.0],
     [
      182.0, 195.0, 208.0]])
    checkDyad2 = numpy.array([
     [[2.0, 3.0, 4.0],
      [
       4.0, 6.0, 8.0],
      [
       6.0, 9.0, 12.0]],
     [
      [
       8.0, 12.0, 16.0],
      [
       10.0, 15.0, 20.0],
      [
       12.0, 18.0, 24.0]],
     [
      [
       14.0, 21.0, 28.0],
      [
       16.0, 24.0, 32.0],
      [
       18.0, 27.0, 36.0]]])
    checkSimpContr1 = numpy.array([[102.0, 126.0, 150.0],
     [
      30.0, 36.0, 42.0],
     [
      66.0, 81.0, 96.0]])
    checkDoubContr = numpy.array([[2895.0, 2940.0, 2985.0],
     [
      29220.0, 29670.0, 30120.0],
     [
      29355.0, 29805.0, 30255.0]])
    checkCrossProd = numpy.array([-3, 6, -3])
    checkTransBas1 = numpy.array([-2.13362382, 2.12589891, 4.46410162])
    checkTransBas2 = numpy.array([[2.17355052, -4.63962191, -6.89623119],
     [
      -1.9075711, 1.49632246, 3.60728478],
     [
      -4.3026066, 6.73921868, 11.33012702]])
    checkTransBas3 = numpy.array([
     [[-8.10947434, 17.45653246, 25.86849009],
      [
       7.55569281, -7.19852174, -15.49539021],
      [
       16.46930138, -26.84538969, -44.36507019]],
     [
      [
       21.9119103, -46.865591, -69.61012692],
      [
       -19.50899151, 16.11054265, 37.65873274],
      [
       -43.63961833, 69.01944568, 115.54943629]],
     [
      [
       30.09812502, -64.50262715, -95.73794769],
      [
       -27.18211907, 23.54641191, 53.51406803],
      [
       -60.308377, 96.29708062, 160.55350227]]])
    checkTransBas4 = numpy.array([
     [
      [[78.15546773, -168.45041641, -249.51033029],
       [
        -73.45498502, 71.72145455, 152.29390452],
       [
        -159.3280813, 261.19356712, 430.60805138]],
      [
       [
        -211.19836053, 452.31446378, 671.50772271],
       [
        189.83652638, -161.90848183, -371.32738048],
       [
        422.32889, -672.22352258, -1122.30768534]],
      [
       [
        -290.09234233, 622.50287454, 923.51459592],
       [
        264.42437924, -235.92411653, -527.09790745],
       [
        583.57904796, -937.58740944, -1559.0781627]]],
     [
      [
       [
        -40.65422024, 87.56216465, 129.73027883],
       [
        38.02651359, -36.63453119, -78.37070482],
       [
        82.70447609, -135.15675872, -223.11841421]],
      [
       [
        109.85320762, -235.09589113, -349.11629909],
       [
        -98.22598772, 82.3143763, 190.74677841],
       [
        -219.18125349, 347.64996926, 581.29764596]],
      [
       [
        150.89189006, -323.56257483, -480.14600216],
       [
        -136.84147534, 120.14034252, 270.92373182],
       [
        -302.88559668, 484.97470709, 807.62040819]]],
     [
      [
       [
        -74.65580149, 160.93124625, 238.36023462],
       [
        70.23703987, -68.77235195, -145.8052783],
       [
        152.26127825, -249.77415691, -411.66608191]],
      [
       [
        201.74363912, -432.13278333, -641.51001375],
       [
        -181.53938238, 155.40182221, 355.63882199],
       [
        -403.61355757, 642.91091934, 1073.02563247]],
      [
       [
        277.10476499, -594.72401837, -882.2546905],
       [
        -252.85868358, 226.36637818, 504.76589986],
       [
        -557.71082251, 896.66916373, 1490.57865328]]]])
    checkTransRot1 = numpy.array([-2.3431833, 3.58370855, 3.26596464])
    checkTransRot2 = numpy.array([[2.60333334, -7.43454568, -5.25170955],
     [
      -3.65087521, 6.30171961, 5.4261594],
     [
      -3.47310699, 7.97966889, 6.09494704]])

    def test_UnitVector(self):
        print 'test_UnitVector'
        outVal = UnitVector(self.R1)
        checkVal = self.R1 / numpy.linalg.norm(self.R1)
        self.assertTrue(numpy.allclose(checkVal, outVal, rtol=1e-05, atol=1e-08))

    def test_UnitMatrix(self):
        print 'test_UnitMatrix'
        outVal = UnitMatrix(3)
        checkVal = numpy.eye(3)
        self.assertTrue(numpy.allclose(checkVal, outVal, rtol=1e-05, atol=1e-08))

    def test_Length(self):
        print 'test_Length'
        outVal = Length(self.R1)
        checkVal = numpy.linalg.norm(self.R1)
        self.assertTrue(numpy.allclose(checkVal, outVal, rtol=1e-05, atol=1e-08))

    def test_MatrixNorm(self):
        print 'test_MarixNorm'
        outVal = MatrixNorm(self.S2, 'Frobenius')
        checkVal = numpy.linalg.norm(self.S2)
        self.assertTrue(numpy.allclose(checkVal, outVal, rtol=1e-05, atol=1e-08))
        outVal2 = MatrixNorm(self.S2, 'Euklid')
        checkVal2 = numpy.linalg.norm(self.S2, 2)
        self.assertTrue(numpy.allclose(checkVal2, outVal2, rtol=1e-05, atol=1e-08))

    def test_CrossProduct(self):
        print 'test_CrossProduct'
        outVal = CrossProduct(self.T1, self.U1)
        self.assertTrue(numpy.allclose(self.checkCrossProd, outVal, rtol=1e-05, atol=1e-08))
        self.assertTrue(numpy.allclose(numpy.cross(self.T1, self.U1), outVal, rtol=1e-05, atol=1e-08))

    def test_DyadicProduct(self):
        print 'test_DyadicProduct'
        outVal = DyadicProduct(self.R1, self.S1)
        self.assertTrue(numpy.allclose(self.checkDyad1, outVal, rtol=1e-05, atol=1e-08))
        outVal = DyadicProduct(self.T2, self.T1)
        self.assertTrue(numpy.allclose(self.checkDyad2, outVal, rtol=1e-05, atol=1e-08))

    def test_SimpleContraction(self):
        print 'test_SimpleContraction'
        outVal = SimpleContraction(self.R1, self.S1)
        self.assertTrue(numpy.allclose(542.0, outVal, rtol=1e-05, atol=1e-08))
        outVal1 = numpy.dot(self.R1, self.S1)
        self.assertTrue(numpy.allclose(outVal1, outVal, rtol=1e-05, atol=1e-08))
        outVal = SimpleContraction(self.T2, self.T1)
        self.assertTrue(numpy.allclose(numpy.array([20.0, 47.0, 74.0]), outVal, rtol=1e-05, atol=1e-08))
        outVal2 = numpy.dot(self.T2, self.T1)
        self.assertTrue(numpy.allclose(outVal2, outVal, rtol=1e-05, atol=1e-08))
        outVal = SimpleContraction(self.S2, self.T2)
        self.assertTrue(numpy.allclose(self.checkSimpContr1, outVal, rtol=1e-05, atol=1e-08))
        outVal3 = numpy.dot(self.S2, self.T2)
        self.assertTrue(numpy.allclose(outVal3, outVal, rtol=1e-05, atol=1e-08))

    def test_DoubleContraction(self):
        print 'test_DoubleContraction'
        outVal = DoubleContraction(self.S2, self.T2)
        self.assertTrue(numpy.allclose(204.0, outVal, rtol=1e-05, atol=1e-08))
        outVal = DoubleContraction(self.T4, self.T2)
        self.assertTrue(numpy.allclose(self.checkDoubContr, outVal, rtol=1e-05, atol=1e-08))

    def test_TensorTransformBas(self):
        print 'test_TensorTransformBas'
        outVal1 = TensorTransformBas(self.T1, self.BO, self.BN)
        self.assertTrue(numpy.allclose(self.checkTransBas1, outVal1, rtol=1e-05, atol=1e-08))
        outVal2 = TensorTransformBas(self.T2, self.BO, self.BN)
        self.assertTrue(numpy.allclose(self.checkTransBas2, outVal2, rtol=1e-05, atol=1e-08))
        outVal3 = TensorTransformBas(self.T3, self.BO, self.BN)
        self.assertTrue(numpy.allclose(self.checkTransBas3, outVal3, rtol=1e-05, atol=1e-08))
        outVal4 = TensorTransformBas(self.T4, self.BO, self.BN)
        self.assertTrue(numpy.allclose(self.checkTransBas4, outVal4, rtol=1e-05, atol=1e-08))
        outVal1b = TensorTransformBas(outVal1, self.BN, self.BO)
        self.assertTrue(numpy.allclose(self.T1, outVal1b, rtol=1e-05, atol=1e-08))
        outVal2b = TensorTransformBas(outVal2, self.BN, self.BO)
        self.assertTrue(numpy.allclose(self.T2, outVal2b, rtol=1e-05, atol=1e-08))
        outVal3b = TensorTransformBas(outVal3, self.BN, self.BO)
        self.assertTrue(numpy.allclose(self.T3, outVal3b, rtol=1e-05, atol=1e-08))
        outVal4b = TensorTransformBas(outVal4, self.BN, self.BO)
        self.assertTrue(numpy.allclose(self.T4, outVal4b, rtol=1e-05, atol=1e-08))

    def test_TensorTransformRot(self):
        print 'test_TensorTransformRot'
        outVal1 = TensorTransformRot(self.T1, self.RO1, self.RA1)
        self.assertTrue(numpy.allclose(self.checkTransRot1, outVal1, rtol=1e-05, atol=1e-08))
        outVal2 = TensorTransformRot(self.T2, self.RO1, self.RA1)
        self.assertTrue(numpy.allclose(self.checkTransRot2, outVal2, rtol=1e-05, atol=1e-08))
        outVal3 = TensorTransformRot(self.T3, self.RO1, self.RA1)
        outVal4 = TensorTransformRot(self.T4, self.RO1, self.RA1)
        outVal1b = TensorTransformRot(outVal1, self.RO2, self.RA2)
        self.assertTrue(numpy.allclose(self.T1, outVal1b, rtol=1e-05, atol=1e-08))
        outVal2b = TensorTransformRot(outVal2, self.RO2, self.RA2)
        self.assertTrue(numpy.allclose(self.T2, outVal2b, rtol=1e-05, atol=1e-08))
        outVal3b = TensorTransformRot(outVal3, self.RO2, self.RA2)
        self.assertTrue(numpy.allclose(self.T3, outVal3b, rtol=1e-05, atol=1e-08))
        outVal4b = TensorTransformRot(outVal4, self.RO2, self.RA2)
        self.assertTrue(numpy.allclose(self.T4, outVal4b, rtol=1e-05, atol=1e-08))


if __name__ == '__main__':
    unittest.main()
    exit(0)
