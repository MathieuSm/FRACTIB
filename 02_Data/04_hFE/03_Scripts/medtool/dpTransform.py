# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/misc/dpTransform.py
# Compiled at: 2016-08-24 21:21:41
"""
#########################################################################
File   : dpTransform.py
Author : D.H.Pahr
Date   : 14.2.2015
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage: python dpTensor.py   

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requirements:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Description: 
  This program is a transform library   
  
#########################################################################
"""
import numpy
from sys import stdout, exit
import os
import unittest

class Transformer():
    __version__ = 'V_14.02.2015'
    __author__ = 'D. H. Pahr'

    def __init__(self):
        """ Transformer class  """
        self.transformMatrix = numpy.eye(4)
        self.transformMatrixInv = numpy.eye(4)
        self.transformDirection = 1

    def setTransformDirection(self, direction, verbose=0):
        """
        Set the transform direction. Forward (+1) or backward/inverse (-1).
        
        direction : int
                    plus (forward) or minus (backward) 
        """
        if verbose > 0:
            print ' .. set transform direction to %i' % direction
        if abs(direction) == 1:
            self.transformDirection = int(direction)
        else:
            stdout.write('\n **ERROR**: dpTransform.setTransformDirection(): direction has to be -1 or 1!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

    def transformVector(self, vecIn, verbose=0):
        """
        Transform a 3x1 vector.
        
        Parameters
        ----------
        vecIn : (3,1) ndarray
                Vector to be transformed.
                
        Returns
        -------
        (3,1) ndarray
                Transformed vector.
        """
        if verbose > 0:
            print ' .. transform vector'
        vecIn4 = numpy.ones(4)
        vecIn4[:3] = vecIn
        if self.transformDirection > 0:
            M = self.transformMatrix
        elif self.transformDirection < 0:
            M = self.transformMatrixInv
        vecOut4 = numpy.dot(M, vecIn4)
        if verbose > 1:
            print 'Transformed Vector=\n', vecOut4, '\n'
        return vecOut4[:3]

    def transformTensor2(self, matIn, verbose=0):
        """
        Transform a 3x3 matrix (tensor of second order)
        
        Parameters
        ----------
        matIn : (3,3) ndarray
                2nd oder tensor (matrix) to be transformed.
                
        Returns
        -------
        (3,3) ndarray
                Transformed 2nd order tensor.
        """
        if verbose > 0:
            print ' .. transform tensor 2 (matrix)'
        M = self.transformMatrix[:3, :3]
        Minv = self.transformMatrixInv[:3, :3]
        if self.transformDirection > 0:
            matOut = numpy.dot(M, numpy.dot(matIn, Minv))
        elif self.transformDirection < 0:
            matOut = numpy.dot(Minv, numpy.dot(matIn, M))
        if verbose > 1:
            print 'Transformed matrix tensor 2 (matrix)\n', matOut, '\n'
        return matOut

    def readTransformMatrix(self, inFileName, verbose=0):
        """ 
        Read and sets the transformation matrix form a file (generic reader).
        
        Parameters
        ----------
        inFileName : string
                     Name of the file where transform information is stored.      
        
        """
        name, ext = os.path.splitext(inFileName)
        ext = ext.lower()
        if ext == '.mhd':
            self.readTransformMatrixMhd(inFileName, verbose=verbose)
        elif ext == '.nhdr':
            self.readTransformMatrixNhdr(inFileName, verbose=verbose)
        elif ext == '.tfm':
            self.readTransformMatrixTfm(inFileName, verbose=verbose)
        else:
            stdout.write('\n **ERROR**: dpTransform.readTransformMatrix() file extension of file: "%s" not known!\n\n' % inFileName)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

    def readTransformMatrixMhd(self, inFileName, verbose=0):
        """
        Read and sets the transformation matrix form a file (mhd reader).
        
        Transformation matrix inside this file maps points from input space 
        (measurement space) to output space (world space). This is known as 
        forward transform. The backward transform is obtained by appling the 
        inverse. 
        
        
        Parameters
        ----------
        inFileName : string
                     Name of the file where transform information is stored.  
        
        """
        if verbose > 0:
            print ' .. read transform matrix from mhd file'
        f = open(inFileName, 'r')
        lines = f.readlines()
        f.close()
        parDict = {}
        for line in lines:
            if len(line) > 1:
                line = line.replace('\n', '')
                key, data = line.split('=')
                key = key.replace(' ', '')
                parDict[key] = data

        if 'TransformMatrix' in parDict:
            stm = parDict['TransformMatrix'].split()
        else:
            stdout.write('\n **ERROR**: dpTransform.readTransformMatrixMhd() keyword  "TransformMatrix" not found!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if 'Offset' in parDict:
            sof = parDict['Offset'].split()
        else:
            stdout.write('\n **ERROR**: dpTransform.readTransformMatrixMhd() keyword  "Offset" not found!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if 'CenterOfRotation' in parDict:
            scor = parDict['CenterOfRotation'].split()
            if abs(float(scor[0])) > 1e-07 or abs(float(scor[1])) > 1e-07 or abs(float(scor[2])) > 1e-07:
                stdout.write('\n **ERROR**: dpTransform.readTransformMatrixMhd() keyword  "CenterOfRotation" not at \'0 0 0\'!\n\n')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        self.transformMatrix = numpy.array([[float(stm[0]), float(stm[3]), float(stm[6]), float(sof[0])],
         [
          float(stm[1]), float(stm[4]), float(stm[7]), float(sof[1])],
         [
          float(stm[2]), float(stm[5]), float(stm[8]), float(sof[2])],
         [
          0.0, 0.0, 0.0, 1.0]])
        self.transformMatrixInv = numpy.linalg.inv(self.transformMatrix)
        if verbose > 1:
            print 'Transform Matrix=\n', self.transformMatrix, '\n'

    def readTransformMatrixNhdr(self, inFileName, verbose=0):
        """
        Read and sets the transformation matrix form a file (nhdr reader).
        
        The transformation matrix inside this file maps points from input 
        space (measurement or old space) to output space (world or new 
        space). This is known as forward transform. The backward transform 
        is obtained by appling the inverse. 
        
        
        Parameters
        ----------
        inFileName : string
                     Name of the file where transform information is stored.  
        """
        if verbose > 0:
            print ' .. read transform matrix from nhdr file'
        f = open(inFileName, 'r')
        lines = f.readlines()
        f.close()
        parDict = {}
        lineNo = 0
        for line in lines:
            if lineNo > 0 and line[0] != '#' and len(line) > 1:
                line = line.replace('\n', '')
                key, data = line.split(':')
                key = key.replace(' ', '')
                parDict[key] = data
            lineNo += 1

        nvec = []
        length = []
        if 'spacedirections' in parDict:
            svec = [
             ' ', ' ', ' ']
            svec[0], svec[1], svec[2] = parDict['spacedirections'].replace(' ', '').split(')(')
            svec[0] = svec[0].replace('(', '')
            svec[2] = svec[2].replace(')', '')
            vec = []
            for ii in range(3):
                sx1, sx2, sx3 = svec[ii].split(',')
                vec.append(numpy.array([float(sx1), float(sx2), float(sx3)]))
                length.append(numpy.linalg.norm(vec[ii]))
                nvec.append(vec[ii] / length[ii])

        else:
            stdout.write("\n **ERROR**: dpTransform.readTransformMatrixNhdr() keyword 'space directions' not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        OF = [0.0, 0.0, 0.0]
        if 'spaceorigin' in parDict:
            svec = parDict['spaceorigin'].replace(' ', '').replace('(', '').replace(')', '')
            ox, oy, oz = svec.split(',')
            OF = [float(ox), float(oy), float(oz)]
        vec[0] = vec[0] / length[0]
        vec[1] = vec[1] / length[1]
        vec[2] = vec[2] / length[2]
        self.transformMatrix = numpy.array([[vec[0][0], vec[1][0], vec[2][0], OF[0]],
         [
          vec[0][1], vec[1][1], vec[2][1], OF[1]],
         [
          vec[0][2], vec[1][2], vec[2][2], OF[2]],
         [
          0.0, 0.0, 0.0, 1.0]])
        self.transformMatrixInv = numpy.linalg.inv(self.transformMatrix)
        if verbose > 1:
            print 'Transform Matrix=\n', self.transformMatrix, '\n'

    def readTransformMatrixTfm(self, inFileName, verbose=0):
        """
        Read and sets the transformation matrix form a file (tfm reader).
        
        The information is stored in ITK style an could be directly used in 
        an ITK resampling filter. Compared to a classical forward transform 
        (from input space to output space) the given transformation matrix 
        inside this file maps points from the output space (world or new 
        space) back to the input space (measurement or old space). This is 
        known as backward transform. The forward transform is obtained by 
        appling the inverse of the transformation matrix. 
        
        Parameters
        ----------
        inFileName : string
                     Name of the file where transform information is stored.  
        """
        if verbose > 0:
            print ' .. read transform matrix from tfm file'
        f = open(inFileName, 'r')
        foundFlag = False
        for line in f:
            if line.find('Parameters') == 0:
                foundFlag = True
                name, data = line.split(':')
                sld = data.split()
                self.transformMatrix = numpy.array([[float(sld[0]), float(sld[1]), float(sld[2]), float(sld[9])],
                 [
                  float(sld[3]), float(sld[4]), float(sld[5]), float(sld[10])],
                 [
                  float(sld[6]), float(sld[7]), float(sld[8]), float(sld[11])],
                 [
                  0.0, 0.0, 0.0, 1.0]])

        f.close()
        if not foundFlag:
            stdout.write('\n **ERROR**: dpTransform.readTransformMatrixTfm() keyword  "Parameters" not found!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        self.transformMatrixInv = numpy.linalg.inv(self.transformMatrix)
        if verbose > 1:
            print 'Transform Matrix=\n', self.transformMatrix, '\n'

    def getTransformMatrix4(self):
        """
        Returns 4x4 Transform Matrix.
        
        Returns
        -------
        
        (4,4) ndarray 
              The transform matrix in homogeneous coordinates. 
               
        """
        return self.transformMatrix

    def getTransformMatrixInv4(self):
        """
        Returns 4x4 Inverse Transform Matrix.
        
        Returns
        -------
        
        (4,4) ndarray 
              The transform matrix in homogeneous coordinates. 
               
        """
        return self.transformMatrixInv

    def setTransformMatrix4(self, matrix):
        """
        Sets 4x4 Transform Matrix directly. 
        
        Parameters
        ----------
        matrix : (4,4) ndarray
                 The transform matrix in homogeneous coordinates. 
        """
        if not matrix.shape == (4, 4):
            stdout.write('\n **ERROR**: dpTransform.setTransformMatrix4: Matrix dimension has to be 4x4 !\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        self.transformMatrix = matrix
        self.transformMatrixInv = numpy.linalg.inv(self.transformMatrix)

    def setTransformMatrixBasis(self, MO, BN):
        """
        Sets 4x4 transformation matrix from basis vectors .
        
        Parameters
        ----------
        MO : (3,3) ndarray
             The old or measurement basis where data are stored. 
             The data are stored as MO = [m1, m2, m3] where m1, m2, m3 are 
             the three eigenvectors of this basis measured in the 
             new world-space bases B. MO[0] gives the first eigenvector m1, 
             MO[1] the second eigenvector m2 etc.  
             
        BN : (3,3) ndarray
             The new or world-space basis to which data should be transformed.
             The data are stored as BN = [b1, b2, b3] where b1, b2, b3 are 
             the three eigenvectors and are usually b1=[1,0,0], b2=[0,1,0], etc.
             BO[0] gives the first eigenvector b1, BO[1] the second ie. b2, and 
             BO[2] the third eigenvector b3. 
        
        """
        if not BN.shape == (3, 3) and not MO.shape == (3, 3):
            stdout.write('\n **ERROR**: dpTransform.setTransformMatrixBasis: Matrix dimensions have to be 3x3 !\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        b1 = BN[0]
        b2 = BN[1]
        b3 = BN[2]
        m1 = MO[0]
        m2 = MO[1]
        m3 = MO[2]
        self.transformMatrix = numpy.array([[numpy.dot(b1, m1), numpy.dot(b1, m2), numpy.dot(b1, m3), 0.0],
         [
          numpy.dot(b2, m1), numpy.dot(b2, m2), numpy.dot(b2, m3), 0.0],
         [
          numpy.dot(b3, m1), numpy.dot(b3, m2), numpy.dot(b3, m3), 0.0],
         [
          0.0, 0.0, 0.0, 1.0]])
        self.transformMatrixInv = numpy.linalg.inv(self.transformMatrix)


class TransformerTest(unittest.TestCase):
    verbose = 1
    stfm = []
    stfm.append('#Insight Transform File V1.0')
    stfm.append('#Transform 0')
    stfm.append('Transform: AffineTransform_double_3_3')
    stfm.append('Parameters: 0.866025 0.5 0 -0.5  0.866025 0 0 0 1.0 -1.79903811 -0.1160254 -0.5')
    stfm.append('FixedParameters: 0 0 0')
    smhd = []
    smhd.append('ObjectType = Image')
    smhd.append('NDims = 3')
    smhd.append('BinaryData = True')
    smhd.append('BinaryDataByteOrderMSB = False')
    smhd.append('CompressedData = False')
    smhd.append('TransformMatrix = 0.866025 0.5 0.0 -0.5  0.866025 0.0 0.0 0.0 1.0')
    smhd.append('Offset = 1.5 1.0 0.5')
    smhd.append('CenterOfRotation = 0 0 0')
    smhd.append('AnatomicalOrientation = LPS')
    smhd.append('ElementSpacing = 0.1 0.1 0.1')
    smhd.append('DimSize = 50 40 30')
    smhd.append('ElementType = MET_UCHAR')
    smhd.append('ElementDataFile = coordinates_50_40_30.raw')
    snhdr = []
    snhdr.append('NRRD0005')
    snhdr.append('# Complete NRRD file format specification at:')
    snhdr.append('# http://teem.sourceforge.net/nrrd/format.html')
    snhdr.append('type: unsigned char')
    snhdr.append('dimension: 3')
    snhdr.append('space: LPS')
    snhdr.append('sizes: 50 40 30')
    snhdr.append('space directions: (0.4330125, 0.25, 0.0) (-0.25, 0.4330125, 0.0) (0.0, 0.0, 0.5) ')
    snhdr.append('kinds: domain domain domain')
    snhdr.append('endian: little')
    snhdr.append('encoding: raw')
    snhdr.append('space origin: ( 1.5, 1.0, 0.5 )')
    snhdr.append('data file: coordinates_50_40_30.raw')
    inTM = numpy.array([[0.866025, -0.5, 0.0, 1.5],
     [
      0.5, 0.866025, 0.0, 1.0],
     [
      0.0, 0.0, 1.0, 0.5],
     [
      0.0, 0.0, 0.0, 1.0]])
    inVec0 = numpy.array([0.0, 0.0, 0.0])
    inVec1 = numpy.array([1.0, 0.0, 0.0])
    inVec2 = numpy.array([0.0, 1.0, 0.0])
    inVec3 = numpy.array([0.0, 0.0, 1.0])
    outVec0 = numpy.array([1.5, 1.0, 0.5])
    outVec1 = numpy.array([2.366025, 1.5, 0.5])
    outVec2 = numpy.array([1.0, 1.866025, 0.5])
    outVec3 = numpy.array([1.5, 1.0, 1.5])
    inTM2 = numpy.array([[1.0, 0.0, 0.0, 0.0],
     [
      0.0, 0.866, 0.5, 0.0],
     [
      0.0, -0.5, 0.866, 0.0],
     [
      0.0, 0.0, 0.0, 1.0]])
    inMat = numpy.array([[1.0, 2.0, 3.0],
     [
      2.0, 4.0, 5.0],
     [
      3.0, 5.0, 6.0]])
    outMat = numpy.array([[1.0, 3.232142, 1.59807],
     [
      3.232, 8.830212, 3.365928],
     [
      1.598, 3.3659281, 1.169787]])
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
    T1 = numpy.array([2, 3, 4])
    T2 = numpy.array([[1, 2, 3],
     [
      4, 5, 6],
     [
      7, 8, 9]])
    checkTransBas1 = numpy.array([-2.13362382, 2.12589891, 4.46410162])
    checkTransBas2 = numpy.array([[2.17355052, -4.63962191, -6.89623119],
     [
      -1.9075711, 1.49632246, 3.60728478],
     [
      -4.3026066, 6.73921868, 11.33012702]])
    inTM3 = numpy.array([[1.0, 0.0, 0.0, 10.0],
     [
      0.0, 1.0, 0.0, 20.0],
     [
      0.0, 0.0, 1.0, 30.0],
     [
      0.0, 0.0, 0.0, 1.0]])

    def test_readTransformMatrixTfm(self):
        print 'test_readTransformMatrixTfm'
        f = open('test.tfm', 'w')
        for line in self.stfm:
            f.write('%s\n' % line)

        f.close()
        myTransformer = Transformer()
        myTransformer.readTransformMatrix('test.tfm', self.verbose)
        outTM = myTransformer.getTransformMatrixInv4()
        self.assertTrue(numpy.allclose(self.inTM, outTM, rtol=1e-05, atol=1e-08))
        os.remove('test.tfm')

    def test_readTransformMatrixMhd(self):
        print 'test_readTransformMatrixMhd'
        f = open('test.mhd', 'w')
        for line in self.smhd:
            f.write('%s\n' % line)

        f.close()
        myTransformer = Transformer()
        myTransformer.readTransformMatrix('test.mhd', self.verbose)
        outTM = myTransformer.getTransformMatrix4()
        self.assertTrue(numpy.allclose(self.inTM, outTM, rtol=1e-05, atol=1e-08))
        os.remove('test.mhd')

    def test_readTransformMatrixNhdr(self):
        print 'test_readTransformMatrixNhdr'
        f = open('test.nhdr', 'w')
        for line in self.snhdr:
            f.write('%s\n' % line)

        f.close()
        myTransformer = Transformer()
        myTransformer.readTransformMatrix('test.nhdr', self.verbose)
        outTM = myTransformer.getTransformMatrix4()
        self.assertTrue(numpy.allclose(self.inTM, outTM, rtol=1e-05, atol=1e-08))
        os.remove('test.nhdr')

    def test_transformVector(self):
        print 'test_transformVector'
        myTransformer = Transformer()
        myTransformer.setTransformMatrix4(self.inTM)
        outVec = myTransformer.transformVector(self.inVec0, self.verbose)
        self.assertTrue(numpy.allclose(self.outVec0, outVec, rtol=1e-05, atol=1e-08))
        outVec = myTransformer.transformVector(self.inVec1, self.verbose)
        self.assertTrue(numpy.allclose(self.outVec1, outVec, rtol=1e-05, atol=1e-08))
        outVec = myTransformer.transformVector(self.inVec2, self.verbose)
        self.assertTrue(numpy.allclose(self.outVec2, outVec, rtol=1e-05, atol=1e-08))
        outVec = myTransformer.transformVector(self.inVec3, self.verbose)
        self.assertTrue(numpy.allclose(self.outVec3, outVec, rtol=1e-05, atol=1e-08))

    def test_transformVectorBackward(self):
        print 'test_transformVectorBackward'
        myTransformer = Transformer()
        myTransformer.setTransformMatrix4(self.inTM)
        myTransformer.setTransformDirection(1, self.verbose)
        outVec1 = myTransformer.transformVector(self.inVec0, self.verbose)
        myTransformer.setTransformDirection(-1, self.verbose)
        outVec2 = myTransformer.transformVector(outVec1, self.verbose)
        self.assertTrue(numpy.allclose(self.inVec0, outVec2, rtol=1e-05, atol=1e-08))
        myTransformer.setTransformDirection(1, self.verbose)
        outVec1 = myTransformer.transformVector(self.inVec1, self.verbose)
        myTransformer.setTransformDirection(-1, self.verbose)
        outVec2 = myTransformer.transformVector(outVec1, self.verbose)
        self.assertTrue(numpy.allclose(self.inVec1, outVec2, rtol=1e-05, atol=1e-08))

    def test_transformTensor2(self):
        print 'test_transformTensor2'
        myTransformer = Transformer()
        myTransformer.setTransformMatrix4(self.inTM2)
        outMat = myTransformer.transformTensor2(self.inMat, self.verbose)
        self.assertTrue(numpy.allclose(self.outMat, outMat, rtol=1e-05, atol=1e-08))

    def test_transformTensorBasis(self):
        print 'test_transformTensorBasis'
        myTransformer = Transformer()
        myTransformer.setTransformMatrixBasis(self.BO, self.BN)
        outVec = myTransformer.transformVector(self.T1, self.verbose)
        self.assertTrue(numpy.allclose(self.checkTransBas1, outVec, rtol=1e-05, atol=1e-08))
        outMat = myTransformer.transformTensor2(self.T2, self.verbose)
        self.assertTrue(numpy.allclose(self.checkTransBas2, outMat, rtol=1e-05, atol=1e-08))
        myTransformer.setTransformDirection(-1, self.verbose)
        outVec2 = myTransformer.transformVector(outVec, self.verbose)
        self.assertTrue(numpy.allclose(self.T1, outVec2, rtol=1e-05, atol=1e-08))
        outMat2 = myTransformer.transformTensor2(outMat, self.verbose)
        self.assertTrue(numpy.allclose(self.T2, outMat2, rtol=1e-05, atol=1e-08))

    def test_tfm_transform(self):
        print 'test_tfm_transform'
        myTransformer = Transformer()
        myTransformer.setTransformMatrix4(self.inTM3)
        print myTransformer.getTransformMatrixInv4()


if __name__ == '__main__':
    unittest.main()
# okay decompiling dpTransform.pyc
