# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/dpMesher/dpMesher.py
# Compiled at: 2016-01-16 23:35:21
"""
Script dpMesher
---------------


Info
~~~~

- File:   dpMesher.py
- Author: D. H. Pahr

"""
import time
import os
import sys
#if 'MEDTOOLPATH' in os.environ.keys():
#    myPath = os.environ['MEDTOOLPATH']
#    if myPath not in sys.path:
#        sys.path.append(os.environ['MEDTOOLPATH'])
#else:
#    sys.stdout.write('\n **ERROR** : MEDTOOLPATH not set on your system!\n')
#    sys.stdout.write('\n E N D E D  with  ERRORS oups \n\n')
#    sys.stdout.flush()
#    sys.exit(1)
import dpUtils
import mic
import dpMesherf77
import numpy
import scipy

def computeImageToSmoothMesh2D(curVoxelModel, dimList, smooth=None, fast=False, echo=1):
    """
    Compute smooth 2D surface mesh file from image file. 
    """
    fecho = 0
    if echo == 1:
        sys.stdout.write(' ... compute smooth 2D mesh from image\n')
        sys.stdout.write("     -> recast model from '%s' to 'i'\n" % curVoxelModel.dtype.char)
        fecho = 2
        sys.stdout.flush()
    if echo == 2:
        sys.stdout.write('     -> compute smooth 2D mesh from image\n')
        sys.stdout.flush()
    curVoxelModel = mic.mic().castType(curVoxelModel, 'i')
    if dimList == None:
        print '\n **ERROR** computeSmoothMesh(): Voxel size not optional for this function!\n'
        sys.stdout.write('\n E N D E D  with ERRORS \n\n')
        sys.stdout.flush()
        exit(1)
    if smooth == None:
        smooth = [7, 0.6307, 0.1, 0, 1, 0.05, 0.6]
    activeNodes, nElem, nNode, nIntElem, nIntFaces, nIntNode = dpMesherf77.check_voxmesh2d(curVoxelModel, smooth[4], fecho)
    if fast == True:
        nodeCoord, nodeCoordInt, noidF77PY, noidIntVoxF77 = dpMesherf77.get_voxmesh2d_nodes_fast(curVoxelModel, dimList, activeNodes, nNode, nIntNode, fecho)
    else:
        nodeCoord, nodeCoordInt, noidF77PY, noidIntVoxF77 = dpMesherf77.get_voxmesh2d_nodes(curVoxelModel, dimList, smooth, activeNodes, nNode, nIntElem, nIntNode, fecho)
    faces, facesSets = dpMesherf77.get_voxmesh2d_elements(curVoxelModel, nodeCoordInt, noidIntVoxF77, nIntFaces, fecho)
    return (
     nodeCoordInt, faces, facesSets)


if __name__ == '__main__':
    guiInit = "*modulName              'dpMesher'                                                      \n" + "*fileEntryIn     -in    'Binary image file name'     image.mhd    no   mhd:nhdr:raw:png:jpg:gif:bmp:tif:aim:bin:AIM:bst:isq:ISQ \n" + "*fileEntryOut    -out   'Output File Name'           testOut.txt  yes  txt;dat     \n" + "*subWindowStart         'Trabecular Bone'                                          \n" + "*splitEntry      -bs    'Bone Surface BS'            1:5          yes  int         \n" + "*splitEntry      -bv    'Bone Volume BV'             1:5          yes  int         \n" + "*entry           -tv    'Total Volume TV'            1            yes  int         \n" + "*entry           -bvtv  'Bone to Total Volume BV/TV' 1:5          yes  int         \n" + "*entry           -tbth  'Thickness Tb.Th'            None         yes  str         \n" + "*entry           -tbsp  'Spacing Tb.Sp'              None         yes  str         \n" + "*entry           -tbn   'Number Tb.N'                None         yes  str         \n" + "*splitEntry      -smi   'Struct Model Index SMI'     1:5          yes  int         \n" + "*splitEntry      -da    'Degree of Anisotropy DA'    1:5          yes  int         \n" + "*subWindowEnd           'Trabecular Bone'                                          \n" + "*subWindowStart         'Cortical Bone'                                            \n" + "*splitEntry      -pos   'Pore Surface Po.S'          1:5          yes  int         \n" + "*splitEntry      -pov   'Pore Volume Po.V'           1:5          yes  int         \n" + "*splitEntry      -por   'Porosity Po'                1:5          yes  int         \n" + "*splitEntry      -pon   'Pores Number Po.N'          1:5          yes  int         \n" + "*subWindowEnd           'Cortical Bone'                                            \n"
    args = dpUtils.getDataDict(guiInit)
    defaultDocStr = dpUtils.getDocString(guiInit, __file__, 'D. H. Pahr')
    dpUtils.initializeArguments(sys.argv, args, guiInit, __doc__, doc2=defaultDocStr)
    sys.stdout.write('\n S T A R T  %s %s - %s \n\n' % (
     __file__, 'V_07.04.2015', 'D. H. Pahr'))
    sys.stdout.flush()
    startTime = time.clock()
    ctime1 = time.time()
    endTime = time.clock()
    ctime2 = time.time()
    sys.stdout.write('\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (
     endTime - startTime, ctime2 - ctime1))
    sys.stdout.flush()
# okay decompiling dpMesher.pyc
