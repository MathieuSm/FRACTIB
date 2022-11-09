# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/fec/fec.py
# Compiled at: 2016-12-18 19:39:57
"""
Converter
---------

Script converts different unstructured grid (FEM) file formats. 
Implemented grid elements which can be read and write are: 
  
  - bar2
  - bar3
  - tria3
  - tria6
  - quad4
  - quad8
  - hexa8
  - hexa20
  - pyra5
  - penta6
  - penta15
  - tetra4
  - tetra10
                                                                                                                                  
Meshes can be refined, smoothed, cropped, translated and elset can be                                                             
extracted.                                                                                                                        
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
Usage                                                                                                                             
~~~~~                                                                                                                             
                                                                                                                                  
Module: ::                                                                                                                        
                                                                                                                                  
  import fec                                                                                                                      
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
Command line: ::                                                                                                                  
                                                                                                                                  
  python fec.py ...                                                                                                               
                                                                                                                                  
    -in         filename                                                                                                          
    -out        filename                                                                                                          
    [-bbcrop]   start1;start2;start3;end1;end2;end3    (optional)                                                                 
    [-extr]     setname                                (optional)                                                                 
    [-refin]    flag                                   (optional)                                                                 
    [-order]    order                                  (optional)
    [-smooth]   type;niter;[param2];[param3]           (optional)
    [-tra]      tra1;tra2;tra3                         (optional)
    [-tfil]     filename                               (optional)
    [-texp]     exponent                               (optional)
    [-histo]    filename;errorVal;normalize;nIntervals (optional)



Parameters
~~~~~~~~~~

-in   :  filename

         Name of the file which is to be converted. The file type is
         a recognized by the extension: Possible extensions are: 

           - 'inp': ABAQUS input file. Implemented keywords: ::
           
               *HEADING
               *NODE
               *ELEMENT
               *NSET
               *ELSET 
               
           - '.raw': ASCII file from dualcontour program
           - '.off': Geomview (only triangles implemented)
           - '.msh': Gmsh mesher
           - '.vtk': vtk ASCII file
           - '.stl': stereolithography file
           - '.obj': wavefront obj file
           - '.fem': Nastran/Optistruct file
         

-out  :  filename

         Name of the result files, converted file. Possible types are:

           - '.geom':  TAHOE II format 
           - '.geo':   PARAVIEW : Ensight gold format (same as '.case')
           - '.case':  PARAVIEW : Ensight gold format (same as '.geo')
           - '.inp':   ABAQUS input file
           - '.off':   Geomview file
           - '.msh':   GMSH file
           - '.stl':   stereolithography file
           - '.surf':  netgen's surface mesh file
           - '.smesh': tetgen's surface mesh file           
           - '.vtk':   vtk ASCII file


-bbcrop  : start1;start2;start3;end1;end2;end3

           Output only a part of the model which is inside the given bounding
           box given by 'start1', 'start2', 'start3', 'end1', 'end2', 'end3'. 
           The COG of the elements are used to define what is inside. 
           Parameters are the start and end points of the bounding box in 
           physical dimensions. 

           EXAMPLE :
           
           .. image:: images/fec-2.png


-extr  :  setname 

          Extracts a mesh from a given elset name 'setname'. 

          EXAMPLE :
           
          .. image:: images/fec-3.png


-refine  : flag

           If the 'flag' parameter is set to "ON" the mesh will be refined. 
           Only linear meshes will be refined.
           
           EXAMPLE :
           
           .. image:: images/fec-1.png
          

-order  : order

          Order of the mesh. This function converts linear to quadratic meshes.
          or quadratic to linear meshes. Allowed values for 'order' are 1 or 2. 
          Currently only 1st to 2nd order change is implemented. The nsets are not 
          modified by this function. Quadradic output is only fully supported 
          for ABAQUS and ENSIGHT (Paraview).


-smooth  : type;niter;[param2];[param3]

           Smooth the mesh. Currently only triangulated meshes (tria3) can be 
           smoothed. Laplacian and Taubin smoothing are possible. 
          
           Possibel 'type' values are: 
          
            - 'Laplace' for Laplacian. where 'niter'is the the number of 
              iterations. 'param2' and 'param3' are not needed. 
 
            - 'Taubin' for Taubin smoothing. where :

               - 'niter' =  is the number of iterations
               - 'param2' =  lam which is the scaling factor lambda (0<lambda<1)
               - 'param3' =  kPB which is the pass-band frequency kPB 

              Note that kPB = 1/lambda + 1/mu > 0).
              The implementation follows Taubin's algorithm (see, Taubin, 
              Eurographics 2000, Boyd and Mueller 2006). 


-tra    : tra1;tra2;tra3

          Translation of nodes in x,y,z (1,2,3) direction. This is the 
          same as 'space origin' or 'offset' in image files. This parameter 
          is ignored on case of 'tfil' option because there the movement
          comes directly from the transformation file. 


-tfil   : filename

          Transform nodal coordinates by using the information given in 
          a transformation 'filename'. Implemented file formats are:
          
            - 'tfm': ITK transformation file. For example, exported  
              from 3D Slicer. The important informations in this file 
              are 'Parameters' which gives the transformation matrix and 
              'FixedParameters' which gives the offset. 
              The ordering of 'Parameters' is ::
              
                11 12 13 21 22 23 31 32 33                
              
              and differes from the 'nhdr' and 'mhd' order scheme. 
              The information is stored in ITK style an could be directly used in 
              an ITK resampling filter. Compared to a classical forward transform 
              (from input space to output space) the given transformation matrix 
              inside this file maps points from the output space (world or new 
              space) back to the input space (measurement or old space). This is 
              known as backward transform. The forward transform is obtained by 
              appling the inverse of the transformation matrix. 
                        
            - 'nhdr': NRRD image header file with transformation 
              information in it. The important informations in this file 
              are 'space directions' which gives indirectly the transfromation 
              matrix and 'space origin' which gives the offset. The 
              transformation matrix is obtained from 'space directions' by 
              normalizing it because 'space directions' include the voxel length!
              The ordering of 'space directions' is ::
              
                11 21 31 12 22 32 13 23 33 
              
              The transformation matrix inside this file maps points from input 
              space (measurement or old space) to output space (world or new 
              space). This is known as forward transform. The backward transform 
              is obtained by appling the inverse. 
              
            - 'mhd': ITK meta image header file with transformation 
              information in it. The important informations in this file 
              are 'TransformMatrix' which gives the transfromation matrix and 
              'Offset' which gives the offset. The 'CenterOfRotation' has to be 
              (0,0,0). Non zero center of rotations are not implemented. The 
              ordering of 'TransformMatrix' is ::
              
                11 21 31 12 22 32 13 23 33 
                
              The transformation matrix inside this file maps points from input 
              space (measurement or old space) to output space (world or new 
              space). This is known as forward transform. The backward transform 
              is obtained by appling the inverse. 

          The implemented transformation are based on following theoretical
          background. Measurements are given in a basis M and should be
          transformed to a basis B. These are give by unit vectors namely 
          
            - b_i: unit vectors of world space basis B (new basis). This
              is usually b_1= (1,0,0), b_2=(0,1,0), b_3=(0,0,1).
                    
            - m_i: unit vectors of measurement basis M (old basis) 
              which are given in cooradinates measured in basis B. 
              For example, a counter-clockwise rotation of 30 deg
              around the 3-axis gives m_1 = (0.866, 0.5, 0.0), 
              m_2 = (-0.5,0.866,0.0), m_3 = (0, 0, 1). 


          The corresponding transformation matrix is given by the dot 
          product indicated by a point '.' ::
          
                   | b_1.m_1  b_1.m_2  b_1.m_3 |   | T_11  T_12  T_13 |
            T_BM = | b_2.m_1  b_2.m_2  b_2.m_3 | = | T_21  T_22  T_23 |
                   | b_3.m_1  b_3.m_2  b_3.m_3 |   | T_31  T_32  T_33 |

          and to transform a vector v measured in basis M denoted as [v]_M  
          into [v]_B ie. vector v measured in B it holds :: 

            [v]_B = T_BM . [v]_M
            
          the back-transform is given by the inverse ::
          
            [v]_M = T_BM^(-1) . [v]_B     

          These relations are written in homogenous coordinates as ::

            | v_1 |     |  T_11  T_12  T_13  t_1  | | v_1 |
            | v_2 |   = |  T_21  T_22  T_23  t_2  | | v_2 |
            | v_3 |     |  T_31  T_32  T_33  t_3  | | v_3 |
            |  1  |_B   |_  0     0     0     1  _| |  1  |_M
            
          In this way translations (tra can be realized. Note that the 
          translation vector is given by the file and not taken from the
          'tra' option. 


-texp   : exponent

          Transformation exponent. Possible choices are 

            - '1'  forward transform    
            - '-1' backward transform
           
          in case of vectors this means (see option 'tfil') ::
           
             [v]_B = T_BM^( 1).[v]_M     (forward)
             [v]_M = T_BM^(-1).[v]_B     (backward)   


-histo  : filename;errorVal;normalize;nIntervals

          Output histogram for the specified "errorValue". 


            - 'filename' where the computed info is written. 
            - 'errorVal' error value. Implemented values are:
            
               - 'ar' : edge aspect ratio of min/max edge length
            
            - 'normalize' option. Can only have the values 'y' or 'n'. 
              If 'y' the percentage of volume within the interval instead 
              of the number of voxel will be given.  
            - 'nIntervals' number of intervals, if not given, 
              the interval size will be "1". 


-help   : Print usage


Info
~~~~

- File:   fec.py
- Author: D. H. Pahr

"""

from sys import argv, exit, stdout
from time import *
from string import split, replace, atof
from dpFem import *
import dpMesh
import os
import struct
import dpTensor
import dpTransform
import numpy
import dpUtils
#import sys

class fec():
    __version__ = 'V_02.11.2016'
    __author__ = 'D. H. Pahr'

    def __init__(self, modelName='default'):
        self.modelName = modelName

    def read(self, inFileName, isBinary=False):
        filename, ext = os.path.splitext(inFileName)
        if ext.upper().find('INP') > -1:
            title, nodes, nsets, elems, elsets = self.readAbaqus(inFileName)
        elif ext.upper().find('RAW') > -1:
            title, nodes, nsets, elems, elsets = self.readRaw(inFileName)
        elif ext.upper().find('OFF') > -1:
            title, nodes, nsets, elems, elsets = self.readGeomview(inFileName)
        elif ext.upper().find('MSH') > -1:
            title, nodes, nsets, elems, elsets = self.readGmsh(inFileName)
        elif ext.upper().find('VTK') > -1:
            title, nodes, nsets, elems, elsets = self.readVTK(inFileName)
        elif ext.upper().find('STL') > -1:
            if isBinary:
                title, nodes, nsets, elems, elsets = self.readSTLbin(inFileName)
            else:
                title, nodes, nsets, elems, elsets = self.readSTL(inFileName)
        elif ext.upper().find('OBJ') > -1:
            title, nodes, nsets, elems, elsets = self.readOBJ(inFileName)
        elif ext.upper().find('FEM') > -1:
            title, nodes, nsets, elems, elsets = self.readFem(inFileName)
        else:
            print ' **ERROR** fec(): intput file extension of file: "%s" not known!\n\n' % inFileName
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return (title, nodes, nsets, elems, elsets)

    def write(self, outFileName, title, nodes, nsets, elems, elsets, NscaResults=None, EscaResults=None, vecResults=None, EvecResults=None, isBinary=False):
        filename, ext = os.path.splitext(outFileName)
        if ext.upper().find('GEOM') > -1:
            self.writeTahoeII(outFileName, title, nodes, nsets, elems, elsets)
        elif ext.upper().find('INP') > -1:
            self.writeAbaqus(outFileName, title, nodes, nsets, elems, elsets, NscaResults=NscaResults)
        elif ext.upper().find('GEO') > -1 or ext.upper().find('CASE') > -1:
            self.writeEnsight(outFileName, title, nodes, nsets, elems, elsets, NscaResults, EscaResults, vecResults, EvecResults)
        elif ext.upper().find('OFF') > -1:
            self.writeOFF(outFileName, title, nodes, elems)
        elif ext.upper().find('STL') > -1:
            if isBinary:
                self.writeSTLbin(outFileName, title, nodes, elems)
            else:
                self.writeSTL(outFileName, title, nodes, elems)
        elif ext.upper().find('MSH') > -1:
            self.writeGmsh(outFileName, title, nodes, elems)
        elif ext.upper().find('SURF') > -1:
            self.writeSurf(outFileName, title, nodes, elems)
        elif ext.upper().find('SMESH') > -1:
            self.writeSmesh(outFileName, title, nodes, elems)
        elif ext.upper().find('VTK') > -1:
            self.writeVTK(outFileName, title, nodes, elems)
        else:
            print ' **ERROR** fec(): output file extension of file: "%s" not known!\n\n' % outFileName
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

    def readAbaqus(self, inFileName, props = False):
        print ' ... read Abaqus file       : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        lines = []
        line = ' '
        while line != '':
            line = inStream.readline()
            if line != '':
                lines.append(line)
            LINE = line.upper()
            if LINE.find('*INCLUDE') == 0:
                inputfile = self.getAbaqusArgument(line, 'INPUT')
                inStream2 = open(inputfile, 'r')
                line2 = ' '
                while line2 != '':
                    line2 = inStream2.readline()
                    if line2 != '':
                        lines.append(line2)
        read = False
        title = None
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        properties = {}
        unknownElems = []
        lineNo = 0
        while lineNo < len(lines):
            lineNo = lineNo + 1
            line = lines[lineNo - 1]
            LINE = line.upper()
            if LINE.find('*HEADING') == 0:
                lineNo = lineNo + 1
                line = lines[lineNo - 1]
                LINE = line.upper()
                title = lines[lineNo - 1]
                title = title.replace('\n', '')
            if LINE.find('*NODE') == 0 and LINE.upper().find('*NODE PRINT') == -1 and LINE.upper().find('*NODE FILE') == -1 and LINE.upper().find('*NODE OUTPUT') == -1:
                nsetName = None
                if LINE.upper().find('NSET') > 0:
                    nsetName = self.getAbaqusArgument(line, 'NSET')
                while lineNo < len(lines):
                    lineNo = lineNo + 1
                    line = lines[lineNo - 1]
                    if line.find('*') == 0 and line.find('**') == -1:
                        lineNo = lineNo - 1
                        break
                    if line.find('**') == 0 or line.find('\n') == 0:
                        pass
                    else:
                        vList = split(line, ',')
                        nNo = int(vList[0])
                        curNode = node(nNo, float(vList[1]))
                        if len(vList) > 2:
                            curNode.set_y(float(vList[2]))
                        if len(vList) > 3:
                            curNode.set_z(float(vList[3]))
                        nodes[nNo] = curNode
                        if nsetName != None:
                            if nsets.has_key(nsetName):
                                nsets[nsetName].append(nNo)
                            else:
                                nsets[nsetName] = [nNo]

                continue
            LINE = line.upper()
            if LINE.find('*NSET') == 0:
                print '  -> found *NSET    at Line %s' % repr(lineNo)
                nsetName = self.getAbaqusArgument(line, 'NSET')
                while lineNo < len(lines):
                    lineNo = lineNo + 1
                    line = lines[lineNo - 1]
                    if line.find('*') == 0 and line.find('**') == -1:
                        lineNo = lineNo - 1
                        break
                    if line.find('**') == 0 or line.find('\n') == 0:
                        pass
                    else:
                        line = line.replace('\n', '')
                        line = line.replace(' ', '')
                        vList = split(line, ',')
                        for Id in vList:
                            if len(Id) > 0:
                                if nsets.has_key(nsetName):
                                    nsets[nsetName].append(int(Id))
                                else:
                                    nsets[nsetName] = [int(Id)]

                continue
            LINE = line.upper()
            if LINE.find('*ELEMENT') == 0 and LINE.upper().find('*ELEMENT OUTPUT') == -1:
                elType = ''
                aElType = self.getAbaqusArgument(line, 'TYPE')
                aElType = aElType.upper()
                nExpNo = 0
                if aElType.find('B32') == 0:
                    elType = 'bar3'
                    noExpNo = 3
                elif aElType.find('B3') == 0 or aElType.find('T3') == 0:
                    elType = 'bar2'
                    noExpNo = 2
                elif aElType.find('CPS3') == 0 or aElType.find('CPE3') == 0 or aElType.find('S3') == 0 or aElType.find('STRI3') == 0:
                    elType = 'tria3'
                    noExpNo = 3
                elif aElType.find('STRI65') == 0:
                    elType = 'tria6'
                    noExpNo = 6
                elif aElType.find('CPS4') == 0 or aElType.find('CPE4') == 0 or aElType.find('S4') == 0:
                    elType = 'quad4'
                    noExpNo = 4
                elif aElType.find('CPS8') == 0 or aElType.find('CPE8') == 0 or aElType.find('S8') == 0:
                    elType = 'quad8'
                    noExpNo = 8
                elif aElType.find('C3D4') == 0:
                    elType = 'tetra4'
                    noExpNo = 4
                elif aElType.find('C3D5') == 0:
                    elType = 'pyra5'
                    noExpNo = 5
                elif aElType.find('C3D8') == 0 or aElType.find('SC8') == 0:
                    elType = 'hexa8'
                    noExpNo = 8
                elif aElType.find('C3D6') == 0 or aElType.find('SC6') == 0:
                    elType = 'penta6'
                    noExpNo = 6
                elif aElType.find('C3D10') == 0:
                    elType = 'tetra10'
                    noExpNo = 10
                elif aElType.find('C3D15') == 0:
                    elType = 'penta15'
                    noExpNo = 15
                elif aElType.find('C3D20') == 0:
                    elType = 'hexa20'
                    noExpNo = 20
                else:
                    if aElType not in unknownElems:
                        unknownElems.append(aElType)
                    continue
                elsetName = ''
                if LINE.find('ELSET') > 0:
                    elsetName = self.getAbaqusArgument(line, 'ELSET')
                while lineNo < len(lines):
                    lineNo += 1
                    line = lines[lineNo - 1]
                    vList = []
                    if line.find('*') == 0 and line.find('**') == -1:
                        lineNo = lineNo - 1
                        break
                    if line.find('**') == 0 or line.find('\n') == 0:
                        pass
                    else:
                        line = line.replace('\n', '')
                        line = line.replace(',', ' ')
                        vList1 = split(line)
                        if len(vList1) - 1 != noExpNo:
                            lineNo += 1
                            line = lines[lineNo - 1]
                            line = line.replace('\n', '')
                            line = line.replace(',', ' ')
                            vList2 = split(line)
                            if len(vList1) + len(vList2) - 1 != noExpNo:
                                lineNo += 1
                                line = lines[lineNo - 1]
                                line = line.replace('\n', '')
                                line = line.replace(',', ' ')
                                vList3 = split(line)
                                if len(vList1) + len(vList2) + len(vList3) - 1 != noExpNo:
                                    stdout.write('\n **ERROR**: fec.readAbaqus(): Line %i ff: Number of nodes for this' % (lineNo - 2))
                                    stdout.write('\n            element and expected nodes to not coincide !\n\n')
                                    stdout.write('\n E N D E D  with ERRORS \n\n')
                                    stdout.flush()
                                    exit(1)
                                else:
                                    vList = vList1 + vList2 + vList3
                            else:
                                vList = vList1 + vList2
                        else:
                            vList = vList1
                        eNo = int(vList[0])
                        nList = []
                        for nNo in range(1, len(vList)):
                            nList.append(int(vList[nNo]))

                        curElem = element(eNo, nList, elType)
                        elems[eNo] = curElem
                        if elsets.has_key(elsetName) > 0:
                            elsets[elsetName].append(eNo)
                        else:
                            elsets[elsetName] = [eNo]

                continue
            if LINE.find('*ELSET') == 0:
                print '\n ** WARNING ** :  *ELSET keyword not supported\n '
            if LINE.find('*BEAM SECTION') == 0:
                elsetName = self.getAbaqusArgument(line, 'ELSET')
                sectName = self.getAbaqusArgument(line, 'SECTION')
                matName = self.getAbaqusArgument(line, 'MATERIAL')
                if sectName.find('CIRC') == 0:
                    lineNo += 1
                    line = lines[lineNo - 1]
                    data = line.split(',')
                    if len(data) == 1:
                        radius = [float(data[0])]
                    else:
                        radius = [float(data[0]), float(data[1])]
                else:
                    stdout.write('\n ** WARNING ** :  *BEAM SECTION, SECTION=%s not implemented\n ' % sectName)
                properties[elsetName] = {'type': 'BEAM',
                 'material': matName,
                 'section': sectName,
                 'geometry': radius}
                continue
            if LINE.find('*SHELL SECTION') == 0:
                elsetName = self.getAbaqusArgument(line, 'ELSET')
                matName = self.getAbaqusArgument(line, 'MATERIAL')
                lineNo += 1
                line = lines[lineNo - 1]
                thickness = float(line)
                properties[elsetName] = {'type': 'SHELL',
                 'material': matName,
                 'thickness': thickness}
                continue
            if LINE.find('*SOLID SECTION') == 0:
                elsetName = self.getAbaqusArgument(line, 'ELSET')
                matName = self.getAbaqusArgument(line, 'MATERIAL')
                properties[elsetName] = {'type': 'SOLID',
                 'material': matName}
                lineNo += 1
                continue

        if len(unknownElems) > 0:
            stdout.write("\n **WARNING**: fec.readAbaqus() Element Types '%s' not implemented!\n" % str(unknownElems))
            stdout.flush()
        if props == True:
            return (title,
             nodes,
             nsets,
             elems,
             elsets,
             properties)
        else:
            return (title,
             nodes,
             nsets,
             elems,
             elsets)

    def readRaw(self, inFileName):
        print ' ... read Raw file           : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        NoVertices, NoElements = split(inStream.readline())
        NoVertices = int(NoVertices)
        NoElements = int(NoElements)
        uniqueNodes = {}
        oldNewNode = {}
        oldprogress = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for nNo in range(NoVertices):
            progress = float(nNo) / float(NoVertices) * 10.0
            if progress > oldprogress:
                dpUtils.progressNext(progress)
            oldprogress = progress
            sCoord = inStream.readline()
            if uniqueNodes.has_key(sCoord) > 0:
                oldNewNode[nNo + 1] = uniqueNodes[sCoord]
            else:
                uniqueNodes[sCoord] = nNo + 1
                oldNewNode[nNo + 1] = nNo + 1

        dpUtils.progressEnd()
        maxNodes = len(uniqueNodes)
        count = 0
        dpUtils.progressStart('     -> Processed Nodes      : ')
        for sCoord in uniqueNodes.keys():
            count += 1
            nNo = uniqueNodes[sCoord]
            if count % 1000 == 0:
                progress = float(count) / float(maxNodes) * 100.0
                dpUtils.progressNext(progress)
            nArray = split(sCoord)
            curNode = node(nNo + 1)
            x = float(nArray[0])
            y = float(nArray[1])
            z = float(nArray[2])
            curNode = node(nNo + 1, x, y, z)
            nodes[nNo] = curNode

        dpUtils.progressEnd()
        dpUtils.progressStart('     -> Processed Elements   : ')
        for eNo in range(NoElements):
            elsetName = 'Part-1'
            if eNo % int(NoElements / 100.0) == 0:
                progress = float(eNo) / float(NoElements) * 10.0
                dpUtils.progressNext(progress)
            eArray = split(inStream.readline())
            nList = []
            for nNo in range(len(eArray)):
                cornerNodeNo = int(eArray[nNo]) + 1
                nList.append(oldNewNode[cornerNodeNo])
                NodesPerElement = len(nList)
                elType = ''
                if NodesPerElement == 8:
                    elType = 'hexa8'
                elif NodesPerElement == 4:
                    stdout.write('\n **ERROR** writeRaw() : Element %s has %s nodes: Not clear if quad4 or tetra4!\n' % (key, NodesPerElement))
                    stdout.flush()
                    stdout.write('           Use elType parameter.\n\n')
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                elif NodesPerElement == 3:
                    elType = 'tria3'
                    noNodes = 3
                elif NodesPerElement == 2:
                    elType = 'bar2'
                    noNodes = 2
                else:
                    stdout.write('\n **ERROR** writeRaw() : Element %s has %s nodes: Not implemented!\n\n' % (key, NodesPerElement))
                    stdout.flush()
                    stdout.write('           Use elType parameter.\n\n')
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

            curElem = element(eNo + 1, nList, elType)
            elems[eNo + 1] = curElem
            if elsets.has_key(elsetName) > 0:
                elsets[elsetName].append(eNo + 1)
            else:
                elsets[elsetName] = [
                 eNo + 1]

        dpUtils.progressEnd()
        return (
         title, nodes, nsets, elems, elsets)

    def readGeomview(self, inFileName):
        print ' ... read Geomview file      : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        allLines = inStream.readlines()
        cleanLines = []
        lineNo = 0
        for line in allLines:
            if line.find('#') == 0 or line.find('\n') == 0:
                continue
                cleanLines.append(line)

        firstLine = cleanLines[lineNo]
        lineNo += 1
        if not (firstLine.find('NOFF') == 0 or firstLine.find('OFF') == 0):
            print '\n **ERROR** readGeomview():'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        NoVertices, NoElements, dummy = split(cleanLines[lineNo])
        lineNo += 1
        NoVertices = int(NoVertices)
        NoElements = int(NoElements)
        oldprogress = 0
        dpUtils.progressStart('     -> Processed Nodes      : ')
        for nNo in range(NoVertices):
            nArray = split(cleanLines[lineNo])
            lineNo += 1
            progress = int(float(nNo + 1) / float(NoVertices) * 10.0)
            if progress > oldprogress:
                dpUtils.progressNext(progress)
            oldprogress = progress
            x = float(nArray[0])
            y = float(nArray[1])
            z = float(nArray[2])
            curNode = node(nNo + 1, x, y, z)
            nodes[nNo + 1] = curNode

        dpUtils.progressEnd()
        dpUtils.progressStart('     -> Processed Elements   : ')
        for eNo in range(NoElements):
            elsetName = 'Part-1'
            if eNo % int(NoElements / 100.0) == 0:
                progress = float(eNo) / float(NoElements) * 10.0
                dpUtils.progressNext(progress)
            eArray = split(cleanLines[lineNo])
            lineNo += 1
            elType = ''
            if int(eArray[0]) == 2:
                elType = 'bar2'
            elif int(eArray[0]) == 3:
                elType = 'tria3'
            else:
                print " **ERROR** readGeomview(): '%s' number of nodes per element not implemented!" % repr(eArray[0])
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            nList = []
            for nNo in range(1, len(eArray)):
                cornerNodeNo = int(eArray[nNo]) + 1
                nList.append(cornerNodeNo)

            curElem = element(eNo + 1, nList, elType)
            elems[eNo + 1] = curElem
            if elsets.has_key(elsetName) > 0:
                elsets[elsetName].append(eNo + 1)
            else:
                elsets[elsetName] = [
                 eNo + 1]

        dpUtils.progressEnd()
        return (
         title, nodes, nsets, elems, elsets)

    def readGmsh(self, inFileName):
        print ' ... read Gmsh file          :', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        lines = inStream.readlines()
        lineNo = 0
        for line in lines:
            if line.upper().find('$MESHFORMAT') == 0:
                version_number, file_type, data_size = split(lines[lineNo + 1])
                if float(version_number) < 2.0:
                    stdout.write('\n **ERROR**: readGmsh() version number < 2.0!')
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                if int(file_type) != 0:
                    stdout.write('\n **ERROR**: readGmsh() file type has to be 0!')
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                break
            lineNo += 1

        if lineNo == len(lines):
            stdout.write('\n **ERROR**: readGmsh() $MESHFORMAT not found!')
            stdout.write('\n            Maybe you use the old gmsh format!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        for line in lines:
            if line.upper().find('$NOD') == 0:
                break
            lineNo += 1

        if lineNo == len(lines):
            stdout.write('\n **ERROR**: readGmsh() $NOD not found!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lineNo += 1
        NoVertices = int(lines[lineNo].replace('\n', ''))
        dpUtils.progressStart('     -> Processed Nodes      : ')
        for nNo in range(NoVertices):
            lineNo += 1
            nArray = split(lines[lineNo])
            if nNo % int(NoVertices / 100.0) == 0:
                progress = float(nNo + 1) / float(NoVertices) * 10.0
                dpUtils.progressNext(progress)
            x = float(nArray[1])
            y = float(nArray[2])
            z = float(nArray[3])
            curNode = node(int(nArray[0]), x, y, z)
            nodes[int(nArray[0])] = curNode

        dpUtils.progressEnd()
        lineNo = 0
        for line in lines:
            if line.upper().find('$ELM') == 0 or line.upper().find('$ELEMENTS') == 0:
                break
            lineNo += 1

        if lineNo == len(lines):
            stdout.write('\n **ERROR**: readGmsh() $NOD not found!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lineNo += 1
        NoElements = int(lines[lineNo].replace('\n', ''))
        dpUtils.progressStart('     -> Processed Elements   : ')
        for eNo in range(NoElements):
            if eNo % int(NoElements / 100.0) == 0:
                progress = float(eNo + 1) / float(NoElements) * 10.0
                dpUtils.progressNext(progress)
            lineNo += 1
            eArray = split(lines[lineNo])
            gElType = int(eArray[1])
            elType = ''
            if gElType == 1:
                elType = 'bar2'
            elif gElType == 2:
                elType = 'tria3'
            elif gElType == 9:
                elType = 'tria6'
            elif gElType == 3:
                elType = 'quad4'
            elif gElType == 6:
                elType = 'penta6'
            elif gElType == 5:
                elType = 'hexa8'
            elif gElType == 4:
                elType = 'tetra4'
            else:
                stdout.write("\n **ERROR** readGmsh() : Element Type '%s' not implemented!\n\n" % repr(gElType))
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            elsetName = 'elset_' + eArray[3]
            noOfNodes = len(eArray) - (3 + int(eArray[2]))
            startNode = 3 + int(eArray[2])
            nList = []
            for nNo in range(noOfNodes):
                cornerNodeNo = int(eArray[nNo + startNode])
                nList.append(cornerNodeNo)

            curElem = element(int(eArray[0]), nList, elType)
            elems[int(eArray[0])] = curElem
            if elsets.has_key(elsetName) > 0:
                elsets[elsetName].append(int(eArray[0]))
            else:
                elsets[elsetName] = [
                 int(eArray[0])]

        dpUtils.progressEnd()
        return (
         title, nodes, nsets, elems, elsets)

    def chunkstring(self, string, length):
        return (string[0 + i:length + i] for i in range(0, len(string), length))

    def userfloat(self, number):
        try:
            return float(number)
        except ValueError:
            if number[1:].find('-') > 1 and number.find('e') < 0 and number.find('E') < 0 and number.find('d') < 0 and number.find('D') < 0:
                number = number[0] + number[1:].replace('-', 'E-')
            if number.find('.') < 0:
                number = '0.' + number
            if number[1:].find('+') > 1 and number.find('e') < 0 and number.find('E') < 0 and number.find('d') < 0 and number.find('D') < 0:
                number = number[0] + number[1:].replace('+', 'E+')
            return float(number)

    def readFem(self, inFileName, props=False):
        print ' ... read Altair Fem file          :', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        properties = {}
        isPBEAM = False
        lines = inStream.readlines()
        lineNo = 0
        for line in lines:
            lineNo += 1
            if line.upper().find('GRID') == 0:
                if line.find(',') > -1:
                    data = line.split(',')
                else:
                    data = list(self.chunkstring(line, 8))
                nid = int(data[1])
                x = self.userfloat(data[3])
                y = self.userfloat(data[4])
                z = self.userfloat(data[5])
                curNode = node(nid, x, y, z)
                nodes[nid] = curNode
            if line.upper().find('CBEAM') == 0:
                data = list(self.chunkstring(line, 8))
                elType = 'bar2'
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nList = [nid1, nid2]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CTRIA3') == 0:
                data = list(self.chunkstring(line, 8))
                elType = 'tria3'
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nList = [nid1, nid2, nid3]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CTRIA6') == 0:
                data = list(self.chunkstring(line, 8))
                elType = 'tria6'
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nid5 = int(data[7])
                nid6 = int(data[8])
                nList = [nid1, nid2, nid3, nid4, nid5, nid6]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CQUAD4') == 0:
                elType = 'quad4'
                data = list(self.chunkstring(line, 8))
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nList = [nid1, nid2, nid3, nid4]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CQUAD8') == 0:
                elType = 'quad8'
                data = list(self.chunkstring(line, 8))
                data2 = list(self.chunkstring(lines[lineNo], 8))
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nid5 = int(data[7])
                nid6 = int(data[8])
                nid7 = int(data2[1])
                nid8 = int(data2[2])
                nList = [
                 nid1, nid2, nid3, nid4, nid5, nid6, nid7, nid8]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CTETRA') == 0:
                elType = 'tetra4'
                data = list(self.chunkstring(line, 8))
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nList = [nid1, nid2, nid3, nid4]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CHEXA') == 0:
                elType = 'hexa8'
                data = list(self.chunkstring(line, 8))
                data2 = list(self.chunkstring(lines[lineNo], 8))
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nid5 = int(data[7])
                nid6 = int(data[8])
                nid7 = int(data2[1])
                nid8 = int(data2[2])
                nList = [
                 nid1, nid2, nid3, nid4, nid5, nid6, nid7, nid8]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CPENTA') == 0:
                data = list(self.chunkstring(line, 8))
                elType = 'penta6'
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nid5 = int(data[7])
                nid6 = int(data[8])
                nList = [nid1, nid2, nid3, nid4, nid5, nid6]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('CPYRA') == 0:
                data = list(self.chunkstring(line, 8))
                elType = 'pyra5'
                elid = int(data[1])
                pid = int(data[2])
                nid1 = int(data[3])
                nid2 = int(data[4])
                nid3 = int(data[5])
                nid4 = int(data[6])
                nid5 = int(data[7])
                nList = [nid1, nid2, nid3, nid4, nid5]
                curElem = element(elid, nList, elType)
                elems[elid] = curElem
                elsetName = str(pid)
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]
            if line.upper().find('PBEAML') == 0:
                data = list(self.chunkstring(line.strip(), 8))
                data2 = list(self.chunkstring(lines[lineNo].strip(), 8))
                pid = int(data[1])
                mid = int(data[2])
                sectName = ''
                if data[4].strip().upper().find('ROD') == 0 and numpy.allclose(self.userfloat(data2[4]), 1.0):
                    sectName = 'CIRC'
                    radius = [self.userfloat(data2[1]), self.userfloat(data2[5])]
                elif data[4].strip().upper().find('I1') == 0:
                    sectName = 'CIRC'
                    radius = [self.userfloat(data2[1]), self.userfloat(data2[1])]
                else:
                    stdout.write("\n **ERROR** readFem() : PBEAML Type '%s' not implemented!\n\n" % repr(data[4]))
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                elsetName = str(pid)
                properties[elsetName] = {'type': 'BEAM','material': mid,'section': sectName,'geometry': radius}
            if line.upper().find('PBEAM') == 0 and line.upper().find('PBEAML') < 0 and line.upper().find('PBEAMX') < 0:
                isPBEAM = True
                data = list(self.chunkstring(line.strip(), 8))
                pid = int(data[1])
                mid = int(data[2])
                area = self.userfloat(data[3])
                sectName = 'CIRC'
                radius = [numpy.sqrt(4 * area / numpy.pi)]
                elsetName = str(pid)
                properties[elsetName] = {'type': 'BEAM','material': mid,'section': sectName,'geometry': radius}
            if line.upper().find('PSHELL') == 0 and line.upper().find('PSHELLX') < 0:
                data = list(self.chunkstring(line.strip(), 8))
                pid = int(data[1])
                mid = int(data[2])
                thickness = float(data[3])
                elsetName = str(pid)
                properties[elsetName] = {'type': 'SOLID','material': mid,'thickness': thickness}
            if line.upper().find('PSOLID') == 0 and line.upper().find('PSOLIDX') < 0:
                data = list(self.chunkstring(line.strip(), 8))
                pid = int(data[1])
                mid = int(data[2])
                elsetName = str(pid)
                properties[elsetName] = {'type': 'SOLID','material': mid}

        if isPBEAM:
            stdout.write('\n **WARNING** readFem() : PBEAM Type is interpreted as circluar section!\n\n')
        for elsetName in elsets.keys():
            if elsetName not in properties:
                stdout.write(" **WARNING** readFem() : PID '%s' used for elements but no property card found!\n" % elsetName)
                stdout.write('                         Concerned elements are ignored.\n')
                del elsets[elsetName]

        stdout.write('\n')
        if props == True:
            return (title, nodes, nsets, elems, elsets, properties)
        else:
            return (
             title, nodes, nsets, elems, elsets)

    def readVTKTrias(self, inFileName):
        print ' ... read vtk file          : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        line1 = inStream.readline()
        line2 = inStream.readline()
        line3 = inStream.readline()
        line4 = inStream.readline()
        if not line3.find('ASCII') == 0:
            print '\n **ERROR** readVtk(): ASCII Format required'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if not line4.find('DATASET UNSTRUCTURED_GRID') == 0:
            print '\n **ERROR** readVtk(): UNSTRUCTURED_GRID format required'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        while 1:
            line = inStream.readline()
            if not line:
                break
            nArray = split(line)
            if len(nArray) > 0:
                CellType = nArray[0]
            else:
                CellType = None
            if CellType == 'POINTS':
                NoVertices = int(nArray[1])
                ptNo = 1
                oldprogress = 0
                dpUtils.progressStart('     -> Processed Data      : ')
                while ptNo <= NoVertices:
                    nArray = split(inStream.readline())
                    progress = int(float(ptNo) / float(NoVertices) * 10.0)
                    if progress > oldprogress:
                        dpUtils.progressNext(progress)
                    oldprogress = progress
                    ptPerLine = len(nArray) / 3
                    for nNo in range(ptPerLine):
                        x = float(nArray[0 + nNo * 3])
                        y = float(nArray[1 + nNo * 3])
                        z = float(nArray[2 + nNo * 3])
                        curNode = node(ptNo, x, y, z)
                        nodes[ptNo] = curNode
                        ptNo += 1

                dpUtils.progressEnd()
            if CellType == 'CELLS':
                NoElements = int(nArray[1])
                dpUtils.progressStart('     -> Processed Elements  : ')
                for eNo in range(NoElements):
                    elsetName = 'Part-1'
                    if eNo % int(NoElements / 100.0) == 0:
                        progress = float(eNo) / float(NoElements) * 10.0
                        dpUtils.progressNext(progress)
                    eArray = split(inStream.readline())
                    elType = ''
                    if int(eArray[0]) == 3:
                        elType = 'tria3'
                    else:
                        print " **ERROR** readVTK(): '%s' number of nodes per element not implemented!" % repr(eArray[0])
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                    nList = []
                    for nNo in range(1, len(eArray)):
                        cornerNodeNo = int(eArray[nNo]) + 1
                        nList.append(cornerNodeNo)

                    curElem = element(eNo + 1, nList, elType)
                    elems[eNo + 1] = curElem
                    if elsets.has_key(elsetName) > 0:
                        elsets[elsetName].append(eNo + 1)
                    else:
                        elsets[elsetName] = [
                         eNo + 1]

                dpUtils.progressEnd()

        return (title, nodes, nsets, elems, elsets)

    def auxGetVtkDataLine(self, fileLineList, keyword):
        matches = [ s for s in fileLineList if keyword in str(s) ]
        index = fileLineList.index(matches[0])
        data = matches[0].split()
        nEntries = int(data[1])
        return (
         index, nEntries)

    def readVTK(self, inFileName, props=False):
        print ' ... read vtk file          : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        properties = {}
        elTypes = {}
        elTypes[3] = 'bar2'
        elTypes[21] = 'bar3'
        elTypes[5] = 'tria3'
        elTypes[22] = 'tria6'
        elTypes[9] = 'quad4'
        elTypes[23] = 'quad8'
        elTypes[10] = 'tetra4'
        elTypes[24] = 'tetra10'
        elTypes[12] = 'hexa8'
        elTypes[25] = 'hexa20'
        elTypes[14] = 'pyra5'
        elTypes[13] = 'penta6'
        elTypes[26] = 'penta15'
        fileLineList = inStream.readlines()
        inStream.close()
        matching = [ s for s in fileLineList if 'ASCII' in s ]
        if not len(matching) == 1:
            dpUtils.throwError('fec.readVtk(): ASCII Format required!')
        matching = [ s for s in fileLineList if 'DATASET UNSTRUCTURED_GRID' in s ]
        if not len(matching) == 1:
            dpUtils.throwError('fec.readVtk(): UNSTRUCTURED_GRID format required!')
        lineNo, NoVertices = self.auxGetVtkDataLine(fileLineList, 'POINTS')
        ptNo = 1
        oldprogress = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        while ptNo <= NoVertices:
            lineNo += 1
            nArray = split(fileLineList[lineNo])
            progress = int(float(ptNo) / float(NoVertices) * 10.0)
            if progress > oldprogress:
                dpUtils.progressNext(progress)
            oldprogress = progress
            ptPerLine = len(nArray) / 3
            for nNo in range(ptPerLine):
                x = float(nArray[0 + nNo * 3])
                y = float(nArray[1 + nNo * 3])
                z = float(nArray[2 + nNo * 3])
                curNode = node(ptNo, x, y, z)
                nodes[ptNo] = curNode
                ptNo += 1

        dpUtils.progressEnd()
        lineNoCells, NoElements = self.auxGetVtkDataLine(fileLineList, 'CELLS')
        lineNoCellType, NoCellTypes = self.auxGetVtkDataLine(fileLineList, 'CELL_TYPE')
        if NoElements != NoCellTypes:
            dpUtils.throwError('fec.readVtk(): Number of CELLS and CELL_TYPE not equivalent!')
        dpUtils.progressStart('     -> Processed Elements  : ')
        oldprogress = 0
        notImplementedCellTypes = []
        for eNo in range(NoElements):
            elsetName = 'Part-1'
            progress = float(eNo) / float(NoElements) * 10.0
            if progress > oldprogress:
                dpUtils.progressNext(progress)
            oldprogress = progress
            lineNoCells += 1
            eArray = split(fileLineList[lineNoCells])
            lineNoCellType += 1
            cellType = int(fileLineList[lineNoCellType])
            if cellType in elTypes:
                nList = []
                for nNo in range(1, len(eArray)):
                    cornerNodeNo = int(eArray[nNo]) + 1
                    nList.append(cornerNodeNo)

                curElem = element(eNo + 1, nList, elTypes[cellType])
                elems[eNo + 1] = curElem
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(eNo + 1)
                else:
                    elsets[elsetName] = [
                     eNo + 1]
            elif cellType not in notImplementedCellTypes:
                notImplementedCellTypes.append(cellType)

        dpUtils.progressEnd()
        if len(notImplementedCellTypes) > 0:
            dpUtils.throwWarning('fec.readVtk() : CELL_TYPEs=%s not implemented and not read!\n' % repr(notImplementedCellTypes))
        properties[elsetName] = {'type': 'VTK','material': 'noMaterial','geometry': 1.0}
        if props == True:
            return (title, nodes, nsets, elems, elsets, properties)
        else:
            return (
             title, nodes, nsets, elems, elsets)

    def readTfm(inFileName):
        print ' ... read tfm file          : ', inFileName
        try:
            f = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**:  fec.readTfm - file '%s' not found!\n\n" % inFileName)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        MR = None
        vo = None
        vc = None
        typ = None
        for line in f:
            if line.find('Parameters') == 0:
                name, data = line.split(':')
                sld = data.split()
                MR = numpy.array([[float(sld[0]), float(sld[1]), float(sld[2])],
                 [
                  float(sld[3]), float(sld[4]), float(sld[5])],
                 [
                  float(sld[6]), float(sld[7]), float(sld[8])]])
                vt = numpy.array([float(sld[9]), float(sld[10]), float(sld[11])])
            elif line.find('FixedParameters') == 0:
                name, data = line.split(':')
                sld = data.split()
                vc = numpy.array([float(sld[0]), float(sld[1]), float(sld[2])])
            elif line.find('Transform') == 0:
                name, data = line.split(':')
                typ = name.replace(' ', '')

        if MR == None:
            stdout.write("\n **ERROR**:  fec.readTfm - 'Parameters' keyword not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if vc == None:
            stdout.write("\n **ERROR**:  fec.readTfm - 'FixedParameters' keyword not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if typ == None:
            stdout.write("\n **ERROR**:  fec.readTfm - 'Transform' keyword not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if typ.fine('AffineTransform_double_3_3'):
            stdout.write("\n **ERROR**:  fec.readTfm - Transform '%s' not implemented!\n" % typ)
            stdout.write("             Use 'AffineTransform_double_3_3' transform.\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return (
         MR, vt, vc)

    def readSTL(self, inFileName):
        print ' ... read stl file          : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        line1 = inStream.readline()
        header = line1.split()
        noid = 0
        elid = 0
        elsetName = 'elset_' + header[1]
        while 1:
            line = inStream.readline()
            if line.upper().find('ENDSOLID') >= 0:
                break
            if line.upper().find('VERTEX') >= 0:
                dummy, posx1, posy1, posz1 = line.split()
                noid += 1
                curNode = node(noid, float(posx1), float(posy1), float(posz1))
                nodes[noid] = curNode
                line = inStream.readline()
                dummy, posx2, posy2, posz2 = line.split()
                noid += 1
                curNode = node(noid, float(posx2), float(posy2), float(posz2))
                nodes[noid] = curNode
                line = inStream.readline()
                dummy, posx3, posy3, posz3 = line.split()
                noid += 1
                curNode = node(noid, float(posx3), float(posy3), float(posz3))
                nodes[noid] = curNode
                elid += 1
                curElem = element(elid, [noid - 2, noid - 1, noid], 'tria3')
                elems[elid] = curElem
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]

        myTriaMesh = dpMesh.dpMesh(nodes=nodes, elems=elems)
        minDist = myTriaMesh.getMinNodalDistance()
        myTriaMesh.equivalenceNodes(tolerance=0.01 * minDist, echo=True)
        elems = myTriaMesh.getElems()
        nodes = myTriaMesh.getNodes()
        return (
         title, nodes, nsets, elems, elsets)

    def readSTLbin(self, inFileName):
        print ' ... read stl file binary   : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'rb')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        header = inStream.read(80)
        title = struct.unpack('80s', header)[0].strip()
        nn = inStream.read(4)
        nElems = struct.unpack('i', nn)[0]
        record_dtype = numpy.dtype([
         (
          'normals', numpy.float32, (3, )),
         (
          'Vertex1', numpy.float32, (3, )),
         (
          'Vertex2', numpy.float32, (3, )),
         (
          'Vertex3', numpy.float32, (3, )),
         (
          'atttr', '<i2', (1, ))])
        data = numpy.fromfile(inStream, dtype=record_dtype, count=nElems)
        inStream.close()
        normals = data['normals']
        vertex1 = data['Vertex1']
        vertex2 = data['Vertex2']
        vertex3 = data['Vertex3']
        noid = 0
        elid = 0
        elsetName = 'elset_tria'
        for i in range(len(normals)):
            noid += 1
            curNode = node(noid, vertex1[i][0], vertex1[i][1], vertex1[i][2])
            nodes[noid] = curNode
            noid += 1
            curNode = node(noid, vertex2[i][0], vertex2[i][1], vertex2[i][2])
            nodes[noid] = curNode
            noid += 1
            curNode = node(noid, vertex3[i][0], vertex3[i][1], vertex3[i][2])
            nodes[noid] = curNode
            elid += 1
            curElem = element(elid, [noid - 2, noid - 1, noid], 'tria3')
            elems[elid] = curElem
            if elsets.has_key(elsetName) > 0:
                elsets[elsetName].append(elid)
            else:
                elsets[elsetName] = [
                 elid]

        myTriaMesh = dpMesh.dpMesh(nodes=nodes, elems=elems)
        minDist = myTriaMesh.getMinNodalDistance()
        myTriaMesh.equivalenceNodes(tolerance=0.01 * minDist, echo=True)
        elems = myTriaMesh.getElems()
        nodes = myTriaMesh.getNodes()
        return (
         title, nodes, nsets, elems, elsets)

    def readOBJ(self, inFileName):
        print ' ... read obj file          : ', inFileName
        stdout.flush()
        try:
            inStream = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        title = inFileName
        nodes = {}
        elems = {}
        nsets = {}
        elsets = {}
        noid = 0
        elid = 0
        elsetName = 'elset_surface'
        while 1:
            line = inStream.readline()
            if not line:
                break
            if line.upper().find('V') == 0:
                dummy, posx1, posy1, posz1 = line.split()
                noid += 1
                curNode = node(noid, float(posx1), float(posy1), float(posz1))
                nodes[noid] = curNode
            if line.upper().find('F') == 0:
                dummy, noid1, noid2, noid3 = line.split()
                elid += 1
                curElem = element(elid, [int(noid1), int(noid2), int(noid3)], 'tria3')
                elems[elid] = curElem
                if elsets.has_key(elsetName) > 0:
                    elsets[elsetName].append(elid)
                else:
                    elsets[elsetName] = [
                     elid]

        return (
         title, nodes, nsets, elems, elsets)

    def writeTahoeII(self, outFileName, title, nodes, nsets, elems, elsets):
        print ' ... write TahoeII file      : ', outFileName
        stdout.flush()
        keys = nodes.keys()
        nkey1 = keys[1]
        del keys
        os = open(outFileName, 'w')
        os.write('*version\n')
        os.write('1.0\n')
        os.write('*title\n')
        os.write('%s\n' % title)
        os.write('############################################################\n')
        os.write('*dimensions\n')
        os.write('%s   \t\t# number of nodes\n' % repr(len(nodes)))
        noSpatDim = nodes[nkey1].get_dimension()
        os.write('%s   \t\t# number of spatial dimensions\n' % repr(3))
        os.write('# ... element sets \n')
        os.write('%s   \t\t# number of element sets\n' % repr(len(elsets)))
        os.write('# [ID] [nEl] [nNo/El] \n')
        elsetId = 10 ** len(repr(len(elsets)))
        for setName in elsets:
            if len(elsets) > 0:
                elsetId = elsetId + 1
                nEl = len(elsets[setName])
                nNoEl = len(elsets[setName][0].get_nodes())
                os.write('%s %s %s   \t# set name : %s\n' % (repr(elsetId), repr(nEl), repr(nNoEl), setName))

        os.write('# ... node sets \n')
        os.write('%s   \t\t# number of node sets\n' % repr(len(nsets)))
        if len(nsets) > 0:
            os.write('# [ID]  [nNodes]\n')
            nsetId = 10 ** len(repr(len(nsets))) + 10 ** len(repr(len(elsets)))
            for setName in nsets:
                nsetId = nsetId + 1
                nNo = len(nsets[setName])
                os.write('%s %s   \t# set name : %s\n' % (repr(nsetId), repr(nNo), setName))

        os.write('# ... side sets \n')
        os.write('%s   \t\t# number of side sets\n' % '0')
        os.write('############################################################\n')
        os.write('*nodesets\n')
        if len(nsets) > 0:
            nsetId = 10 ** len(repr(len(nsets))) + 10 ** len(repr(len(elsets)))
            for setName in nsets:
                os.write('*set\n')
                nsetId = nsetId + 1
                nNo = len(nsets[setName])
                os.write('%s   \t# number of nodes in set : %s = %s\n' % (repr(nNo), repr(nsetId), setName))
                count = 0
                for nodeId in nsets[setName]:
                    count += 1
                    os.write('%s ' % nodeId)
                    if count > 9:
                        os.write('\n')
                        count = 0

                os.write('\n')

        os.write('############################################################\n')
        os.write('*sidesets\n')
        os.write('############################################################\n')
        os.write('*elements\n')
        if len(elsets) > 0:
            elsetId = 10 ** len(repr(len(elsets)))
            for setName in elsets:
                os.write('*set\n')
                elsetId = elsetId + 1
                elNo = len(elsets[setName])
                os.write('%s   \t# number of elements in set : %s = %s\n' % (repr(elNo), repr(elsetId), setName))
                count = 0
                for elId in elsets[setName]:
                    count += 1
                    if count == 1:
                        nNoEl = len(elems[elId].get_nodes())
                        os.write('%s   \t# number of element nodes\n' % repr(nNoEl))
                    os.write('%s ' % repr(elId))
                    for node in elems[elId].get_nodes():
                        os.write('%s ' % repr(node))

                    os.write('\n')

        os.write('############################################################\n')
        os.write('*nodes\n')
        os.write('%s   \t# number of nodes\n' % repr(len(nodes)))
        os.write('%s   \t# number of spatial dimensions\n' % repr(noSpatDim))
        for nodeId in nodes:
            os.write('%s ' % repr(nodeId))
            os.write('%13.7e ' % nodes[nodeId].get_x())
            if noSpatDim > 1:
                os.write('%13.7e ' % nodes[nodeId].get_y())
            if noSpatDim > 2:
                os.write('%13.7e ' % nodes[nodeId].get_z())
            os.write('\n')

        os.write('############################################################\n')

    def writeEnsight(self, outFileName, title, nodes, nsets, elems, elsets, NscaResults=None, EscaResults=None, vecResults=None, EvecResults=None):
        AllEResults = {}
        if elsets and len(elsets) > 0:
            for elset in elsets:
                AllEResults[elems[elsets[elset][0]].get_type()] = {}

            matid = {}
            i = 1
            for setName in elsets:
                matid[setName] = i
                i += 1

            for setName in elsets:
                for elid in elsets[setName]:
                    AllEResults[elems[elid].get_type()][elid] = float(matid[setName])

            for elid in AllEResults:
                if EscaResults:
                    EscaResults.append(AllEResults[elid])
                else:
                    EscaResults = [
                     AllEResults[elid]]

        fileList = outFileName.split('/')
        ffilename = fileList.pop()
        filename, ext = self.getFilenameAndExtension(ffilename)
        pathFilename, ext = self.getFilenameAndExtension(outFileName)
        cfilename = filename + '.case'
        print ' ... write Ensight case file :', cfilename
        caseOS = open(pathFilename + '.case', 'w')
        caseOS.write('FORMAT\n')
        caseOS.write('type: ensight gold\n\n')
        caseOS.write('GEOMETRY\n')
        caseOS.write('model:                           %s\n\n' % (filename + '.geo'))
        if NscaResults != None or EscaResults != None or vecResults != None or EvecResults != None:
            caseOS.write('VARIABLE\n')
            if vecResults != None:
                vres = 0
                for vecRes in vecResults:
                    vres += 1
                    caseOS.write('vector per node:  %10i Vector  %s\n' % (1, filename + '.vec' + repr(vres)))

            if EvecResults != None:
                vres = 0
                for vecRes in EvecResults:
                    vres += 1
                    caseOS.write('vector per element:  %10i Vector  %s\n' % (1, filename + '.Evec' + repr(vres)))

            if NscaResults != None:
                sres = 0
                for scaRes in NscaResults:
                    sres += 1
                    caseOS.write('scalar per node:  %10i Scalar  %s\n' % (1, filename + '.sca' + repr(sres)))

            if EscaResults != None:
                sres = 0
                for scaRes in EscaResults:
                    sres += 1
                    caseOS.write('scalar per element:  ESca%s%s  %s\n' % (repr(sres), filename, filename + '.Esca' + repr(sres)))

        caseOS.close()
        gfilename = filename + '.geo'
        print ' ... write Ensight Data File :', gfilename
        geoOS = open(pathFilename + '.geo', 'w')
        geoOS.write('Title: %s\n' % title)
        geoOS.write('Description 2\n')
        geoOS.write('node id given\n')
        geoOS.write('element id given\n')
        geoOS.write('part \n')
        geoOS.write('         1\n')
        geoOS.write('Description PART \n')
        geoOS.write('coordinates \n')
        geoOS.write('%10i\n' % len(nodes))
        nnode = 0
        ensightNodeDict = {}
        for key in nodes.keys():
            nnode += 1
            geoOS.write('%10i\n' % key)
            ensightNodeDict[key] = nnode

        for key in nodes.keys():
            geoOS.write('%12.5e\n' % nodes[key].get_x())

        for key in nodes.keys():
            geoOS.write('%12.5e\n' % nodes[key].get_y())

        for key in nodes.keys():
            geoOS.write('%12.5e\n' % nodes[key].get_z())

        geoOS.write('\n')
        elTypes = {}
        for key in elems.keys():
            elType = elems[key].get_type()
            if not elTypes.has_key(elType):
                elTypes[elType] = 1
            else:
                elTypes[elType] += 1

        for elType in elTypes:
            geoOS.write('%s\n' % elType)
            geoOS.write('%10i\n' % elTypes[elType])
            for key in elems.keys():
                if elems[key].get_type() == elType:
                    geoOS.write('%10i\n' % key)

            for key in elems.keys():
                if elems[key].get_type() == elType:
                    nodeList = elems[key].get_nodes()
                    for noid in nodeList:
                        geoOS.write('%10i' % ensightNodeDict[noid])

                    geoOS.write('\n')

        geoOS.write('\n')
        geoOS.write('\n')
        geoOS.close()
        if NscaResults != None or vecResults != None or EscaResults != None or EvecResults != None:
            if vecResults != None:
                vres = 0
                for vecRes in vecResults:
                    vres += 1
                    filename2 = pathFilename + '.vec' + repr(vres)
                    disOS = open(filename2, 'w')
                    disOS.write('Disp\n')
                    disOS.write('part \n')
                    disOS.write('         1\n')
                    disOS.write('coordinates\n')
                    for dir in range(3):
                        for resId in vecRes:
                            disOS.write('%12.5e\n' % vecRes[resId][dir])

                    disOS.close()

            if EvecResults != None:
                vres = 0
                for vecRes in EvecResults:
                    vres += 1
                    filename2 = pathFilename + '.Evec' + repr(vres)
                    disOS = open(filename2, 'w')
                    disOS.write('Disp\n')
                    disOS.write('part \n')
                    disOS.write('         1\n')
                    disOS.write('%s\n' % elType)
                    for dir in range(3):
                        for resId in vecRes:
                            disOS.write('%12.5e\n' % vecRes[resId][dir])

                    disOS.close()

            if NscaResults != None:
                sres = 0
                for scaRes in NscaResults:
                    sres += 1
                    filename2 = pathFilename + '.sca' + repr(sres)
                    strOS = open(filename2, 'w')
                    strOS.write('%s \n' % ('Scalar' + repr(sres)))
                    strOS.write('part \n')
                    strOS.write('         1\n')
                    strOS.write('coordinates\n')
                    for resId in scaRes:
                        strOS.write('%12.5e\n' % scaRes[resId])

                    strOS.close()

            if EscaResults != None:
                sres = 0
                for scaRes in EscaResults:
                    sres += 1
                    filename2 = pathFilename + '.Esca' + repr(sres)
                    strOS = open(filename2, 'w')
                    strOS.write('%s \n' % ('ESca' + repr(sres) + filename))
                    strOS.write('part \n')
                    strOS.write('         1\n')
                    for resId in scaRes:
                        elType = elems[resId].get_type()
                        break

                    strOS.write('%s\n' % elType)
                    for resId in scaRes:
                        strOS.write('%12.5e\n' % scaRes[resId])

                    strOS.close()

        return

    def writeAbaqus(self, outFileName, title, nodes, nsets, elems, elsets, NscaResults=None):
        time1 = clock()
        print ' ... write ABAQUS file       : ', outFileName
        stdout.flush()
        keys = nodes.keys()
        nkey1 = keys[1]
        del keys
        noSpatDim = nodes[nkey1].get_dimension()
        os = open(outFileName, 'w')
        if not title == None:
            os.write('*HEADING\n')
            os.write('%s\n' % title)
        os.write('***********************************************************\n')
        os.write('*NODE\n')
        for nodeId in nodes:
            os.write('%s, ' % repr(nodeId))
            os.write('%13.7e, ' % nodes[nodeId].get_x())
            if noSpatDim > 1:
                os.write('%13.7e, ' % nodes[nodeId].get_y())
            else:
                os.write(', ')
            if noSpatDim > 2:
                os.write('%13.7e ' % nodes[nodeId].get_z())
            os.write('\n')

        if NscaResults != None:
            os.write('***********************************************************\n')
            os.write('*NODAL THICKNESS\n')
            nodeThick = NscaResults[0]
            for nodeId in nodes:
                os.write('%s, ' % repr(nodeId))
                os.write('%13.7e\n' % nodeThick[nodeId])

        os.write('***********************************************************\n')
        if nsets != None:
            if len(nsets) > 0:
                for setName in nsets:
                    if setName != '':
                        os.write('*NSET, NSET=%s\n' % setName)
                        count = 0
                        for nodeId in nsets[setName]:
                            count += 1
                            if count == 16:
                                os.write('%s' % nodeId)
                                os.write('\n')
                                count = 0
                            else:
                                os.write('%s, ' % nodeId)

                        if count != 0:
                            os.write('\n')

        else:
            os.write('** no NSET written\n')
        os.write('***********************************************************\n')
        if elsets != None:
            if len(elsets) > 0:
                for setName in elsets:
                    aElType = ''
                    elType = elems[elsets[setName][0]].get_type()
                    if elType == 'bar2':
                        aElType = 'B3'
                    elif elType == 'tria3':
                        aElType = 'S3'
                    elif elType == 'quad4':
                        aElType = 'S4'
                    elif elType == 'penta6':
                        aElType = 'C3D6'
                    elif elType == 'hexa8':
                        aElType = 'C3D8'
                    elif elType == 'tetra4':
                        aElType = 'C3D4'
                    elif elType == 'pyra5':
                        aElType = 'C3D5'
                    elif elType == 'bar3':
                        aElType = 'B32'
                    elif elType == 'tria6':
                        aElType = 'STRI65'
                    elif elType == 'quad8':
                        aElType = 'S8'
                    elif elType == 'penta15':
                        aElType = 'C3D15'
                    elif elType == 'hexa20':
                        aElType = 'C3D20'
                    elif elType == 'tetra10':
                        aElType = 'C3D10'
                    else:
                        stdout.write("\n **ERROR** writeAbaqus() : Element Type '%s' not implemented!\n\n" % repr(elType))
                        stdout.flush()
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                    os.write('*ELEMENT, TYPE=%s, ELSET=%s\n' % (aElType, setName))
                    for elId in elsets[setName]:
                        os.write('%s' % repr(elId))
                        count = 1
                        for node in elems[elId].get_nodes():
                            count += 1
                            if count == 8:
                                os.write(', %s,\n' % repr(node))
                                count = 0
                            elif count == 1:
                                os.write('%s' % repr(node))
                            else:
                                os.write(', %s' % repr(node))

                        if count != 0:
                            os.write('\n')

        else:
            os.write('** no ELEMENTS and ELSET written\n')
        time2 = clock()
        print '     -> write finished in    :   %8.1f sec' % (time2 - time1)
        os.close
        return

    def getAbaqusArgument(self, string, argName):
        argument = ''
        string = string.replace('\n', '')
        STRING = string.upper()
        ARGNAME = argName.upper()
        if STRING.find(ARGNAME) > 0:
            string1 = split(string, ',')
            for stringpart in string1:
                stringpart = stringpart.replace(' ', '')
                if stringpart.upper().find(ARGNAME) == 0:
                    command, argument = split(stringpart, '=')

        else:
            print " **ERROR** getAbaqusArgument(): Argument '%s' not found!" % argName
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            stdout.flush()
            exit(1)
        return argument

    def writeOFF(self, outFileName, title, nodes, elems):
        print ' ... write OFF file         : ', outFileName
        stdout.flush()
        os = open(outFileName, 'w')
        os.write('OFF\n')
        os.write('%s %s 0\n' % (repr(len(nodes)), repr(len(elems))))
        nodeIdTable = {}
        count = 0
        for nodeId in nodes:
            nodeIdTable[nodeId] = count
            count += 1
            os.write('%13.7e %13.7e %13.7e\n' % (nodes[nodeId].get_x(), nodes[nodeId].get_y(), nodes[nodeId].get_z()))

        for elId in elems:
            nodes = elems[elId].get_nodes()
            noNodes = len(nodes)
            if noNodes < 2 or noNodes > 3:
                print " **ERROR** writeOFF(): '%s' number of nodes per element not implemented!" % repr(noNodes)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            os.write('%13i ' % noNodes)
            for node in nodes:
                os.write('%13i ' % nodeIdTable[node])

            os.write('\n')

    def writeVTK(self, outFileName, title, nodes, elems):
        print ' ... write VTK file         : ', outFileName
        stdout.flush()
        os = open(outFileName, 'w')
        os.write('# vtk DataFile Version 2.0\n')
        os.write('%s\n' % title)
        os.write('ASCII\n')
        os.write('DATASET UNSTRUCTURED_GRID\n')
        nodeIdTable = {}
        count = 0
        os.write('POINTS %i float\n' % len(nodes))
        for nodeId in nodes:
            nodeIdTable[nodeId] = count
            count += 1
            x, y, z = nodes[nodeId].get_coord()
            os.write('%13.7e %13.7e %13.7e\n' % (x, y, z))

        cell_types = []
        elem_active = {}
        noNodes = 0
        for elId in elems:
            foundType = True
            aElType = None
            elType = elems[elId].get_type()
            if elType == 'bar2':
                aElType = 3
                noNodes += 2
            elif elType == 'tria3':
                aElType = 5
                noNodes += 3
            elif elType == 'quad4':
                aElType = 8
                noNodes += 4
            elif elType == 'pyra5':
                aElType = 14
                noNodes += 5
            elif elType == 'penta6':
                aElType = 13
                noNodes += 6
            elif elType == 'hexa8':
                aElType = 12
                noNodes += 8
            elif elType == 'tetra4':
                aElType = 10
                noNodes += 4
            elif elType == 'bar3':
                aElType = 21
                noNodes += 3
            elif elType == 'tria6':
                aElType = 22
                noNodes += 6
            elif elType == 'quad8':
                aElType = 23
                noNodes += 8
            elif elType == 'hexa20':
                aElType = 25
                noNodes += 20
            elif elType == 'tetra10':
                aElType = 24
                noNodes += 10
            elif elType == 'penta15':
                aElType = 26
                noNodes += 15
            else:
                stdout.write(" **WARNING** writeVTK(): Element type '%s' is not implemented!\n\n" % repr(elType))
                foundType = False
            if foundType:
                cell_types.append(aElType)
                elem_active[elId] = elems[elId]

        noElems = len(elem_active)
        os.write('\nCELLS %i %i\n' % (noElems, noNodes + noElems))
        for elId in elem_active:
            nodes = elem_active[elId].get_nodes()
            if elems[elId].get_type() == 'quad4':
                n3 = nodes[2]
                n4 = nodes[3]
                nodes[2] = n4
                nodes[3] = n3
            noNodes = len(nodes)
            os.write('%i ' % noNodes)
            for node in nodes:
                os.write('%i ' % nodeIdTable[node])

            os.write('\n')

        os.write('\nCELL_TYPES %i\n' % len(cell_types))
        for cell_type in cell_types:
            os.write('%i\n' % cell_type)

        os.close()
        return

    def writeSTL(self, outFileName, title, nodes, elems):
        print ' ... write STL file         : ', outFileName
        stdout.flush()
        os = open(outFileName, 'w')
        os.write('solid ascii\n')
        for elId in elems:
            elNodes = elems[elId].get_nodes()
            noNodes = len(elNodes)
            elType = elems[elId].get_type()
            if elType == 'tria3':
                n1 = numpy.array(nodes[elNodes[0]].get_coord())
                n2 = numpy.array(nodes[elNodes[1]].get_coord())
                n3 = numpy.array(nodes[elNodes[2]].get_coord())
                vec1 = n2 - n1
                vec2 = n3 - n1
                nVec = dpTensor.UnitVector(dpTensor.CrossProduct(vec1, vec2))
                os.write('  facet normal %13.7e %13.7e %13.7e\n' % (nVec[0], nVec[1], nVec[2]))
                os.write('    outer loop\n')
                count = 0
                for nodeId in elNodes:
                    count += 1
                    os.write('      vertex %13.7e %13.7e %13.7e\n' % (nodes[nodeId].get_x(), nodes[nodeId].get_y(), nodes[nodeId].get_z()))

                os.write('    end loop\n')
                os.write('  endfacet\n')

        os.write('endsolid ascii\n')

    def writeSTLbin(self, outFileName, title, nodes, elems):
        print ' ... write STL file binary  : ', outFileName
        stdout.flush()
        outFile = open(outFileName, 'wb')
        for i in range(80):
            if i < len(title):
                outFile.write(struct.pack('c', title[i]))
            else:
                outFile.write(struct.pack('c', ' '))

        outFile.write(struct.pack('I', len(elems)))
        for elId in elems:
            elNodes = elems[elId].get_nodes()
            noNodes = len(elNodes)
            elType = elems[elId].get_type()
            if elType == 'tria3':
                n1 = numpy.array(nodes[elNodes[0]].get_coord())
                n2 = numpy.array(nodes[elNodes[1]].get_coord())
                n3 = numpy.array(nodes[elNodes[2]].get_coord())
                a = n2 - n1
                b = n3 - n1
                nVec = numpy.array([a[1] * b[2] - a[2] * b[1],
                 a[2] * b[0] - a[0] * b[2],
                 a[0] * b[1] - a[1] * b[0]])
                outFile.write(struct.pack('12fH', nVec[0], nVec[1], nVec[2], n1[0], n1[1], n1[2], n2[0], n2[1], n2[2], n3[0], n3[1], n3[2], 0))

        outFile.close()

    def writeSmesh(self, outFileName, title, nodes, elems):
        print ' ... write Tetgen smesh file   : ', outFileName
        stdout.flush()
        os = open(outFileName, 'w')
        os.write('# Part 1 - node list\n')
        os.write('# node count, 3 dim, no attribute, no boundary marker\n')
        os.write('%i %i %i %i\n' % (len(nodes), 3, 0, 0))
        for nodeId in nodes:
            os.write('%i %13.7e %13.7e %13.7e\n' % (nodeId, nodes[nodeId].get_x(), nodes[nodeId].get_y(), nodes[nodeId].get_z()))

        os.write('\n')
        os.write('# Part 2 - facet list\n')
        os.write('# facet count, no boundary marker\n')
        os.write('%i %i\n' % (len(elems), 0))
        os.write('# facets\n')
        for elId in elems:
            nodeList = elems[elId].get_nodes()
            noNodes = len(nodeList)
            elType = elems[elId].get_type()
            if elType == 'tria3':
                os.write('3 %i %i %i\n' % (nodeList[0], nodeList[1], nodeList[2]))
            elif elType == 'quad4':
                os.write('4 %i %i %i\n' % (nodeList[0], nodeList[1], nodeList[2], nodeList[2]))
            else:
                print " **ERROR** writeTetgen(): Element type '%s' for faces not implemented!" % elType
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)

        os.write('\n')
        os.write('# Part 3 - hole list\n')
        os.write('0 # no hole\n')
        os.write('# Part 4 - region list\n')
        os.write('0 # no region\n')
        os.write('\n')
        os.close()
        os = open(outFileName.replace('smesh', 'mtr'), 'w')
        os.write('%i %i\n' % (len(nodes), 1))
        for nodeId in nodes:
            os.write(' 2.0000000e+01\n')

        os.close()

    def writeSurf(self, outFileName, title, nodes, elems):
        print ' ... write netgen surf file : ', outFileName
        stdout.flush()
        os = open(outFileName, 'w')
        os.write('surfacemesh\n')
        nodeIdTable = {}
        count = 0
        os.write('%i\n' % len(nodes))
        for nodeId in nodes:
            count += 1
            nodeIdTable[nodeId] = count
            os.write('%13.7e %13.7e %13.7e\n' % (nodes[nodeId].get_x(), nodes[nodeId].get_y(), nodes[nodeId].get_z()))

        os.write('%i\n' % len(elems))
        for elId in elems:
            nodes2 = elems[elId].get_nodes()
            noNodes = len(nodes2)
            if noNodes != 3:
                print " **ERROR** writeSurf(): '%s' number of nodes per element not implemented!" % repr(noNodes)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            for node in nodes2:
                os.write('%13i ' % nodeIdTable[node])

            os.write('\n')

    def writeGmsh(self, outFileName, title, nodes, elems):
        print ' ... write Gmsh file         :', outFileName
        stdout.flush()
        os = open(outFileName, 'w')
        os.write('$MeshFormat\n')
        os.write('2.1 0 8\n')
        os.write('$EndMeshFormat\n')
        os.write('$Nodes\n')
        os.write('%s\n' % repr(len(nodes)))
        for nodeId in nodes:
            os.write('%s %13.7e %13.7e %13.7e\n' % (repr(nodeId), nodes[nodeId].get_x(), nodes[nodeId].get_y(), nodes[nodeId].get_z()))

        os.write('$EndNodes\n')
        os.write('$Elements\n')
        os.write('%s\n' % repr(len(elems)))
        for elId in elems:
            os.write('%s ' % repr(elId))
            NodesPerElement = len(elems[elId].get_nodes())
            elType = elems[elId].get_type()
            if elType == None:
                if NodesPerElement == 8:
                    elType = 5
                elif NodesPerElement == 3:
                    elType = 2
                elif NodesPerElement == 2:
                    elType = 1
                else:
                    stdout.write('\n **ERROR** writeGmsh() : Element %s has %s nodes: Not implemented!\n\n' % (elId, NodesPerElement))
                    stdout.flush()
                    stdout.write('           Use elType parameter.\n\n')
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
            elif elType == 'bar2':
                elType = 1
            elif elType == 'tria3':
                elType = 2
            elif elType == 'quad4':
                elType = 3
            elif elType == 'penta6':
                elType = 6
            elif elType == 'hexa8':
                elType = 5
            elif elType == 'tetra4':
                elType = 4
            elif elType == 'bar3':
                elType = 8
            elif elType == 'tria6':
                elType = 9
            else:
                stdout.write("\n **ERROR** writeGmsh() : Element Type '%s' not implemented!\n\n" % repr(elType))
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            reg_phys = 1
            reg_elem = 1
            os.write('%s 3 %s %s 0 ' % (repr(elType), repr(reg_phys), repr(reg_elem)))
            for node in elems[elId].get_nodes():
                os.write('%s ' % repr(node))

            os.write('\n')

        os.write('$EndElements\n')
        return

    def getFilenameAndExtension(self, fileName):
        """ Function returns file extension and file name. """
        parts = fileName.split('.')
        nparts = len(parts)
        ext = parts[nparts - 1]
        filename = ''
        for part in range(nparts - 1):
            if part < nparts - 2:
                filename = filename + parts[part] + '.'
            else:
                filename = filename + parts[part]

        return (
         filename, ext)

    def writeHistogram(self, outFileName, errorType, normalize, noInterval, nodes, elems):
        """
        Writes histogram data into an output file. 
        """
        stdout.write(' ... write histogram\n')
        values = []
        minVal = 0.0
        maxVal = 0.0
        if errorType.upper() == 'AR':
            minVal = 0.0
            maxVal = 1.0
            for elId in elems:
                nodeids = elems[elId].get_nodes()
                armax = 0.0
                armin = 1000000000000000.0
                ar = 0.0
                for lnoid in range(len(nodeids)):
                    a = nodes[nodeids[lnoid]].get_coord()
                    l1 = numpy.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
                    if lnoid == len(nodeids) - 1:
                        a = nodes[nodeids[0]].get_coord()
                    else:
                        a = nodes[nodeids[lnoid + 1]].get_coord()
                    l2 = numpy.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
                    if l1 > l2:
                        ar = l2 / l1
                    else:
                        ar = l1 / l2
                    if ar > armax:
                        armax = ar
                    if ar < armin:
                        armin = ar

                values.append(armin / armax)

        else:
            stdout.write("\n **ERROR**: fec.writeHistogram() Error type '%s' not supported !\n\n" % errorType)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        norm = False
        if normalize == 'y' or normalize == 'Y':
            norm = True
        noInterval = int(noInterval)
        histogram = {}
        dInt = (maxVal - minVal) / float(noInterval - 1)
        noVal = len(values)
        for i in range(noInterval):
            histogram[i] = 0

        for value in values:
            i = int((value - minVal) / dInt)
            histogram[i] += 1

        os = open(outFileName, 'w')
        sum = 0.0
        for i in range(noInterval):
            if norm:
                sum += histogram[i] / float(noVal) * 100.0
                os.write('%12.5f %13.6g  %13.6g\n' % (i * dInt + minVal,
                 histogram[i] / float(noVal) * 100.0, sum))
            else:
                sum += histogram[i]
                os.write('%12.5f %12i %13i\n' % (i * dInt + minVal, histogram[i], sum))

        os.close()

    def changeOrder(self, nodes, elems, order):
        """
        Change the element order
        """
        stdout.write(' ... change element order to %s\n' % repr(order))
        if int(order) == 2:
            midNodeMap = {}
            midNodeMap['bar2'] = [(1, 2)]
            midNodeMap['tria3'] = [(1, 2), (2, 3), (1, 3)]
            midNodeMap['quad4'] = [(1, 2), (2, 3), (3, 4), (1, 4)]
            midNodeMap['penta6'] = [(1, 2), (2, 3), (1, 3), (4, 5), (5, 6), (4, 6), (1, 4), (2, 5), (3, 6)]
            midNodeMap['hexa8'] = [(1, 2), (2, 3), (3, 4), (1, 4), (5, 6), (6, 7), (7, 8), (5, 8), (1, 5), (2, 6), (3, 7), (4, 8)]
            midNodeMap['tetra4'] = [(1, 2), (2, 3), (1, 3), (1, 4), (2, 4), (3, 4)]
            elType = {}
            elType['bar2'] = 'bar3'
            elType['tria3'] = 'tria6'
            elType['quad4'] = 'quad8'
            elType['penta6'] = 'penta15'
            elType['hexa8'] = 'hexa20'
            elType['tetra4'] = 'tetra10'
            maxNodeId = 0
            for noid in nodes:
                if noid > maxNodeId:
                    maxNodeId = noid

            nodePairMap = {}
            for elid in elems:
                noids = elems[elid].get_nodes()
                midNodes = midNodeMap[elems[elid].get_type()]
                for noPair in midNodes:
                    noid1 = noids[noPair[0] - 1]
                    noid2 = noids[noPair[1] - 1]
                    if not nodePairMap.has_key((noid1, noid2)):
                        maxNodeId += 1
                        nodePairMap[noid1, noid2] = maxNodeId
                        nodePairMap[noid2, noid1] = maxNodeId
                        x = (nodes[noid1].get_x() + nodes[noid2].get_x()) / 2.0
                        y = (nodes[noid1].get_y() + nodes[noid2].get_y()) / 2.0
                        z = (nodes[noid1].get_z() + nodes[noid2].get_z()) / 2.0
                        newNode = node(maxNodeId, x, y, z)
                        nodes[maxNodeId] = newNode
                    midNodeId = nodePairMap[noid1, noid2]
                    noids.append(midNodeId)

                elems[elid].set_nodes(noids)
                elems[elid].set_type(elType[elems[elid].get_type()])

        else:
            stdout.write('\n **ERROR** changeOrder() : Order change from 2->1 not implemented!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return (
         nodes, elems)

    def refine(self, nodes, elems, elsets):
        """
        Refines the mesh 
        """
        stdout.write(' ... refine mesh \n')
        stdout.flush()
        conNodeMap = {}
        conNodeMap['bar2'] = [
         (1, 2), (2, 3)]
        conNodeMap['tria3'] = [(1, 4, 6), (4, 2, 5), (4, 5, 6), (6, 5, 3)]
        conNodeMap['quad4'] = [(1, 5, 9, 8), (5, 2, 6, 9), (9, 6, 3, 7), (8, 9, 7, 4)]
        conNodeMap['penta6'] = [(1, 7, 9, 13, 16, 18), (7, 2, 8, 16, 14, 17), (9, 8, 3, 18, 17, 15), (7, 8, 9, 16, 17, 18),
         (13, 16, 18, 4, 10, 12), (16, 14, 17, 10, 5, 11), (18, 17, 15, 12, 11, 6), (16, 17, 18, 10, 11, 12)]
        conNodeMap['hexa8'] = [(1, 9, 21, 12, 17, 23, 27, 26), (9, 2, 10, 21, 23, 18, 24, 27), (12, 21, 11, 4, 26, 27, 25, 20), (21, 10, 3, 11, 27, 24, 19, 25),
         (17, 23, 27, 26, 5, 13, 22, 16), (23, 18, 24, 27, 13, 6, 14, 22), (26, 27, 25, 20, 16, 22, 15, 8), (27, 24, 19, 25, 22, 14, 7, 15)]
        conNodeMap['tetra4'] = [
         (1, 5, 7, 8), (5, 2, 6, 9), (6, 3, 7, 10), (8, 9, 10, 4), (5, 7, 8, 10), (5, 6, 7, 10), (5, 8, 9, 10), (6, 5, 9, 10)]
        addNodeMap = {}
        addNodeMap['bar2'] = []
        addNodeMap['tria3'] = []
        addNodeMap['quad4'] = [
         (5, 7, 8, 6)]
        addNodeMap['penta6'] = [(13, 14, 7, 10), (14, 15, 8, 11), (15, 13, 9, 12)]
        addNodeMap['hexa8'] = [(10, 12, 9, 11), (14, 16, 13, 15), (17, 18, 9, 13), (18, 19, 10, 14), (19, 20, 11, 15), (20, 17, 12, 16), (21, 22, 26, 24)]
        addNodeMap['tetra4'] = []
        elType = {}
        elType['bar3'] = 'bar2'
        elType['tria6'] = 'tria3'
        elType['quad8'] = 'quad4'
        elType['penta15'] = 'penta6'
        elType['hexa20'] = 'hexa8'
        elType['tetra10'] = 'tetra4'
        nodes, elems = self.changeOrder(nodes, elems, 2)
        maxNodeId = 0
        for noid in nodes:
            if noid > maxNodeId:
                maxNodeId = noid

        nodePairMap = {}
        newElems = {}
        newElemsId = 0
        mapElset = {}
        newElsets = {}
        for elset in elsets:
            newElsets[elset] = []
            for elid in elsets[elset]:
                mapElset[elid] = elset

        for elid in elems:
            noids = elems[elid].get_nodes()
            curElType = elType[elems[elid].get_type()]
            for noPair in addNodeMap[curElType]:
                noid1 = noids[noPair[0] - 1]
                noid2 = noids[noPair[1] - 1]
                noid3 = noids[noPair[2] - 1]
                noid4 = noids[noPair[3] - 1]
                if nodePairMap.has_key((noid1, noid2)) == False and nodePairMap.has_key((noid3, noid4)) == False:
                    maxNodeId += 1
                    nodePairMap[noid1, noid2] = maxNodeId
                    nodePairMap[noid2, noid1] = maxNodeId
                    nodePairMap[noid3, noid4] = maxNodeId
                    nodePairMap[noid4, noid3] = maxNodeId
                    x = (nodes[noid1].get_x() + nodes[noid2].get_x()) / 2.0
                    y = (nodes[noid1].get_y() + nodes[noid2].get_y()) / 2.0
                    z = (nodes[noid1].get_z() + nodes[noid2].get_z()) / 2.0
                    newNode = node(maxNodeId, x, y, z)
                    nodes[maxNodeId] = newNode
                midNodeId = nodePairMap[noid1, noid2]
                noids.append(midNodeId)

            for noPair in conNodeMap[curElType]:
                nList = []
                for noid in noPair:
                    nList.append(noids[noid - 1])

                newElemsId += 1
                newElems[newElemsId] = element(newElemsId, nList, curElType)
                newElsets[mapElset[elid]].append(newElemsId)

        del elems
        return (
         nodes, newElems, newElsets)

    def bboxcrop(self, bbox, nodes, elems, nsets, elsets):
        """
        Crop part of the model 
        """
        stdout.write(' ... crop mesh \n')
        stdout.flush()
        x1 = float(bbox[0])
        y1 = float(bbox[1])
        z1 = float(bbox[2])
        x2 = float(bbox[3])
        y2 = float(bbox[4])
        z2 = float(bbox[5])
        keepNid = {}
        keepElid = {}
        allElids = elems.keys()
        for elid in allElids:
            nodeList = elems[elid].get_nodes()
            xs = 0
            ys = 0
            zs = 0
            nn = float(len(nodeList))
            for nd in nodeList:
                xs += nodes[nd].get_x() / nn
                ys += nodes[nd].get_y() / nn
                zs += nodes[nd].get_z() / nn

            if xs >= x1 and xs <= x2 and ys >= y1 and ys <= y2 and zs >= z1 and zs <= z2:
                for nid in nodeList:
                    keepNid[nid] = True

            else:
                del elems[elid]
                for elsetName in elsets:
                    if elsets[elsetName].count(elid) > 0:
                        elsets[elsetName].remove(elid)

        allNodes = nodes.keys()
        for nid in allNodes:
            if not keepNid.has_key(nid):
                del nodes[nid]
                for nsetName in nsets:
                    if nsets[nsetName].count(nid) > 0:
                        nsets[nsetName].remove(nid)

        return (
         nodes, elems, nsets, elsets)

    def extract(self, setname, nodes, elems, nsets, elsets):
        """
        Crop part of the model 
        """
        stdout.write(' ... extract mesh \n')
        stdout.flush()
        keepNid = {}
        keepElid = {}
        if elsets.has_key(setname):
            for elid in elsets[setname]:
                keepElid[elid] = True
                for nid in elems[elid].get_nodes():
                    keepNid[nid] = True

            allElsetsNames = elsets.keys()
            for elsetName in allElsetsNames:
                if elsetName != setname:
                    del elsets[elsetName]

            allNodes = nodes.keys()
            for nid in allNodes:
                if not keepNid.has_key(nid):
                    del nodes[nid]
                    for nsetName in nsets:
                        if nsets[nsetName].count(nid) > 0:
                            nsets[nsetName].remove(nid)

            allNsetsNames = nsets.keys()
            for nsetName in allNsetsNames:
                if len(nsets[nsetName]) == 0:
                    del nsets[nsetName]

            newElems = {}
            for elid in elsets[setname]:
                newElems[elid] = elems[elid]

        else:
            stdout.write('\n **ERROR** extract() : elset=%s not found! Available Sets:\n' % setname)
            stdout.write('\n %s \n\n' % elsets.keys())
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return (
         nodes, newElems, nsets, elsets)

    def translate(self, transXYZ, nodes):
        """
        Crop part of the model 
        """
        stdout.write(' ... translate mesh \n')
        stdout.flush()
        x1 = float(transXYZ[0])
        y1 = float(transXYZ[1])
        z1 = float(transXYZ[2])
        for noid in nodes:
            nodes[noid].set_x(nodes[noid].get_x() + x1)
            nodes[noid].set_y(nodes[noid].get_y() + y1)
            nodes[noid].set_z(nodes[noid].get_z() + z1)

        return nodes

    def transform(self, nodes, tFilename, exponent=1):
        stdout.write(' ... transform mesh\n')
        stdout.flush()
        myTransformer = dpTransform.Transformer()
        myTransformer.readTransformMatrix(tFilename)
        myTransformer.setTransformDirection(exponent)
        for nid in nodes:
            node = nodes[nid]
            pointOld = node.get_coord_numpy()
            pointNew = myTransformer.transformVector(pointOld)
            node.set_coord_numpy(pointNew)

        return nodes

    def getElementsAroundNode(self, nodes, elems):
        nodeList = {}
        for node in nodes:
            nodeList[node] = []

        for elemId in elems:
            elemNodes = elems[elemId].get_nodes()
            for elNode in elemNodes:
                nodeList[elNode].append(elemId)

        return nodeList

    def getNodesAroundNode(self, nodes, elems):
        elemsAroundNode = self.getElementsAroundNode(nodes, elems)
        nodeList = {}
        for node in nodes:
            nodeList[node] = []
            elemList = elemsAroundNode[node]
            for elemId in elemList:
                elemNodes = elems[elemId].get_nodes()
                for elNode in elemNodes:
                    if nodeList[node].count(elNode) == 0 and elNode != node:
                        nodeList[node].append(elNode)

        return nodeList

    def smoothLaplacian(self, nodes, elems, itMax, noidsToSmooth=None):
        stdout.write(' ... Laplace smoothing \n')
        stdout.flush()
        nodesAroundNode = self.getNodesAroundNode(nodes, elems)
        if not noidsToSmooth:
            noidsToSmooth = nodes
        for it in range(itMax):
            for noid in noidsToSmooth:
                xi = 0.0
                yi = 0.0
                zi = 0.0
                for nnoid in nodesAroundNode[noid]:
                    x, y, z = nodes[nnoid].get_coord()
                    xi += x
                    yi += y
                    zi += z

                nNodes = len(nodesAroundNode[noid])
                if nNodes > 0:
                    nk = float(len(nodesAroundNode[noid]))
                    nodes[noid].set_x(xi / nk)
                    nodes[noid].set_y(yi / nk)
                    nodes[noid].set_z(zi / nk)

        return nodes

    def smoothTaubin(self, nodes, elems, itMax, lam, kPB, noidsToSmooth=None):
        stdout.write(' ... Taubin smoothing \n')
        stdout.flush()
        weight = 1.0
        mu = 1.0 / (kPB - 1.0 / lam)
        nodesAroundNode = self.getNodesAroundNode(nodes, elems)
        if not noidsToSmooth:
            noidsToSmooth = nodes
        for it in range(itMax):
            for noid in noidsToSmooth:
                xi, yi, zi = nodes[noid].get_coord()
                dxi = 0.0
                dyi = 0.0
                dzi = 0.0
                nNodes = len(nodesAroundNode[noid])
                if nNodes > 0:
                    wij = weight / float(nNodes)
                    for nnoid in nodesAroundNode[noid]:
                        xj, yj, zj = nodes[nnoid].get_coord()
                        dxi += wij * (xi - xj)
                        dyi += wij * (yi - yj)
                        dzi += wij * (zi - zj)

                    if it % 2 == 0 or it == 0:
                        nodes[noid].set_x(xi + lam * dxi)
                        nodes[noid].set_y(yi + lam * dyi)
                        nodes[noid].set_z(zi + lam * dzi)
                    else:
                        nodes[noid].set_x(xi + mu * dxi)
                        nodes[noid].set_y(yi + mu * dyi)
                        nodes[noid].set_z(zi + mu * dzi)

        return nodes

    def get_mesh(self, nodeCoords, faceConn, faceSets, mtype):
        nodes = {}
        i = 0
        for nodeCoord in nodeCoords:
            i += 1
            nodes[i] = node(i, nodeCoord[0], nodeCoord[1], nodeCoord[2])

        elems = {}
        i = 0
        for face in faceConn:
            i += 1
            nList = list(face + 1)
            elems[i] = element(i, nList, mtype)

        elsets = {}
        i = 0
        for faceSet in faceSets:
            i += 1
            gv1 = faceSet[0]
            gv2 = faceSet[1]
            name = 'SET'
            if gv1 == -1:
                gv1 = 0
                name = 'BSET'
            key = name + '_' + repr(gv1) + '_' + repr(gv2)
            if key not in elsets:
                elsets[key] = []
            elsets[key].append(i)

        return (
         nodes, elems, elsets)


if __name__ == '__main__':
    inName = None
    outName = None
    histo = None
    order = None
    refin = None
    smooth = None
    bbcrop = None
    extr = None
    trans = None
    tfil = None
    texp = None
    binFlag = None
    guiInit = "*fileEntryIn    -in     'Input File Name       | filename '                              test.raw             no   inp;raw;off;msh;vtk;stl;obj;fem              \n" + "*fileEntryOut   -out    'Output File Name      | filename '                              test.inp             no   geom;inp;geo;case;vtk;off;msh;stl;surf;smesh \n" + "*radioEntry     -bin    'Binary File Flag      | flag'                                   ON                   yes  -        \n*subWindowStart 'Crop'                                                                                                                                          \n" + "*splitEntry     -bbcrop 'Crop part of model    | start1;start2;start3;end1;end2;end3'    0.;0.;0.;1.;2.;3.    yes  6                                            \n" + "*entry          -extr   'Extract Elset         | setname'                                SET123               yes  1                                            \n" + "*subWindowEnd   'Crop'                                                                                                                                          \n" + "*subWindowStart 'Modify'                                                                                                                                        \n" + "*combo          -refin  'Refine Mesh           | flag'                                   OFF                  yes  ON;OFF                                       \n" + "*combo          -order  'Change mesh order     | order'    1                             yes  1;2                                                               \n" + "*entry          -smooth 'Smooth tria3 mesh     | type;niter;[param2];[param3]'           Taubin;31;0.1;0.4    yes  3                                            \n" + "*subWindowEnd   'Modify'                                                                                                                                        \n" + "*subWindowStart 'Transform'                                                                                                                                     \n" + "*splitEntry     -tra    'Translate             | tra1;tra2;tra3'                         0.0;0.0;0.0          yes  3                                            \n" + "*fileEntryIn    -tfil   'Transform Matrix File | filename'                               trans.mhd            yes  mhd;nhdr;tfm                                 \n" + "*combo          -texp   'Direction Exponent    | exponent'                               1                    yes  1;-1                                         \n" + "*subWindowEnd   'Transform'                                                                                                                                     \n" + "*subWindowStart 'Info'                                                                                                                                          \n" + "*splitEntry     -histo  'Histogram output      | filename;errorVal;normalize;nIntervals' test.txt;ar;y;100    yes  4                                            \n" + "*subWindowEnd   'Info'                                                                                                                                          \n"
    argList = argv
    argc = len(argList)
    i = 0
    while i < argc:
        if argList[i][:2] == '-i':
            i += 1
            inName = argList[i]
        elif argList[i][:4] == '-out':
            i += 1
            outName = argList[i]
        elif argList[i][:4] == '-bin':
            i += 1
            binFlag = argList[i]
        elif argList[i][:6] == '-histo':
            i += 1
            histo = argList[i]
        elif argList[i][:6] == '-order':
            i += 1
            order = argList[i]
        elif argList[i][:6] == '-refin':
            i += 1
            refin = argList[i]
        elif argList[i][:7] == '-smooth':
            i += 1
            smooth = argList[i]
        elif argList[i][:7] == '-bbcrop':
            i += 1
            bbcrop = argList[i]
        elif argList[i][:5] == '-extr':
            i += 1
            extr = argList[i]
        elif argList[i][:4] == '-tra':
            i += 1
            trans = argList[i]
        elif argList[i][:5] == '-tfil':
            i += 1
            tfil = argList[i]
        elif argList[i][:5] == '-texp':
            i += 1
            texp = argList[i]
        elif argList[i][:4] == '-gui':
            print guiInit
            exit(0)
        elif argList[i][:5] == '-help':
            print __doc__
            exit(0)
        i += 1

    if not inName:
        print __doc__
        print ' **ERROR** input file name not given\n'
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if not outName:
        print __doc__
        print ' **ERROR** output file name not given\n'
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    stdout.write('\n S T A R T  %s %s - %s \n\n' % (fec.__name__, fec.__version__, fec.__author__))
    stdout.flush()
    myModel = fec()
    if not binFlag:
        binFlag = 'OFF'
    if binFlag.upper() == 'ON':
        binFlag = True
    else:
        binFlag = False
    title, nodes, nsets, elems, elsets = myModel.read(inName, isBinary=binFlag)
    stdout.write(' ... read summary\n')
    stdout.write('     -> Processed Nodes      : %10i        \n' % len(nodes))
    stdout.write('     -> Processed Elements   : %10i        \n' % len(elems))
    stdout.write('     -> Processed Node Sets  : %10i        \n' % len(nsets))
    stdout.write('     -> Processed Elems Set  : %10i        \n' % len(elsets))
    stdout.flush()
    if histo != None:
        filename, errorType, normalize, noInterval = dpUtils.userSplit(histo)
        myModel.writeHistogram(filename, errorType, normalize, noInterval, nodes, elems)
    if refin != None:
        if refin.upper() == 'ON':
            nodes, elems, elsets = myModel.refine(nodes, elems, elsets)
    if order != None:
        nodes, elems = myModel.changeOrder(nodes, elems, int(order))
    if smooth != None:
        smoothPar = dpUtils.userSplit(smooth)
        if smoothPar[0].upper() == 'LAPLACE':
            nodes = myModel.smoothLaplacian(nodes, elems, int(smoothPar[1]))
        elif smoothPar[0].upper() == 'TAUBIN':
            nodes = myModel.smoothTaubin(nodes, elems, int(smoothPar[1]), float(smoothPar[2]), float(smoothPar[3]))
        else:
            stdout.write(" **ERROR** : fec() option '-smooth': type = %s elements!\n" % smoothPar[0])
            stdout.write("             It has to be 'Laplace' or 'Taubin'!\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    if bbcrop != None:
        bbcropList = dpUtils.userSplit(bbcrop)
        if len(bbcropList) != 6:
            stdout.write(" **ERROR** : fec() option '-bbcrop' requires six parameters!\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        bbox = [
         float(bbcropList[0]), float(bbcropList[1]), float(bbcropList[2]),
         float(bbcropList[3]), float(bbcropList[4]), float(bbcropList[5])]
        nodes, elems, nsets, elsets = myModel.bboxcrop(bbox, nodes, elems, nsets, elsets)
    if trans != None:
        transList = dpUtils.userSplit(trans)
        if len(transList) != 3:
            stdout.write(" **ERROR** : fec() option '-trans' requires three parameters!\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        transXYZ = [
         float(transList[0]), float(transList[1]), float(transList[2])]
        nodes = myModel.translate(transXYZ, nodes)
    if extr != None:
        nodes, elems, nsets, elsets = myModel.extract(extr, nodes, elems, nsets, elsets)
    if tfil != None:
        if texp != None:
            nodes = myModel.transform(nodes, tfil, exponent=int(texp))
        else:
            stdout.write(" **ERROR** : fec() option '-tfil' requires alos '-texp' option!\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    myModel.write(outName, title, nodes, nsets, elems, elsets, isBinary=binFlag)
    stdout.write('\n E N D E D  SUCCESSFULLY\n\n')
    stdout.flush()

# file fec.pyc
# Deparsing stopped due to parse error
