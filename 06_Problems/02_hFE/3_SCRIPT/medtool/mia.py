

"""
Fabric
------

Script performs morpholocial analyses of medical voxelized CT images 
based on MIL, SLD, SVD, and GST. 

.. image:: images/mia.png

For MIL, SLD, SVD, morphological distributions are computed by a "ray 
field" method (ray field through RVE). The original distribution can 
be approximated by an ellipsoid, general 2nd or general 4th order tensor.
Eigenvalues and eigenvectors are computed from 2nd order tensors. Graphical 
output of original and approximated distribution, as well as for eigen
vectors is done in Ensight Gold Format. The main information is written 
in a text file (extension '.fab'). 

This class also computes a fabric tensor based on the gradient structure 
tensor (GST). Beside the required parameters only 'in', 'out', 'dtype', 
'gstpow', 'norm', 'sort', 'grap', and 'mod' are implemented. 

Furthermore, the class mia can be used as library.    

Usage
~~~~~

Module ::

  import mia

Command line: ::

  python mia.py ...

    -in       infileName
    -out      outfileName(:mode)
    -dtype    distribType
    -ftype    fabricType     
    [-thres]  threshold        [optional]
    [-step]   voxelStep        [optional]
    [-power]  starDirectParam  [optional]
    [-valid]  minIntersectLen  [optional]
    [-code]   code             [optional]
    [-gstpow] gstpow           [optional]
    [-norm]   normalizeType    [optional]
    [-mod]    outputMode       [optional]
    [-graph]  graphicOutFlag   [optional]
    [-distr]  textOutFlag      [optional]
    [-sort]   sortFlag         [optional]
    [-fout]   fabricOutfile    [optional]
    [-help]                    [optional]


Parameters
~~~~~~~~~~
 
-in    :  filename

          'filename' is a thresholded voxel file. 
          A segmented image (one threshold) with an isotropic resolution 
          is needed. The file is allowed to contain two values: 
          min=0, max=threshold (>0). If the file is not segmented 
          it is necessary to give the option 'thres'. 
          
          
-out   :  filename

          Text output file where general informations like eigenvectors 
          and eigenvalues are written. 2 possible ways to use it: 
          PARAM: filename -> '.fab' file is written
          PARAM: filename 'mod' option -> Write general information in more 
          compact form (one line) to a file.
          The "mode" can be "w" (create new file, write header + data) and
          "a" (append to file, write data).


-dtype  : type

          Distribution type. 'type' can be 'MIL', 'SLD', 'SVD', 'GST'


-ftype  : type

          Type of computed fabric tensor. An approach which is discribed  
          by Zysset 1995 et. al. is used. There, the original distribution 
          computed with the MIL, SLD or SVD method is fitted with a 
          general tensor of second or fourth order. ('type'=2 or 'type'=4). 
          Furthermore, an ellipsoid can be used for the fit ('type' = 1).
          For 'GST' only 'type' = 1 is implemented. 
 
 
-thres  : threshold

          Threshold used for the computation of BVTV in case of a 
          given gray-value file. The images is segmented and interally 
          scaled to 0...'threshold'.
          NOTE: not required for 'GST' or in case of a thresholded file.


-step  :  distance

          'distance' is the number of voxels between individual rays in 
          ray field for the fabric computations. 
          Default = 5;
          NOTE: not required for 'GST'


-power  : power

          Parameter for the computation of the number of star directions   
          on a unit sphere. The number of directions is given by ::
          
             nDirections = 8*4^'power' 
          
          The direction algorithm projects triangles on a unit sphere 
          surface and refines this mesh. The refinement is done by dividing 
          each triangle into four triangles. Default = 2 which gives in sum 
          128 directions.
          NOTE: not required for 'GST'


-valid  : validLength

          Smallest valid intersection length. Intersections are only 
          counted if ::
          
             isLen > validLength 
          
          The intersection length 'isLen' is computed from mid-to-mid point 
          of the bone voxels in voxel dimensions (i.e. voxel length = 1). 
          
          For example, 'isLen' = sqrt(3) means that more then 2 voxel 
          in each direction build a valid intersection. 
          If only one voxel builds the intersections than 'isLen'=0 and 
          the voxel is not counted if 'validLength' = 0.
          MIL is sensitive to this parameter. 
          Default = 0.0;   
          NOTE: not required for 'GST'


-code  :  flag

          Type of the used code for the core routines. Possible choices 
          for 'flag' are "py2 or "f77".
          NOTE: not required for 'GST'


-gstpow  : power

           Power coefficient only needed for GST method where the original
           eigenvalues from the GST evaluation (say l_i) are taken 
           to compute the eigenvalue (m_i) of the fabric M as ::
            
             m_i = 1.0/l_i**'power'
           
           Default = 0.5 


-norm   : type

          Type how the fabric should be normalize. Currently implemented 
          'type' parameters are:

            - "det"  : in that case ``m1*m2*m3 = 1`` (det(M)=1)
            - "rank"/"trace" : in that case ``m1+m2+m3 = 3`` 

          Default is 'trace'


-mod   :  mode
 
          The output 'mode' can be 'w' (create new file, write header
          + data) and 'a' (append to file, write data).
          Default = 'w'


-graph  : flag

          Switch on the writing of the original and approximated 
          distributions are written in Ensight Gold Files ('.case' and 
          '.geo'). 'flag' can be "ON" or "OFF". The results and can be 
          viewed with Paraview
          Default = "OFF".


-distr  : filename

          Direction and intercept lengths of original distributions 
          are written to a text file. 
          NOTE: not required for 'GST'


-sort   : flag

          Switch sorting of eigen values and corresponding eigen 
          vectors. 'flag' can be "ON" or "OFF". 
          Default = "OFF".
          
          
-fout   : filename

          Output quantities to compute the fabric tensors by hand. 
          The parameter specifies the file name. 


-help   : Print usage

 
Info
~~~~

- File:   mia.py 
- Author: D. H. Pahr
  
"""
from sys import argv, exit, stdout
from time import clock, time, localtime
from string import split, replace, atof
import numpy
import numpy.linalg as linalg
import dpTensor
import dpFem
import mic
import fec
import dpStiffness
import dpMesh
import os
import miaf77

class mia():
    __version__ = 'V_02.05.2015'
    __author__ = 'D. H. Pahr'

    def __init__(self, modelName='default'):
        """ (M)edical (I)mage (A)nalysis Class  """
        self.modelName = modelName

    def computeBVTV(self, voxelModel, threshold, rveType='hex'):
        """      
        Computes density/BVTV of a given voxel model. Only one threhold is 
        allowed. 
        
        @param  voxelModel: segmented voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue 
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.       
          - int grayValue  ... gray value of voxel, 0..255  
                                       
        @param  threshold:  current threshold for segmentation 
          - TYPE: int
          
        @param  rveType:  type of rve, default = hexahedral  
          - TYPE: string  ... 'sphere', 'hex' 
                                     
        @return: 
          density: BVTV of the model                                       
        """
        stdout.write(' ... compute BV/TV          = ')
        stdout.flush()
        nx, ny, nz = self.get_Shape(voxelModel)
        flatArray = numpy.reshape(voxelModel, (nz * ny * nx,))
        BV = numpy.sum(numpy.greater(flatArray, threshold - 1))
        TV = nz * ny * nx
        del flatArray
        rho = 0.0
        if rveType == 'sphere':
            rho = BV / float(TV) / numpy.pi * 6.0
        elif rveType == 'hex':
            rho = BV / float(TV)
        else:
            stdout.write("\n **ERROR**: mia.computeBVTV(): rveType %s' not supported!\n\n" % rveType)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        stdout.write('%10.3f %s\n' % (rho * 100.0, '%'))
        stdout.flush()
        return rho

    def computeOrigDistribution_STAR(self, thresVoxelModel, step, power):
        """
        Computes MIL/SLD/SVD distributions for directions vectors "n" 
        using a star in bone voxels. Normals n = (nix,niy,niz) 
        are the directions from the midpoint of a unit sphere to the COG of triangles 
        which build the surface of the sphere.  
        A segement voxel model with isotropic resolution is needed.   
        
        @param thresVoxelModel: segmented voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue 
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.       
          - int grayValue  ... gray value of voxel, 0..255                  
        @param step: Step in the considered voxel 
          - TYPE: int > 0 
        @param power: Power for the number of star directions = 8*4^power
          - TYPE: int > 1
                             
        @return: 
             MIL             : Mean Intercept Length 
               - TYPE: dict[ (nix, niy, niz) ] = value 
               - float nix, niy, niz ... components of normal vectors 
               - float value         ... MIL for that direction          
             SLD             : Star Length Distribution 
               - TYPE: same as for MIL 
             SVD             : Star Volume Distribution 
               - TYPE: same as for MIL       
             Area            : Weight (Area of triange on unit sphere) for each direction 
               - TYPE: same as for MIL                         
        """
        stdout.write(' ... compute original distribution\n')
        stdout.flush()
        time1 = clock()
        nx, ny, nz = self.get_Shape(thresVoxelModel)
        self.nx = nx
        self.ny = ny
        self.nz = nz
        faceInfo = self.getFaceInfo(thresVoxelModel)
        aArray, Area_n = self.computeNormalAndArea(power)
        randPointArray = []
        jumper = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if thresVoxelModel[k, j, i] > 0:
                        jumper += 1
                        if jumper == step:
                            randPointArray.append((i, j, k))
                            jumper = 0

        c = 0
        sumL = {}
        sumL2 = {}
        sumL4 = {}
        for direction in aArray:
            if direction[2] > 0.0:
                progress = float(c) / float(len(aArray)) * 100.0
                stdout.write('     -> Processed Data      : %10i %s\r' % (progress, '%'))
                stdout.flush()
                minLen = 0.0
                direct1 = direction
                direct2 = (-direct1[0], -direct1[1], -direct1[2])
                for point in randPointArray:
                    minLen1 = self.computeLength(faceInfo, point, direct1, nx, ny, nz)
                    minLen2 = self.computeLength(faceInfo, point, direct2, nx, ny, nz)
                    sum = minLen1 + minLen2
                    if sumL.has_key(direct1):
                        sumL[direct1] += sum
                        sumL2[direct1] += sum * sum
                        sumL4[direct1] += sum * sum * sum * sum
                    else:
                        sumL[direction] = sum
                        sumL2[direction] = sum * sum
                        sumL4[direction] = sum * sum * sum * sum

                c += 2

        n = len(randPointArray)
        MIL = {}
        SVD = {}
        SLD = {}
        for direc in sumL:
            direct1 = direc
            direct2 = (-direct1[0], -direct1[1], -direct1[2])
            MIL[direct1] = sumL[direct1] / float(n)
            MIL[direct2] = sumL[direct1] / float(n)
            SLD[direct1] = sumL2[direct1] / sumL[direct1]
            SLD[direct2] = sumL2[direct1] / sumL[direct1]
            SVD[direct1] = numpy.pi / 3.0 * sumL4[direct1] / sumL[direct1]
            SVD[direct2] = numpy.pi / 3.0 * sumL4[direct1] / sumL[direct1]

        stdout.write('     -> Processed Data      : %10i   \n' % (c * len(randPointArray)))
        stdout.flush()
        time2 = clock()
        stdout.write('     -> Compute in          :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        return (
         MIL, SVD, SLD, Area_n)

    def computeOrigDistribution_MIL(self, thresVoxelModel, step, power, valid, echo=True):
        """
        Function computes MIL/SLD/SVD distributions for direction vectors "n" 
        using a voxel ray field going trought the RVE. Normals n = tuple(nix,niy,niz) 
        are the directions from the midpoint of a unit sphere to the COG of triangles 
        which build the surface of the sphere. Very similar to 
        self.computeOrigDistribution_STAR(). 
        A segement voxel model with isotropic resolution is needed.   
        
        @param thresVoxelModel: segmented voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue 
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.       
          - int grayValue  ... gray value of voxel, 0..255                  
        @param step: Step in the considered voxel 
          - TYPE: int > 0 
        @param power: Power for the number of star directions = 8*4^power
          - TYPE: int > 1
        @param valid: Smallest valid intersection length.
          - TYPE: float
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
                                     
        @return: 
             MIL             : Mean Intercept Length 
               - TYPE: dict[ (nix, niy, niz) ] = value 
               - float nix, niy, niz ... components of normal vectors 
               - float value         ... MIL for that direction          
             SLD             : Star Length Distribution 
               - TYPE: same as for MIL 
             SVD             : Star Volume Distribution 
               - TYPE: same as for MIL       
             Area            : Weight (Area of triange on unit sphere) for each direction 
               - TYPE: same as for MIL          
        """
        if echo == True:
            stdout.write(' ... compute original distribution\n')
            stdout.flush()
        time1 = clock()
        ctime1 = time()
        self.nx, self.ny, self.nz = self.get_Shape(thresVoxelModel)
        valid2 = valid * valid
        MIL = {}
        SVD = {}
        SLD = {}
        sumL = {}
        sumL2 = {}
        sumL4 = {}
        corners = {}
        corners['swb'] = (0.0, 0.0, 0.0)
        corners['seb'] = (float(self.nx), 0.0, 0.0)
        corners['neb'] = (float(self.nx), float(self.ny), 0.0)
        corners['nwb'] = (0.0, float(self.ny), 0.0)
        corners['swt'] = (0.0, 0.0, float(self.nz))
        corners['set'] = (float(self.nx), 0.0, float(self.nz))
        corners['net'] = (float(self.nx), float(self.ny), float(self.nz))
        corners['nwt'] = (0.0, float(self.ny), float(self.nz))
        modelPlanes = {}
        modelPlanes['s'] = dict([('r-dir', (1.0, 0.0, 0.0)), ('s-dir', (0.0, 0.0, 1.0)), ('base', 'swb')])
        modelPlanes['e'] = dict([('r-dir', (0.0, 1.0, 0.0)), ('s-dir', (0.0, 0.0, 1.0)), ('base', 'seb')])
        modelPlanes['n'] = dict([('r-dir', (1.0, 0.0, 0.0)), ('s-dir', (0.0, 0.0, 1.0)), ('base', 'nwb')])
        modelPlanes['w'] = dict([('r-dir', (0.0, 1.0, 0.0)), ('s-dir', (0.0, 0.0, 1.0)), ('base', 'swb')])
        modelPlanes['b'] = dict([('r-dir', (1.0, 0.0, 0.0)), ('s-dir', (0.0, 1.0, 0.0)), ('base', 'swb')])
        modelPlanes['t'] = dict([('r-dir', (1.0, 0.0, 0.0)), ('s-dir', (0.0, 1.0, 0.0)), ('base', 'swt')])
        viewerAt = {}
        viewerAt['swb'] = (1.0, 1.0, 1.0)
        viewerAt['seb'] = (-1.0, 1.0, 1.0)
        viewerAt['neb'] = (-1.0, -1.0, 1.0)
        viewerAt['nwb'] = (1.0, -1.0, 1.0)
        viewerTo = {}
        viewerTo['swb'] = 'net'
        viewerTo['seb'] = 'nwt'
        viewerTo['neb'] = 'swt'
        viewerTo['nwb'] = 'set'
        normals, area = self.computeNormalAndArea(power)
        lookUpTable = {}
        direct = ''
        for n in normals:
            nX = n[0]
            nY = n[1]
            nZ = n[2]
            if nX >= 0.0 and nY >= 0.0 and nZ >= 0.0:
                voxX = 1
                voxY = 1
                voxZ = 1
                stepX = 1
                stepY = 1
                stepZ = 1
                voxelRay = []
                voxelRay.append((voxX, voxY, voxZ))
                if abs(nX) > abs(nY):
                    if abs(nZ) > abs(nX):
                        direct = 'z'
                    else:
                        direct = 'x'
                elif abs(nZ) > abs(nY):
                    direct = 'z'
                else:
                    direct = 'y'
                preVoxX = 1
                preVoxY = 1
                preVoxZ = 1
                while abs(voxX) <= self.nx and abs(voxY) <= self.ny and abs(voxZ) <= self.nz:
                    tMaxX = voxX / nX
                    tMaxY = voxY / nY
                    tMaxZ = voxZ / nZ
                    if abs(tMaxX) < abs(tMaxY):
                        if abs(tMaxX) < abs(tMaxZ):
                            voxX = voxX + stepX
                        else:
                            voxZ = voxZ + stepZ
                    elif abs(tMaxY) < abs(tMaxZ):
                        voxY = voxY + stepY
                    else:
                        voxZ = voxZ + stepZ
                    if abs(voxX) <= self.nx and abs(voxY) <= self.ny and abs(voxZ) <= self.nz:
                        if direct == 'x':
                            if voxX > preVoxX:
                                voxelRay.append((voxX, voxY, voxZ))
                        if direct == 'y':
                            if voxY > preVoxY:
                                voxelRay.append((voxX, voxY, voxZ))
                        if direct == 'z':
                            if voxZ > preVoxZ:
                                voxelRay.append((voxX, voxY, voxZ))
                    preVoxX = voxX
                    preVoxY = voxY
                    preVoxZ = voxZ

                lookUpTable[n] = voxelRay

        i = 0
        sum = len(viewerAt) * 4.0 ** float(power)
        for vpt in viewerAt:
            curLookUpTable = {}
            cornVoxX = int(corners[vpt][0])
            cornVoxY = int(corners[vpt][1])
            cornVoxZ = int(corners[vpt][2])
            if cornVoxX == 0:
                cornVoxX = 1
            if cornVoxY == 0:
                cornVoxY = 1
            if cornVoxZ == 0:
                cornVoxZ = 1
            stepX = int(viewerAt[vpt][0])
            stepY = int(viewerAt[vpt][1])
            stepZ = int(viewerAt[vpt][2])
            for n in lookUpTable:
                voxelray = lookUpTable[n]
                newVoxelRay = []
                for voxel in voxelray:
                    VoxelX = cornVoxX + stepX * voxel[0] - stepX
                    VoxelY = cornVoxY + stepY * voxel[1] - stepY
                    VoxelZ = cornVoxZ + stepZ * voxel[2] - stepZ
                    newVoxelRay.append((VoxelX, VoxelY, VoxelZ))

                curLookUpTable[n[0] * viewerAt[vpt][0], n[1] * viewerAt[vpt][1], n[2] * viewerAt[vpt][2]] = newVoxelRay

            entryPlanes = [
             vpt[0], vpt[1], vpt[2]]
            for n in curLookUpTable:
                nL = 0
                nNotValid0 = 0
                nValid0 = 0
                nNotValid1 = 0
                nValid1 = 0
                nNotValid3 = 0
                nValid3 = 0
                sumL[n] = 0.0
                sumL2[n] = 0.0
                sumL4[n] = 0.0
                i += 1
                newVoxelRay = curLookUpTable[n]
                stdout.write('     -> Setup Data          : %10i %s\r' % (int(float(i) / float(sum) * 100), '%'))
                stdout.flush()
                nn = numpy.array([n[0], n[1], n[2]])
                nb = numpy.array((0.0, 0.0, 1.0))
                ng = dpTensor.UnitVector(dpTensor.CrossProduct(nn, nb))
                ns = dpTensor.UnitVector(dpTensor.CrossProduct(ng, nn))
                nr = dpTensor.UnitVector(dpTensor.CrossProduct(ns, nn))
                rmax = 0.0
                rmin = 0.0
                smax = 0.0
                smin = 0.0
                r1c = numpy.array(corners[vpt])
                for c in corners:
                    r0c = numpy.array(corners[c])
                    b = r0c - r1c
                    a11 = nr[0]
                    a12 = ns[0]
                    a13 = -nn[0]
                    a21 = nr[1]
                    a22 = ns[1]
                    a23 = -nn[1]
                    a31 = nr[2]
                    a32 = ns[2]
                    a33 = -nn[2]
                    DET = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13)
                    x = [0.0, 0.0, 0.0]
                    x[0] = 1.0 / DET * ((a33 * a22 - a32 * a23) * b[0] - (a33 * a12 - a32 * a13) * b[1] + (a23 * a12 - a22 * a13) * b[2])
                    x[1] = 1.0 / DET * (-(a33 * a21 - a31 * a23) * b[0] + (a33 * a11 - a31 * a13) * b[1] - (a23 * a11 - a21 * a13) * b[2])
                    if x[0] > rmax:
                        rmax = x[0]
                    if x[0] < rmin:
                        rmin = x[0]
                    if x[1] > smax:
                        smax = x[1]
                    if x[1] < smin:
                        smin = x[1]

                for curR in range(int(rmin), int(rmax + 1), step):
                    for curS in range(int(smin), int(smax + 1), step):
                        curPlanes = entryPlanes
                        check = 0
                        for surf in curPlanes:
                            cutPlane = modelPlanes[surf]
                            r1 = numpy.array(corners[cutPlane['base']])
                            r0 = curR * nr + curS * ns + r1c
                            at = nn
                            br = numpy.array(cutPlane['r-dir'])
                            cs = numpy.array(cutPlane['s-dir'])
                            b = r0 - r1
                            a11 = br[0]
                            a12 = cs[0]
                            a13 = -at[0]
                            a21 = br[1]
                            a22 = cs[1]
                            a23 = -at[1]
                            a31 = br[2]
                            a32 = cs[2]
                            a33 = -at[2]
                            DET = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13)
                            x = [0.0, 0.0, 0.0]
                            x[0] = 1.0 / DET * ((a33 * a22 - a32 * a23) * b[0] - (a33 * a12 - a32 * a13) * b[1] + (a23 * a12 - a22 * a13) * b[2])
                            x[1] = 1.0 / DET * (-(a33 * a21 - a31 * a23) * b[0] + (a33 * a11 - a31 * a13) * b[1] - (a23 * a11 - a21 * a13) * b[2])
                            ipt = x[0] * br + x[1] * cs + r1
                            if ipt[0] >= 0.0 and ipt[1] >= 0.0 and ipt[2] >= 0.0 and ipt[0] <= float(self.nx) and ipt[1] <= float(self.ny) and ipt[2] <= float(self.nz):
                                check += 1
                                entryVoxX = int(ipt[0] + 0.5)
                                entryVoxY = int(ipt[1] + 0.5)
                                entryVoxZ = int(ipt[2] + 0.5)
                                if entryVoxX == 0:
                                    entryVoxX = 1
                                if entryVoxY == 0:
                                    entryVoxY = 1
                                if entryVoxZ == 0:
                                    entryVoxZ = 1
                                startBone = (1, 1, 1)
                                endBone = (1, 1, 1)
                                preVox = (1, 1, 1)
                                count = 0
                                startFlag = False
                                for startRayVox in newVoxelRay:
                                    voxX = startRayVox[0] - (cornVoxX - entryVoxX)
                                    voxY = startRayVox[1] - (cornVoxY - entryVoxY)
                                    voxZ = startRayVox[2] - (cornVoxZ - entryVoxZ)
                                    count += 1
                                    if voxX < 1 or voxY < 1 or voxZ < 1 or voxX > self.nx or voxY > self.ny or voxZ > self.nz:
                                        if startFlag == True:
                                            if voxX > self.nx or voxY > self.ny or voxZ > self.nz:
                                                startFlag = False
                                                endBone = (prevVox[0], prevVox[1], prevVox[2])
                                                lx = startBone[0] - endBone[0]
                                                ly = startBone[1] - endBone[1]
                                                lz = startBone[2] - endBone[2]
                                                L2 = lx * lx + ly * ly + lz * lz
                                                if L2 > valid2:
                                                    nL += 1
                                                    sumL[n] += L2 ** 0.5
                                                    sumL2[n] += L2
                                                    sumL4[n] += L2 * L2
                                        if count == 1:
                                            print '**WARNING: computeOrigDistribution_MIL(): Point not inside Voxelmodel!'
                                            print 'StartRayVox=', startRayVox[0], startRayVox[1], startRayVox[2]
                                            print 'CornerVox  =', cornVoxX, cornVoxY, cornVoxZ
                                            print 'EntryVox   =', entryVoxX, entryVoxY, entryVoxZ,
                                            print 'curVox     =', voxX, voxY, voxZ
                                            print '-------------'
                                        break
                                    elif thresVoxelModel[voxZ - 1, voxY - 1, voxX - 1] == 0:
                                        if startFlag == True:
                                            startFlag = False
                                            endBone = (prevVox[0], prevVox[1], prevVox[2])
                                            lx = startBone[0] - endBone[0]
                                            ly = startBone[1] - endBone[1]
                                            lz = startBone[2] - endBone[2]
                                            L2 = lx * lx + ly * ly + lz * lz
                                            if L2 > valid2:
                                                nL += 1
                                                sumL[n] += L2 ** 0.5
                                                sumL2[n] += L2
                                                sumL4[n] += L2 * L2
                                    elif startFlag == False:
                                        startBone = (voxX, voxY, voxZ)
                                        startFlag = True
                                    prevVox = (
                                     voxX, voxY, voxZ)

                                break

                        if check > 1:
                            print '**WARNING: computeOrigDistribution_MIL(): Intersection Problems: '

                if nL > 0:
                    n2 = (
                     -n[0], -n[1], -n[2])
                    MIL[n] = sumL[n] / float(nL)
                    MIL[n2] = sumL[n] / float(nL)
                    SLD[n] = sumL2[n] / sumL[n]
                    SLD[n2] = sumL2[n] / sumL[n]
                    SVD[n] = numpy.pi / 3.0 * sumL4[n] / sumL[n]
                    SVD[n2] = numpy.pi / 3.0 * sumL4[n] / sumL[n]
                else:
                    print '**WARNING: computeOrigDistribution_MIL(): number of intersection nL=0'

        if echo == True:
            time2 = clock()
            ctime2 = time()
            stdout.write('     -> Compute in CPU/TOT  :   %8.1f/%8.1f sec\n' % (time2 - time1, ctime2 - ctime1))
            stdout.flush()
        return (MIL, SVD, SLD, area)

    def set_nx_ny_nz(self, voxelModel):
        nx, ny, nz = self.get_Shape(voxelModel)
        self.nx = nx
        self.ny = ny
        self.nz = nz

    def computeOrigDistribution_F77(self, thresVoxelModel, step, power, valid, echo=True):
        """
        Function computes MIL/SLD/SVD distributions for direction vectors "n" 
        using a voxel ray field going trought the RVE. Normals n = tuple(nix,niy,niz) 
        are the directions from the midpoint of a unit sphere to the COG of triangles 
        which build the surface of the sphere. Very similar to 
        self.computeOrigDistribution_STAR(). 
        A segement voxel model with isotropic resolution is needed.   
        
        @param thresVoxelModel: segmented voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue 
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.       
          - int grayValue  ... gray value of voxel, 0..255                  
        @param step: Step in the considered voxel 
          - TYPE: int > 0 
        @param power: Power for the number of star directions = 8*4^power
          - TYPE: int > 1
        @param valid: Smallest valid intersection length.
          - TYPE: float
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
                                     
        @return: 
             MIL             : Mean Intercept Length 
               - TYPE: dict[ (nix, niy, niz) ] = value 
               - float nix, niy, niz ... components of normal vectors 
               - float value         ... MIL for that direction          
             SLD             : Star Length Distribution 
               - TYPE: same as for MIL 
             SVD             : Star Volume Distribution 
               - TYPE: same as for MIL       
             Area            : Weight (Area of triange on unit sphere) for each direction 
               - TYPE: same as for MIL          
        """
        if echo == True:
            stdout.write(' ... compute original distribution F77\n')
            stdout.flush()
        time1 = clock()
        ctime1 = time()
        nx, ny, nz = self.get_Shape(thresVoxelModel)
        self.nx = nx
        self.ny = ny
        self.nz = nz
        normals, area = self.computeNormalAndArea(power)
        posNormals = len(normals) / 8
        arrNormals = numpy.zeros((3, posNormals), numpy.float)
        j = 0
        for n in normals:
            nX = n[0]
            nY = n[1]
            nZ = n[2]
            if nX >= 0.0 and nY >= 0.0 and nZ >= 0.0:
                arrNormals[0, j] = nX
                arrNormals[1, j] = nY
                arrNormals[2, j] = nZ
                j += 1

        viewerAt, orgF77MIL, orgF77SVD, orgF77SLD = miaf77.computeorigdistribution_f77(thresVoxelModel, arrNormals, step, valid)
        MIL = {}
        SLD = {}
        SVD = {}
        for i in range(4):
            for j in range(posNormals):
                nX = arrNormals[0, j] * viewerAt[0, i]
                nY = arrNormals[1, j] * viewerAt[1, i]
                nZ = arrNormals[2, j] * viewerAt[2, i]
                n = (nX, nY, nZ)
                n2 = (-n[0], -n[1], -n[2])
                MIL[n] = orgF77MIL[j, i]
                MIL[n2] = orgF77MIL[j, i]
                SLD[n] = orgF77SLD[j, i]
                SLD[n2] = orgF77SLD[j, i]
                SVD[n] = orgF77SVD[j, i]
                SVD[n2] = orgF77SVD[j, i]

        if echo == True:
            time2 = clock()
            ctime2 = time()
            stdout.write('     -> Compute in CPU/TOT  :   %8.1f/%8.1f sec\n' % (time2 - time1, ctime2 - ctime1))
            stdout.flush()
        return (MIL, SVD, SLD, area)

    def readOrigDistribution(self, infile):
        """
        Reads previous computed MIL/SLD/SVD distributions from a file. Return 
        values a same as in function L{self.computeOrigDistribution_MIL()} and 
        self.computeOrigDistribution_STAR()
              
        @param infile :  name of input file ::
        
            file format :
                  n1x   n1y   n1z   value1
                  ... 
                  nix   niy   niz   value  
                  ...           
                  nNx   nNy   nNz   valueN         
                                    
        @return: 
             MIL             : Mean Intercept Length 
               - TYPE: dict[ (nix, niy, niz) ] = value 
               - float nix, niy, niz ... components of normal vectors 
               - float value         ... MIL for that direction          
             SLD             : Star Length Distribution 
               - TYPE: same as for MIL 
             SVD             : Star Volume Distribution 
               - TYPE: same as for MIL       
             Area            : Weight (Area of triange on unit sphere) for each direction 
               - TYPE: same as for MIL          
        """
        stdout.write(' ... read original distribution\n')
        stdout.flush()
        DIST_n = {}
        area_n = {}
        try:
            data = open(infile)
        except IOError:
            stdout.write("\n **ERROR**: '-in' intput file '%s' not found!\n\n" % infile)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        dataList = ['None']
        ndir = 0
        while len(dataList) > 0:
            line = data.readline()
            dataList = line.split()
            if len(dataList) > 0:
                ndir += 2
                nx = float(dataList[0])
                ny = float(dataList[1])
                nz = float(dataList[2])
                length = numpy.sqrt(nx * nx + ny * ny + nz * nz)
                nx = nx / length
                ny = ny / length
                nz = nz / length
                value = float(dataList[3])
                DIST_n[nx, ny, nz] = value
                area_n[nx, ny, nz] = 1.0
                DIST_n[-nx, -ny, -nz] = value

        AKugelSeg = 4.0 * numpy.pi / float(ndir)
        for direct in DIST_n:
            area_n[direct] = AKugelSeg

        return (
         DIST_n, area_n)

    def computeEigenValueAndVector(self, origDistr, ftype, Area_n=None, normType='det', echo=True):
        """
        Function computes the eigenvalue and eigenvectors by fitting an ellipsoid 
        or a general second order tensor through a given distribution. 
        If ftype = 4 only the second order part of the tensor will be used.  
        
        @param  origDistr: Original distribution 
          - TYPE: dict[ (nix, niy, niz) ] = value 
          - float nix, niy, niz ... components of normal vector 
          - float value         ... value for that direction             
        @param  Area_n: Weight (Area of triange on unit sphere) for each direction. 
          Not needed for ftype = 1. 
            - TYPE: same as for origDistr
        @param  ftype: type of approximation: 
          - TYPE: int = 1,2,4     
          - ftype = 1 : Ellipsoid
          - ftype = 2 : general tensor 2nd order 
          - ftype = 4 : general tensor 4th order
        @param  normType: Type of fabric normalization 'det' or 'rank' / 'trace'
            - TYPE: string
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
               
        @return:      
          evalue          : Eigenvalues of fabric tensor 
            - TYPE: numpy.array[evalID] = eval
            - int evalID ... eigenvalues ID. 0,1, or 2
            - float eval ... current eigenvalue 
          evector         : Eigenvectors of fabric tensor
            - TYPE: numpy.array[evectID, evalID] = evect
            - int evectID ... eigenvector ID. 0,1, or 2                         
            - int evalID  ... eigenvalue ID. 0,1, or 2
            - float evect ... component of eigenvectors, e.g. evector[0,2] = ev1_z
        """
        if echo == True:
            stdout.write(' ... compute eigenvalues and vectors\n')
            stdout.flush()
        time1 = clock()
        ftype = int(ftype)
        M = numpy.array(numpy.zeros((3, 3), numpy.float))
        if ftype == 2 or ftype == 4:
            if ftype == 4:
                stdout.write('\n **WARNING** : mia.computeEigenValueAndVector() only second order tensor will be used!\n')
                stdout.flush()
            M = self.computeFabricTensorsSpectral(origDistr, Area_n)
            evalue, evector = linalg.eig(M)
        elif ftype == 1:
            M = self.computeFabricTensorsEllipsoid(origDistr)
            evalue, evector = linalg.eig(M)
            evalue[0] = 1.0 / numpy.sqrt(evalue[0])
            evalue[1] = 1.0 / numpy.sqrt(evalue[1])
            evalue[2] = 1.0 / numpy.sqrt(evalue[2])
        else:
            stdout.write('\n **ERROR** : mia.computeEigenValueAndVector() ftype = %2i not known!' % ftype)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        normVal = 0
        if normType.lower() == 'det':
            normVal = (evalue[0] * evalue[1] * evalue[2]) ** (1.0 / 3.0)
        elif normType.lower() == 'rank' or normType.lower() == 'trace':
            normVal = (evalue[0] + evalue[1] + evalue[2]) / 3.0
        else:
            stdout.write('\n **ERROR** : mia.computeEigenValueAndVector() normType = %s not known!' % normType)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        evalue[0] = evalue[0] / normVal
        evalue[1] = evalue[1] / normVal
        evalue[2] = evalue[2] / normVal
        return (
         evalue, evector)

    def computeApproxDistribution(self, origDistr, triaArray, ftype, Area_n=None, nameOfDist=None, normType='det'):
        """
        Function computes an approximated distribution at triangle nodes 
        (see L{self.computeOrigDistribution_MIL()} for more details)
        by using an original MIL/SLD/SVD distribution at COG of that triangles.
        Function used for graphical output in PARAVIEW. 
        
        @param  origDistr: Original distribution 
          - TYPE: dict[ (nix, niy, niz) ] = value 
          - float nix, niy, niz ... components of normal vector 
          - float value         ... value for that direction             
        @param  Area_n: Weight (Area of triange on unit sphere) for each direction. 
          Not needed for ftype = 1. 
            - TYPE: same as for origDistr
        @param  triaArray: List of triangles 
          - TYPE: list[ ntriangles ] = ( (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ) 
          - ntriangles ... number of triangles 
          - float xi, yi, zi ... x,y,z coordinates of point i    
        @param  ftype: type of approximation: 
          - TYPE: int = 1,2,4     
          - ftype = 1 : Ellipsoid
          - ftype = 2 : general tensor 2nd order 
          - ftype = 4 : general tensor 4th order                 
        @param nameOfDist: name of original distribution, needed for tdout output echo 
          - TYPE: string  
        @param  normType: Type of fabric normalization 'det' or 'rank' / 'trace'
          - TYPE: string
          
                     
        @return:      
             approxDistr: approximated unscaled distribution value at triangle nodes (not COG)
               - TYPE: dict[ (nix,niy,niz) ] = value 
               - float nix, niy, niz ... components of normal vectors at nodal points.
               - float value         ... MIL/SLD/SVD value for that nodal point.
             scale: scale factor such that det(M) = 1 (M..fabric tensor)
               - TYPE: float    
        """
        if nameOfDist == None:
            stdout.write(' ... compute distribution, type=%2i\n, norm=%s' % (ftype, normType))
            stdout.flush()
        else:
            stdout.write(' ... compute %s distribution, type=%2i, norm=%s\n' % (nameOfDist, ftype, normType))
            stdout.flush()
        time1 = clock()
        ftype = int(ftype)
        if ftype == 2 or ftype == 4:
            g, G, G4 = self.computeSpectralDecomposition(origDistr, Area_n)
            M = self.computeFabricTensorsSpectral(origDistr, Area_n)
            evalue, evector = linalg.eig(M)
        if ftype == 1:
            H = self.computeFabricTensorsEllipsoid(origDistr)
            evalue, evector = linalg.eig(H)
            evalue[0] = 1.0 / numpy.sqrt(evalue[0])
            evalue[1] = 1.0 / numpy.sqrt(evalue[1])
            evalue[2] = 1.0 / numpy.sqrt(evalue[2])
        normVal = 0
        if normType.lower() == 'det':
            normVal = (evalue[0] * evalue[1] * evalue[2]) ** (1.0 / 3.0)
        elif normType.lower() == 'rank' or normType.lower() == 'trace':
            normVal = (evalue[0] + evalue[1] + evalue[2]) / 3.0
        else:
            stdout.write('\n **ERROR** : mia.computeEigenValueAndVector() normType = %s not known!' % normType)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        scale = 1.0 / normVal
        F4 = numpy.array(numpy.zeros((3, 3, 3, 3), numpy.float))
        fabricNodes = {}
        I = dpTensor.UnitMatrix(3)
        approxDistr = {}
        checkdouble = {}
        for tria in triaArray:
            for point in tria:
                if not checkdouble.has_key(point):
                    checkdouble[point] = True
                    n = numpy.array(point)
                    N = dpTensor.DyadicProduct(n, n)
                    if ftype == 1:
                        L = 1.0 / numpy.sqrt(dpTensor.DoubleContraction(N, H))
                    else:
                        F = N - 1.0 / 3.0 * I
                        if ftype == 4:
                            F4 = dpTensor.DyadicProduct(N, N) - 1.0 / 7.0 * (dpTensor.DyadicProduct(I, N) + dpTensor.DyadicProduct(N, I)) - 2.0 / 7.0 * (dpTensor.SymmetricProduct(I, N) + dpTensor.SymmetricProduct(N, I)) + 1.0 / 35.0 * dpTensor.DyadicProduct(I, I) + 2.0 / 35.0 * dpTensor.SymmetricProduct(I, I)
                        L = 0.0
                        if int(ftype) == 2:
                            L = g + dpTensor.DoubleContraction(G, F)
                        else:
                            L = g + dpTensor.DoubleContraction(G, F) + dpTensor.DoubleContraction(G4, F4)
                    approxDistr[point] = L

        time2 = clock()
        stdout.write('     -> Compute in          :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        return (
         approxDistr, scale)

    def computeApproxDistributionGST(self, M, triaArray, ftype, normalize):
        """
        Same as computeApproxDistribution but by directly given the Fabric tensor M
        """
        stdout.write(' ... compute distribution GST\n')
        stdout.flush()
        time1 = clock()
        ftype = int(ftype)
        EVAL, EVEC = linalg.eig(M)
        normVal = 1.0
        if normalize.lower() == 'det':
            normVal = (EVAL[0] * EVAL[1] * EVAL[2]) ** (1.0 / 3.0)
        elif normalize == 'rank' or normalize == 'trace':
            normVal = (EVAL[0] + EVAL[1] + EVAL[2]) / 3.0
        scale = 1.0 / normVal
        EVAL[0] = EVAL[0] ** (-2)
        EVAL[1] = EVAL[1] ** (-2)
        EVAL[2] = EVAL[2] ** (-2)
        EVALM = numpy.zeros((3, 3), numpy.float)
        EVALM[(0, 0)] = EVAL[0]
        EVALM[(1, 1)] = EVAL[1]
        EVALM[(2, 2)] = EVAL[2]
        H = numpy.dot(numpy.dot(EVEC, EVALM), numpy.transpose(EVEC))
        I = dpTensor.UnitMatrix(3)
        approxDistr = {}
        checkdouble = {}
        for tria in triaArray:
            for point in tria:
                if not checkdouble.has_key(point):
                    checkdouble[point] = True
                    n = numpy.array(point)
                    N = dpTensor.DyadicProduct(n, n)
                    L = 1.0 / numpy.sqrt(dpTensor.DoubleContraction(N, H))
                    approxDistr[point] = L

        time2 = clock()
        stdout.write('     -> Compute in          :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        return (
         approxDistr, scale)

    def computeSpectralDecomposition(self, fN, A):
        """ 
        This function is called by self.computeApproxDistribution() and 
        L{self.computeFabricTensorsSpectral()} in order to compute 
        the spectral decomposition of a 4th order tensor which is fitted through fN.  
        
         @param  fN: Original distribution 
           - TYPE: dict[ (nix, niy, niz) ] = value 
           - float nix, niy, niz ... components of normal vector 
           - float value         ... value for that direction             
         @param  A: Weight (Area of triange on unit sphere) for each direction. 
           Not needed for ftype = 1. 
             - TYPE: same as for origDistr     
                      
         @return:      
           g: first order tensor from spectral decomposition 
             - TYPE: float 
           G: second order tensor from spectral decomposition 
             - TYPE: float numpy.array[3,3]   
           G4: fourth order tensor from spectral decomposition  
             - TYPE: float numpy.array[3,3,3,3]       
        """
        g = 0.0
        G = numpy.zeros([3, 3], numpy.float)
        G4 = numpy.zeros([3, 3, 3, 3], numpy.float)
        F4 = numpy.array(numpy.zeros((3, 3, 3, 3), numpy.float))
        I = dpTensor.UnitMatrix(3)
        for di in fN:
            n = numpy.array(di)
            g += fN[di] * A[di]
            N = dpTensor.DyadicProduct(n, n)
            F = N - 1.0 / 3.0 * I
            G += fN[di] * F * A[di]
            F4 = dpTensor.DyadicProduct(N, N) - 1.0 / 7.0 * (dpTensor.DyadicProduct(I, N) + dpTensor.DyadicProduct(N, I)) - 2.0 / 7.0 * (dpTensor.SymmetricProduct(I, N) + dpTensor.SymmetricProduct(N, I)) + 1.0 / 35.0 * dpTensor.DyadicProduct(I, I) + 2.0 / 35.0 * dpTensor.SymmetricProduct(I, I)
            G4 += fN[di] * F4 * A[di]

        G4 = G4 * (315.0 / (32.0 * numpy.pi))
        G = G * (15.0 / (8.0 * numpy.pi))
        g = g / (4.0 * numpy.pi)
        return (
         g, G, G4)

    def computeFabricTensorsSpectral(self, fN, A):
        """ 
        This function is called by L{self.computeEigenValueAndVector()} in order to 
        compute the fabric tensors using spectral decomposition. M = I*g + G 
        
         @param  fN: Original distribution 
           - TYPE: dict[ (nix, niy, niz) ] = value 
           - float nix, niy, niz ... components of normal vector 
           - float value         ... value for that direction             
         @param  A: Weight (Area of triange on unit sphere) for each direction. 
           Not needed for ftype = 1. 
             - TYPE: same as for origDistr              
                      
         @return:      
           M: fabric tensor from ellipsoidal fit 
             - TYPE: float numpy.array[3,3]   
        """
        g, G, G4 = self.computeSpectralDecomposition(fN, A)
        I = dpTensor.UnitMatrix(3)
        M = I * g + G
        return M

    def computeFabricTensorsEllipsoid(self, fN):
        """ 
        This function is called by self.computeApproxDistribution() and 
        self.computeEigenValueAndVector() in order to compute the fabric tensors 
        using an ellipsoidal fit. 
        
         @param  fN: Original distribution 
           - TYPE: dict[ (nix, niy, niz) ] = value 
           - float nix, niy, niz ... components of normal vector 
           - float value         ... value for that direction             
                      
         @return:      
           M: fabric tensor from ellipsoidal fit 
             - TYPE: float numpy.array[3,3]            
        """
        ndir = len(fN)
        nhat = numpy.array(numpy.zeros((ndir, 6), numpy.float))
        An = numpy.array(numpy.zeros(ndir, numpy.float))
        H = numpy.array(numpy.zeros((3, 3), numpy.float))
        d = 0
        for n in fN:
            nhat[d, 0] = n[0] * n[0]
            nhat[d, 1] = n[1] * n[1]
            nhat[d, 2] = n[2] * n[2]
            nhat[d, 3] = numpy.sqrt(2.0) * n[1] * n[2]
            nhat[d, 4] = numpy.sqrt(2.0) * n[2] * n[0]
            nhat[d, 5] = numpy.sqrt(2.0) * n[0] * n[1]
            An[d] = 1.0 / fN[n] * (1.0 / fN[n])
            d += 1

        N1 = numpy.dot(numpy.transpose(nhat), nhat)
        N2 = numpy.dot(numpy.transpose(nhat), An)
        VM = numpy.dot(linalg.inv(N1), N2)
        H[(0, 0)] = VM[0]
        H[(1, 1)] = VM[1]
        H[(2, 2)] = VM[2]
        H[(1, 2)] = VM[3] / numpy.sqrt(2.0)
        H[(2, 0)] = VM[4] / numpy.sqrt(2.0)
        H[(0, 1)] = VM[5] / numpy.sqrt(2.0)
        H[(2, 1)] = VM[3] / numpy.sqrt(2.0)
        H[(0, 2)] = VM[4] / numpy.sqrt(2.0)
        H[(1, 0)] = VM[5] / numpy.sqrt(2.0)
        return H

    def computeFabricTensorGST(self, voxelModel, power, normalize, threshold=None, gst=None, echo=None):
        """ 
         This function computes the fabric tensor based on the 
         Gradient structure tensor 
        """
        if echo == 0:
            pass
        else:
            stdout.write(' ... compute GST based fabric\n')
            stdout.flush()
        if normalize.lower() == 'det' or normalize == 'rank' or normalize == 'trace':
            pass
        else:
            stdout.write('\n **ERROR** : mia.computeFabricTensorGST() -norm = %s not known!' % normalize)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if gst == None:
            if threshold:
                gst = miaf77.compute_gst_thres(voxelModel, threshold, 2)
            else:
                gst = miaf77.compute_gst(voxelModel, 2)
        if numpy.abs(linalg.det(gst)) < 1e-07:
            gst = numpy.eye(3)
            EVAL, EVEC = linalg.eig(gst)
        else:
            EVAL, EVEC = linalg.eig(gst)
        EVAL[0] = EVAL[0] ** (-0.5)
        EVAL[1] = EVAL[1] ** (-0.5)
        EVAL[2] = EVAL[2] ** (-0.5)
        normVal = 1.0
        if normalize.lower() == 'det':
            normVal = (EVAL[0] * EVAL[1] * EVAL[2]) ** (1.0 / 3.0)
        elif normalize == 'rank' or normalize == 'trace':
            normVal = (EVAL[0] + EVAL[1] + EVAL[2]) / 3.0
        EVAL[0] = EVAL[0] / normVal
        EVAL[1] = EVAL[1] / normVal
        EVAL[2] = EVAL[2] / normVal
        EVAL[0] = EVAL[0] ** power
        EVAL[1] = EVAL[1] ** power
        EVAL[2] = EVAL[2] ** power
        normVal = 1.0
        if normalize.lower() == 'det':
            normVal = (EVAL[0] * EVAL[1] * EVAL[2]) ** (1.0 / 3.0)
        elif normalize == 'rank' or normalize == 'trace':
            normVal = (EVAL[0] + EVAL[1] + EVAL[2]) / 3.0
        EVAL[0] = EVAL[0] / normVal
        EVAL[1] = EVAL[1] / normVal
        EVAL[2] = EVAL[2] / normVal
        EVALM = numpy.zeros((3, 3), numpy.float)
        EVALM[(0, 0)] = EVAL[0]
        EVALM[(1, 1)] = EVAL[1]
        EVALM[(2, 2)] = EVAL[2]
        M_hat = numpy.dot(numpy.dot(EVEC, EVALM), numpy.transpose(EVEC))
        return M_hat

    def computeNormalAndArea(self, power):
        """
        Computes the normals at COG and area (weight) of a triangulated unit 
        sphere. Called from L{self.computeOrigDistribution_STAR()} and  
        L{self.computeOrigDistribution_MIL()}.        
        
        @param power: Parameter for number of triangles on unit sphere 
                      (No of triangles = 8*4^power). 
                        - TYPE: int
                     
        @return:      
             normals: normals from COG with unit length 
               - TYPE: list[ (nix, niy, niz) ]  
               - float nix, niy, niz ... components of normal vectors 
                                  
             Area_n: area of triangles which build the surface of sphere    
               - TYPE: dict[ (nix, niy, niz) ] = value    
               - float nix, niy, niz ... components of normal vectors 
               - float value         ... Area for that direction                         
        """
        triaArray = self.setupSphereTriangles(power, True)
        normals = []
        Area_n = {}
        Asum = 0.0
        for tria in triaArray:
            A, COG = self.computeAreaAndCOG(tria)
            normals.append(COG)
            Area_n[COG] = A
            Asum += A

        k = 4.0 * numpy.pi / Asum
        for n in Area_n:
            Area_n[n] = Area_n[n] * k

        return (
         normals, Area_n)

    def computeLength(self, faceInfo, point, direc, nx, ny, nz):
        """ 
        Called from L{self.computeOrigDistribution_STAR()} in order to compute 
        the lengths of a "star".        
        """
        xg = point[0] + 0.5
        yg = point[1] + 0.5
        zg = point[2] + 0.5
        ax = direc[0]
        ay = direc[1]
        az = direc[2]
        l1arr = []
        bigL = numpy.sqrt(nx * nx + ny * ny + nz * nz)
        xF = 0.0
        yF = 0.0
        zF = 0.0
        if ax > 0.0 or ax < 0.0:
            start = -1
            end = -1
            if ax > 0.0:
                start = point[0]
                end = nx
            elif ax < 0.0:
                start = 0
                end = point[0] + 1
            for ii in range(start, end):
                if ax > 0.0:
                    i = ii
                elif ax < 0.0:
                    i = point[0] - ii - 1
                t = -(-(i + 1) + xg) / ax
                fy = int(t * ay + yg)
                fz = int(t * az + zg)
                if fz < nz and fz >= 0 and fy < ny and fy >= 0:
                    if faceInfo['x'][fz, fy, i] == 0 or i == -1:
                        xF = t * ax
                        yF = t * ay
                        zF = t * az
                        break
                else:
                    xF = bigL + 10.0
                    yF = bigL + 10.0
                    zF = bigL + 10.0
                    break

        else:
            xF = bigL + 10.0
            yF = bigL + 10.0
            zF = bigL + 10.0
        xlength1 = numpy.sqrt(xF * xF + yF * yF + zF * zF)
        xF = 0.0
        yF = 0.0
        zF = 0.0
        if ay > 0.0 or ay < 0.0:
            start = -1
            end = -1
            if ay > 0.0:
                start = point[1]
                end = ny
            elif ay < 0.0:
                start = 0
                end = point[1] + 1
            for ii in range(start, end):
                if ay > 0.0:
                    i = ii
                elif ay < 0.0:
                    i = point[1] - ii - 1
                t = -(-(i + 1) + yg) / ay
                fx = int(t * ax + xg)
                fz = int(t * az + zg)
                if fz < nz and fz >= 0 and fx < nx and fx >= 0:
                    if faceInfo['y'][fz, i, fx] == 0 or i == -1:
                        xF = t * ax
                        yF = t * ay
                        zF = t * az
                        break
                else:
                    xF = bigL + 10.0
                    yF = bigL + 10.0
                    zF = bigL + 10.0
                    break

        else:
            xF = bigL + 10.0
            yF = bigL + 10.0
            zF = bigL + 10.0
        ylength1 = numpy.sqrt(xF * xF + yF * yF + zF * zF)
        xF = 0.0
        yF = 0.0
        zF = 0.0
        if az > 0.0 or az < 0.0:
            start = -1
            end = -1
            if az > 0.0:
                start = point[2]
                end = nz
            elif az < 0.0:
                start = 0
                end = point[2] + 1
            for ii in range(start, end):
                if az > 0.0:
                    i = ii
                elif az < 0.0:
                    i = point[2] - ii - 1
                t = -(-(i + 1) + zg) / az
                fx = int(t * ax + xg)
                fy = int(t * ay + yg)
                if fy < ny and fy >= 0 and fx < nx and fx >= 0:
                    if faceInfo['z'][i, fy, fx] == 0 or i == -1:
                        xF = t * ax
                        yF = t * ay
                        zF = t * az
                        break
                else:
                    xF = bigL + 10.0
                    yF = bigL + 10.0
                    zF = bigL + 10.0
                    break

        else:
            xF = bigL + 10.0
            yF = bigL + 10.0
            zF = bigL + 10.0
        zlength1 = numpy.sqrt(xF * xF + yF * yF + zF * zF)
        l1arr = []
        l1arr.append(xlength1)
        l1arr.append(ylength1)
        l1arr.append(zlength1)
        minL = min(l1arr)
        if minL == 0.0 or minL > bigL:
            stdout.write('\n **ERROR** computeLength(): SVD/SLD algorithm found no valid length!\n\n')
            stdout.flush()
            stdout.write('   Point     %s: \n' % repr(point))
            stdout.flush()
            stdout.write('   Direction %s: \n' % repr(direc))
            stdout.flush()
            stdout.write('   Length    %s: \n' % repr(minL))
            stdout.flush()
            stdout.write('   xL1=%7.3f   \n' % l1arr[0])
            stdout.flush()
            stdout.write('   yL1=%7.3f   \n' % l1arr[1])
            stdout.flush()
            stdout.write('   zL1=%7.3f   \n' % l1arr[2])
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return minL

    def getFaceInfo(self, thresVoxelModel):
        """ 
        Called from L{self.computeOrigDistribution_STAR()} in order to setup 
        the information which face between voxels is a bone-marrow interface. 
        Not used for L{self.computeOrigDistribution_MIL()}.      
        """
        stdout.write('     -> set up face info    :')
        stdout.flush()
        stdout.flush()
        nx, ny, nz = self.get_Shape(thresVoxelModel)
        xFace = numpy.zeros((nz, ny, nx), numpy.uint8)
        yFace = numpy.zeros((nz, ny, nx), numpy.uint8)
        zFace = numpy.zeros((nz, ny, nx), numpy.uint8)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    if i + 1 < nx:
                        if thresVoxelModel[k, j, i] == thresVoxelModel[k, j, i + 1]:
                            xFace[k, j, i] = 1
                    if j + 1 < ny:
                        if thresVoxelModel[k, j, i] == thresVoxelModel[k, j + 1, i]:
                            yFace[k, j, i] = 1
                    if k + 1 < nz:
                        if thresVoxelModel[k, j, i] == thresVoxelModel[k + 1, j, i]:
                            zFace[k, j, i] = 1

        faceInfo = {}
        faceInfo['x'] = xFace
        faceInfo['y'] = yFace
        faceInfo['z'] = zFace
        stdout.write('       done  \n')
        stdout.flush()
        return faceInfo

    def setupSphereTriangles(self, nDirs, correct=False):
        """ 
         Called from L{self.computeNormalAndArea()} and from the main program in order to 
         setup a mesh for a unit sphere. 
         
         @param nDirs: Parameter for number of triangles on unit sphere 
                       (No of triangles = 8*4^power). 
                         - TYPE: int          
         @param correct: Flag in order to prevent that the direction of 
                         n has angles of 45 deg. If this is the case the 
                         direction will be disturbed slightly. 
                           - TYPE: bool 
                      
         @return:      
           triaArray: List of triangles 
             - TYPE: list[ (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ]
             - float xi, yi, zi ... x,y,z coordinates of rtriangle corner point i                          
        """
        triaArray = []
        triaArray.append(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
        for cDir in range(int(nDirs)):
            newTriaArray = []
            for tria in triaArray:
                nTri1, nTri2, nTri3, nTri4 = self.splitTriangle(tria)
                newTriaArray.append(nTri1)
                newTriaArray.append(nTri2)
                newTriaArray.append(nTri3)
                newTriaArray.append(nTri4)

            triaArray = newTriaArray

        newTriaArray = []
        for tria in triaArray:
            nP1 = self.projectionToUnitSphereM2(tria[0], correct)
            nP2 = self.projectionToUnitSphereM2(tria[1], correct)
            nP3 = self.projectionToUnitSphereM2(tria[2], correct)
            nTr = (nP1, nP2, nP3)
            newTriaArray.append(nTr)

        triaArray = newTriaArray
        newTriaArray2 = []
        for tria in triaArray:
            newTriaArray2.append(tria)
            newTriaArray2.append(((tria[0][0], -tria[0][1], tria[0][2]), (tria[1][0], -tria[1][1], tria[1][2]), (tria[2][0], -tria[2][1], tria[2][2])))

        triaArray = newTriaArray2
        newTriaArray3 = []
        for tria in triaArray:
            newTriaArray3.append(tria)
            newTriaArray3.append(((-tria[0][0], tria[0][1], tria[0][2]), (-tria[1][0], tria[1][1], tria[1][2]), (-tria[2][0], tria[2][1], tria[2][2])))

        triaArray = newTriaArray3
        newTriaArray4 = []
        for tria in triaArray:
            newTriaArray4.append(tria)
            newTriaArray4.append(((tria[0][0], tria[0][1], -tria[0][2]), (tria[1][0], tria[1][1], -tria[1][2]), (tria[2][0], tria[2][1], -tria[2][2])))

        triaArray = newTriaArray4
        return triaArray

    def setupSphereFem(self, triaArray, distribut=None):
        """ 
         Setup unit sphere FEM mesh (delete double nodes) from a "triagle Array" 
         and map nodal scalar values if required. 
         Called from L{self.writeParaview()}. 
         
         @param triaArray: List of triangles 
           - TYPE: list[ (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ]
           - float xi, yi, zi ... x,y,z coordinates of rtriangle corner point i       
         @param distribut: Flag for mapping scalar values on nodes. The value is 
                           the distance from Point(0,0,0) to the current nodal point.  
                             - TYPE: bool 
                      
         @return:      
              sphereNodes: Nodes of the mesh 
                - TYPE: dict[ nodeId ] = node 
                - int nodeId      : node identifier  
                - dpFem.node node : FEM node 
              sphereElems: Elements of the mesh 
                - TYPE: dict[ elemId ] = element 
                - int elemId            : element identifier  
                - dpFem.element element : FEM element object 
              NscaResults: Scalar results on the nodes 
                - TYPE: dict[ nodeId ] = value 
                - int nodeId  : node identifier  
                - float value : nodal scalar value                                         
        """
        uniqueNodes = {}
        sphereNodes = {}
        sphereElems = {}
        NscaResults = {}
        nid = 0
        elid = 0
        for tria in triaArray:
            elid += 1
            cn = [0, 0, 0]
            for lnoid in range(3):
                if uniqueNodes.has_key(tria[lnoid]):
                    cn[lnoid] = uniqueNodes[tria[lnoid]]
                else:
                    nid += 1
                    uniqueNodes[tria[lnoid]] = nid
                    mul = 1.0
                    if not distribut == None:
                        mul = distribut[tria[lnoid]]
                    xc = float(tria[lnoid][0]) * mul
                    yc = float(tria[lnoid][1]) * mul
                    zc = float(tria[lnoid][2]) * mul
                    curNode = dpFem.node(nid, xc, yc, zc)
                    sphereNodes[nid] = curNode
                    cn[lnoid] = nid
                    NscaResults[nid] = mul

            nList = [cn[0], cn[1], cn[2]]
            curElem = dpFem.element(elid, nList, 'tria3')
            sphereElems[elid] = curElem

        return (
         sphereNodes, sphereElems, NscaResults)

    def projectionToUnitSphereM2(self, PointRS, correct=False):
        """
        Projects an equally sided triangle patch to a unit sphere
        """
        sc45 = numpy.sin(numpy.pi / 4.0)
        XYZ = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.0, 0.0), (sc45, sc45, 0.0), (0.0, 0.5, 0.0), (sc45, 0.0, sc45), (0.0, sc45, sc45), (0.0, 0.0, 0.5)]
        r = PointRS[0]
        s = PointRS[1]
        t = PointRS[2]
        N5 = 4.0 * r * (1.0 - r - s - t)
        N6 = 4.0 * r * s
        N7 = 4.0 * s * (1.0 - r - s - t)
        N8 = 4.0 * r * t
        N9 = 4.0 * s * t
        N10 = 4.0 * t * (1.0 - r - s - t)
        N1 = 1.0 - r - s - t - 0.5 * N5 - 0.5 * N7 - 0.5 * N10
        N2 = r - 0.5 * N5 - 0.5 * N6 - 0.5 * N8
        N3 = s - 0.5 * N6 - 0.5 * N7 - 0.5 * N9
        N4 = t - 0.5 * N8 - 0.5 * N9 - 0.5 * N10
        aN = [
         N1, N2, N3, N4, N5, N6, N7, N8, N9, N10]
        x = 0.0
        y = 0.0
        z = 0.0
        for node in range(10):
            x += XYZ[node][0] * aN[node]
            y += XYZ[node][1] * aN[node]
            z += XYZ[node][2] * aN[node]

        if correct == True:
            x, y, z = self.correctValues(x, y, z)
        fac = 1.0 / numpy.sqrt(x * x + y * y + z * z)
        return (
         fac * x, fac * y, fac * z)

    def correctValues(self, x, y, z, prec=1e-06):
        """
        Ensure that current direction do not go through corner or edge  
        i.e. has an angle of 45 deg.
        """
        prec = 1e-06
        if abs(int(x / prec)) == abs(int(y / prec)) and abs(int(x / prec)) == abs(int(z / prec)):
            x += prec
            z += 2.0 * prec
        elif abs(int(x / prec)) == abs(int(y / prec)):
            x += prec
        elif abs(int(x / prec)) == abs(int(z / prec)):
            x += prec
        elif abs(int(z / prec)) == abs(int(y / prec)):
            z += prec
        arr = []
        arr.append(self.nx)
        arr.append(self.ny)
        arr.append(self.nz)
        nmax = max(arr)
        for n in range(nmax):
            if abs(int(x / prec)) == abs((n + 1) * int(y / prec)):
                x += prec
                stdout.write('\n **WARNING** projectionToUnitSphereM2(): line goes through edge - check!\n\n')
                stdout.flush()
            elif abs(int(x / prec)) == abs((n + 1) * int(z / prec)):
                x += prec
                stdout.write('\n **WARNING** projectionToUnitSphereM2(): line goes through edge - check!\n\n')
                stdout.flush()
            elif abs(int(z / prec)) == abs((n + 1) * int(y / prec)):
                z += prec
                stdout.write('\n **WARNING** projectionToUnitSphereM2(): line goes through edge - check!\n\n')
                stdout.flush()
            elif abs((n + 1) * int(x / prec)) == abs(int(y / prec)):
                x += prec
                stdout.write('\n **WARNING** projectionToUnitSphereM2(): line goes through edge - check!\n\n')
                stdout.flush()
            elif abs((n + 1) * int(x / prec)) == abs(int(z / prec)):
                x += prec
                stdout.write('\n **WARNING** projectionToUnitSphereM2(): line goes through edge - check!\n\n')
                stdout.flush()
            elif abs((n + 1) * int(z / prec)) == abs(int(y / prec)):
                z += prec
                stdout.write('\n **WARNING** projectionToUnitSphereM2(): line goes through edge - check!\n\n')
                stdout.flush()

        return (
         x, y, z)

    def splitTriangle(self, Tri):
        """ 
        Splits one triange into four triangles. 
        """
        P1 = Tri[0]
        P2 = Tri[1]
        P3 = Tri[2]
        P1x = P1[0]
        P1y = P1[1]
        P1z = P1[2]
        P2x = P2[0]
        P2y = P2[1]
        P2z = P2[2]
        P3x = P3[0]
        P3y = P3[1]
        P3z = P3[2]
        P4 = ((P1x + P2x) / 2, (P1y + P2y) / 2, (P1z + P2z) / 2)
        P5 = ((P3x + P2x) / 2, (P3y + P2y) / 2, (P3z + P2z) / 2)
        P6 = ((P1x + P3x) / 2, (P1y + P3y) / 2, (P1z + P3z) / 2)
        nTri1 = (
         P1, P4, P6)
        nTri2 = (P4, P2, P5)
        nTri3 = (P4, P5, P6)
        nTri4 = (P5, P3, P6)
        return (
         nTri1, nTri2, nTri3, nTri4)

    def computeAreaAndCOG(self, tria):
        """ 
        Computes area and center of gravity of a triangle. Called
        from L{self.computeNormalAndArea()}. The length of the normal is "1". 
        """
        p1 = numpy.array(tria[0])
        p2 = numpy.array(tria[1])
        p3 = numpy.array(tria[2])
        p21 = p2 - p1
        p31 = p3 - p1
        A = 0.5 * dpTensor.Length(dpTensor.CrossProduct(p21, p31))
        xs = (p1[0] + p2[0] + p3[0]) / 3.0
        ys = (p1[1] + p2[1] + p3[1]) / 3.0
        zs = (p1[2] + p2[2] + p3[2]) / 3.0
        fac = 1.0 / numpy.sqrt(xs * xs + ys * ys + zs * zs)
        return (
         A, (fac * xs, fac * ys, fac * zs))

    def setupStarFem(self, dist, basePts=None, bothDir=False, starLength=None):
        """
        Function creates a "star" FEM mesh. Line mesh goes from a base point 
        (default = (0.,0.,0.)) to a given point in space. 
        
        @param dist: List of points in space. 
          - TYPE: list[ point1, point2, ... ]
          - TYPE point : tuple(xcoord, ycoord, zcoord) 
          - float xcoord, ycoord, zcoord ... coordiantes of point 
        @param basePts: base point of line mesh. Default = (0.,0.,0.)
          - TYPE point : dict[ (nX,nY,nZ) ] = (xcoord, ycoord, zcoord) 
          - float nX, nY, nZ ... coordinates of the corresponding normal vector 
          - float xcoord, ycoord, zcoord ... coordiantes of point                                  
        @param bothDir: Flag if line mesh should also be created in opposite direction 
          - TYPE: bool
        @param starLength: List of values for a given length of line mesh. 
          The vector from base point to points in space will be scale such that 
          it length is equal to starLength value. 
            - TYPE point : dict[ (nX,nY,nZ) ] = value 
            - float nX, nY, nZ ... coordinates of the corresponding normal vector           
            - float value ... length value              
                     
        @return:                 
             Nodes: Nodes of the mesh 
               - TYPE: dict[ nodeId ] = node 
               - int nodeId      : node identifier  
               - dpFem.node node : FEM node 
             voxElems: Elements of the mesh 
               - TYPE: dict[ elemId ] = element 
               - int elemId            : element identifier  
               - dpFem.element element : FEM element object 
             EscaResults: Element Scalar results  
               - TYPE: dict[ elemId ] = value 
               - int elemId  : element identifier  
               - float value : element scalar value                                                  
        """
        nid = 0
        elid = 0
        Nodes = {}
        starElems = {}
        EscaResults = {}
        newDist = {}
        if basePts != None:
            newBase = {}
        if starLength != None:
            newLen = {}
        if bothDir == True:
            for di in dist:
                newDist[di] = dist[di]
                newDist[-di[0], -di[1], -di[2]] = dist[di]
                if basePts != None:
                    newBase[di] = basePts[di]
                    newBase[-di[0], -di[1], -di[2]] = basePts[di]
                if starLength != None:
                    newLen[di] = starLength[di]
                    newLen[-di[0], -di[1], -di[2]] = starLength[di]

            dist = newDist
            if basePts != None:
                basePts = newBase
            if starLength != None:
                starLength = newLen
        baseId = len(dist) + 1
        if basePts == None:
            baseNode = dpFem.node(baseId, 0.0, 0.0, 0.0)
            Nodes[baseId] = baseNode
        for di in dist:
            if basePts != None:
                baseId += 1
                baseNode = dpFem.node(baseId, basePts[di][0], basePts[di][1], basePts[di][2])
                Nodes[baseId] = baseNode
            elid += 1
            nid += 1
            Node = dpFem.node(nid, 0.0)
            if starLength == None:
                length = dist[di]
            else:
                length = float(starLength[di])
            if basePts == None:
                Node.set_x(float(length * di[0]))
                Node.set_y(float(length * di[1]))
                Node.set_z(float(length * di[2]))
            else:
                Node.set_x(float(length * di[0] + basePts[di][0]))
                Node.set_y(float(length * di[1] + basePts[di][1]))
                Node.set_z(float(length * di[2] + basePts[di][2]))
            Nodes[nid] = Node
            nList = [baseId, nid]
            curElem = dpFem.element(nid, nList, 'bar2')
            starElems[elid] = curElem
            EscaResults[elid] = float(dist[di])

        return (
         Nodes, starElems, EscaResults)

    def setupRayFieldFEM(self, rayfield):
        """ 
        Setup an FEM model for a ray field. For test purposes.  
        """
        nid = 0
        elid = 0
        Nodes = {}
        starElems = {}
        for n in rayfield:
            rayPts = rayfield[n]
            for ray in rayPts:
                for npos in range(2):
                    nid += 1
                    Node = dpFem.node(nid, float(ray[npos][0]), float(ray[npos][1]), float(ray[npos][2]))
                    Nodes[nid] = Node

                elid += 1
                nList = [nid - 1, nid]
                curElem = dpFem.element(nid, nList, 'bar2', 'rays')
                starElems[elid] = curElem

        return (Nodes, starElems)

    def setupVoxelStarFem(self, normalsAndVoxels):
        """
        Setup an FEM model for voxelized "stars"
        """
        nid = 0
        elid = 0
        Nodes = {}
        voxelElems = {}
        corners = {}
        corners[1] = (0.0, 0.0, 0.0)
        corners[2] = (1.0, 0.0, 0.0)
        corners[3] = (1.0, 1.0, 0.0)
        corners[4] = (0.0, 1.0, 0.0)
        corners[5] = (0.0, 0.0, 1.0)
        corners[6] = (1.0, 0.0, 1.0)
        corners[7] = (1.0, 1.0, 1.0)
        corners[8] = (0.0, 1.0, 1.0)
        for n in normalsAndVoxels:
            voxelArray = normalsAndVoxels[n]
            for voxel in voxelArray:
                nList = []
                for corn in corners:
                    pt = corners[corn]
                    nid += 1
                    Node = dpFem.node(nid, float(voxel[0] - 1.0 + pt[0]), float(voxel[1] - 1.0 + pt[1]), float(voxel[2] - 1.0 + pt[2]))
                    Nodes[nid] = Node
                    nList.append(nid)

                elid += 1
                curElem = dpFem.element(elid, nList, 'hexa8', 'voxels')
                voxelElems[elid] = curElem

        return (Nodes, voxelElems)

    def setupVoxelFem(self, bboxList, EscaList=None):
        """
         Setup a voxel FEM model from a list of bounding boxes. 
         
         @param bboxList: Bounding box list 
            - TYPE: list[ anaId ] = (cenX, cenY, cenZ, lenX, lenY, lenZ) 
            - int anaId            ... analysis identifier                              
            - int cenX, cenY, cenZ ... coordinates x,y,z of center of RVE 
              (number of voxels in x,y,z) 
            - int lenX, lenY, lenZ ... length of RVE (number of voxels in x,y,z)
         @param EscaList: List of element scalar values
            - TYPE: list[ anaId ] = value  
            - int anaId   ... analysis identifier  
            - float value ... scalar value
                                    
         @return:      
              Nodes: Nodes of the mesh 
                - TYPE: dict[ nodeId ] = node 
                - int nodeId      : node identifier  
                - dpFem.node node : FEM node 
              voxElems: Elements of the mesh 
                - TYPE: dict[ elemId ] = element 
                - int elemId            : element identifier  
                - dpFem.element element : FEM element object 
              EscaResults: Element Scalar results  
                - TYPE: dict[ elemId ] = value 
                - int elemId  : element identifier  
                - float value : element scalar value       
        """
        nid = 0
        elid = 0
        Nodes = {}
        voxelElems = {}
        EscaResults = {}
        corners = {}
        corners[1] = (0.0, 0.0, 0.0)
        corners[2] = (1.0, 0.0, 0.0)
        corners[3] = (1.0, 1.0, 0.0)
        corners[4] = (0.0, 1.0, 0.0)
        corners[5] = (0.0, 0.0, 1.0)
        corners[6] = (1.0, 0.0, 1.0)
        corners[7] = (1.0, 1.0, 1.0)
        corners[8] = (0.0, 1.0, 1.0)
        for aId in range(len(bboxList)):
            bbox = bboxList[aId]
            staX = bbox[0]
            staY = bbox[1]
            staZ = bbox[2]
            endX = staX + bbox[3] - 1
            endY = staY + bbox[4] - 1
            endZ = staZ + bbox[5] - 1
            corners[1] = (
             staX, staY, staZ)
            corners[2] = (endX, staY, staZ)
            corners[3] = (endX, endY, staZ)
            corners[4] = (staX, endY, staZ)
            corners[5] = (staX, staY, endZ)
            corners[6] = (endX, staY, endZ)
            corners[7] = (endX, endY, endZ)
            corners[8] = (staX, endY, endZ)
            nList = []
            for corn in corners:
                pt = corners[corn]
                nid += 1
                Node = dpFem.node(nid, float(pt[0]), float(pt[1]), float(pt[2]))
                Nodes[nid] = Node
                nList.append(nid)

            elid += 1
            curElem = dpFem.element(elid, nList, 'hexa8', 'Elset' + repr(aId + 1))
            voxelElems[elid] = curElem
            if EscaList == None:
                EscaResults[elid] = None
            else:
                EscaResults[elid] = float(EscaList[aId])

        return (Nodes, voxelElems, EscaResults)

    def writeHeader(self, OS, step=None, power=None, ftype=None, valid=None, dtype=None):
        """
        Write a header to an output file. 
        
        @param step: Step in the considered voxel.   
           - TYPE: int > 0 
        @param   power: Power for the number of star directions = 8*4^power
           - TYPE: int > 1  
        @param   ftype: type of approximation: 
           - TYPE: int = 1,2,4
           - ftype = 1 : Ellipsoid
           - ftype = 2 : general tensor 2nd order 
           - ftype = 4 : general tensor 4th order      
                                       
        @param valid: Smallest valid intersection length.
           - TYPE: float
                                 
        @return:  
               no return value
        """
        OS.write('\n=====================================================================')
        OS.write('\n     Medical Image Analysis  by  D.H. Pahr ')
        OS.write('\n=====================================================================')
        dt = localtime()
        date = str(dt[0]) + '/' + str(dt[1]) + '/' + str(dt[2]) + ' - ' + str(dt[3]) + ':' + str(dt[4])
        OS.write('\n Date - Time   = %s' % date)
        if not dtype.upper() == 'GST':
            if step != None:
                OS.write('\n Param: -step  = %6i' % step)
            if power != None:
                OS.write('\n Param: -power = %6i   (%5i Directions)' % (power, int(8.0 * 4.0 ** float(power))))
            if ftype != None:
                OS.write('\n Param: -ftype = %6i   (Ellip..1, 2nd Order..2)' % ftype)
            if valid != None:
                OS.write('\n Param: -valid = %6.1f' % valid)
        return

    def writeEigenValueAndVectors(self, fabricName, Eval, Evect, normalize, OS):
        """
        Write eigenvalues and eigenvectors to an output file.  
        """
        OS.write('\n\n\n--------------------------------------------------------')
        OS.write('\nRESULTS: Morphological analysis normalized by : %s' % normalize)
        OS.write('\n--------------------------------------------------------')
        OS.write('\n %s - Eigen values   : %9.5f  %9.5f  %9.5f' % (fabricName, Eval[0], Eval[1], Eval[2]))
        OS.write('\n %s - Eigen vector 1 : %9.5f  %9.5f  %9.5f' % (fabricName, Evect[0][0], Evect[1][0], Evect[2][0]))
        OS.write('\n %s - Eigen vector 2 : %9.5f  %9.5f  %9.5f' % (fabricName, Evect[0][1], Evect[1][1], Evect[2][1]))
        OS.write('\n %s - Eigen vector 3 : %9.5f  %9.5f  %9.5f' % (fabricName, Evect[0][2], Evect[1][2], Evect[2][2]))
        OS.write('\n')

    def writeEigenValueAndVectors2(self, Name, Eval, Evect, OS):
        """
        Write eigenvalues and eigenvectors to an output file.  
        """
        OS.write('%20s %9.5f  %9.5f  %9.5f ' % (Name, Eval[0], Eval[1], Eval[2]))
        OS.write('%9.5f  %9.5f  %9.5f ' % (Evect[0][0], Evect[0][1], Evect[0][2]))
        OS.write('%9.5f  %9.5f  %9.5f ' % (Evect[1][0], Evect[1][1], Evect[1][2]))
        OS.write('%9.5f  %9.5f  %9.5f\n' % (Evect[2][0], Evect[2][1], Evect[2][2]))

    def writeParaview(self, name, typeDIST, orgDIST, approxDIST, triaArray, scale, EvectDIST, EvalDIST):
        """
        Write a PARAVIEW FEM file for eigenvectors, approximated and origial 
        distribution of one morphology analysis. 
        
        @param name: filename of the output file 
                          - TYPE: string      
        @param typeDIST: name/type of distribution e.g. MIL, SLD, SVD
                          - TYPE: string
        @param orgDIST: Original distribution written as "star", line elements 
                        from point(0.,0.,0.), direction = n, length = value 
                          - TYPE: dict[ (nix, niy, niz) ] = value 
                          - float nix, niy, niz ... components of normal vectors 
                          - float value         ... distrubtion value for that direction    
        @param approxDIST: Approximated distribution written as triangulated surface 
                          - TYPE: dict[ (nix, niy, niz) ] = value 
                          - float nix, niy, niz ... components of normal vectors 
                          - float value         ... distrubtion value for that direction 
        @param triaArray: List of triangles which build the surface of approxDIST
                          - TYPE: list[ ntriangles ] = ( (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ) 
                          - ntriangles ... number of triangles 
                          - float xi, yi, zi ... x,y,z coordinates of point i  
        @param scale: scale factor for the orgDIST / approxDIST data 
                          - TYPE: float      
        @param EvectDIST: Eigenvectors of fabric tensor
                          - TYPE: numpy.array[evectID, evalID] = evect
                          - int evectID ... eigenvector ID. 0,1, or 2                          
                          - int evalID  ... eigenvalue ID. 0,1, or 2
                          - float evect ... component of eigenvectors 
        @param EvalDIST: Eigenvalues of fabric tensor 
                          - TYPE: numpy.array[evalID] = eval
                          - int evalID ... eigenvalues ID. 0,1, or 2
                          - float eval ... current eigenvalue 
                                               
        @return:      
             no return value                                   
        """
        fecModel = fec.fec()
        mul = 1.0
        DISTStar = {}
        DISTStar[EvectDIST[0][0], EvectDIST[1][0], EvectDIST[2][0]] = mul * EvalDIST[0]
        DISTStar[EvectDIST[0][1], EvectDIST[1][1], EvectDIST[2][1]] = mul * EvalDIST[1]
        DISTStar[EvectDIST[0][2], EvectDIST[1][2], EvectDIST[2][2]] = mul * EvalDIST[2]
        Nodes, starElems, Eresults = self.setupStarFem(DISTStar, bothDir=True)
        sEResults = [Eresults]
        fecModel.write(name + '_eigen' + typeDIST + '.geo', typeDIST, Nodes, None, starElems, None, EscaResults=sEResults)
        for keyId in approxDIST:
            approxDIST[keyId] = approxDIST[keyId] * scale

        approxNodes, approxElems, results = self.setupSphereFem(triaArray, approxDIST)
        sResults = [results]
        fecModel.write(name + '_approx' + typeDIST + '.geo', typeDIST, approxNodes, None, approxElems, None, NscaResults=sResults)
        for keyId in orgDIST:
            orgDIST[keyId] = orgDIST[keyId] * scale

        Nodes, starElems, Eresults = self.setupStarFem(orgDIST)
        sEResults = [Eresults]
        fecModel.write(name + '_orig' + typeDIST + '.geo', typeDIST, Nodes, None, starElems, None, EscaResults=sEResults)
        return

    def writeParaview2(self, name, typeDIST, approxDIST, triaArray, scale, EvectDIST, EvalDIST):
        """
        Same as writeParaview but without writing orgDIST
        """
        fecModel = fec.fec()
        mul = 1.0
        DISTStar = {}
        DISTStar[EvectDIST[0][0], EvectDIST[1][0], EvectDIST[2][0]] = mul * EvalDIST[0]
        DISTStar[EvectDIST[0][1], EvectDIST[1][1], EvectDIST[2][1]] = mul * EvalDIST[1]
        DISTStar[EvectDIST[0][2], EvectDIST[1][2], EvectDIST[2][2]] = mul * EvalDIST[2]
        Nodes, starElems, Eresults = self.setupStarFem(DISTStar, bothDir=True)
        sEResults = [Eresults]
        fecModel.write(name + '_eigen' + typeDIST + '.geo', typeDIST, Nodes, None, starElems, None, EscaResults=sEResults)
        for keyId in approxDIST:
            approxDIST[keyId] = approxDIST[keyId] * scale

        approxNodes, approxElems, results = self.setupSphereFem(triaArray, approxDIST)
        sResults = [results]
        fecModel.write(name + '_approx' + typeDIST + '.geo', typeDIST, approxNodes, None, approxElems, None, NscaResults=sResults)
        return

    def writeParaviewAverage(self, name, typeDIST, mul, offset, EvectDIST, EvalDIST):
        """
        Write a PARAVIEW FEM file for eigenvectors.
        """
        DISTStar = {}
        DISTStar[EvectDIST[0][0], EvectDIST[0][1], EvectDIST[0][2]] = mul * EvalDIST[0]
        DISTStar[EvectDIST[1][0], EvectDIST[1][1], EvectDIST[1][2]] = mul * EvalDIST[1]
        DISTStar[EvectDIST[2][0], EvectDIST[2][1], EvectDIST[2][2]] = mul * EvalDIST[2]
        Nodes, starElems, Eresults = self.setupStarFem(DISTStar, bothDir=True)
        for noid in Nodes:
            no = Nodes[noid]
            no.set_x(no.get_x() + offset[0])
            no.set_y(no.get_y() + offset[1])
            no.set_z(no.get_z() + offset[2])

        sEResults = [Eresults]
        fecModel = fec.fec()
        fecModel.write(name + '_Average' + typeDIST + '.geo', typeDIST, Nodes, None, starElems, None, EscaResults=sEResults)
        return

    def writeParaviewMulti(self, name, stiffList, bboxList, EvectList, rhoList):
        """
        Write a PARAVIEW FEM file with orthotropic axes and RVE sizes of multiple 
        morphology analyses.  
        
        @param name: filename of the output file 
           - TYPE: string                        
        @param stiffList: List of elasticities for each analysis 
           - TYPE list[ anaId ] = stiffness 
           - int anaId             ... analysis identifier 
           - dpStiffness stiffness ... stiffness information 
        @param bboxList: Bounding box list 
           - TYPE: list[ anaId ] = (cenX, cenY, cenZ, lenX, lenY, lenZ, type) 
           - int anaId            ... analysis identifier                              
           - int cenX, cenY, cenZ ... coordinates x,y,z of center of RVE 
             (number of voxels in x,y,z) 
           - int lenX, lenY, lenZ ... length of RVE (number of voxels in x,y,z)   
           - string type ... RVE type: 'sphere' or 'hex'  
        @param EvectList: List of Eigenvectors 
           - TYPE: list[ anaId ] = evector
           - int anaId ... analysis identifier 
           - evector   ... eigenvectors see self.computeEigenValueAndVector() for details
        @param rhoList: List of densities (BVTV) 
           - TYPE: list[ anaId ] = BVTV
           - int anaId  ... analysis identifier 
           - float BVTV ... density           
                       
        @return:      
             no return value 
        """
        fecModel = fec.fec()
        aMax = len(stiffList)
        basePtsMin = {}
        basePtsMid = {}
        basePtsMax = {}
        directsMin = {}
        directsMid = {}
        directsMax = {}
        starLenMin = {}
        starLenMid = {}
        starLenMax = {}
        EList = []
        lList = []
        for aNo in range(aMax):
            stiff = stiffList[aNo]
            engConst = stiff.getEngConstants()
            EList.append(engConst['E_1'])
            EList.append(engConst['E_2'])
            EList.append(engConst['E_3'])
            bbox = bboxList[aNo]
            lenX = bbox[3]
            lenY = bbox[4]
            lenZ = bbox[5]
            lList.append(lenX)
            lList.append(lenY)
            lList.append(lenZ)

        Emax = max(EList)
        Emin = min(EList)
        Lmax = max(lList)
        scale = 1.0
        L3 = 0.3 * Lmax
        L1 = 1.2 * Lmax
        for aNo in range(aMax):
            stiff = stiffList[aNo]
            engConst = stiff.getEngConstants()
            E1 = engConst['E_1']
            E2 = engConst['E_2']
            E3 = engConst['E_3']
            bbox = bboxList[aNo]
            lenX = bbox[3]
            lenY = bbox[4]
            lenZ = bbox[5]
            cenX = bbox[0] + lenX / 2
            cenY = bbox[1] + lenY / 2
            cenZ = bbox[2] + lenZ / 2
            basPt = (
             cenX, cenY, cenZ)
            Evect = EvectList[aNo]
            ev1 = Evect[0]
            ev2 = Evect[1]
            ev3 = Evect[2]
            Evects = {}
            Evects[E1] = (
             ev1[0], ev1[1], ev1[2])
            Evects[E2] = (ev2[0], ev2[1], ev2[2])
            Evects[E3] = (ev3[0], ev3[1], ev3[2])
            Emods = [
             E1, E2, E3]
            Emods.sort()
            Emods.reverse()
            if directsMin.has_key(Evects[Emods[2]]):
                print '**WARNING: writeParaviewMulti(): Two similar directions. First will be overwritten!'
            else:
                directsMin[Evects[Emods[2]]] = scale * Emods[2]
                basePtsMin[Evects[Emods[2]]] = basPt
                starLenMin[Evects[Emods[2]]] = (L3 - L1) / (Emin - Emax) * Emods[2] + (Emax * L3 - Emin * L1) / (Emax - Emin)
            if directsMid.has_key(Evects[Emods[1]]):
                print '**WARNING: writeParaviewMulti(): Two similar directions. First will be overwritten!'
            else:
                directsMid[Evects[Emods[1]]] = scale * Emods[1]
                basePtsMid[Evects[Emods[1]]] = basPt
                starLenMid[Evects[Emods[1]]] = (L3 - L1) / (Emin - Emax) * Emods[1] + (Emax * L3 - Emin * L1) / (Emax - Emin)
            if directsMax.has_key(Evects[Emods[0]]):
                print '**WARNING: writeParaviewMulti(): Two similar directions. First will be overwritten!'
            else:
                directsMax[Evects[Emods[0]]] = scale * Emods[0]
                basePtsMax[Evects[Emods[0]]] = basPt
                starLenMax[Evects[Emods[0]]] = (L3 - L1) / (Emin - Emax) * Emods[0] + (Emax * L3 - Emin * L1) / (Emax - Emin)

        Nodes, starElems, Eresults = self.setupStarFem(directsMin, basePts=basePtsMin, bothDir=True, starLength=starLenMin)
        for elid in Eresults:
            Eresults[elid] = Eresults[elid] / scale

        sEResults = [
         Eresults]
        fecModel.write(name + '_Orient' + 'Min' + '.geo', 'Min', Nodes, None, starElems, None, EscaResults=sEResults)
        Nodes, starElems, Eresults = self.setupStarFem(directsMid, basePts=basePtsMid, bothDir=True, starLength=starLenMid)
        for elid in Eresults:
            Eresults[elid] = Eresults[elid] / scale

        sEResults = [
         Eresults]
        fecModel.write(name + '_Orient' + 'Mid' + '.geo', 'Mid', Nodes, None, starElems, None, EscaResults=sEResults)
        Nodes, starElems, Eresults = self.setupStarFem(directsMax, basePts=basePtsMax, bothDir=True, starLength=starLenMax)
        for elid in Eresults:
            Eresults[elid] = Eresults[elid] / scale

        sEResults = [
         Eresults]
        fecModel.write(name + '_Orient' + 'Max' + '.geo', 'Max', Nodes, None, starElems, None, EscaResults=sEResults)
        Nodes, starElems, Eresults = self.setupVoxelFem(bboxList, rhoList)
        sEResults = [Eresults]
        fecModel.write(name + '_rve' + '.geo', 'RVE', Nodes, None, starElems, None, EscaResults=sEResults)
        return

    def writeParaviewMulti2(self, name, stiffList, EvectList, rhoList, isoFlagList, femNodes, femElems):
        """
        Write a PARAVIEW FEM file with orthotropic axes and RVE sizes of multiple 
        morphology analyses.  
        
        @param name: filename of the output file 
           - TYPE: string                        
        @param stiffList: List of elasticities for each analysis 
           - TYPE list[ anaId ] = stiffness 
           - int anaId             ... analysis identifier 
           - dpStiffness stiffness ... stiffness information 
        @param EvectList: List of Eigenvectors 
           - TYPE: list[ anaId ] = evector
           - int anaId ... analysis identifier 
           - evector   ... eigenvectors see self.computeEigenValueAndVector() for details
        @param rhoList: List of densities (BVTV) 
           - TYPE: list[ anaId ] = BVTV
           - int anaId  ... analysis identifier 
           - float BVTV ... density           
        @param femNodes: Nodes of the mesh 
           - TYPE: dict[ nodeId ] = node 
           - int nodeId      : node identifier  
           - dpFem.node node : FEM node 
        @param femElems: Elements of the mesh 
           - TYPE: dict[ elemId ] = element 
           - int elemId            : element identifier  
           - dpFem.element element : FEM element object                
                       
        @return:      
             no return value 
        """
        fecModel = fec.fec()
        basePtsMin = {}
        basePtsMid = {}
        basePtsMax = {}
        directsMin = {}
        directsMid = {}
        directsMax = {}
        starLenMin = {}
        starLenMid = {}
        starLenMax = {}
        EList = []
        for elid in stiffList:
            stiff = stiffList[elid]
            engConst = stiff.getEngConstants()
            EList.append(engConst['E_1'])
            EList.append(engConst['E_2'])
            EList.append(engConst['E_3'])

        Emax = max(EList)
        Emin = min(EList)
        myMesh = dpMesh.dpMesh(nodes=femNodes, elems=femElems)
        Lmean = myMesh.getMeanNodalDistance()
        scale = 1.0
        L3 = 0.3 * Lmean
        L1 = 1.2 * Lmean
        warnFlag = False
        for elid in stiffList:
            stiff = stiffList[elid]
            engConst = stiff.getEngConstants()
            E1 = engConst['E_1']
            E2 = engConst['E_2']
            E3 = engConst['E_3']
            sumNo = 0
            xNo = 0.0
            yNo = 0.0
            zNo = 0.0
            for noid in femElems[elid].get_nodes():
                sumNo += 1
                xNo += femNodes[noid].get_x()
                yNo += femNodes[noid].get_y()
                zNo += femNodes[noid].get_z()

            basPt = (
             xNo / float(sumNo), yNo / float(sumNo), zNo / float(sumNo))
            Evect = EvectList[elid]
            ev1 = Evect[0]
            ev2 = Evect[1]
            ev3 = Evect[2]
            Evects = {}
            Evects[E1] = (
             ev1[0], ev1[1], ev1[2])
            Evects[E2] = (ev2[0], ev2[1], ev2[2])
            Evects[E3] = (ev3[0], ev3[1], ev3[2])
            Emods = [
             E1, E2, E3]
            Emods.sort()
            Emods.reverse()
            if directsMin.has_key(Evects[Emods[2]]):
                warnFlag = True
            else:
                directsMin[Evects[Emods[2]]] = scale * Emods[2]
                basePtsMin[Evects[Emods[2]]] = basPt
                starLenMin[Evects[Emods[2]]] = (L3 - L1) / (Emin - Emax) * Emods[2] + (Emax * L3 - Emin * L1) / (Emax - Emin)
            if directsMid.has_key(Evects[Emods[1]]):
                warnFlag = True
            else:
                directsMid[Evects[Emods[1]]] = scale * Emods[1]
                basePtsMid[Evects[Emods[1]]] = basPt
                starLenMid[Evects[Emods[1]]] = (L3 - L1) / (Emin - Emax) * Emods[1] + (Emax * L3 - Emin * L1) / (Emax - Emin)
            if directsMax.has_key(Evects[Emods[0]]):
                warnFlag = True
            else:
                directsMax[Evects[Emods[0]]] = scale * Emods[0]
                basePtsMax[Evects[Emods[0]]] = basPt
                starLenMax[Evects[Emods[0]]] = (L3 - L1) / (Emin - Emax) * Emods[0] + (Emax * L3 - Emin * L1) / (Emax - Emin)

        if warnFlag == True:
            print '\n **WARNING: writeParaviewMulti(): Two similar directions dedected.'
            print '            Usually in case of isotropy. First will be overwritten! \n '
        Nodes, starElems, Eresults = self.setupStarFem(directsMin, basePts=basePtsMin, bothDir=True, starLength=starLenMin)
        for elid in Eresults:
            Eresults[elid] = Eresults[elid] / scale

        sEResults = [
         Eresults]
        fecModel.write(name + '_Orient' + 'Min' + '.geo', 'Min', Nodes, None, starElems, None, EscaResults=sEResults)
        Nodes, starElems, Eresults = self.setupStarFem(directsMid, basePts=basePtsMid, bothDir=True, starLength=starLenMid)
        for elid in Eresults:
            Eresults[elid] = Eresults[elid] / scale

        sEResults = [
         Eresults]
        fecModel.write(name + '_Orient' + 'Mid' + '.geo', 'Mid', Nodes, None, starElems, None, EscaResults=sEResults)
        Nodes, starElems, Eresults = self.setupStarFem(directsMax, basePts=basePtsMax, bothDir=True, starLength=starLenMax)
        for elid in Eresults:
            Eresults[elid] = Eresults[elid] / scale

        sEResults = [
         Eresults]
        fecModel.write(name + '_Orient' + 'Max' + '.geo', 'Max', Nodes, None, starElems, None, EscaResults=sEResults)
        fecModel.write(name + '_femRho' + '.geo', 'RVE', femNodes, None, femElems, None, EscaResults=[rhoList])
        if isoFlagList != None:
            fecModel.write(name + '_femIsoFlag' + '.geo', 'RVE', femNodes, None, femElems, None, EscaResults=[isoFlagList])
        return

    def writeAbaqusRVE(self, name, bboxList):
        """
        Write a PARAVIEW FEM file with orthotropic axes and RVE sizes of multiple 
        morphology analyses.  
        
        @param name: filename of the output file 
           - TYPE: string         
        @param bboxList: Bounding box list 
           - TYPE: list[ anaId ] = (cenX, cenY, cenZ, lenX, lenY, lenZ) 
           - int anaId            ... analysis identifier                              
           - int cenX, cenY, cenZ ... coordinates x,y,z of center of RVE 
             (number of voxels in x,y,z) 
           - int lenX, lenY, lenZ ... length of RVE (number of voxels in x,y,z)        
                       
        @return:      
             no return value 
        """
        fecModel = fec.fec()
        Nodes, voxElems, Eresults = self.setupVoxelFem(bboxList)
        NodeSets = []
        ElemsSets = []
        for i in range(len(bboxList)):
            aId = i + 1
            ElemsSets.append('Elset' + repr(aId))

        fecModel.write(name + '_rve' + '.inp', None, Nodes, NodeSets, voxElems, ElemsSets)
        return

    def writeDistribution(self, filename, DISTtype, origDIST):
        """
        Write a morphology distribution to an output file. 
        """
        stdout.write(' ... write original distribution\n')
        stdout.flush()
        OS = open(filename + '-' + DISTtype + '.dis', 'w')
        for n in origDIST:
            OS.write('%13.7f %13.7f %13.7f %13.7g\n' % (n[0], n[1], n[2], origDIST[n]))

    def initValues(self, sVal, valName='unknown', minVal=0, valType=''):
        """
        Initialize argument values. 
        """
        if sVal != None:
            sValList = []
            sVal2 = sVal.replace(':', ' ')
            sValStr = sVal2.split()
            if minVal != 0 and len(sValStr) != minVal:
                stdout.write("\n **ERROR** : initVales(): Option '%s' needs '%s' values!\n\n" % (valName, repr(minVal)))
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
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

    def get_Shape(self, voxelModel):
        """
        Return the shape of a voxel model. 
        """
        shapes = voxelModel.shape
        return (
         shapes[2], shapes[1], shapes[0])

    def orderEvalEvec(self, iEval, iEvec, prec=1e-06):
        isEval = numpy.zeros(3, numpy.float) + 1.0
        isEvec = dpTensor.UnitMatrix(3)
        eigen_tuples = [
         (
          0, iEval[0]), (1, iEval[1]), (2, iEval[2])]
        b = sorted(eigen_tuples, key=lambda eigen: eigen[1], reverse=True)
        isEval[0] = iEval[b[0][0]]
        isEval[1] = iEval[b[1][0]]
        isEval[2] = iEval[b[2][0]]
        isEvec[0] = iEvec[b[0][0]]
        isEvec[1] = iEvec[b[1][0]]
        isEvec[2] = iEvec[b[2][0]]
        return (
         isEval, isEvec)


if __name__ == '__main__':
    inName = None
    outName = None
    step = None
    power = None
    ftype = None
    dtype = None
    valid = None
    code = None
    gstpow = None
    thres = None
    graph = None
    distr = None
    mod = None
    norm = None
    sort = None
    fout = None
    guiInit = "*modulName           'MIA - Medical Image Analyzer'                                          \n" + "*fileEntryIn  -in    'Input File Name'                       test.mhd   no   mhd;rawdp       \n" + "*entry        -out   'Analysis Output File Name'             test.fab   no   fab;txt;dat     \n" + "*combo        -dtype 'Distribution Type'                     MIL        no   MIL;SLD;SVD;GST \n" + "*combo        -ftype 'Fabric Approximation Type'             1          no   1;2;4           \n" + "*entry        -thres 'Threshold (value)'                     75         yes  1               \n" + "*subWindowStart      'Solution Parameter'                                                    \n" + "*entry        -step  'Ray Distance in # of Voxels (>1)'      5          yes  1               \n" + "*entry        -power 'Power of Star Directions = 8*4^n'      2          yes  1               \n" + "*entry        -valid 'Min. Valid Intersection Lenght (>=0)'  0.0        yes  1               \n" + "*combo        -code  'Code of the core routines'             f77        yes  f77;py          \n" + "*entry        -gstpow 'Power for GST method'                 0.5        yes  1               \n" + "*combo        -norm  'Fabric Normalization'                  det        yes  det;rank;trace  \n" + "*subWindowEnd        'Solution Parameter'                                                    \n" + "*subWindowStart      'Output Options'                                                        \n" + "*combo        -mod   'Output Mode'                           w          yes  w;a             \n" + "*combo        -graph 'Graphical Distribution Output'         ON         yes  ON;OFF          \n" + "*combo        -distr 'Original Distribution Output'          ON         yes  ON;OFF          \n" + "*combo        -sort  'Sort EigVal and EigVec'                OFF        yes  ON;OFF          \n" + "*fileEntryOut -fout  'Fabric Output'                         fab.txt    yes  txt;dat         \n" + "*subWindowEnd        'Output Options'                                                        \n"
    argList = argv
    argc = len(argList)
    i = 0
    while i < argc:
        if argList[i][:3] == '-in':
            i += 1
            inName = argList[i]
        elif argList[i][:3] == '-ou':
            i += 1
            outName = argList[i]
        elif argList[i][:4] == '-ste':
            i += 1
            step = argList[i]
        elif argList[i][:3] == '-po':
            i += 1
            power = argList[i]
        elif argList[i][:6] == '-ftype':
            i += 1
            ftype = argList[i]
        elif argList[i][:6] == '-dtype':
            i += 1
            dtype = argList[i]
        elif argList[i][:3] == '-gu':
            print guiInit
            exit(0)
        elif argList[i][:3] == '-va':
            i += 1
            valid = argList[i]
        elif argList[i][:3] == '-co':
            i += 1
            code = argList[i]
        elif argList[i][:7] == '-gstpow':
            i += 1
            gstpow = argList[i]
        elif argList[i][:6] == '-thres':
            i += 1
            thres = argList[i]
        elif argList[i][:6] == '-graph':
            i += 1
            graph = argList[i]
        elif argList[i][:6] == '-distr':
            i += 1
            distr = argList[i]
        elif argList[i][:4] == '-mod':
            i += 1
            mod = argList[i]
        elif argList[i][:5] == '-norm':
            i += 1
            norm = argList[i]
        elif argList[i][:5] == '-sort':
            i += 1
            sort = argList[i]
        elif argList[i][:5] == '-fout':
            i += 1
            fout = argList[i]
        elif argList[i][:2] == '-h':
            print __doc__
            exit(0)
        i += 1

    if not inName:
        stdout.write(__doc__)
        stdout.flush()
        stdout.write('\n **ERROR** input file name not given\n\n')
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if not outName:
        stdout.write(__doc__)
        stdout.flush()
        stdout.write('\n **ERROR** output file name not given\n\n')
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if not dtype:
        stdout.write(__doc__)
        stdout.flush()
        stdout.write('\n **ERROR** -dtype option not given\n\n')
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if not ftype:
        stdout.write(__doc__)
        stdout.flush()
        stdout.write('\n **ERROR** -ftype option not given\n\n')
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if dtype.upper() == 'GST':
        if int(ftype) == 1:
            pass
        else:
            stdout.write(__doc__)
            stdout.flush()
            stdout.write('\n **ERROR** -ftype = %s not implemented for GST\n\n' % ftype)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    stdout.write('\n S T A R T  %s %s - %s \n\n' % (mia.__name__, mia.__version__, mia.__author__))
    stdout.flush()
    startTime = clock()
    ctime1 = time()
    myMiaModel = mia()
    myMicModel = mic.mic()
    dtype = dtype.upper()
    if step == None:
        step = 5
    else:
        step = int(step)
    if power == None:
        power = 2
    else:
        power = int(power)
    if ftype == None:
        ftype = 2
    else:
        ftype = int(ftype)
    if valid == None:
        valid = 0.0
    else:
        valid = float(valid)
    if code == None:
        code = 'f77'
    if gstpow == None:
        gstpow = 1.0
    else:
        gstpow = float(gstpow)
    if norm == None:
        norm = 'trace'
    elif not (norm.lower() == 'det' or norm.lower() == 'rank' or norm.lower() == 'trace'):
        stdout.write('\n **ERROR** : mia.computeEigenValueAndVector() -norm = %s not known!' % norm)
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if sort == None:
        sort = 'OFF'
    elif sort == 'ON' or sort == 'OFF':
        pass
    else:
        stdout.write("\n **ERROR** : The option -sort can only be ON or OFF and not '%s'" % sort)
        stdout.flush()
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    threshold = None
    ldim = None
    ndim = None
    thresVoxelModel, additData = myMicModel.read(inName)
    myMiaModel.set_nx_ny_nz(thresVoxelModel)
    if additData.has_key('-thres'):
        threshold = additData['-thres']
        if len(threshold) > 1:
            stdout.write('\n **ERROR** : Only one threshold (bw image) allowed!')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    elif dtype.upper() == 'GST':
        pass
    elif thres == None:
        unique = numpy.unique(thresVoxelModel)
        if len(unique) > 2:
            stdout.write('\n **ERROR** : It seems the given voxel image file is not thresholded!')
            stdout.write('\n             More than two gray-values are in the file!\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        minVal, maxVal = (
         numpy.min(thresVoxelModel), numpy.max(thresVoxelModel))
        if minVal != 0:
            stdout.write('\n **ERROR** : Minimal value of thresholded file needs to be zero!')
            stdout.write('\n             Current min=%s!\n' % repr(minVal))
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        additData['-thres'] = [
         int(maxVal)]
        threshold = additData['-thres']
    else:
        additData['-thres'] = [
         int(thres)]
        threshold = additData['-thres']
        thresVoxelModel = myMicModel.threshold(thresVoxelModel, [int(thres)], echo=True)
        thresVoxelModel = myMicModel.scaleModel(thresVoxelModel, 0, int(thres))
    if additData.has_key('-ldim'):
        ldim = additData['-ldim']
        if ldim[0] != ldim[1] or ldim[0] != ldim[2] or ldim[1] != ldim[2]:
            stdout.write('\n **ERROR** : Non-isotropic resolution not implemented!\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    else:
        stdout.write('\n **ERROR** : Voxel length (-ldim) not given in the input file')
        stdout.flush()
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if additData.has_key('-ndim'):
        ndim = additData['-ndim']
    else:
        stdout.write('\n **ERROR** : Voxel numbers (-ndim) not given in the input file')
        stdout.flush()
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    power2 = 4
    triaArray = myMiaModel.setupSphereTriangles(power2, True)
    if dtype.upper() == 'GST':
        M = myMiaModel.computeFabricTensorGST(thresVoxelModel, gstpow, norm, thres)
        EvalGST, EvectGST = linalg.eig(M)
        if sort == 'ON':
            EvectGST = numpy.transpose(EvectGST)
            EvalGST, EvectGST = myMiaModel.orderEvalEvec(EvalGST, EvectGST)
            EvectGST = numpy.transpose(EvectGST)
        approxGST, scaleGST = myMiaModel.computeApproxDistributionGST(M, triaArray, ftype, norm)
    else:
        BVTV = myMiaModel.computeBVTV(thresVoxelModel, threshold[0])
        if BVTV < 0.01:
            stdout.write('\n **ERROR** : BVTV < 1% - no morphological analysis done!\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if code == 'py':
            orgMIL, orgSVD, orgSLD, area = myMiaModel.computeOrigDistribution_MIL(thresVoxelModel, step, power, valid)
        elif code == 'f77':
            orgMIL, orgSVD, orgSLD, area = myMiaModel.computeOrigDistribution_F77(thresVoxelModel, step, power, valid)
        else:
            stdout.write('\n **ERROR** : Code type %s not implemented' % code)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dtype.find('MIL') > -1:
            approxMIL, scaleMIL = myMiaModel.computeApproxDistribution(orgMIL, triaArray, int(ftype), area, nameOfDist='MIL', normType=norm)
        if dtype.find('SLD') > -1:
            approxSLD, scaleSLD = myMiaModel.computeApproxDistribution(orgSLD, triaArray, int(ftype), area, nameOfDist='SLD', normType=norm)
        if dtype.find('SVD') > -1:
            approxSVD, scaleSVD = myMiaModel.computeApproxDistribution(orgSVD, triaArray, int(ftype), area, nameOfDist='SVD', normType=norm)
        if dtype.find('MIL') > -1:
            EvalMIL, EvectMIL = myMiaModel.computeEigenValueAndVector(orgMIL, int(ftype), area, normType=norm)
            if sort == 'ON':
                EvectMIL = numpy.transpose(EvectMIL)
                EvalMIL, EvectMIL = myMiaModel.orderEvalEvec(EvalMIL, EvectMIL)
                EvectMIL = numpy.transpose(EvectMIL)
        if dtype.find('SLD') > -1:
            EvalSLD, EvectSLD = myMiaModel.computeEigenValueAndVector(orgSLD, int(ftype), area, normType=norm)
            if sort == 'ON':
                EvectSLD = numpy.transpose(EvectSLD)
                EvalSLD, EvectSLD = myMiaModel.orderEvalEvec(EvalSLD, EvectSLD)
                EvectSLD = numpy.transpose(EvectSLD)
        if dtype.find('SVD') > -1:
            EvalSVD, EvectSVD = myMiaModel.computeEigenValueAndVector(orgSVD, int(ftype), area, normType=norm)
            if sort == 'ON':
                EvectSVD = numpy.transpose(EvectSVD)
                EvalSVD, EvectSVD = myMiaModel.orderEvalEvec(EvalSVD, EvectSVD)
                EvectSVD = numpy.transpose(EvectSVD)
    parts = inName.split('/')
    nParts = len(parts)
    inName2 = parts[nParts - 1]
    if mod == None:
        stdout.write(' ... write output \n')
        stdout.flush()
        OS = open(outName, 'w')
        myMiaModel.writeHeader(OS, step, power, ftype, valid, dtype=dtype)
        OS.write('\n Filename      = %s' % inName2)
        OS.write('\n nx,ny,nz      = %s' % repr(ndim))
        if not dtype.upper() == 'GST':
            OS.write('\n Threshold     = %s' % repr(threshold))
            OS.write('\n BV/TV         = %6.4f' % BVTV)
        if dtype.find('MIL') > -1:
            myMiaModel.writeEigenValueAndVectors('MIL', EvalMIL, EvectMIL, norm, OS)
        if dtype.find('SLD') > -1:
            myMiaModel.writeEigenValueAndVectors('SLD', EvalSLD, EvectSLD, norm, OS)
        if dtype.find('SVD') > -1:
            myMiaModel.writeEigenValueAndVectors('SVD', EvalSVD, EvectSVD, norm, OS)
        if dtype.find('GST') > -1:
            myMiaModel.writeEigenValueAndVectors('GST', EvalGST, EvectGST, norm, OS)
        OS.close()
    elif mod.lower() == 'a' or mod.lower() == 'w':
        Eval = numpy.array(numpy.zeros((3, 1), numpy.float))
        Evec = numpy.array(numpy.zeros((3, 3), numpy.float))
        if dtype.find('MIL') > -1:
            Eval = EvalMIL
            Evec = EvectMIL
        if dtype.find('SLD') > -1:
            Eval = EvalSLD
            Evec = EvectSLD
        if dtype.find('SVD') > -1:
            Eval = EvalSVD
            Evec = EvectSVD
        if dtype.find('GST') > -1:
            Eval = EvalGST
            Evec = EvectGST
        if mod == 'a':
            osfil = open(outName, 'a')
        elif mod == 'w':
            osfil = open(outName, 'w')
            osfil.write('#filename;step;power;ftype;valid;dtype;code;normType;BVTV[%];lam1;lam2;lam3;m1x;m1y;m1z;m2x;m2y;m2z;m3x;m3y;m3z\n')
        if dtype.find('GST') > -1:
            osfil.write('%s;%s;%g;%s;%s;%s;%s;%s;%s;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n' % (
             inName2, '-', gstpow, '-', '-', dtype, '-', norm, '-',
             Eval[0], Eval[1], Eval[2],
             Evec[0][0], Evec[1][0], Evec[2][0],
             Evec[0][1], Evec[1][1], Evec[2][1],
             Evec[0][2], Evec[1][2], Evec[2][2]))
        else:
            osfil.write('%s;%i;%i;%i;%g;%s;%s;%s;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g\n' % (
             inName2, step, power, ftype, valid, dtype, code, norm, BVTV * 100.0,
             Eval[0], Eval[1], Eval[2],
             Evec[0][0], Evec[1][0], Evec[2][0],
             Evec[0][1], Evec[1][1], Evec[2][1],
             Evec[0][2], Evec[1][2], Evec[2][2]))
    else:
        stdout.write('\n ** ERROR **: -mod has to be None, w, a arguments! \n ')
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if fout:
        OS = open(fout, 'w')
        if ftype == 2 or ftype == 4:
            if dtype.find('MIL') > -1:
                OS.write('** MIL - Tensors,  M = I*g + G **\n')
                g, G, G4 = myMiaModel.computeSpectralDecomposition(orgMIL, area)
            if dtype.find('SLD') > -1:
                OS.write('** SLD - Tensors,  M = I*g + G **\n')
                g, G, G4 = myMiaModel.computeSpectralDecomposition(orgSLD, area)
            if dtype.find('SVD') > -1:
                OS.write('** SVD - Tensors,  M = I*g + G **\n')
                g, G, G4 = myMiaModel.computeSpectralDecomposition(orgSVD, area)
            OS.write('g = %s\n\n' % repr(g))
            for i in range(3):
                for j in range(3):
                    OS.write('G2_%s%s = %s\n' % (repr(i), repr(j), repr(G[i, j])))

            OS.write('\n')
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            OS.write('G4_%s%s%s%s = %s\n' % (repr(i), repr(j), repr(k), repr(l), repr(G4[i, j, k, l])))

        if ftype == 1:
            if dtype.upper() == 'GST':
                M = myMiaModel.computeFabricTensorGST(thresVoxelModel, gstpow, norm, thres)
                OS.write('** GST - Fabric Tensor **\n')
                for i in range(3):
                    for j in range(3):
                        OS.write('M_%s%s = %s\n' % (repr(i), repr(j), repr(M[i, j])))

            else:
                if dtype.find('MIL') > -1:
                    H = myMiaModel.computeFabricTensorsEllipsoid(orgMIL)
                    OS.write('** MIL - H Tensor,  M = 1/sqrt(H) **\n')
                if dtype.find('SLD') > -1:
                    H = myMiaModel.computeFabricTensorsEllipsoid(orgSLD)
                    OS.write('** SLD - H Tensor,  M = 1/sqrt(H) **\n')
                if dtype.find('SVD') > -1:
                    H = myMiaModel.computeFabricTensorsEllipsoid(orgSVD)
                    OS.write('** SVD - H Tensor,  M = 1/sqrt(H) **\n')
                for i in range(3):
                    for j in range(3):
                        OS.write('H_%s%s = %s\n' % (repr(i), repr(j), repr(H[i, j])))

        OS.close()
    if distr:
        if distr.upper() == 'ON':
            filename, ext = outName.split('.')
            if dtype.find('MIL') > -1:
                myMiaModel.writeDistribution(filename, 'MIL', orgMIL)
            if dtype.find('SLD') > -1:
                myMiaModel.writeDistribution(filename, 'SLD', orgSLD)
            if dtype.find('SVD') > -1:
                myMiaModel.writeDistribution(filename, 'SVD', orgSVD)
    if graph:
        if graph.upper() == 'ON':
            name, ext = outName.split('.')
            if dtype.find('MIL') > -1:
                myMiaModel.writeParaview(name, 'MIL', orgMIL, approxMIL, triaArray, scaleMIL, EvectMIL, EvalMIL)
            if dtype.find('SLD') > -1:
                myMiaModel.writeParaview(name, 'SLD', orgSLD, approxSLD, triaArray, scaleSLD, EvectSLD, EvalSLD)
            if dtype.find('SVD') > -1:
                myMiaModel.writeParaview(name, 'SVD', orgSVD, approxSVD, triaArray, scaleSVD, EvectSVD, EvalSVD)
            if dtype.find('GST') > -1:
                myMiaModel.writeParaview2(name, 'GST', approxGST, triaArray, scaleGST, EvectGST, EvalGST)
    endTime = clock()
    ctime2 = time()
    stdout.write('\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime - startTime, ctime2 - ctime1))
    stdout.flush()
