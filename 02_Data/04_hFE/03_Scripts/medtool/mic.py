# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 19 2016, 06:48:10) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/mic/mic.py
# Compiled at: 2017-02-19 17:10:14
"""
Processor
---------

This module contains a comprehensive set of 3d image processing algorithms
(conversion, editing, segmentation, ...). It takes 3d medical CT image data sets
(voxel models) as input. Such a voxel data set (3d image) is modified by selected
filters. The modified voxel data set can be written to a new 3D image data file.
Also triangulated surface and Finite Element Model output can be generated
and stored in files. 

The software is mainly based on numpy and Fortran routines (f2py). The internal
data type is a float32 (4 bytes). Depending on the filter the memory requirement
is up to 20 times of the size of the input image. 

The processing order of the filters is the same as the order of the optional
parameters given in the GUI. The -mfil options allows to give a user defined 
order of the filters.

The software is designed to provide/get data from VTK / ITK. Thus, the
input/output format is usually '.mhd' and surfaces/grids are written in a format
which can be read with ParaView. 

All filters which are indictated by a ``*`` at the end of the description 
modify or uses the physical offset of an image.


Usage 
~~~~~

Module: ::

  import mic

Command line: :: 

  python mic.py ... 

    -in       filename
    [-out]    filename                                          [optional]
    [-form]   format                                            [optional]
    [-mid]    filename                                          [optional]
    [-imr]    startId;step;endId                                [optional]
    [-raw]    binFormat;headerSize;endianess;fastestdirection   [optional]
    [-ldim]   len1;len2;len3                                    [optional]
    [-ndim]   nVox1;nVox2;nVox3                                 [optional]
    [-temp]   filename                                          [optional]
    [-muscal] scaleFactor                                       [optional]
    [-smooth]  niter;lambda;kPB;interface;boundary;jacobian     [optional]
    [-re2d3] filename                                         [optional]
    [-mesh2d4] filename                                         [optional]
    [-sing]   filename;direction;sliceId                        [optional]
    [-geom]   type;grayValue;voxelSize;filename;diameter;height [optional]
    [-flip]   axis1;axis2;axis3;...                             [optional]
    [-mask]   filename                                          [optional]
    [-bool]   operator;filename                                 [optional]
    [-repl]   filename;elset;grayvalue                          [optional]
    [-arith]  operation1;operation2;...                         [optional]
    [-scale]  newMin;newMax;[oldMin];[oldMax];[format]          [optional]
    [-rota]   axi1;axi2;axi3;angle1;angle2;angle3;interpolate   [optional]
    [-rotm]   R11;R12;R13;R21;R22;R23;R31;R32;R33;interpolate   [optional]
    [-rotf]   filename;interpolate                              [optional]
    [-fill]   threshold;valid;type;kernel;[minThick];...        [optional]
    [-cont]   threshold                                         [optional]
    [-avg]    filename;phyWinSize;thres;maskValue               [optional]
    [-laham]  weight;cutoff;amplitude                           [optional]
    [-sobel]  None                                              [optional]
    [-mean]   kernel;thres1;thres2                              [optional]
    [-median] kernel;thres1;thres2                              [optional]
    [-gauss]  radius;sigma;thres1;thres2                        [optional]
    [-morph]  radius;type;shape;thres                           [optional]
    [-morpk]  kernel;type;thres                                 [optional]
    [-grad]   None                                              [optional]
    [-lapl]   None                                              [optional]
    [-cfill]  nlevels                                           [optional]
    [-cut]    n0Vox1;n0Vox2;n0Vox3;dnVox1;dnVox2;dnVox3         [optional]
    [-bbcut]  threshold;extVox1;extVox2;extVox3                 [optional]
    [-roicut] x0;y0;z0;x1;y1;z1                                 [optional]
    [-mcut]   dnVox1;dnVox2;dnVox3                              [optional]
    [-res2]   res1;res2;res3                                    [optional]
    [-res3]   len1;len2;len3                                    [optional]
    [-resf]   factor                                            [optional]
    [-refi]   direction;factor                                  [optional]
    [-mirr]   axis1;axis2;axis3;...                             [optional]
    [-mir2]   nVox1;nVox2;nVox3                                 [optional]
    [-autot]  BVTV;error;estimate                               [optional]
    [-fixt]   thresRatio;minimum;maximum                        [optional]
    [-slevt]  threshold                                         [optional]
    [-dlevt]  thres1;thres2                                     [optional]
    [-lthres] type;alpha;LT;UT                                  [optional]    
    [-gthres] type                                              [optional]       
    [-thres]  thres1;thres2;thres3;...                          [optional]
    [-clean]  type                                              [optional]
    [-extend] direction;thickVox;[newGrayvalue]                 [optional]
    [-close]  direction;threshold;kernel;newGrayValue           [optional]
    [-embed]  direction;thickVoxIn;thickVoxOut;newGrayValue     [optional]
    [-cap]    direction;thickVox;newGrayValue                   [optional]
    [-cover]  thickVox;newGrayValue                             [optional]
    [-block]  nVox1;nVox2;nVox3;newGrayValue                    [optional]
    [-mfil]   filter1;;par1;par2;;;filter2;;par3                [optional]    
    [-cen]    threshold;filename;mode                           [optional]
    [-bbox]   type;fileName;mode;i;j;k                          [optional]
    [-cog]    newGrayValue                                      [optional]
    [-ifo]    fileName;mode                                     [optional]
    [-histo]  filename;normalize;nIntervals                     [optional]
    [-shist]  fitNo;filename;mode                               [optional]

Parameters
~~~~~~~~~~

-in     : filename 

          Read a image file (voxel data) into the memory. Note that internally the voxel model 
          it is converted to 4 byte float. The file type is recognized by the file extension. 
          Possible extensions are:
          
          * '.mhd' : MetaImage file(ITK). Consist of a meta data file '.mhd' and a
            binary '.raw' file. The mhd file contains all important informations. It is
            readable by a text editor. The data itself are stored in a  separate raw 
            file. This or the '.nhdr' is the recommended file format.
          * '.nhdr' : NRRD image file composed of a header '.nhdr' and a binary 
            '.raw' file like the mhd format. This is the recommended file format if 
            transformations are important and Slicer3D is used. 
          * '.isq' : SCANCO 'isq' file format (primary format of SCANCO scanners) 
          * '.aim' : SCANCO 'aim' file (common for als SCANCO scanners). Only for small 
             files < 2 GB. 
          * '.png/.jpg/.gif/.tif/.bmp' : Single image or stack of image files. 
            In case of a stack (multiple images) a wildcard '#' have to be used in 
            the file name for every counter digit. For example if files with three 
            counter digits are to be read a file name could look like: test-###.png              
            Images are dedected automatically and simply sorted by filenames by 
            using Python's sort function. Note that the image content nor the image 
            numbers are checked. 
            A more save way is to give the 'imr' option which specifies the 
            start;step;end number of the image. 
            If no wildcards are given, a single image file is read. 
          * '.bin' : ISOSURF file format. It has a 0 Byte Header and contains short integer 
            data. See http://mi.eng.cam.ac.uk/~gmt11/software/isosurf/isosurf.html
            x is the fastest direction, than y and z. 'ndim' is required for this option!
          * '.raw' : Standard binary raw file. 'ndim' option is required for this option! 
            '.raw' reading option is required for this output format e.g. header, voxel 
            size, voxel number, and endianess have to be given! 
            Other formats like BST are simple raw formats which are used in the way
            like this option.


-out    : filename 

          Write image file (voxel data) to a new file. Use the 'form' option to control the 
          size of the data type of a singe voxel value. The file type is recognized 
          by the file extension. Possible file formats are (for more detailed 
          explanation see 'in' option):
                   
          * '.mhd' : ITK meta image file composed of a header and raw data file.
          * '.nhdr' : NRRD image file composed of a header and raw data file.
          * '.aim' : SCANCO 'aim' file.
          * '.png/.jpg/.gif/.tif/.bmp' : Single or Stack of image files.
          * '.bin' : ISOSURF file format. 
          * '.raw' : Standard binary raw file. 
          * '.svx' : Shapeways file format for 3D printing
          * '.inr.gz' : Inria image format. Binary file with header in gzip format. 
            A '.inr' format is working in Linux but not Windows. The following 
            parameters are fixed: ::  
            
              TYPE=float 
              PIXSIZE=32 bits
              SCALE=2**0 
              CPU=pc
         
          * '.inp' : Abaqus/Calculix input file format. Note: 
              - Needs as input a gray-values from 0..255. 
              - GV=0 are per default pores
              - GV=1..250 are used for bone tissue
              - GV=251..255 are used for other materials
                (polymer, implants, ...)           
                
            If the 'temp' option is not specified, all voxels (except those with GV=0) 
            are converted to FE elements and the same  material card is written for 
            all of them. Use '.inp' file format in combination with 'temp' option to 
            control which voxel should be converted and which element should get 
            which material card (gray level dependent material). For more details on 
            the format see http://www.calculix.de/
             
            EXAMPLES:
  
            .. image:: images/test_out_inp_1.png

          * '.in' : OOFEM input file format. For details with the 'temp' option check 
            the '.inp' file format (above), it uses the same algorithms. For more 
            details on the formats see http://www.oofem.org 
          * '.fem' : Altair Optistruct format. All thresholds and bulk material data 
            are written into the fem file. Additionally a density file ('.sh') is 
            written. For more details on the formats see http://www.altair.de

          The general idea of FE output is shown in the following Figure: 
          
          .. image:: images/mic-fe.png


-form   : format

          Byte format of image file (voxel data) for RAW. Currently implemented
          B...uint8, H...unit16,  h...int16, i...int32, f...float32
          Note that the supported format depend on the platform 32 or 
          64 bit. Depending on the format the file has to be converted 
          to with 'scale' option first e.g. for uint8 gray-values from 
          0-255 are required. 


-mid   :  filename
 
          Write mid-plane images of image file (voxel data) in the selected 
          file format. The routine automatically appends "_XM","_YM","_ZM" at 
          the end of the output file. 
          Available file formats: '.png', '.jpg', '.gif', '.tif', '.bmp', 
          '.case' (Ensight Gold) file format. 
         
          EXAMPLES:
         
          Trabecular bone:
           
          .. image:: images/test_mic_mid_1-XM.png
             :scale: 200 %
          .. image:: images/test_mic_mid_1-YM.png
             :scale: 200 %
          .. image:: images/test_mic_mid_1-ZM.png
             :scale: 200 %
                
          Femur (clinical resolution):
           
          .. image:: images/test_mic_mid_3-XM.png
          .. image:: images/test_mic_mid_3-YM.png
          .. image:: images/test_mic_mid_3-ZM.png
              

-mfil    : filter1;;par1;par2;;;filter2;;par3
           
           This option allows to run most of the filters in an arbitrary 
           order. Note that first the active options (filters) are 
           applied on the voxel model. After this run the given filters in 
           this option are applied on the previously modified voxel model.
           It is recommended do not mix other filters with this option.
           This means if 'mfil' is active, no other options (except 'in')
           should be active). The option names and default values are 
           the same as for all options implemented in this module. 
           Note: Not all filters are supported with this option. 
           
           If you not using the GUI, the different filters are separated by ";;;". 
           The filter name is separated from the filter parameters with ";;". 
           The current voxel model is passed to the next filter without 
           storing the results. Use the 'out', 'mid', etc option to store 
           the results during a multi run in between.   


-imr     : startId;step;endId

           Image file read parameters. Only needed if image stacks are to 
           be read by given the start, step, and end image number.

           The 'start id', 'step', 'stop id' (first, step in, last image number) 
           have to be given.
           Example: If file0001.png file0003.png file0005.png should be
           read the format is 1;2;5 and the file name in 'in' option should
           be file####.png (using # as wild card). If single images are 
           read (no wildcards), this option is not needed. 


-raw  :  binFormat;headerSize;endianess;fastestDirection

         Raw specifiers for reading raw image data ('.raw' files).
         
           - 'binFormat': The supported binary format (based on python) like : ::
           
               B   unsigned integer  1 Byte
               h   short             2 Byte
               H   unsigned short    2 Byte            
               i   integer           4 Byte
               f   float             4 Byte
               1   bit               1 Bit (0/1)
             
           - 'headerSize' ... size of file header which should be overwritten
           - 'endianess'   ... endian type: 'big' or 'little'
           - 'fastestDirection'  ... fastest varying direction: 'x' or 'z'


-ldim  : len1;len2;len3

         Physical sizes of voxels (length in mm) in 1, 2, 3 direction.
         Not optional for Abaqus, Vista & Optistruct output! Note that 
         this option influences the physicals units of the FE input 
         files (written by this script).


-ndim  : nVox1;nVox2;nVox3

         Number of voxels in the corresponding 1, 2, 3 direction.
         Not optional for some input file readers (see 'in' option)!


-temp    : filename

           Name of an ABAQUS/Calculix/OOFEM Template file. This file
           allows you to control which gray values are converted to 
           elements, gray-level dependent material cards, user specific
           solver commands etc. Internal keywords (see below) are replaced
           by the script. All other text is simply copied to the newly 
           generated input deck. 
           
           Internal supported keywords are: ::
           
             *USER NODE
               -generate nodes, requires USER PROPERTY keyword
             *USER ELEMENT 
               -generate elements, requires USER PROPERTY keyword
             *USER NSET, type=point, location=arbitrary
               -generate NSET: ARB_NODE_S, ARB_NODE_N, ARB_NODE_E, 
                               ARB_NODE_W, ARB_NODE_T, ARB_NODE_B 
             *USER NSET, type=point, location=addcorner
               -generate NSET: ACOR_NODE_SWB, ACOR_NODE_SEB, ACOR_NODE_NEB, 
                               ACOR_NODE_NWB, ACOR_NODE_SWT, ACOR_NODE_SET, 
                               ACOR_NODE_NET, ACOR_NODE_NWT
             *USER NSET, type=face 
               -generate NSET: ALL_NODE_S, ALL_NODE_N, ALL_NODE_E, 
                               ALL_NODE_W, ALL_NODE_T, ALL_NODE_B 
             *USER ELSET, type=face
               -generate ELSET: ALL_ELEM_S, ALL_ELEM_N, ALL_ELEM_E, 
                                ALL_ELEM_W, ALL_ELEM_T, ALL_ELEM_B      
             *USER PROPERTY, file=property_temp.inp, [range=minVal:maxVal]
               -This line can be repeated multiple times (see example
                below) to give individal material properties for different
                gray-value ranges. 
               -file= is required and contains the property template file
                (see example below). Maybe it is needed to give also the b
                path to this file.
               -range= is optional. If not given everything (including the
                GV=0) is processed i.e. range=0:255. 
                If given gray-value from minVal till maxVal are processed 
                i.e. nodes, elements and material cards are genereated
                accordingly.  
               -internal variables are for the template file are 
                SetName, CardName, GrayValue
               -GV >= minVal and GV <= maxValue use the material card 
                given with file=... In case of more gray values ranges 
                they look e.g. like 1:69, 70:250, 251:251, 252:255

             Note: _S, _N, _E, _W, _T, _B standsfor South, North, East, 
             West, Top, Bottom). Usually not all *USER NSET and *USER 
             ELSET keywords are used. 
         
           The strings $filename and $pathfilename can be used inside 
           the template file to write the current inp filename without or 
           with the path to the inp file (no extension is written). 
           Only implemented in the case of ABAQUS and CALCULIX.

           Example files for CalculiX (similar for Abaqus) are given
           below. Store them in separate files like main_temp.inp, 
           property_temp_bone.inp, and property_temp_embed.inp to use them. ::
           
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ** main_temp.inp 
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             *HEADING
             main input $filename
             ** Nodal data from voxel model 
             *USER NODE
             ** Elements + elset from voxel model
             ** type=C3D8, elset="SetName" (see USER PROPERTY)
             *USER ELEMENT 
             ** Additional nset from voxel model. New generated nsets are:
             ** ARB_NODE_S, ARB_NODE_N, ARB_NODE_E, ARB_NODE_W, ARB_NODE_T, 
             ** ARB_NODE_B 
             *USER NSET, type=point, location=arbitrary
             ** Additional nset from voxel model. New generated nsets are:
             ** ACOR_NODE_SWB, ACOR_NODE_SEB, ACOR_NODE_NEB, ACOR_NODE_NWB,
             ** ACOR_NODE_SWT, ACOR_NODE_SET, ACOR_NODE_NET, ACOR_NODE_NWT
             ** Note: the nodes are not connected to the model! 
             *USER NSET, type=point, location=addcorner
             ** Additional nset from voxel model. New generated nsets are:
             ** ALL_NODE_S, ALL_NODE_N, ALL_NODE_E, ALL_NODE_W, ALL_NODE_T, 
             ** ALL_NODE_B
             *USER NSET, type=face
             ** Additional elset from voxel model. New generated nsets are:
             ** ALL_ELEM_S, ALL_ELEM_N, ALL_ELEM_E, ALL_ELEM_W, ALL_ELEM_T, 
             ** ALL_ELEM_B
             *USER ELSET, type=face
             ** User property (section & material) from  template file 
             ** internal variables are: "SetName", "CardName", "GrayValue"
             *USER PROPERTY, file=property_temp_bone.inp, range=1:250
             *USER PROPERTY, file=property_temp_embed.inp, range=251:255
             ***********************************************************
             *STEP
             *STATIC
             *BOUNDARY, TYPE=DISPLACEMENT
             ALL_NODE_B, 1, 3, 0
             *BOUNDARY, TYPE=DISPLACEMENT
             ALL_NODE_T, 1, 3, 0.1
             *NODE FILE 
             U
             *EL FILE 
             S
             *NODE PRINT, NSET=ALL_NODE_T
             RF
             *END STEP
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
           ::
           
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ** property_temp_bone.inp
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             *SOLID SECTION, ELSET=SetName, MATERIAL=CardName
             1.
             *MATERIAL,NAME=CardName
             *ELASTIC
             5400.*(GrayValue/250.)**2.5, 0.3
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
           ::
           
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             ** property_temp_embed.inp
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             *SOLID SECTION, ELSET=SetName, MATERIAL=CardName
             1.
             *MATERIAL, NAME=CardName
             *ELASTIC
             2000., 0.3
             ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


-muscal   : scaleFactor

            Scaling factor between the linear attenuation coefficients and 
            the gray-values of the image in case of writing SCANCO '.aim' 
            files. Examples :
         
              - Scanco Xtreme CT (AKH Vienna) 8192
              - Scanco muCT 40   (TU Vienna) 4096


-smooth   : niter;lambda;kPB;interface;boundary;jacobian

            Smooth mesh before output. This filter is only active if 
            Abaqus/Calculix/ParFE grids (meshes) are written. 
            The implementation follows Taubin's algorithm (see, Taubin, 
            Eurographics 2000, Boyd and Mueller 2006).
         
             - 'niter' is the number of iterations. One iteration means forward 
               (with lambda) and backward smoothing (with mu) 
               smoothing i.e. 2 smoothing steps.
             - 'lambda' scaling factor (0<lambda<1)
             - 'kPB' pass-band frequency kPB  (note kPB=1/lambda+1/mu > 0) => mu
             - 'interface' on/off switch for including near interface node smoothing.  
               This option is in the current Fortran code not implemented!!
             - 'boundary' boundary smoothing id, have to be provided with
             
                - boundary=0 ... no smoothing of boundary elements/nodes
                - boundary=1 ... smoothing of boundary elements/nodes but preserving 
                  cutting planes
                  
             - 'jacobian' minimal allowed scaled Jacobian. E.g. 0.02 means the 
               smallest allowed volume change is 2% of the initial volume. 
               
            EXAMPLES:
           
            .. image:: images/test_out_inp_3.png


-mesh2d3  : filename

            Write 2D triangular meshes of the interface between gray-values. 
            Only meaningful for segmented voxel models.  
            Use the the 'smooth' option to smooth the surface mesh. 
            Possible output formats are:
            
             * geom       : TAHOE II format 
             * geo        : PARAVIEW : Ensight gold format (same as case)
             * case       : PARAVIEW : Ensight gold format (same as geo)
             * inp        : ABAQUS input file
             * off        : Geomview file
             * msh        : GMSH file

            EXAMPLES:
            
            .. image:: images/test_mesh2d3_1.png
            
            .. image:: images/test_mic_mid_2-YM.png
               :scale: 200 %


-mesh2d4  : filename

            Write 2D quadrilateral meshes of the interface between gray-values. 
            Only meaningful for segmented voxel models.  
            Use the the 'smooth' option to smooth the surface mesh. 
            Possible output formats are:
            
             * geom       : TAHOE II format 
             * geo        : PARAVIEW : Ensight gold format (same as case)
             * case       : PARAVIEW : Ensight gold format (same as geo)
             * inp        : ABAQUS input file
             * off        : Geomview file
             * msh        : GMSH file

            EXAMPLES:
           
            .. image:: images/test_mesh2d4_1.png
            
            .. image:: images/test_mic_mid_2-YM.png
               :scale: 200 %


-sing  : filename;direction;sliceId

         Write single images in the selected file format (known by filename 
         extension). 

         - 'filename' The filename where the routine  automatically appends 
           the direction and the slice number (e.g. filname_1+_21.png) at the 
           end of the output file. Available file formats: png, jpg, gif. 
         - 'direction' The direction (1+,2+,3+,1-,2-,3-) and the slice number 
           in this direction. Positive direction mean one is looking in the 
           direction of the axis.
         - 'sliceId' The number of the required slices starting with '1'.


-geom  : type;grayValue;voxelSize;filename;diameter;height

         Writes a geometrical object as voxel model to the output. This 
         option ignores the 'in' options. 
         
         Currently implemented objects are:

          - sphere   : type='sph', 'diameter' is the diameter, height is not given 
          - cylinder : type='cyl', 'diameter' is the diameter, height the 
            height of the cylinder. 
          
         The 'grayValue' is the voxel gray scale value of the geometrical object. 
         The 'voxelSize' is the physical voxel length. The remaining dimensions are
         given as number of voxels. 


-flip  : axis1;axis2;axis3;...

         Flip the model around the given flip axes: The size of the model 
         is not changed. 
         
         EXAMPLES:
             
         .. image:: images/test_mic_mid_1-XM.png 
            :scale: 300 %
         .. image:: images/test_flip_1-XM.png
            :scale: 300 %


-mask  : filename

         Mask voxel model with the given second voxel model. File have to
         be segmented and scaled to 0/1 i.e. exactly two gray values (0/1)!
         Instead the 'arith' option can be use to obtain the same results.

         EXAMPLES:
             
         .. image:: images/test_mic_mid_1-XM.png 
            :scale: 300 %
         .. image:: images/test_mask_1-XM.png
            :scale: 300 %


-bool  : operator;filename

         Add a second voxel model to the existing one. Compared to 'arith'
         and 'mask' the two voxel models do not need to have the same dimensions.
         The 'operator' says how the models are combined. Currently 
         implemented ``+,-,*,r`` (r..replace if gray value larger than 0). 
         If the voxel model is given as 'mhd' file, the offset will be taken 
         into account.

         EXAMPLES:

         .. image:: images/test_mic_mid_3-XM.png 
            :scale: 150 %
         .. image:: images/test_bool_1-XM.png
            :scale: 150 %


-repl  : filename;elset;grayvalue;resolution

         Replace the gray-values of a voxel model based on a meshed domain 
         given by 'filename'. The 'elset' is a string of the corresponding
         element set (name or id) or it can be 'ALL' if the whole FE mesh
         should be used. The 'grayvalue' is used for the new gray-value 
         inside the voxel model. If 'grayvalue'='None' than the gray-value 
         is generated from the material IDs starting after the highest 
         grayvalue in the model and steping by one. 
         This function can be used to insert an implant into a bone.

         EXAMPLES:
             
         .. image:: images/test_mic_mid_3-XM.png 
            :scale: 150 %
         .. image:: images/test_repl_1-XM.png
            :scale: 150 %

         If the 'resolution' is given, a new voxel model is generated 
         with zero background and the only the meshed domain is marked. 
         The voxel model given by -out is ignored. It can be used to 
         create a 3D print model from an FE mesh. 

         EXAMPLES:
             
         .. image:: images/test_repl_2.png
            :scale: 80 %

         The meshed domain have to be given in a format which can be 
         read by the finite elemente converter 'fec' (e.g. inp, fem, ...).

         Linear solid elements namely:

           - tetrahedrons (tet4), 
           - hexahedrons (hex8), 
           - pentahedron (penta6)  
           
         lineare structural elements are implemented:      
        
           - tapered beams (bar2)
           - triangles (tria3)
           - quadriliterals (quad4) 

         In case of structural elements the thickness information is read from 
         the property entry. The beam type can be 'CIRC' in case of Abaqus 
         and 'ROD' in case of Optistruct. If elements are outside the given 
         voxel model a warning will be thrown. 


-arith  : operation1;operation2;...

          Apply an equation/formula on the data array. Implemented operators 
          are: :: 
         
             +A  ... add value A to loaded voxelmodel
             -A  ... subtract value A from loaded voxelmodel
             *A  ... multiply with A
             /A  ... divide by A
             <B  ... store current value in B
             >B  ... load B in current value
             ^I  ... power I (only integers!) of current value
             >>C ... load file with name C into memory (mhd, nhdr, or isq)
             ?   ... logical conditions like ?<0=0 or ?==3=2 or ?>7=12
             &sin ... compute sinus of array in memory
             &cos ... compute cosinus of array in memory
             &arcsin ... compute arcus sinus of array in memory
             &arccos ... compute arcus cosinus of array in memory
             &tan    ... compute tangens of array in memory
             &arctan ... compute arcus tangens of array in memory
          
          EXAMPLES:
             
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 300 %
          .. image:: images/test_arith_1-XM.png
             :scale: 300 %


-scale  : newMin;newMax;oldMin;oldMax;format

          Scale the min/max values of a voxel model. If old values (PARAM
          "oldMin;oldMax") are not given the min and max values of the original
          voxel model will taken instead of given values. The string 'oldMin' or 
          'oldMax' can be used in old or new values. 
          (Example oldmin;255;70;oldmax). If the format specifier is given the 
          return format of the filter can be controlled. The default format is 
          the input format of the voxel array.


-rota  : axi1;axi2;axi3;angle1;angle2;angle3;interpolate

         Rotate the voxel model by rotation type 'angle':

             - axi1;axi2;axi3 = Rotation order around axis. For example 3;1;2
                            means first the model is rotated around 3-axis, 
                            second around 1 axis and third around 2 axis 
             - angle1;angle2;angle3 = Rotation angles in Deg around 1,2,3 axis. That
                            means ang1 will be applied in the previous example
                            as second rotation. 
             - interpolate    = YES/NO if YES data will be interpolated (preferable)
             - Note: if the angles are taken from sop.py (Stiffness Optimizer) the 
               model have to be rotated such that rotation order = 1;2;3 and 
               rotation angle = x;y;z (from sop)

         Note: The image is increased in size and shifted such that the rotated 
         model fits fully inside. The offset is updated accordingly such that 
         the new image is still in the correct position in the world coordinate 
         system. This can be checked e.g. with Slicer. 

         EXAMPLES:
               
         .. image:: images/test_mic_mid_1-XM.png 
            :scale: 300 %
         .. image:: images/test_rot_1-XM.png
            :scale: 300 %


-rotm  : R11;R12;R13;R21;R22;R23;R31;R32;R33;interpolate 

         Rotate the voxel model by rotation type 'matrix':

             - Rij         = coefficients of the rotation matrix
             - interpolate = YES/NO if YES data will be interpolated (preferable)

         Note: The image is increased in size and shifted such that the rotated 
         model fits fully inside. The offset is updated accordingly such that 
         the new image is still in the correct position in the world coordinate 
         system. This can be checked e.g. with Slicer. 


-rotf  : filename;interpolate 

         Rotate the voxel model by using transformation file. Implemented 
         'filename' can be : 

            - 'tfm': ITK transformation file. For example, exported  
              from 3D Slicer. 

            - 'nhdr': NRRD image header file with transformation 
              information in it. The important informations in this file 
              are 'space directions' which gives indirectly the transfromation 
              matrix and 'space origin' which gives the offset. 
              
            - 'mhd': ITK meta image header file with transformation 
              information in it. The important informations in this file 
              are 'TransformMatrix' which gives the transfromation matrix and 
              'Offset' which gives the offset. The 'CenterOfRotation' has to be 
              (0,0,0). 

         For example, this can the the file given in the 'in' option if 
         the applied transformation matrices are inside. More details about 
         the transformation is given in the *Mesh - Converter* module (fec). 

         Note: The image is increased in size and shifted such that the rotated 
         model fits fully inside. The offset is updated accordingly such that 
         the new image is still in the correct position in the world coordinate 
         system. This can be checked e.g. with Slicer. 


-fill  : threshold;valid;type;kernel;[minThick];[niter];[kernel2];[bx0];[bx1];[by0];[by1];[bz0];[bz1]

         Fill the pores of a segmented image. Details can be found in :
         Pahr and Zysset, Comput Methods Biomech Biomed Engin, 2009, 12, 45-57
         Parameters are described in this paper and are:

          - 'threshold' the threshold value which is used to segment the image
          - 'valid' is also some kind of a threshold. The algorithm counts 
            how often a voxel was marked as a fill voxel. Maximum marks are 
            seven. E.g. 'valid'=5 means that also voxels which are five times 
            marked are treated as fill voxels. 
          - 'type' gives the type of the algorithm. Implemented options are:

            - out ... outer surface (param1 = smooth)
            - in  ... inner surface (param1 = smooth)
            - t   ... cortical shell (param1 = smooth)
            - v   ... value from filling (for test purposes!)
                    (param1 = smooth)
            - c   ... contour of filling (thickness is given by min_t)
            - Add '_2d?' with ?=1,2,3 if 2D filling should be selected. 
              The "?" gives the direction of the normal of the 2D plane. 
              E.g. 'out_2d3' means outer surface fill in 2D where the 3-axis
              is normal to this plane. If '_2d?' is not given the 3d 
              search option is used as default.
            
          - 'kernel' : size (radius) of smooth kernel in voxels dimensions 
          - 'minthick': optional - minimal cortex thickness, only used in 'type'=c
          - 'niter' : optional - maximum number of fill error correction iterations. 
             if not given or set to 0 no correction will be done. 
          - 'kernel2': optional -size (radius) of kernel in voxels for the fill
            error correction steps. If not given it is the same as "kernel".
          - 'bx0;bx1;by0;by1;bz0;bz1' : optional - bounding box which is used to 
            apply the fill error correction.

         EXAMPLES:
             
         Thickness:
           
         .. image:: images/test_mic_mid_3-XM.png 
            :scale: 150 %
         .. image:: images/test_fill_1-XM.png
            :scale: 150 %
              
         Outside (with 2d3, without):
           
         .. image:: images/test_fill_2-XM.png
            :scale: 150 %
         .. image:: images/test_fill_3-XM.png
            :scale: 150 %
                
         Outside with correction (full, bounding box):
           
         .. image:: images/test_fill_4-XM.png
            :scale: 150 %
         .. image:: images/test_fill_5-XM.png
            :scale: 150 %
            
            
-cont   :  threshold 

           Extracts the contour based on the given threshold. All voxels >=
           threshold will be included in the evaluation. Edges and corners are 
           not implemented in the contour search. 
           
           EXAMPLES:
             
           .. image:: images/test_mic_mid_1-XM.png 
              :scale: 300 %
           .. image:: images/test_cont_1-XM.png
              :scale: 300 %
        

-avg    :  filename;phyWinSize;thres;maskValue

           Averages BVTV or gray values (selected by 'threshold')
           over a 'phyWindowSize' and writes the output to a text file 
           given by 'filename'. Output variables if 'threshold'='None' are :: 

              nx0;ny0;nz0;dnx;dny;dnz;lx;ly;lz;ZV;TV;grayValueAverage

           otherwise 'threshold'=somevalue ::
         
              nx0;ny0;nz0;dnx;dny;dnz;lx;ly;lz;ZV;TV;BVTV
 
           The algorithm starts at the physical 0/0/0 location and steps with 
           'phyWindowSize' (cube) over the image. It selects all voxels which 
           have there COG inside the window.

           If the 'maskValue' is given (not ='None') it will be considered in 
           the computations e.g BVTV = BV/(TV-ZV) where ZV = maskValue.


-laham  :  weight;cutoff;amplitude

           Laplace-hamming filter: Transfers image to the frequency domain and  

             - applies the laplace operator
             - weights the laplacian filtered image with parameter 'weight'
             - multiplies  the original image with 1-'weight'
             - sum up both images 
             - multiplies the resulting image with a low pass filter. 
             
           The low pass filter follows the formula ampl='amplitude' ::
           
             A(k) = (1-ampl/2) - ampl/2*cos(2*pi*k/(M*cutoff-1))
             
           Some characteristics are:
           
             - ampl = 0 no low pass filtering is done
             - ampl = 1 and cutoff=0 ... Hanning filter
             - ampl = 0.92 and cutoff=0 ... Hamming filter
             
           where ::

             cutoff*Nyquist_frequency = cut_off_frequency

           Note thatboth parameters are ratios and have to be between 0..1.

           EXAMPLES:
             
           .. image:: images/test_mic_mid_1-XM.png 
              :scale: 300 %
           .. image:: images/test_laham_1-XM.png
              :scale: 300 %

-sobel  :  None 

         Sobel filter (only a dummy argument is required but not used). This 
         is a classical edge detection filter

         EXAMPLES:
             
         .. image:: images/test_mic_mid_1-XM.png 
            :scale: 300 %
         .. image:: images/test_sobel_1-XM.png
            :scale: 300 %


-mean   : kernel;thres1;thres2

          Mean filter: 'kernel' gives the type of the kernel window. 
          Only voxels  with a gray value >= 'thres1' and value <= 'thres2' 
          are considered in the computations. 

          Possible values for 'kernel' are : 

             - 'kernel' is a specifier of following type: 
                
                - "k3x3", "k5x5", "k7x7", ... which is a squared kernel 
                  window of size 3x3, 5x5, ... where the numbers have to 
                  be odd numbers which are seperated by an 'x'. 

                - "k6": This is a 3x3 kernel with 6 non-zero values. 
                  This is the same as "k1.0". 

                - "k18": This is a 3x3 kernel with 18 non-zero values. 
                  This is the same as "k1.732" (<sqrt(3)). 

                - "k?.??" where ?.?? is an abitrary float number > 1.0. 
                  This number is the radius of a sphere in voxel dimensions 
                  which will set up the non-zeros in the kernel window.
                  It is important to write '.' as decimal point. 
                  Examples are 'kernel'="k1.0", 'kernel'="k2.34", ... 
 
             - 'kernel' = number e.g. 'kernel'=3. In this case the 
               size of the kernel window in voxel dimensions is 
               n = int(2*('kernel'+0.5)). 

          The kernel window is cropped by the bounding box of the image.

          EXAMPLES:
             
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 300 %
          .. image:: images/test_mean_1-XM.png
             :scale: 300 %


-median : kernel;thres1;thres2

          Median filter: 'kernel' gives the type of the kernel window. 
          Only voxels  with a gray value >= 'thres1' and value <= 'thres2' 
          are considered in the computations. 

          Possible values for 'kernel' are : 

             - 'kernel' is a specifier of following type: 
                
                - "k3x3", "k5x5", "k7x7", ... which is a squared kernel 
                  window of size 3x3, 5x5, ... where the numbers have to 
                  be odd numbers which are seperated by an 'x'. 

                - "k6": This is a 3x3 kernel with 6 non-zero values. 
                  This is the same as "k1.0". 

                - "k18": This is a 3x3 kernel with 18 non-zero values. 
                  This is the same as "k1.732" (<sqrt(3)). 

                - "k?.??" where ?.?? is an abitrary float number > 1.0. 
                  This number is the radius of a sphere in voxel dimensions 
                  which will set up the non-zeros in the kernel window.
                  It is important to write '.' as decimal point. 
                  Examples are 'kernel'="k1.0", 'kernel'="k2.34", ... 
 
          The kernel window is cropped by the bounding box of the image.

          EXAMPLES:

          .. image:: images/test_mic_median_sp_1-XM.png
             :scale: 300 %
          .. image:: images/test_mic_median_1-XM.png
             :scale: 300 %


-gauss  : radius;sigma;thres1;thres2

          Gaussian filter: Filter computes the Gaussian smoothing of all 
          voxels around a considered voxel within the 'radius' by using 
          'sigma'. Only voxel with a value >= 'thres1' and value <= 'thres2' 
          are considered.

          EXAMPLES:
             
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 300 %
          .. image:: images/test_gauss_1-XM.png
             :scale: 300 %


-morph  : radius;type;shape;thres
 
          Morhological filters: compute morphological opening ("o"), 
          closing ("c"), erosion ("e"), dilation ("d"). Similar to 
          'morpk' but with different kernel definitions.  Parameters are:

            - 'radius' is a radius of specifying the kernel (in number of
              voxels). The size of the kernel is nk = int(2*(r+0.5)). 
              This means that a radius of "1" gives a kernel of 3x3x3 voxels.  
              In case of a spherical kernel everything which is inside 
              rc <= nk/2 is taken as kernel. 
            - 'type' specifies if opening (dilation + erosion) or closing 
              (erosion+dilation), or erosion or dilation will be performed. 
            - 'shape' specifies if a rectangle ("1") or spherical kernel 
              shape ("2") should be used.
            - 'threshold' is a value which is used to segment the image 
              before applying the filters 

          EXAMPLES:
             
          Original Images:

          .. image:: images/test_mic_mid_3-XM.png 
             :scale: 150 %
          .. image:: images/test_mic_mid_4-XM.png 

          Opening / Closing:
           
          .. image:: images/test_morph_1-XM.png
             :scale: 150 %
          .. image:: images/test_morph_2-XM.png
             :scale: 150 %
                
          Erosion / Dilation 
           
          .. image:: images/test_morph_3-XM.png
          .. image:: images/test_morph_4-XM.png


-morpk  : kernel;type;thres
 
          Morhological filters: compute morphological opening ("o"), 
          closing ("c"), erosion ("e"), dilation ("d"). Similar to 
          'morph' but with different kernel definitions. Parameters are:

             - 'kernel' is a specifier of following type: 
                
                - "k3x3", "k5x5", "k7x7", ... which is a squared kernel 
                  window of size 3x3, 5x5, ... where the numbers have to 
                  be odd numbers which are seperated by an 'x'. 

                - "k6": This is a 3x3 kernel with 6 non-zero values. 
                  This is the same as "k1.0". 

                - "k18": This is a 3x3 kernel with 18 non-zero values. 
                  This is the same as "k1.732" (<sqrt(3)). 

                - "k?.??" where ?.?? is an abitrary float number > 1.0. 
                  This number is the radius of a sphere in voxel dimensions 
                  which will set up the non-zeros in the kernel window.
                  It is important to write '.' as decimal point. 
                  Examples are 'kernel'="k1.0", 'kernel'="k2.34", ... 

            - 'type' is the type of the morphological operation. It can 
               be 'o', 'c', 'e', or 'd'.

            - 'thres' is a value which is used to segment the image 
              before applying the filters 


-grad   : None

          Gradient filter (only dummy argument required - not used).

          EXAMPLES:
             
          .. image:: images/test_mic_mid_4-XM.png 
          .. image:: images/test_grad_1-XM.png


-lapl   : None

          Laplacian filter (only dummy argument required - not used).

          EXAMPLES:
             
          .. image:: images/test_mic_mid_4-XM.png 
          .. image:: images/test_lapl_1-XM.png


-cfill  : nlevels

          Octree based filling of a segmented voxel model. A voxel model with 
          a single threshold (e.g. 0 and 75) has to be provided. The filter
          colors depending on the given 'nlevels' cubes of size 2*1, 2*2, ... 
          2*nlevels. The lower level has the gray value nlevels+1, the highest
          has a grayvalue of 1.  
         
          EXAMPLES:
             
          .. image:: images/test_mic_mid_4-XM.png 
          .. image:: images/test_cfill_1-XM.png


-cut    : n0Vox1;n0Vox2;n0Vox3;dnVox1;dnVox2;dnVox3

          Cut/crop a small part of a bigger model by using start voxel ids 
          ('n0Vox1', 'n0Vox2', 'n0Vox3') and the size of a ROI 
          ('dnVox1', 'dnVox2', 'dnVox3'). 
          This size will be the new size of the voxel model.
          i.e. three positions of the new origin (old origin at 0;0;0)
          and three extensions (number of voxels) = ROI.

          EXAMPLES:
             
          .. image:: images/test_mic_mid_4-XM.png 
          .. image:: images/test_cut_1-XM.png


-bbcut  : threshold;extVox1;extVox2;extVox3

          Cut/crop a small part of a bigger model by using an automatic
          bounding box search algorithm. The bounding box is detected
          using the 'threshold'. The three optional values 'extVox1', 
          'extVox2', 'extVox3' are the number of voxels in 1, 2, 3 which 
          should increase the bounding in all directions 

          EXAMPLES:
             
          .. image:: images/test_mic_mid_3-XM.png 
             :scale: 150 %
          .. image:: images/test_bbcut_1-XM.png
             :scale: 150 %


-roicut : x0;y0;z0;x1;y1;z1 

          Cut/crop a small part of a bigger model by using a ROI given 
          in physical dimensions. 'x0', y0', 'z0' are the start coordinates 
          and  'x1', y1', 'z1' are the end coordinates the the bounding 
          box of the ROI. Note the "one voxel" inaccuracy due to this 
          method. 

          If the image data contain an offset ('mhd' or 'nhdr') it is taken 
          into account. 


-mcut  :  dnVox1;dnVox2;dnVox3

          Cut/crop a ROI from the middle of a bigger voxel model. 
          Teh size of the model is given by the three extensions 
          'dnVox1', 'dnVox2', 'dnVox3' (number of voxels in 1,2,3).

          EXAMPLES:
             
          .. image:: images/test_mic_mid_3-XM.png 
             :scale: 200 %
          .. image:: images/test_mcut_1-XM.png
             :scale: 200 %


-res2  :  res1;res2;res3

          Change the resolution of the voxel data. 'res1;res2;res3' are 
          the new resolutions in 1,2,3 axis. 
          
          Avoid the usage of this filter because it is slow and interpolates 
          voxel data. 'resf' is recommended. 


-res3  :  len1;len2;len3

          Change the resolution of the voxel data by given voxel lengths. 
          'len1;len2;len3' are approximations of the new voxel lengths. They 
          are adopted such that the voxel model size will be unchanged.

          Avoid the usage of this filter because it is slow and interpolates 
          voxel data. 'resf' is recommended. 


-resf  :  factor

          Change the resolution (recoarse/refine) of the voxel data by the 
          given 'factor'. Note that the overall size of the voxel model might 
          change slightly in case of recoarsening because no interpolations 
          is done. 

          Positive 'factor' values recoarse the image, negative 'factor' 
          values refine the image by this factor.

          EXAMPLES:
             
          .. image:: images/test_mic_mid_4-XM.png 
          .. image:: images/test_resf_1-XM.png
             :scale: 500 %


-refi  :  direction;factor 

          Refine the image in one direction (given by 'direction' either 
          1,2, or 3) by a given factor. 
          This filter can be used to convert anisotropic QCT images to isotropic once.


          EXAMPLES:
             
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 200 %
          .. image:: images/test_refi_1-XM.png
             :scale: 200 %


-mirr  :  axis1;axis2;axis3;...

          Mirror the model around the given mirror axes. The size of the 
          model is doubled in the given direction. 

          EXAMPLES:
             
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 200 %
          .. image:: images/test_mirr_1-XM.png
             :scale: 200 %


-mir2  :  nVox1;nVox2;nVox3

          Mirror the model such that it fills a new array with dimension 
          'nVox1', 'nVox2', 'nVox3'. This filter is used to create arrays for FFT.

          EXAMPLES:
             
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 200 %
          .. image:: images/test_mir2_1-XM.png
             :scale: 200 %
              

-autot  :  BVTV;error;estimate

           Find the threshold for a given 'BVTV' value within a given 'error'.
           In order to increase speed an threshold 'estimate' should also be
           given. The required error


-fixt  :  thresRatio;minimum;maximum

          Fixed threshold applied on the image. The given parameter 'thresRatio' 
          is a ratio for the computations (see example below). The second 
          and third parameter specifies the 'minimum' and 'maximum' range 
          of the threshold. Allowed values numbers as well as strings like 
          'minVal' and 'maxVal'. In this case these values are computed from 
          the image. 
          Example (thresRatio=0.4, minVal = -100, maxVal = 2000) ::
          
            0.4;minVal;maxVal -> threshold = (2000-(-100))*0.4
            0.4;0.0;maxVal    -> threshold = (2000-0.0)*0.4
            0.4;0.0;200       -> threshold = (200-0.0)*0.4 
           
          This method is suggested in the literature to segment HR-pQCT images.

-slevt  :  threshold

           Threshold using the iterative selection method of Ridler
           (single level thresholding) after Scanco code. The given 
           'threshold' is used to decide which values should be considered
           in the computation. E.g. if a value of 1 is given all values 
           lower than 1 are not considered.
           This is the recommended segmentation option for high-resolution
           micro-CT images. 

           EXAMPLES:
           
           .. image:: images/test_mic_mid_4-XM.png 
           .. image:: images/test_slevt_1-XM.png

-dlevt  :  thres1;thres2

           Threshold using double level thresholding. Similar to 'slevt' 
           but two threshold values are computed.
           The given threshold should be chosen properly and are used in 
           the switching function.

           EXAMPLES:
           
           .. image:: images/test_mic_mid_3-XM.png 
              :scale: 150 %
           .. image:: images/test_dlevt_1-XM.png
              :scale: 150 %
                             
                                         
-lthres  : type;alpha;LT;UT

           Local adaptive threshold computed from global measures of the 
           image. The kernel is a 3x3x3 cube. Parameters are:
           
            - 'type' : "stddev": similar to the method proposed by Kang and 
              Engelke 2003.  The criterion is ::

               voxValue < (mean-alpha*stddev) --> marrow voxel 

              where voxValue is the grayvalue of the cube's center, mean and 
              stddev (n-1) is computed from the kernel. 
                  
            - 'alpha' : This is the scaling for the weighting 
           
            - 'LT' : lower threshold. Everything below is marrow. If not 
               given the threshold is automatically estimated by the 
               mean grayvalue (same as -mthres). 
                  
            - 'UT'   : upper threshold. Everything above is bone. If not 
              given the threshold is automatically estimated by the 
              histogram (same as -slevt). 
                  
           EXAMPLES:
           
           .. image:: images/test_mic_mid_5-YM.png 
              :scale: 200 %
           .. image:: images/test_gthres_1-YM.png
              :scale: 200 %


-gthres  : type

           Global threshold computed from global measures of the image. 
           'type' can be either "mean", "median" or "histo" (is the same as 
           'slevt' where "threshold" is the minimum grayvalue). 
         
           EXAMPLES :
         
           .. image:: images/test_mic_mid_5-YM.png 
              :scale: 200 %
           .. image:: images/test_gthres_1-YM.png
              :scale: 200 %


-thres  : thres1;thres2;...

          Classical multi thresholding of the image data by the given values 
          'thres1','thres2', ...
          The first threshold has to be bigger or equal to the lowest 
          gray value in the model.
          For example, for a 25;75 threshold all values from lowestValue-24 
          will be set to lowestValue, all from 25..74 will be set to 25 and 
          all from 75 to biggestValue will be set to 75. 
         
          Note: Only grayvalue>0 will produce an FEM element output i.e.:
         
           - if data from 0..255, thres 25;75 -> two elsets (25, 75) are 
             written because values from 0..24 are set to O!
           - if data from 5..255, thres 25;75 -> three elsets (5, 25, 75) 
             are written

          EXAMPLES : Original / 1 Threshold / 2 Thresholds 
         
          .. image:: images/test_mic_mid_1-XM.png 
             :scale: 300 %
          .. image:: images/test_thres_1-XM.png
             :scale: 300 %
          .. image:: images/test_thres_2-XM.png
             :scale: 300 %


-clean  : type 

          This option cleans the segmented (binary) image file which 
          contains '0' and values >'0'. Different variations are implemented: 

           - 'type' = 'FAST': Delete unconnected regions of 'bone' material 
             (gray value > 0). The background gray value (air) has to be '0'. 

           - 'type' = 'FAST1': Delete unconnected regions of 'air' material 
             (gray value = 0). 
          
           - 'type' =  'FAST2': Runs option 'FAST' (delete unconnected 'bone' 
             regions) followed by 'FAST1' (delete unconnected regions of 'air' 
             material with gray value = 0).         

           - 'type' = nodes;island: Similar to first possibility ('FAST') but 
             much slower. The first variable ('nodes') is the minimum number 
             of shared nodes between connected bone region of deletion algorithm.
             The second value ('island') is a minimum number of bone/marrow 
             voxels. If an island is detected with less voxels, it is deleted. 

           - 'type' = 'NUMPY': Delete unconnected regions of foreground voxels.
             Only two different materials are allowed. The bg has not to be 0.

           - 'type' = 'NUMPY1': Delete unconnected regions of background voxels.
             Only two different materials are allowed. The bg has not to be 0.


          EXAMPLES : Original / Clean ('FAST')
           
          .. image:: images/test_clean_0-XM.png 
             :scale: 150 %
          .. image:: images/test_clean_1-XM.png
             :scale: 150 %
              
          EXAMPLES : Original / Clean ('FAST1')
           
          .. image:: images/test_clean_3-XM.png 
             :scale: 150 %
          .. image:: images/test_clean_4-XM.png
             :scale: 150 %


-extend  : direction;thickVox;[newGrayvalue]

           Extend the model. Following parameters are required:
           
             - 'direction' is the direction of embedding. It can be -1 (x=0), 
               1 (x=nx), -2 (y=0), 2 (y=ny), -3 (z=0), 3 (z=nz).
           
             - 'thickVox' is the number of extended voxel layers (thickness). 
             - 'newGrayvalue' (optional) is a new grayvalue for the extended material. 
               If the thisparamater is given the extended region has not the 
               original grayvalue but this one. Only the voxels with a grayvalue 
               bigger than 0 are modified. 

           EXAMPLES : Original / Threshold + Extended
           
           .. image:: images/test_mic_mid_1-XM.png 
              :scale: 300 %
           .. image:: images/test_extend_1-XM.png
              :scale: 300 %


-close  :  direction;threshold;kernel;newGrayValue

           Close open models at the boundary. Parameters are:  
           
             - 'direction' gives the direction of closing. It can 
               be 1 (x), 2 (y), or 3 (z). 
             - 'threshold' and 'kernel' are parameters
               which are used to close the slice (see 'fill' option for details). 
             - 'newGrayValue' gives the gray value which should be used for 
               the new fill voxels. 
               
           The new model will not change in size, only the first and
           last layer will be replaced. 

           EXAMPLES :
         
           .. image:: images/test_close_0-XM.png 
              :scale: 150 %
           .. image:: images/test_close_1-XM.png
              :scale: 150 %


-embed  : direction;thickVoxIn;thickVoxOut;newGrayValue

          Embed model in a new material. Parameters are: 

            - 'direction' gives the direction of embedding (1,2,3).
              You can use + or - to indicate embed on top (+) or bottom (-) 
              only. If no sign is given, both directions will be embedded.
            - 'thickVoxIn' and 'thickVoxOut'  give the thickness of the
              embedding layer inside and outside of a bounding box (number
              of voxels).  Inside voxels with a value of "0"  and all 
              outside voxel will be set to 'newGrayValue'.

          EXAMPLES : 3- / 3+  / 2- / 2+
         
          .. image:: images/test_embed_1-XM.png 
             :scale: 130 %

          .. image:: images/test_embed_2-XM.png
             :scale: 130 %

          .. image:: images/test_embed_3-XM.png 
             :scale: 130 %

          .. image:: images/test_embed_4-XM.png
             :scale: 130 %
              

-cap   : direction;thickVox;newGrayValue  

         Make an end cap (bone material) on both ends in given direction
         of given thickness (number of voxels). The voxels will set to
         the given gray value. Parameters are: 

           - 'direction' gives the direction of embedding (1,2,3).
             You can use + or - to indicate embed on top (+) or bottom (-) 
             only. If no sign is given, both directions will be embedded.
           - 'thickVox' is the thickness of the cap
           - 'newGrayvalue' gray-value of the generated voxels.
                      
         EXAMPLES :

         .. image:: images/test_cap_1-XM.png 
            :scale: 200 %


-cover  : thickVox;newGrayValue

          Surround the whole model with a layer of material with given 
          thickness (number of voxels). Parameters are: 

            - 'thickVox' thickness of the cover voxels. 
            - 'newGrayValue' gray-value of the generated voxels.

          EXAMPLES :

          .. image:: images/test_cover_1-XM.png 
             :scale: 200 %


-block :  nVox1;nVox2;nVox3;newGrayValue

          Place the voxel model inside a bigger block of material. 
          Parameters are: 

            - 'nVox1', 'nVox2', 'nVox3' which are the number of voxels 
              in 1,2,3 direction of the new block 
            - 'newGrayValue' is the gray-value background of the new block

           EXAMPLES :

           .. image:: images/test_block_1-XM.png 
              :scale: 60 %


-cen   : threshold;filename;mode
         Compute the first voxel near center which has a value bigger or 
         equal than the given 'threshold'. It also computes the radius of a 
         sphere located in this point which captured the full voxel model. 
         This values are needed by the CGAL surface mesher.
         The output is in physical dimensions which is written to 'filename'
         using 'mode'. In detail: 
         
           - 'filename' If given the computed info is written into the file.
           - 'mode' of writing (python style) can be:
           
             - 'w' (create new file) write header + data) and 
             - 'a' (append to file, write data).


-bbox  : type;threshold;filename;mode;i;j;k

         Writes a bounding box to a given filename. Parameters are: 
         
           - 'type' can be "in" for inside or "out" for outside. 
           - 'threshold' is the given to indicate the boundary of the box.
           - 'filename' where the computed info is written. 
           - 'mode' of writing (python style) can be:
           
             - 'w' (create new file) write header + data) and 
             - 'a' (append to file, write data).            
           
           - 'i,j,k' are the starting points i,j,k of the search is 
              used only in the "in" bounding box. If it is not given
              the middle of the voxel model will be used. 
         
         Note: For the inner bounding box the model has to be enclosed by 
         matter.
         
         It writes following info to a file: ::
         
           #filename;threshold;bbType;minX;minY;minZ;maxX;maxY;maxZ
           femur-r8.mhd;100;in;0;0;0;106;161;102


-cog   : newGrayValue

         Compute the center of gravity. There will be an output on the 
         screen and if activagted in the information file ('ifo').
         The parameter 'newGrayValue' is a gray value which is be used 
         to plot the 3 orthogonal planes indicating the cog in the voxel model. 
         If 'newGrayValue' = 'off' the planes are not plotted in the
         voxel file. 


-ifo   : filename:mode

         Write general information to a file. This option can be very 
         helpful to generate a CSV file from your dataset. 
         
           - 'filename' where the computed info is written. 
           - 'mode' of writing (python style) can be:
           
             - 'w' (create new file) write header + data) and 
             - 'a' (append to file, write data).   

         It writes following CSV table to a file:

         .. image:: images/test_ifo.png 
            :scale: 150 %

-histo  : filename;normalize;nIntervals

          Output histogram data in a file. 

            - 'filename' where the computed info is written. 
            - 'normalize' option. Can only have the values 'y' or 'n'. 
              If 'y' the percentage of volume within the interval instead 
              of the number of voxel will be given.  
            - 'nIntervals' number of intervals, if not given, 
              the interval size will be "1". 

         It writes data to a file which can be plotted with Excel e.g.:

         .. image:: images/test_histo.png 


-shist  : type;[filename];[mode]

          Fit 1,2, or three materials which show a bell-shape distribution 
          (exp functions). Show the gray-value histogram as well as this 
          distribution. Parameters are: 

           - 'type' type of fit. 1,2, or 3 different functions (=materials) 
             can be selected by this parameter.  
           - 'filename' (optional) where the computed info is written. 
           - 'mode' of writing (python style) can be:
           
             - 'w' (create new file) write header + data) and 
             - 'a' (append to file, write data).  

         This option shows a histogram using matplotlib e.g. fit 1 or 2:
         
         .. image:: images/test_shist_1.png 

         .. image:: images/test_shist_2.png 

-help  :  Print usage


Info
~~~~ 

- File:   mic.py
- Author: D.H. Pahr

"""

from sys import argv, exit, stdout
import platform
import re
import struct
import time
import random
from string import split, replace, atof
import numpy
import numpy.fft
import shutil
import array
aarray = array
import dpFem
import dpTensor
import fec
import dpUtils
import dpTransform
import gzip
import xml.dom
import xml.dom.minidom
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy import misc
import os
osId = 0
if platform.uname()[4].find('64') != -1:
    osId = 64
else:
    osId = 32
import micf77
import miaf77
import dpMesherf77
import dpMesher

class mic():
    __version__ = 'V_18.12.2016'
    __author__ = 'D. H. Pahr'

    def __init__(self, modelName='default'):
        self.modelName = modelName

    def read(self, inFileName, ndim=None, imread=None, rawInfo=None, echo=True, recast=True):
        """
        Function selects which file reader should be called. Possible types
        file types are: bin, vol, rawiv, rawdp, aim, png, jpg, gif, raw
        
        @param  inFileName: input file name
          - TYPE: string
        @param  ndim: number of voxels in x, y, z. Needed for *.vol, *.bin, *.raw files.
          - TYPE: list[nx, ny, nz]
          - int nx,ny,nz ... number of voxels in x,y,z
        @param  imread: Image stack read parameters. Needed for *.png, *.gif, *.jpg files.
          - TYPE: list[firstImNo, step, lastImNo]
          - firstImNo ... number of the first image in the filename.
              E.g. if filename myFile0005.jpg, the firstImNo = 5
          - step      ... step in image numbers
          - lastImNo  ... number of last image
        @param rawInfo: Reading parameters for raw files.
          - TYPE: dict['binFormat':a, 'headerSize':b, 'bigEndian':c, 'fastestDir':d]
          - char a   ... binary format identifier e.g. 'B', 'f', ...
          - int b    ... Size of a preheader which should be overwritten.
          - bool c   ... True if byte order is big endian, otherwise False
          - string d ... Faster variing direction. Possibilieties: "x","z"
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - grayValue  ... gray value of voxel
          additdata: additional data information. Only in the case of an rawdp file.
            - TYPE: dict['-ndim': ndim, '-ldim':ldim, '-thres':thres]
            - list(int)   thres ... List of threshold values (0<...<255);
            - list(float32) ldim  ... x,y,z size of voxels (length = resolution)
            - list(int)   ndim  ... x,y,z number of voxels
        """
        additData = {}
        stdout.write(' ... read file %s\n' % inFileName)
        stdout.flush()
        if inFileName.find('.xml') > -1:
            voxelModel, additData = self.readXmlFile(inFileName)
        elif inFileName.find('.mhd') > -1 or inFileName.find('.mha') > -1:
            voxelModel, additData = self.readMhdFile(inFileName, echo=echo)
        elif inFileName.find('.nhdr') > -1:
            voxelModel, additData = self.readNhdrFile(inFileName, echo=echo)
        elif inFileName.find('.jpg') > -1 or inFileName.find('.png') > -1 or inFileName.find('.gif') > -1 or inFileName.find('.tif') > -1 or inFileName.find('.bmp') > -1:
            voxelModel = self.readImageFile(inFileName, imread, echo=echo)
        elif inFileName.find('.aim') > -1 or inFileName.find('.AIM') > -1:
            voxelModel, additData = self.readAimFile(inFileName, printLog=False)
        elif inFileName.find('.vol') > -1 or inFileName.find('.VOL') > -1:
            voxelModel = self.readGeneralFile1(inFileName, ndim, startByte=210)
        elif inFileName.find('.bin') > -1 or inFileName.find('.BIN') > -1:
            voxelModel = self.readGeneralFile1(inFileName, ndim, startByte=0)
        elif inFileName.find('.rawdp') > -1 or inFileName.find('.RAWDP') > -1:
            voxelModel, additData = self.readRawdpFile(inFileName)
        elif inFileName.find('.rawiv') > -1 or inFileName.find('.RAWIV') > -1:
            voxelModel = self.readRawivFile(inFileName, format='f')
        elif inFileName.upper().find('.ISQ') > -1:
            voxelModel, additData = self.readIsqFile(inFileName)
        elif inFileName.upper().find('.RAW') > -1 or inFileName.upper().find('.BST') > -1:
            if rawInfo == None:
                stdout.write("\n **ERROR**: read() '-raw' option not optional for reading RAW data! \n")
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            else:
                rawFormat = rawInfo['binFormat']
                headerSize = rawInfo['headerSize']
                bigEndian = rawInfo['bigEndian']
                fastestDir = rawInfo['fastestDir']
                if rawFormat == '1':
                    voxelModel = self.readBitFile(inFileName, ndim, fastestDir)
                elif headerSize == 0 and fastestDir == 'x':
                    voxelModel = self.readGeneralFile3(inFileName, ndim, rawFormat=rawFormat, bigEndian=bigEndian, echo=echo)
                else:
                    voxelModel = self.readGeneralFile2(inFileName, ndim, headerSize=headerSize, rawFormat=rawFormat, bigEndian=bigEndian, fastestDir=fastestDir)
        else:
            stdout.write('\n **ERROR**: \'-in\' intput file extension of file: "%s" not known!\n\n' % inFileName)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        additData['ElementType'] = voxelModel.dtype.name
        if 'Offset' not in additData:
            additData['Offset'] = [0.0, 0.0, 0.0]
        if 'TransformMatrix' not in additData:
            additData['TransformMatrix'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        if recast:
            if echo:
                stdout.write('     -> recast data         :     ')
            voxelModel = self.castType(voxelModel, 'f')
            if echo:
                stdout.write('done\n')
            stdout.flush()
        return (
         voxelModel, additData)

    def readAimFile(self, inFileName, printLog=False):
        """
        Scanco *.aim file reader. It can read *.aim files of format
        "B" (compressed) and "h" (standard). The voxel model is scaled automatically.
        
        @param  inFileName: input file name
          - TYPE: string
        @param  printLog: prints a read log to stdout.
          - TYPE: bool
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        """
        try:
            data = open(inFileName).read()
        except IOError:
            stdout.write("\n **ERROR** mic.readAimFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        start = 0
        stop = struct.calcsize('5i')
        sizePreHeader, sizeImStruct, sizeProcLog, sizeImData, sizeImAssocData = struct.unpack('5i', data[start:stop])
        start = stop
        stop = start + struct.calcsize('i')
        aimVer = struct.unpack('i', data[start:stop])[0]
        start = stop
        stop = start + struct.calcsize('4i')
        start = stop
        stop = start + struct.calcsize('i')
        aimTypeNum = struct.unpack('i', data[start:stop])[0]
        if aimTypeNum == 131074:
            aimType = 'short'
        elif aimTypeNum == 65537:
            aimType = 'char'
        elif aimTypeNum == 1376257:
            aimType = 'bin compressed'
        else:
            aimType = 'unknown'
        start = stop
        stop = start + struct.calcsize('3i')
        px, py, pz = struct.unpack('3i', data[start:stop])
        start = stop
        stop = start + struct.calcsize('3i')
        nx, ny, nz = struct.unpack('3i', data[start:stop])
        start = stop
        stop = start + struct.calcsize('3i')
        ox, oy, oz = struct.unpack('3i', data[start:stop])
        start = stop
        stop = start + struct.calcsize('12i')
        start = stop
        stop1 = start + struct.calcsize('B')
        stop2 = start + struct.calcsize('2B')
        stop3 = start + struct.calcsize('3B')
        stop4 = start + struct.calcsize('4B')
        stop = stop4
        bin_swap = struct.pack('cccc', data[stop2:stop3], data[stop3:stop4], data[start:stop1], data[stop1:stop2])
        lx = struct.unpack('f', bin_swap)[0] / 4.0
        start = stop
        stop1 = start + struct.calcsize('B')
        stop2 = start + struct.calcsize('2B')
        stop3 = start + struct.calcsize('3B')
        stop4 = start + struct.calcsize('4B')
        stop = stop4
        bin_swap = struct.pack('cccc', data[stop2:stop3], data[stop3:stop4], data[start:stop1], data[stop1:stop2])
        ly = struct.unpack('f', bin_swap)[0] / 4.0
        start = stop
        stop1 = start + struct.calcsize('B')
        stop2 = start + struct.calcsize('2B')
        stop3 = start + struct.calcsize('3B')
        stop4 = start + struct.calcsize('4B')
        stop = stop4
        bin_swap = struct.pack('cccc', data[stop2:stop3], data[stop3:stop4], data[start:stop1], data[stop1:stop2])
        lz = struct.unpack('f', bin_swap)[0] / 4.0
        start = stop
        stop = start + struct.calcsize('4i')
        start = sizePreHeader + sizeImStruct
        stop = start + sizeProcLog
        format = repr(sizeProcLog) + 's'
        ProcLog = struct.unpack(format, data[start:stop])[0]
        time1 = time.clock()
        stdout.write(' ... read & process image data\n')
        stdout.flush()
        start = sizePreHeader + sizeImStruct + sizeProcLog
        size = sizeImData
        if aimType == 'char':
            form2 = 'B'
            form = 'i1'
        elif aimType == 'short':
            form2 = 'h'
            form = 'i2'
        else:
            stdout.write(' **ERROR** readAimFile() AIM format unknown!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        del data
        voxelModel = numpy.core.records.fromfile(inFileName, formats=form, shape=nx * ny * nz, offset=start)['f0']
        time2 = time.clock()
        stdout.write('     -> read finished in    :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        voxelModel = voxelModel.reshape((nz, ny, nx))
        time3 = time.clock()
        stdout.write('     -> reshape finished in :   %8.1f sec\n' % (time3 - time2))
        stdout.flush()
        if printLog == True:
            print '========================================================================'
            print '\nL O G F I L E\n'
            print '========================================================================'
            print 'Size of pre header           [bytes] :', sizePreHeader
            print 'Size of image structure      [bytes] :', sizeImStruct
            print 'Size of processing log       [bytes] :', sizeProcLog
            print 'Size of image data           [bytes] :', sizeImData
            print 'Size of image associate data [bytes] :', sizeImAssocData
            print '\nMODEL:'
            print '========================================================================'
            print 'Size in x direction          [pixel] :', nx
            print 'Size in y direction          [pixel] :', ny
            print 'Size in z direction          [pixel] :', nz
            print 'Voxel length in x direction     [mm] :', lx
            print 'Voxel length in y direction     [mm] :', ly
            print 'Voxel length in z direction     [mm] :', lz
            print 'Position in x direction      [pixel] :', px
            print 'Position in y direction      [pixel] :', py
            print 'Position in z direction      [pixel] :', pz
            print 'Offset in x direction        [pixel] :', ox
            print 'Offset in y direction        [pixel] :', oy
            print 'Offset in z direction        [pixel] :', oz
            print '\nPROCESSING LOG:'
            print '========================================================================'
            print ProcLog
            print 'IMAGE DATA:'
            print '========================================================================'
            print 'Shape of imported array              :', voxelModel.shape
        ndim = [
         int(nx), int(ny), int(nz)]
        ldim = [float(lx), float(ly), float(lz)]
        additData = {}
        additData['-ldim'] = ldim
        additData['-ndim'] = ndim
        additData['ElementSpacing'] = ldim
        additData['DimSize'] = ndim
        return (
         voxelModel, additData)

    def readAimFile_OLD(self, inFileName, printLog=False):
        """
        Scanco *.aim file reader. It can read *.aim files of format
        "B" (compressed) and "h" (standard). The voxel model is scaled automatically.
        
        @param  inFileName: input file name
          - TYPE: string
        @param  printLog: prints a read log to stdout.
          - TYPE: bool
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        """
        try:
            data = open(inFileName).read()
        except IOError:
            stdout.write("\n **ERROR** mic.readAimFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        start, stop = 0, struct.calcsize('5i')
        sizePreHeader, sizeImStruct, sizeProcLog, sizeImData, sizeImAssocData = struct.unpack('5i', data[start:stop])
        start = 56
        size = struct.calcsize('3i')
        stop = start + size
        nx, ny, nz = struct.unpack('3i', data[start:stop])
        start = 128
        size = struct.calcsize('3f')
        stop = start + size
        start = sizePreHeader + sizeImStruct
        size = sizeProcLog
        stop = start + size
        format = repr(size) + 's'
        text = struct.unpack(format, data[start:stop])
        ProcLog = text[0]
        noOfScannedPixel = nx * ny * nz
        noOfStoredPixel = sizeImData / 2
        time1 = time.clock()
        stdout.write(' ... read & process image data\n')
        stdout.flush()
        start = sizePreHeader + sizeImStruct + sizeProcLog
        size = sizeImData
        stop = start + size
        if nx * ny * nz == sizeImData:
            form = 'i1'
        elif 2 * nx * ny * nz == sizeImData:
            form = 'i2'
        else:
            stdout.write(' **ERROR** readAimFile() Number of scanned pixels and voxel model size is not equal!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        del data
        voxelModel = numpy.core.records.fromfile(inFileName, formats=form, shape=nx * ny * nz, offset=start)['f0']
        time2 = time.clock()
        stdout.write('     -> read finished in    :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        voxelModel = voxelModel.reshape((nz, ny, nx))
        time3 = time.clock()
        stdout.write('     -> reshape finished in :   %8.1f sec\n' % (time3 - time2))
        stdout.flush()
        if printLog == True:
            print '========================================================================'
            print '\nL O G F I L E\n'
            print '========================================================================'
            print 'Size of pre header           [bytes] :', sizePreHeader
            print 'Size of image structure      [bytes] :', sizeImStruct
            print 'Size of processing log       [bytes] :', sizeProcLog
            print 'Size of image data           [bytes] :', sizeImData
            print 'Size of image associate data [bytes] :', sizeImAssocData
            print '\nMODEL:'
            print '========================================================================'
            print 'Size in x direction          [pixel] :', nx
            print 'Size in y direction          [pixel] :', ny
            print 'Size in z direction          [pixel] :', nz
            print '\nPROCESSING LOG:'
            print '========================================================================'
            print ProcLog
            print 'IMAGE DATA:'
            print '========================================================================'
            print 'Shape of imported array              :', voxelModel.shape
        return voxelModel

    def readRawivFile(self, inFileName, format='f'):
        """
        LBIE - mesher file reader.
        
        @param  inFileName: input file name
          - TYPE: string
        @param  format: binary format secifier
          - TYPE: char
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        """
        try:
            data = open(inFileName).read()
        except IOError:
            stdout.write("\n **ERROR** mic.readRawivFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        fsize3 = struct.calcsize('>3f')
        Isize3 = struct.calcsize('>3I')
        Isize = struct.calcsize('>I')
        stdout.write(' ... read & process image data\n')
        stdout.flush()
        start = 0
        stop = start + fsize3
        minX, minY, minZ = struct.unpack('>3f', data[start:stop])
        start = stop
        stop = start + fsize3
        maxX, maxY, maxZ = struct.unpack('>3f', data[start:stop])
        start = stop
        stop = start + Isize
        numVerts = struct.unpack('>I', data[start:stop])
        start = stop
        stop = start + Isize
        numCells = struct.unpack('>I', data[start:stop])
        start = stop
        stop = start + Isize3
        dimX, dimY, dimZ = struct.unpack('>3I', data[start:stop])
        start = stop
        stop = start + fsize3
        originX, originY, originZ = struct.unpack('>3f', data[start:stop])
        start = stop
        stop = start + fsize3
        spanX, spanY, spanZ = struct.unpack('>3f', data[start:stop])
        voxelModel = self.createFlatVoxelModel(dimX * dimY * dimZ, format)
        count = 0
        for i in xrange(dimX * dimY * dimZ):
            count += 1
            if count % 10000 == 0 or count == 0:
                progress = float(count) / float(dimX * dimY * dimZ) * 100.0
                stdout.write('     -> Processed Data      : %10i %s\r' % (progress, '%'))
                stdout.flush()
            start = stop
            stop = start + struct.calcsize('>f')
            voxelModel[i] = abs(struct.unpack('>f', data[start:stop])[0])

        voxelModel = numpy.reshape(voxelModel, (dimZ, dimY, dimX))
        stdout.write('     -> Processed Data      : %10i      \n' % count)
        stdout.flush()
        return voxelModel

    def readRawdpFile(self, inFileName):
        """
        Raw file reader with special header (from Dieter Pahr).
        A typical header lookds like ::
          voxelnumber   [60, 70, 85]
          voxellength   [15., 15., 35.]
          threshold     [75]
          format        B
          endheader
        
        @param  inFileName: input file name
          - TYPE: string
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - Int grayValue  ... gray value of voxel,
          additdata: additional data information. Only in the case of an rawdp file.
            - TYPE: dict['-ndim': ndim, '-ldim':ldim, '-thres':thres]
            - list(int)   thres ... List of threshold values (0<...<255);
            - list(float32) ldim  ... x,y,z size of voxels (length = resolution)
            - list(int)   ndim  ... x,y,z number of voxels
        """
        try:
            data = open(inFileName, 'r')
        except IOError:
            stdout.write("\n **ERROR** mic.readRawdpFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        stdout.write(' ... read rawdp data\n')
        stdout.flush()
        line = ''
        ndim = []
        ldim = []
        format = ''
        threshold = []
        additData = {}
        while line.find('end') < 0:
            line = data.readline()
            line = line.replace('[', ' ')
            line = line.replace(']', ' ')
            line = line.replace(',', ' ')
            if line.find('voxelnumber') > -1:
                try:
                    key, nx, ny, nz = line.split()
                    ndim = [int(nx), int(ny), int(nz)]
                    additData['-ndim'] = ndim
                    additData['DimSize'] = ndim
                except ValueError:
                    stdout.write("\n **ERROR**: readRawdpFile() 'voxelnumber' has wrong size!\n\n")
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

            elif line.find('voxellength') > -1:
                try:
                    key, lx, ly, lz = line.split()
                    ldim = [float(lx), float(ly), float(lz)]
                    additData['-ldim'] = ldim
                    additData['ElementSpacing'] = ldim
                except ValueError:
                    stdout.write("\n **ERROR**: readRawdpFile() 'voxellength' has wrong size!\n\n")
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

            elif line.find('threshold') > -1:
                alist = []
                alist = line.split()
                if alist[1] != 'None':
                    for thresId in range(1, len(alist)):
                        threshold.append(int(alist[thresId]))

                    additData['-thres'] = threshold
            elif line.find('format') > -1:
                try:
                    key, format = line.split()
                except ValueError:
                    stdout.write("\n **ERROR**: readRawdpFile() 'format' has wrong size!\n\n")
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

            elif line.find('endheader') > -1:
                continue

        time1 = time.clock()
        dimX = ndim[0]
        dimY = ndim[1]
        dimZ = ndim[2]
        data = data.read()
        dataLen = len(data)
        start = 0
        if format == 'B':
            mul = 1
        elif format == 'h':
            mul = 2
        elif format == 'H':
            mul = 2
        elif format == 'i':
            mul = 4
        elif format == 'f':
            mul = 4
        else:
            stdout.write('\n **ERROR** readRawdp(). format=%s! not implemented\n' % format)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        stop = dimX * dimY * dimZ * mul
        form = repr(dimX * dimY * dimZ) + format
        voxelModel = struct.unpack(form, data[start:stop])
        voxelModel = numpy.reshape(voxelModel, (dimZ, dimY, dimX))
        time2 = time.clock()
        stdout.write('     -> read finished in    :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        return (
         voxelModel, additData)

    def readGeneralFile1(self, inFileName, ndim, startByte=0):
        """
        Raw file reader for 1 byte data (format='B'). Faster then self.readGeneralFile2().
        
        @param  inFileName: input file name
          - TYPE: string
        @param  ndim: number of voxels in x, y, z. Needed for *.vol, *.bin, *.raw files.
          - TYPE: list[nx, ny, nz]
          - int nx,ny,nz ... number of voxels in x,y,z
        @param startByte: Location where data section starts. The reading
          routine jumps over the header.
            - TYPE: int
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        """
        if ndim == None:
            stdout.write("\n **ERROR** mic.readGeneralFile1(): '-ndim' ... Dimension of voxel array not given!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        format = 'B'
        stdout.write(' ... read & process image data 1\n')
        stdout.flush()
        time1 = time.clock()
        nx = ndim[0]
        ny = ndim[1]
        nz = ndim[2]
        size = nx * ny * nz + startByte
        time1 = time.clock()
        try:
            file = open(inFileName, mode='rb')
        except IOError:
            stdout.write("\n **ERROR**: mic.readGeneralFile1(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        binvalues = aarray.array(format)
        binvalues.read(file, size)
        time2 = time.clock()
        stdout.write('     -> read    data in     : %8.2f sec\n' % (time2 - time1))
        stdout.flush()
        try:
            voxelModel = numpy.array(binvalues[startByte:size], numpy.uint8)
        except MemoryError:
            stdout.write(' **ERROR** read: Cannot allocate enough memory!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        time3 = time.clock()
        stdout.write('     -> convert data in     : %8.2f sec\n' % (time3 - time2))
        stdout.flush()
        voxelModel = numpy.reshape(voxelModel, (nz, ny, nx))
        time4 = time.clock()
        file.close()
        stdout.write('     -> reshape data in     : %8.2f sec\n' % (time4 - time3))
        stdout.flush()
        voxelModel = self.castType(voxelModel, 'B')
        return voxelModel

    def readGeneralFile2(self, inFileName, ndim, headerSize=0, rawFormat='B', bigEndian=False, fastestDir='x'):
        """
        General raw file reader for other data formats.
        
        @param  inFileName: input file name
          - TYPE: string
        @param  ndim: number of voxels in x, y, z. Needed for *.vol, *.bin, *.raw files.
          - TYPE: list[nx, ny, nz]
          - int nx,ny,nz ... number of voxels in x,y,z
        @param headerSize: Location where data section starts. The reading
          routine jumps over the header.
            - TYPE: int
        @param rawFormat: Binary format specifier.
          - TYPE: char
        @param bigEndian: Flag if data are in big endian format.
          - TYPE: bool
        @param fastestDir: Fastest variing direction. Possible values are 'x' or 'z'.
          - TYPE: char
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - Int grayValue  ... gray value of voxel,
        
        """
        if ndim == None:
            stdout.write("\n **ERROR** mic.readGeneralFile2(): '-ndim' ... Dimension of voxel array not given!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        stdout.write(' ... read & process image data 2\n')
        stdout.flush()
        time1 = time.clock()
        nx = ndim[0]
        ny = ndim[1]
        nz = ndim[2]
        time1 = time.clock()
        try:
            data = open(inFileName, mode='rb').read()
            dataLen = len(data)
        except IOError:
            stdout.write("\n **ERROR** mic.readGeneralFile2(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        start = headerSize
        size = struct.calcsize(repr(nx * ny * nz) + rawFormat)
        stop = start + size
        if bigEndian == True:
            pre = '>'
        elif bigEndian == False:
            pre = ''
        else:
            stdout.write('\n **ERROR** mic.readGeneralFile2(): bigEndian have to be True/False not "%s"\n\n!' % bigEndian)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dataLen != stop:
            stdout.write("\n **WARNING readGeneralFile2(): Unmatch sizes - check '-raw' option!")
            stdout.flush()
            stdout.write('\n           Header size                       %10i  ' % headerSize)
            stdout.flush()
            stdout.write('\n           Number of voxels nx*ny*nz         %10i  ' % (nx * ny * nz))
            stdout.flush()
            stdout.write("\n           Number of bytes of format '%s'     %10i  " % (rawFormat, struct.calcsize(rawFormat)))
            stdout.flush()
            stdout.write('\n           --------------------------------------------')
            stdout.flush()
            stdout.write('\n           => Bytes to be read               %10i  ' % stop)
            stdout.flush()
            stdout.write('\n           => Bytes in datafile              %10i  \n\n' % dataLen)
            stdout.flush()
        if fastestDir == 'x':
            voxelModel = self.createFlatVoxelModel(nz * ny * nx, rawFormat)
            voxelModel = struct.unpack(pre + repr(nx * ny * nz) + rawFormat, data[start:stop])
            voxelModel = numpy.reshape(voxelModel, (nz, ny, nx))
        elif fastestDir == 'z':
            auxModel = self.createFlatVoxelModel(nz * ny * nx, rawFormat)
            auxModel = struct.unpack(pre + repr(nx * ny * nz) + rawFormat, data[start:stop])
            auxModel = numpy.reshape(auxModel, (nx, ny, nz))
            voxelModel = self.createVoxelModel(nx, ny, nz, rawFormat)
            stdout.write('     ->  start renumbering \n')
            stdout.flush()
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        voxelModel[k, j, i] = auxModel[i, j, k]

            del auxModel
        else:
            stdout.write('\n **ERROR** readGeneralFile2(): Fastest direction "%s" not known!\n\n' % fastestDir)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        time2 = time.clock()
        stdout.write('     -> read finished in    :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        return voxelModel

    def readGeneralFile3(self, inFileName, ndim, rawFormat='B', bigEndian=False, echo=True):
        """
        General Raw file reader based on  ideas (see www). Faster then
        self.readGeneralFile2().
        
        @param  inFileName: input file name
          - TYPE: string
        @param  ndim: number of voxels in x, y, z. Needed for *.vol, *.bin, *.raw files.
          - TYPE: list[nx, ny, nz]
          - int nx,ny,nz ... number of voxels in x,y,z
        @param rawFormat: Binary format specifier.
          - TYPE: char
        @param bigEndian: Flag if data are in big endian format.
          - TYPE: bool
          
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - Int grayValue  ... gray value of voxel,
        """
        if echo:
            stdout.write(' ... read & process image data 3 \n')
            stdout.flush()
        time1 = time.clock()
        if ndim == None:
            stdout.write("\n **ERROR** mic.readGeneralFile3(): '-ndim' ... Dimension of voxel array not given!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx = ndim[0]
        ny = ndim[1]
        nz = ndim[2]
        size = nx * ny * nz
        mul = 0
        try:
            file = open(inFileName, 'rb')
        except IOError:
            stdout.write("\n **ERROR**: mic.readGeneralFile3(): intput file '%s' not found!\n\n" % inFileName)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        if bigEndian == False:
            swap = 0
        else:
            swap = 1
        dtype, mul = self.getDtypeMul(rawFormat)
        if swap == 1:
            voxelModel = numpy.fromfile(inFileName, dtype=dtype).byteswap()
        else:
            voxelModel = numpy.fromfile(inFileName, dtype=dtype)
        time2 = time.clock()
        if echo:
            stdout.write('     -> read    data in     : %8.2f sec\n' % (time2 - time1))
            stdout.flush()
        time3 = time.clock()
        voxelModel = voxelModel.reshape((nz, ny, nx))
        if echo:
            stdout.write('     -> reshape data in     : %8.2f sec\n' % (time3 - time2))
            stdout.flush()
        file.close()
        return voxelModel

    def readIsqFile(self, inFileName, echo=True, info=False):
        """
        Read Isq File (Scanco)
        """
        if echo:
            stdout.write(' ... read isq file \n')
            stdout.flush()
        time1 = time.clock()
        try:
            f = open(inFileName, 'rb')
        except IOError:
            stdout.write("\n **ERROR**: mic.readIsqFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        f.seek(44)
        header = numpy.zeros(6)
        for i in range(0, 6):
            header[i] = struct.unpack('i', f.read(4))[0]

        elsp = [header[3] / header[0] / 1000,
         header[4] / header[1] / 1000,
         header[5] / header[2] / 1000]
        f.seek(508)
        headersize = 512 * (1 + struct.unpack('i', f.read(4))[0])
        f.seek(headersize)
        if info == True:
            voxelModel = None
        else:
            voxelModel = numpy.fromfile(f, dtype='i2')
        f.close()
        ndim = [
         int(header[0]), int(header[1]), int(header[2])]
        ldim = [float(elsp[0]), float(elsp[1]), float(elsp[2])]
        additData = {}
        additData['-ldim'] = ldim
        additData['-ndim'] = ndim
        additData['ElementSpacing'] = ldim
        additData['DimSize'] = ndim
        additData['HeaderSize'] = headersize
        additData['TransformMatrix'] = [
         1, 0, 0, 0, 1, 0, 0, 0, 1]
        additData['CenterOfRotation'] = [0.0, 0.0, 0.0]
        additData['Offset'] = [0.0, 0.0, 0.0]
        additData['AnatomicalOrientation'] = 'LPS'
        additData['ElementType'] = 'int16'
        additData['ElementDataFile'] = inFileName
        time2 = time.clock()
        if echo:
            stdout.write('     -> read    data in     : %8.2f sec\n' % (time2 - time1))
            stdout.flush()
        time3 = time.clock()
        if info == False:
            voxelModel = voxelModel.reshape((ndim[2], ndim[1], ndim[0]))
            if echo:
                stdout.write('     -> reshape data in     : %8.2f sec\n' % (time3 - time2))
                stdout.flush()
        return (
         voxelModel, additData)

    def readBitFile(self, infile, ndim, fastestDir):
        """
        Read a bit file 
        """
        stdout.write(' ... read bit file (headerSize and Endian ignored!) \n')
        stdout.flush()
        io = open(infile, 'rb')
        rawData = io.read()
        img = Image.fromstring('1', (ndim[2] * ndim[1] * ndim[0], 1), rawData, 'raw')
        dataList = list(img.getdata())
        voxelModel = numpy.array(dataList)
        if fastestDir == 'x':
            voxelModel = voxelModel.reshape(ndim)
        elif fastestDir == 'z':
            voxelModel = voxelModel.reshape((ndim[2], ndim[1], ndim[0]))
        else:
            stdout.write('\n **ERROR** mic.readBitFile(): direction=%s not implemented!' % fastestDir)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return voxelModel

    def readImageFile(self, inFileName, imread, echo=True):
        """
        Image file writer. Function writes *.png, *.jpg, *.gif files.
        The file type is choosen from the file extension.
        Upper left corner of first image (0,0) pixel is the (0,0,0) voxel at
        the "SWB" (south west bottom) voxel of an RVE. The image is located
        such that a viewer looks from outside to the bottom RVE surface.
        
        Image stack file reader for png, jpg, gif, ... files.
        The file type is choosen from the file extension.
        Upper left corner of first image (0,0 pixel) is the (0,0,0) voxel at
        the "SWB" (south west bottom) voxel of an RVE. The image is located
        such that a viewer looks from outside to the bottom RVE surface.
        
        The last pixel (lower,right) the "NWB" point::
        
                          N
                    nwt         net
                     o-----------o
                    /|     T    /|
           W       / o-nwb- - -/-o nwb    E
              swt o-----------o set
                  |/      B   |/
                  o-----------o
                 swb         seb
        
        
        @param  inFileName: input file name with wildcards '#'
          - TYPE: string
        @param  imread: image read paramters. Parameters are need
          to create the file names for the image stack
            - TYPE: list[0] = firstImNo, list[1] = step, list[2] = lastImNo
            - int firstImNo ... number of first image >=0
            - int step      ... step in image numbers >0
            - int lastImNo  ... number of last image
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        """
        time1 = time.clock()
        fileList = []
        if imread == None and inFileName.find('#') > -1:
            if echo:
                stdout.write(' ... read & process image files automatically \n')
                stdout.flush()
            pathFile = os.path.realpath(inFileName)
            path, fileName = os.path.split(pathFile)
            allFiles = os.listdir(path)
            preFile, postExt = fileName.replace('#', ' ').split()
            for curFile in allFiles:
                if curFile.find(preFile) > -1 and curFile.find(postExt) > -1:
                    fileList.append(path + '/' + curFile)

            fileList.sort()
            if len(fileList) == 0:
                stdout.write('\n **ERROR** mic.readImageFile(): No files are found for auto read!\n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        elif imread == None and inFileName.find('#') == -1:
            if echo:
                stdout.write(' ... read & process image 1 file\n')
                stdout.flush()
            image = misc.imread(inFileName)
            ny, nx = image.shape
            nz = 1
            voxelModel = self.createVoxelModel(nx, ny, nz, 'B')
            voxelModel[0, :, :] = image
        elif imread != None:
            if echo:
                stdout.write(' ... read & process image files manually \n')
                stdout.flush()
            infile, ext = self.getFilenameAndExtension(inFileName)
            firstImNo = imread[0]
            step = imread[1]
            lastImNo = imread[2]
            length = infile.count('#')
            infile = infile.replace('#', '')
            sno = repr(firstImNo)
            lno = len(sno)
            dno = length - lno
            if dno < 0:
                stdout.write('\n **ERROR** mic.readImageFile(): At least %i Digits needed!\n\n' % lno)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            zeros = ''
            for l in range(dno):
                zeros += '0'

            for m in range(firstImNo, lastImNo + 1, step):
                sno = repr(m)
                lno = len(sno)
                dno = length - lno
                if dno < 0:
                    stdout.write('\n **ERROR** readImageFile(): At least %i Digits needed!\n\n' % lno)
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                zeros = ''
                for l in range(dno):
                    zeros += '0'

                curInFileName = infile + zeros + sno + '.' + ext
                fileList.append(curInFileName)

        else:
            stdout.write("\n **ERROR** mic.readImageFile(): Wild Card has to be '#'! \n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if len(fileList) > 0:
            curInFileName = fileList[0]
            image = misc.imread(curInFileName)
            ny, nx = image.shape
            nz = len(fileList)
            voxelModel = self.createVoxelModel(nx, ny, nz, 'B')
            k = 0
            sum = len(fileList)
            notRead = []
            if echo:
                dpUtils.progressStart('     -> Processed Images    : ')
            for curInFileName in fileList:
                image = misc.imread(curInFileName)
                nry, nrx = image.shape
                if nry == ny and nrx == nx:
                    voxelModel[k, :, :] = image
                    k += 1
                else:
                    notRead.append(curInFileName)
                progress = k / float(sum) * 10
                if echo:
                    dpUtils.progressNext(progress)

            if echo:
                dpUtils.progressEnd()
            if echo:
                if len(notRead) > 0:
                    stdout.write('     ->  **WARNING** following file not read (different shape): \n')
                    stdout.flush()
                    for fileName in notRead:
                        stdout.write('        %s \n' % fileName)
                        stdout.flush()

        time2 = time.clock()
        if echo:
            stdout.write('     -> read    data in     : %8.2f sec\n' % (time2 - time1))
            stdout.flush()
        return voxelModel

    def readXmlFile(self, inFileName, info=False):
        """
        General file reader. All important information is in an seperate xml 
        header file.
        
        @param  inFileName: input file name
          - TYPE: string
        @param  info:  if 'True' only additdata are return and voxelModel is set
                       to 'None'
          
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z stareadXmlrt a 0, x fastest.
            - Int grayValue  ... gray value of voxel
        additdata: additional data information. Only in the case of an rawdp file.
          - TYPE: dict['-ndim': ndim, '-ldim':ldim, '-thres':thres]
          - list(int)   thres ... List of threshold values (0<...<255);
          - list(float32) ldim  ... x,y,z size of voxels (length = resolution)
          - list(int)   ndim  ... x,y,z number of voxels
        """
        stdout.write(' ... read & process xml data \n')
        stdout.flush()
        time1 = time.clock()
        try:
            file = open(inFileName)
            file.close()
        except IOError:
            stdout.write("\n **ERROR** mic.readXmlFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        xmlDoc = xml.dom.minidom.parse(inFileName)
        volume = xmlDoc.getElementsByTagName('volume')
        dimension = volume[0].getElementsByTagName('dimensions')[0]
        xnode = dimension.getElementsByTagName('x')[0]
        xchild = xnode.firstChild
        nx = int(xchild.nodeValue)
        ny = int(dimension.getElementsByTagName('y')[0].firstChild.nodeValue)
        nz = int(dimension.getElementsByTagName('z')[0].firstChild.nodeValue)
        ndim = [nx, ny, nz]
        spacing = volume[0].getElementsByTagName('spacing')[0]
        xs = float(spacing.getElementsByTagName('x')[0].firstChild.nodeValue)
        ys = float(spacing.getElementsByTagName('y')[0].firstChild.nodeValue)
        zs = float(spacing.getElementsByTagName('z')[0].firstChild.nodeValue)
        ldim = [xs, ys, zs]
        file_info = volume[0].getElementsByTagName('volumeinonefile_info')[0]
        path = file_info.getElementsByTagName('rawfilepath')[0].firstChild.nodeValue
        xmlpathRawFile = os.path.realpath(inFileName)
        xmlpath, RawFile = os.path.split(xmlpathRawFile)
        path = xmlpath
        name = file_info.getElementsByTagName('rawfilename')[0].firstChild.nodeValue
        filename = path + '/' + name
        if not os.path.isfile(filename):
            stdout.write('\n **ERROR**: mic.readXmlFile(): %s is could not be found!\n\n' % filename)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        file_info = volume[0].getElementsByTagName('volumeinslicefiles_info')[0]
        if file_info.getElementsByTagName('numberofslicefiles')[0].firstChild != None:
            stdout.write('\n **ERROR**: mic.readXmlFile(): Multi slice file not implemented!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        bType = volume[0].getElementsByTagName('scalartype')[0].firstChild.nodeValue
        bits = int(volume[0].getElementsByTagName('bitsused')[0].firstChild.nodeValue)
        header = int(volume[0].getElementsByTagName('headersize')[0].firstChild.nodeValue)
        bigendian = volume[0].getElementsByTagName('bigendian')[0].firstChild.nodeValue
        rawInfo = {}
        if bType == 'uchar' and bits == 8:
            rawInfo['binFormat'] = 'B'
        elif bType == 'short' and bits == 16:
            rawInfo['binFormat'] = 'h'
        elif bType == 'integer' and bits == 32:
            rawInfo['binFormat'] = 'i'
        elif bType == 'float' and bits == 32:
            rawInfo['binFormat'] = 'f'
        else:
            stdout.write('\n **ERROR**: mic.readXmlFile(): scalartype = %s with %s bitsused not implemented!\n\n' % (bType, repr(bits)))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        rawInfo['headerSize'] = header
        if bigendian == 'no':
            rawInfo['bigEndian'] = False
            etype = 'little'
        else:
            rawInfo['bigEndian'] = True
            etype = 'big'
        rawInfo['fastestDir'] = 'x'
        if info == True:
            voxelModel = None
        else:
            voxelModel = self.readGeneralFile3(filename, ndim, rawFormat=rawInfo['binFormat'], bigEndian=rawInfo['bigEndian'])
        additData = {}
        additData['-ldim'] = ldim
        additData['-ndim'] = ndim
        additData['ElementSpacing'] = ldim
        additData['DimSize'] = ndim
        return (
         voxelModel, additData)

    def readMhdFile(self, inFileName, info=False, echo=True):
        """
        Mha header file reader. It will call a raw reader to read the data. 
        
        @param  inFileName: input file name
          - TYPE: string
        @param  info:  if 'True' only additdata are return and voxelModel is set
                       to 'None'
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX]        int ret = QMessageBox::warning(this, tr("Item Cut"),
                     tr("Do you really like to cut this item?"),
                     QMessageBox::Yes | QMessageBox::No);
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - Int grayValue  ... gray value of voxel
        additdata: additional data information. Only in the case of an rawdp file.
          - TYPE: dict['-ndim': ndim, '-ldim':ldim, '-thres':thres]
          - list(int)   thres ... List of threshold values (0<...<255);
          - list(float32) ldim  ... x,y,z size of voxels (length = resolution)
          - list(int)   ndim  ... x,y,z number of voxels
        """
        if echo:
            stdout.write(' ... read & process mhd data \n')
            stdout.flush()
        time1 = time.clock()
        try:
            file = open(inFileName)
            file.close()
        except IOError:
            stdout.write("\n **ERROR** mic.readMhdFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        osIn = open(inFileName, 'r')
        lines = osIn.readlines()
        osIn.close()
        parDict = {}
        for line in lines:
            if len(line) > 1:
                line = line.replace('\n', '')
                key, data = line.split('=')
                key = key.replace(' ', '')
                parDict[key] = data

        if 'CompressedData' not in parDict:
            parDict['CompressedData'] = 'False'
        if 'BinaryData' in parDict and 'CompressedData' in parDict:
            errKey = None
            corrVal = None
            if parDict['BinaryData'].find('True') == -1:
                errKey = 'BinaryData'
                corVal = 'True'
            elif parDict['CompressedData'].find('False') == -1:
                errKey = 'CompressedData'
                corVal = 'False'
            if errKey != None:
                stdout.write("\n **ERROR** mic.readMhdFile(): Key '%s' must be '%s'!\n\n" % (errKey, corVal))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        nx, ny, nz = parDict['DimSize'].split()
        ndim = [int(nx), int(ny), int(nz)]
        xs = 0.0
        ys = 0.0
        zs = 0.0
        if 'ElementSpacing' in parDict:
            xs, ys, zs = parDict['ElementSpacing'].split()
        elif 'ElementSize' in parDict:
            xs, ys, zs = parDict['ElementSize'].split()
        else:
            stdout.write('\n **ERROR**: mic.readMhdFile(): ElementSpacing or ElementSize not found!\n\n' % filename)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        ldim = [
         float(xs), float(ys), float(zs)]
        filename = parDict['ElementDataFile']
        filename = filename.replace(' ', '')
        parDict['ElementDataFile'] = filename
        path, RawMhaFile = os.path.split(filename)
        if path == '':
            mhdPathRawFile = os.path.realpath(inFileName)
            path, RawMhaFile = os.path.split(mhdPathRawFile)
            filename = path + '/' + filename
        if not os.path.isfile(filename):
            stdout.write('\n **ERROR**: mic.readMhdFile(): %s is could not be found!\n\n' % filename)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        rawInfo = {}
        format = parDict['ElementType']
        format = format.replace(' ', '')
        if format == 'MET_UCHAR':
            rawInfo['binFormat'] = 'B'
        elif format == 'MET_SHORT':
            rawInfo['binFormat'] = 'h'
        elif format == 'MET_USHORT':
            rawInfo['binFormat'] = 'H'
        elif format == 'MET_INT':
            rawInfo['binFormat'] = 'i'
        elif format == 'MET_FLOAT':
            rawInfo['binFormat'] = 'f'
        else:
            stdout.write("\n **ERROR** mic.readMhd(): 'ElementType = %s' not supported!\n\n" % format)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if 'BinaryDataByteOrderMSB' in parDict:
            byteorder = parDict['BinaryDataByteOrderMSB']
        elif 'ElementByteOrderMSB' in parDict:
            byteorder = parDict['ElementByteOrderMSB']
        else:
            byteorder = 'False'
        rawInfo['headerSize'] = 0
        byteorder = byteorder.replace(' ', '')
        if byteorder == 'False':
            rawInfo['bigEndian'] = False
        elif byteorder == 'True':
            rawInfo['bigEndian'] = True
        else:
            stdout.write("\n **ERROR** mic.readMhd(): '...ByteOrderMSB = %s' not supported!\n\n" % byteorder)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        rawInfo['fastestDir'] = 'x'
        if info == True:
            voxelModel = None
        else:
            voxelModel = self.readGeneralFile3(filename, ndim, rawFormat=rawInfo['binFormat'], bigEndian=rawInfo['bigEndian'], echo=echo)
        additData = {}
        additData['-ldim'] = ldim
        additData['-ndim'] = ndim
        additData['ElementSpacing'] = ldim
        additData['DimSize'] = ndim
        additData['HeaderSize'] = 0
        if inFileName.find('.mhd') > -1:
            TM = [
             1, 0, 0, 0, 1, 0, 0, 0, 1]
            if 'TransformMatrix' in parDict:
                TM = parDict['TransformMatrix'].split()
            OF = [
             0.0, 0.0, 0.0]
            if 'Offset' in parDict:
                OF = parDict['Offset'].split()
            CR = [
             0.0, 0.0, 0.0]
            if 'CenterOfRotation' in parDict:
                CR = parDict['CenterOfRotation'].split()
            AnatomicalOrientation = 'LPS'
            if 'AnatomicalOrientation' in parDict:
                AnatomicalOrientation = parDict['AnatomicalOrientation']
            additData['TransformMatrix'] = [
             float(TM[0]), float(TM[1]), float(TM[2]),
             float(TM[3]), float(TM[4]), float(TM[5]),
             float(TM[6]), float(TM[7]), float(TM[8])]
            additData['Offset'] = [float(OF[0]), float(OF[1]), float(OF[2])]
            additData['CenterOfRotation'] = [float(CR[0]), float(CR[1]), float(CR[2])]
            additData['AnatomicalOrientation'] = AnatomicalOrientation.replace(' ', '')
            rawDict = self.getRawBinForm2()
            additData['ElementType'] = rawDict[rawInfo['binFormat']]
            additData['ElementDataFile'] = parDict['ElementDataFile']
        return (
         voxelModel, additData)

    def readNhdrFile(self, inFileName, info=False, echo=True):
        """
        Nrrd header file reader. It will call a raw reader to read the data. 
        
        @param  inFileName: input file name
          - TYPE: string
        @param  info:  if 'True' only additdata are return and voxelModel is set
                       to 'None'
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX]        int ret = QMessageBox::warning(this, tr("Item Cut"),
                     tr("Do you really like to cut this item?"),
                     QMessageBox::Yes | QMessageBox::No);
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - Int grayValue  ... gray value of voxel
        additdata: additional data information. Only in the case of an rawdp file.
          - TYPE: dict['-ndim': ndim, '-ldim':ldim, '-thres':thres]
          - list(int)   thres ... List of threshold values (0<...<255);
          - list(float32) ldim  ... x,y,z size of voxels (length = resolution)
          - list(int)   ndim  ... x,y,z number of voxels
        """
        if echo:
            stdout.write(' ... read & process nhdr data \n')
            stdout.flush()
        time1 = time.clock()
        try:
            file = open(inFileName)
            file.close()
        except IOError:
            stdout.write("\n **ERROR** mic.readNhdrFile(): intput file '%s' not found!\n\n" % inFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        osIn = open(inFileName, 'r')
        lines = osIn.readlines()
        osIn.close()
        parDict = {}
        lineNo = 0
        checkParDict = {}
        for line in lines:
            if lineNo > 0 and line[0] != '#' and len(line) > 1:
                line = line.replace('\n', '')
                key, data = dpUtils.userSplit(line)
                key = key.replace(' ', '')
                parDict[key] = data
                checkParDict[key] = False
            lineNo += 1

        if 'dimension' in parDict:
            checkParDict['dimension'] = True
            if int(parDict['dimension']) != 3:
                stdout.write("\n **ERROR**: 'dimension' field, only '3' is implemented!\n\n")
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            stdout.write("\n **ERROR**: 'dimension' field not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if 'type' in parDict:
            checkParDict['type'] = True
            rawInfo = {}
            format = parDict['type']
            format = format.replace(' ', '')
            if format == 'unsignedchar' or format == 'uchar' or format == 'uint8' or format == 'uint8_t':
                rawInfo['binFormat'] = 'B'
            elif format == 'short' or format == 'shortint' or format == 'signedshort' or format == 'signedshortint' or format == 'int16' or format == 'int16_t':
                rawInfo['binFormat'] = 'h'
            elif format == 'unsignedshort' or format == 'ushort' or format == 'unsignedshortint' or format == 'uint16' or format == 'uint16_t':
                rawInfo['binFormat'] = 'H'
            elif format == 'int' or format == 'signedint' or format == 'int32' or format == 'int32_t':
                rawInfo['binFormat'] = 'i'
            elif format == 'float':
                rawInfo['binFormat'] = 'f'
            else:
                stdout.write("\n **ERROR** mic.readMhd(): type field '%s' not supported!\n\n" % format)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            stdout.write("\n **ERROR** mic.readMhd(): 'type' field is not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        AnatomicalOrientation = None
        if 'space' in parDict:
            checkParDict['space'] = True
            AnatomicalOrientation = parDict['space'].replace(' ', '')
            if AnatomicalOrientation == 'right-anterior-superior' or AnatomicalOrientation == 'RAS':
                AnatomicalOrientation = 'RAS'
            elif AnatomicalOrientation == 'left-anterior-superior' or AnatomicalOrientation == 'LAS':
                AnatomicalOrientation = 'LAS'
            elif AnatomicalOrientation == 'left-posterior-superior' or AnatomicalOrientation == 'LPS':
                AnatomicalOrientation = 'LPS'
            else:
                stdout.write("\n **ERROR** mic.readMhd(): 'space' field '%s' not supported!\n\n" % AnatomicalOrientation)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if 'sizes' in parDict:
            checkParDict['sizes'] = True
            nlist = parDict['sizes'].split()
            if len(nlist) != 3:
                stdout.write("\n **ERROR**: mic.readNhdrFile(): 'sizes' field has to have 3 entries!\n\n")
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            else:
                ndim = [
                 int(nlist[0]), int(nlist[1]), int(nlist[2])]
        else:
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'sizes' field not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nvec = []
        length = []
        if 'spacedirections' in parDict:
            checkParDict['spacedirections'] = True
            svec = [' ', ' ', ' ']
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
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'space directions' field not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        ldim = [
         length[0], length[1], length[2]]
        if 'spacings' in parDict:
            checkParDict['spacings'] = True
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'spacings' field not supported!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if 'thicknesses' in parDict:
            checkParDict['thicknesses'] = True
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'thicknesses' field not supported!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if 'kinds' in parDict:
            checkParDict['kinds'] = True
            if parDict['kinds'].replace(' ', '') != 'domaindomaindomain':
                stdout.write("\n **ERROR**: mic.readNhdrFile(): 'kinds' has to be 'domain'!\n\n" % parDict['encoding'])
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'kinds' field not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        rawInfo['bigEndian'] = False
        if 'endian' in parDict:
            checkParDict['endian'] = True
            if parDict['endian'].replace(' ', '') == 'big':
                rawInfo['bigEndian'] = True
        rawInfo['fastestDir'] = 'x'
        if 'encoding' in parDict:
            checkParDict['encoding'] = True
            if parDict['encoding'].replace(' ', '') != 'raw':
                stdout.write("\n **ERROR**: mic.readNhdrFile(): 'encoding' is %s but has to be 'raw'!\n\n" % parDict['encoding'])
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'encoding' field not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        OF = [0.0, 0.0, 0.0]
        if 'spaceorigin' in parDict:
            checkParDict['spaceorigin'] = True
            svec = parDict['spaceorigin'].replace(' ', '').replace('(', '').replace(')', '')
            ox, oy, oz = svec.split(',')
            OF = [float(ox), float(oy), float(oz)]
        if 'datafile' in parDict:
            checkParDict['datafile'] = True
            filename = parDict['datafile']
            filename = filename.replace(' ', '')
            parDict['datafile'] = filename
            path, RawFile = os.path.split(filename)
            if path == '':
                nhdrPathRawFile = os.path.realpath(inFileName)
                path, RawNhdrFile = os.path.split(nhdrPathRawFile)
                filename = path + '/' + filename
            if not os.path.isfile(filename):
                stdout.write('\n **ERROR**: mic.readMhdFile(): %s is could not be found!\n\n' % filename)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            stdout.write("\n **ERROR**: mic.readNhdrFile(): 'data file' field not found!\n\n")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        rawInfo['headerSize'] = 0
        if 'byte skip' in parDict:
            checkParDict['byte skip'] = True
            rawInfo['headerSize'] = int(parDict['byte skip'])
        if info == True:
            voxelModel = None
        else:
            voxelModel = self.readGeneralFile3(filename, ndim, rawFormat=rawInfo['binFormat'], bigEndian=rawInfo['bigEndian'], echo=echo)
        additData = {}
        additData['-ldim'] = ldim
        additData['-ndim'] = ndim
        additData['ElementSpacing'] = ldim
        additData['DimSize'] = ndim
        additData['HeaderSize'] = rawInfo['headerSize']
        CR = [
         0.0, 0.0, 0.0]
        additData['TransformMatrix'] = [nvec[0][0], nvec[0][1], nvec[0][2],
         nvec[1][0], nvec[1][1], nvec[1][2],
         nvec[2][0], nvec[2][1], nvec[2][2]]
        additData['Offset'] = [float(OF[0]), float(OF[1]), float(OF[2])]
        additData['CenterOfRotation'] = [0.0, 0.0, 0.0]
        additData['AnatomicalOrientation'] = AnatomicalOrientation
        rawDict = self.getRawBinForm2()
        additData['ElementType'] = rawDict[rawInfo['binFormat']]
        additData['ElementDataFile'] = parDict['datafile']
        for key in checkParDict:
            if checkParDict[key] == False:
                stdout.write("\n **ERROR**: mic.readNhdrFile(): '%s' field not implemented!\n\n" % key)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)

        return (voxelModel, additData)

    def write(self, outFileName, curVoxelModel, thresList=None, dimList=None, format=None, history=None, template=None, smooth=None, geomList=None, muscal=None, echo=True):
        """
        Function selects a writing routine based on the given file extension.
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... value of voxel
        @param  thresList: Threshold list
              - TYPE: list[treshID] = value
              - int thresID ... threshold ID, 0,1,2 .....
              - int value   ... threshold value for this ID, 0...255
        @param  dimList: list of voxel dimension's
              - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
              - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param  format: binary format of output file
              - string value ... format specifier (Python conventions)
        @param  history: history of the performed modifications on the model
              - TYPE: dict[key] = values
              - int key    ... identifier of the applied option (e.g. -in, -thres)
              - int values ... values for this identifier
        @param  template: template file name for ABAQUS output
              - TYPE: string = filename
        @param  smooth: taubin voxel model smoothing parameter 
              - TYPE: list[0] = iter, list[1] = lambda, list[2] = kPB
              - int iter, float lambda, float kPB
        @param  muscal: Scaling parameter for aim files 
              - int factor
              
        @return:
            newVoxelModelvoxel: cut/cropped model
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... value of voxel
        """
        ext = self.getExtension(outFileName)
        upext = ext.upper()
        stdout.write(' ... write file %s\n' % outFileName)
        stdout.flush()
        if upext == 'RAW' or upext == 'RAWDP':
            if format == None:
                stdout.write("\n **ERROR**: write() '-form' option not optional for writing RAW/RAWDP data! \n")
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if upext == 'XML':
            self.writeXml(outFileName, curVoxelModel, dimList, format, history)
        elif upext == 'MHD':
            self.writeMhd(outFileName, curVoxelModel, dimList, format, history, geomList=geomList, echo=echo)
        elif upext == 'NHDR':
            self.writeNhdr(outFileName, curVoxelModel, dimList, format, history, geomList=geomList, echo=echo)
        elif upext == 'SVX':
            self.writeSvx(outFileName, curVoxelModel, dimList, format, history, geomList=geomList, echo=echo)
        elif upext == 'TXT':
            self.writeTxt(outFileName, curVoxelModel, dimList, echo=echo)
        elif upext == 'JPG' or upext == 'PNG' or upext == 'GIF' or upext == 'TIF' or upext == 'BMP':
            self.writeImages(outFileName, curVoxelModel, echo=echo)
        elif upext == 'BIN':
            self.writeGeneralRaw3(outFileName, curVoxelModel, rawFormat='B', bigEndian=False, echo=echo)
        elif upext == 'RAW':
            self.writeGeneralRaw3(outFileName, curVoxelModel, rawFormat=format, bigEndian=False, echo=echo)
        elif upext == 'RAWIV':
            self.writeRawiv(outFileName, curVoxelModel, dimList)
        elif upext == 'V':
            self.writeVista(outFileName, curVoxelModel, thresList, dimList)
        elif upext == 'RAWDP':
            self.writeRawdp(outFileName, curVoxelModel, thresList, dimList, format=format)
        elif upext == 'INP':
            self.writeAbaqusGeneral(outFileName, curVoxelModel, dimList, template, smooth)
        elif upext == 'MESH':
            self.writeParfe(outFileName, curVoxelModel, dimList, template, smooth)
        elif upext == 'FEM':
            self.writeOptistruct(outFileName, curVoxelModel, thresList, dimList)
        elif upext == 'IN':
            self.writeOOFEM(outFileName, curVoxelModel, dimList, template)
        elif upext == 'GZ':
            self.writeInrimage(outFileName, curVoxelModel, thresList, dimList, format)
        elif upext == 'AIM':
            if muscal == None:
                stdout.write('\n **ERROR** writeAim(): option -muscal has to be provided!\n\n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            else:
                if format == None:
                    stdout.write('\n **ERROR** writeAim(): option -form has to be int16!\n\n')
                    stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
                self.writeAim(outFileName, curVoxelModel, dimList, muscal, format=format)
        else:
            stdout.write('\n **ERROR** write(): output file extension of file: "%s" not known!\n\n' % outFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if not dimList == None:
            filename, ext = self.getFilenameAndExtension(outFileName)
        return

    def writeImages(self, outfilename, curVoxelModel, echo=True, zero=False):
        """
        Image file writer. Function writes *.png, *.jpg, *.gif files.
        The file type is choosen from the file extension.
        Upper left corner of first image (0,0) pixel is the (0,0,0) voxel at
        the "SWB" (south west bottom) voxel of an RVE. The image is located
        such that a viewer looks from outside to the bottom RVE surface.
        
        
        @param outfilename: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        
        @return:
          no return value
        """
        if echo:
            stdout.write(' ... write image files\n')
            stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'B')
        time1 = time.clock()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        if echo:
            dpUtils.progressStart('     -> Processed Images    : ')
        for n in range(nz):
            if zero == True:
                no = repr(n)
            else:
                no = repr(n + 1)
            if len(no) == 1:
                no = '000' + no
            if len(no) == 2:
                no = '00' + no
            if len(no) == 3:
                no = '0' + no
            if len(no) > 4:
                stdout.write('\n **ERROR** writeImages(). To many files to write!\n\n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            filename, ext = self.getFilenameAndExtension(outfilename)
            upext = ext.upper()
            if upext == 'JPG':
                upext = 'JPEG'
            name = filename + no + '.' + ext
            misc.toimage(curVoxelModel[n], cmin=0, cmax=255).save(name)
            progress = float(n + 1) / float(nz) * 10.0
            if echo:
                dpUtils.progressNext(progress)

        if echo:
            dpUtils.progressEnd()
        time2 = time.clock()
        if echo:
            stdout.write('     -> Processed Images    : %10i      \n' % (n + 1))
            stdout.flush()
        if echo:
            stdout.write('     -> write finished in   :   %8.1f sec\n' % (time2 - time1))
            stdout.flush()

    def writeRawiv(self, outFileName, curVoxelModel, dimList):
        """
        Datafile for LBIE - mesher. Not tested!
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        @return:
          no return value
        """
        stdout.write(' ... write rawiv data\n')
        stdout.flush()
        time1 = time.clock()
        if dimList == None:
            stdout.write("\n **ERROR** writeRawiv(): Voxel sizes '-ldim' not optional for this function!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        xvox = dimList[0]
        yvox = dimList[1]
        zvox = dimList[2]
        nx, ny, nz = self.get_Shape(curVoxelModel)
        checklist = []
        checklist.append(nz)
        checklist.append(ny)
        checklist.append(nx)
        octriVal = 0
        val = 0
        Octree = False
        for i in range(100):
            val = 2 ** i
            if (checklist[0] - 1) % val == 0 and (checklist[0] - 1) / val == 1 and (checklist[1] - 1) % val == 0 and (checklist[1] - 1) / val == 1 and (checklist[2] - 1) % val == 0 and (checklist[2] - 1) / val == 1:
                octriVal = val + 1
                Octree = True

        checkval = max(checklist)
        if not Octree:
            for i in range(100):
                val = 2 ** i + 1
                if checkval / val == 0 and checkval == checkval % val:
                    octriVal = val
                    break

            stdout.write(' ** WARNING ** writeRawiv(): Resolution %s not an Octree Structure.\n' % repr(checklist))
            stdout.flush()
            stdout.write('               Set Output Resolution to (2^n+1) = %s,%s,%s\n' % (octriVal, octriVal, octriVal))
            stdout.flush()
        minX = 0.0
        minY = 0.0
        minZ = 0.0
        maxX = float(xvox * (octriVal - 1))
        maxY = float(yvox * (octriVal - 1))
        maxZ = float(zvox * (octriVal - 1))
        dimX = octriVal
        dimY = octriVal
        dimZ = octriVal
        originX = 0.0
        originY = 0.0
        originZ = 0.0
        spanX = (maxX - minX) / (dimX - 1.0)
        spanY = (maxY - minY) / (dimY - 1.0)
        spanZ = (maxZ - minZ) / (dimZ - 1.0)
        numVerts = octriVal * octriVal * octriVal
        numCells = (octriVal - 1) * (octriVal - 1) * (octriVal - 1)
        stdout.write(' ... write rawiv data\n')
        stdout.flush()
        os = open(outFileName, 'w')
        os.write(struct.pack('>3f', minX, minY, minZ))
        os.write(struct.pack('>3f', maxX, maxY, maxZ))
        os.write(struct.pack('>I', numVerts))
        os.write(struct.pack('>I', numCells))
        os.write(struct.pack('>3I', dimX, dimY, dimZ))
        os.write(struct.pack('>3f', originX, originY, originZ))
        os.write(struct.pack('>3f', spanX, spanY, spanZ))
        sum = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for k in range(octriVal):
            for j in range(octriVal):
                progress = float(sum) / float(octriVal * octriVal * octriVal) * 10.0
                for i in range(octriVal):
                    sum = sum + 1
                    if k < nz and j < ny and i < nx:
                        os.write(struct.pack('>f', float(curVoxelModel[k, j, i])))
                    else:
                        os.write(struct.pack('>f', 0.0001))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        time2 = time.clock()
        stdout.write('     -> Processed Data      : %10i      \n' % sum)
        stdout.flush()
        stdout.write('     -> write finished in   : %8.1f sec \n' % (time2 - time1))
        stdout.flush()
        os.close
        return

    def writeGeneralRaw(self, outFileName, curVoxelModel, format='B'):
        """
        General raw file writer.
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param format: binary format specifier
            - TYPE: float
        
        @return:
          no return value
        """
        stdout.write(' ... write raw data\n')
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, format)
        time1 = time.clock()
        fileobj = open(outFileName, mode='wb')
        outvalues = aarray.array(format)
        time1 = time.clock()
        outvalues.fromlist(curVoxelModel.flat.tolist())
        time2 = time.clock()
        stdout.write('     -> make flat data in   : %8.2f sec\n' % (time2 - time1))
        stdout.flush()
        outvalues.tofile(fileobj)
        time3 = time.clock()
        stdout.write('     -> write to file in    : %8.2f sec\n' % (time3 - time2))
        stdout.flush()
        fileobj.close()

    def writeGeneralRaw3(self, outFileName, curVoxelModel, rawFormat='B', bigEndian=False, echo=True):
        """
        Fast general raw file writer using numpyIO ideas (see www).
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param rawFormat: Binary format specifier.
          - TYPE: char
        @param bigEndian: Flag if data are in big endian format.
          - TYPE: bool
          
        @return:
          no return value
        """
        if echo:
            stdout.write(' ... write raw data 3\n')
            stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, rawFormat)
        time1 = time.clock()
        if bigEndian == False:
            fp = numpy.memmap(outFileName, dtype=curVoxelModel.dtype, mode='w+', shape=curVoxelModel.shape)
            fp[:] = curVoxelModel[:]
            del fp
        else:
            curVoxelModel.byteswap().tofile(file=outFileName, sep='')
        time3 = time.clock()
        if echo:
            stdout.write('     -> write to file in    : %8.2f sec\n' % (time3 - time1))
            stdout.flush()

    def writeAim(self, outFileName, curVoxelModel, dimList, muscal, format='h'):
        """
        File writer for Scanco's AIM file format.
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param format: binary format specifier
            - TYPE: short
        
        @return:
          no return value
        """
        stdout.write(' ... write AIM data\n')
        stdout.flush()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        if format == 'h' or format == 'int16':
            dataSize = nx * ny * nz * 2
            aimTypeNum = 131074
        elif format == 'B' or format == 'uint8':
            dataSize = nx * ny * nz
            aimTypeNum = 65537
        else:
            stdout.write(' **ERROR** writeAimFile() format %s not supported!\n\n' % format)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        curVoxelModel = self.castType(curVoxelModel, format)
        time1 = time.clock()
        if dimList == None:
            stdout.write("\n **ERROR** writeRawdp(): Voxel dimension '-ldim' not optional for this function!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx, ny, nz = self.get_Shape(curVoxelModel)
        lx = dimList[0]
        ly = dimList[1]
        lz = dimList[2]
        aimVer = 16
        px = 0
        py = 0
        pz = 0
        ox = 0
        oy = 0
        oz = 0
        strMuscal = str(muscal)
        ProcLogText = '\n!\n! Processing Log\n!\n!-------------------------------------------------------------------------------\nCreated by                    Medtool\n!-------------------------------------------------------------------------------\nMu_Scaling                                       ' + strMuscal + '\n'
        ProcLogTextLen = len(ProcLogText)
        ProcLogFormat = repr(ProcLogTextLen) + 's'
        ProcLog = struct.pack(ProcLogFormat, ProcLogText)
        OF = open(outFileName, mode='w')
        sizePreHeader = 20
        sizeImStruct = 140
        sizeProcLog = ProcLogTextLen
        sizeImData = dataSize
        sizeImAssocData = 0
        OF.write(struct.pack('5i', sizePreHeader, sizeImStruct, sizeProcLog, sizeImData, sizeImAssocData))
        OF.write(struct.pack('i', aimVer))
        OF.write(struct.pack('4i', 0, 0, 0, 0))
        OF.write(struct.pack('i', aimTypeNum))
        OF.write(struct.pack('3i', px, py, pz))
        OF.write(struct.pack('3i', nx, ny, nz))
        OF.write(struct.pack('3i', ox, oy, oz))
        OF.write(struct.pack('12i', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        lx = lx * 4.0
        bin = struct.pack('f', lx)
        b1, b2, b3, b4 = struct.unpack('cccc', bin)
        OF.write(struct.pack('cccc', b3, b4, b1, b2))
        ly = ly * 4.0
        bin = struct.pack('f', ly)
        b1, b2, b3, b4 = struct.unpack('cccc', bin)
        OF.write(struct.pack('cccc', b3, b4, b1, b2))
        lz = lz * 4.0
        bin = struct.pack('f', lz)
        b1, b2, b3, b4 = struct.unpack('cccc', bin)
        OF.write(struct.pack('cccc', b3, b4, b1, b2))
        OF.write(struct.pack('5i', 0, 0, 0, 0, 0))
        OF.write(ProcLog)
        time2 = time.clock()
        sum = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for k in range(nz):
            for j in range(ny):
                progress = float(sum) / float(nx * ny * nz) * 10.0
                for i in range(nx):
                    sum = sum + 1
                    OF.write(struct.pack(format, curVoxelModel[k, j, i]))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        OF.close()
        time3 = time.clock()
        stdout.write('     -> write to file in    : %8.2f sec\n' % (time3 - time2))
        stdout.flush()
        return

    def writeRawdp(self, outFileName, curVoxelModel, thresList, dimList, format='f'):
        """
        File writer for Dieter Pahr's raw file format. This is
        a raw file with ASCII header! See self.readRawdp2().
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  thresList: Threshold list
            - TYPE: list[treshID] = value
            - int thresID ... threshold ID, 0,1,2 .....
            - int value   ... threshold value for this ID, 0...255
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param format: binary format specifier
            - TYPE: float
        
        @return:
          no return value
        """
        stdout.write(' ... write rawdp data\n')
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, format)
        time1 = time.clock()
        if dimList == None:
            stdout.write("\n **ERROR** writeRawdp(): Voxel dimension '-ldim' not optional for this function!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        os = open(outFileName, mode='w')
        nx, ny, nz = self.get_Shape(curVoxelModel)
        nList = [nx, ny, nz]
        os.write('voxelnumber   %s\n' % repr(nList))
        os.write('voxellength   %s\n' % repr(dimList))
        os.write('threshold     %s\n' % repr(thresList))
        os.write('format        %s\n' % format)
        os.write('endheader       \n')
        time2 = time.clock()
        sum = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for k in range(nz):
            for j in range(ny):
                progress = float(sum) / float(nx * ny * nz) * 10.0
                for i in range(nx):
                    sum = sum + 1
                    os.write(struct.pack(format, curVoxelModel[k, j, i]))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        time3 = time.clock()
        stdout.write('     -> write to file in    : %8.2f sec\n' % (time3 - time2))
        stdout.flush()
        return

    def writeInrimage(self, outFileName, curVoxelModel, thresList, dimList, format):
        """
        File writer for Inria raw file format. This is
        a raw file with ASCII header! 
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  thresList: Threshold list
            - TYPE: list[treshID] = value
            - int thresID ... threshold ID, 0,1,2 .....
            - int value   ... threshold value for this ID, 0...255
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        
        @return:
          no return value
        """
        stdout.write(' ... write inr data\n')
        stdout.flush()
        if format == None:
            format = 'f'
        curVoxelModel = self.castType(curVoxelModel, format)
        time1 = time.clock()
        if dimList == None:
            stdout.write("\n **ERROR** writeInrimage(): Voxel dimension '-ldim' not optional for this function!\n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        fos = gzip.open(outFileName, 'wb')
        vType = ''
        bits = ''
        if format == 'B':
            vType = 'unsigned fixed'
            bits = '8'
        elif format == 'H':
            vType = 'unsigned fixed'
            bits = '16'
        elif format == 'h':
            vType = 'signed fixed'
            bits = '16'
        elif format == 'i':
            vType = 'signed fixed'
            bits = '32'
        elif format == 'f':
            vType = 'float'
            bits = '32'
        else:
            stdout.write("\n **ERROR** mic.writeInriaImage(): Format '%s' not supported!\n\n" % format)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx, ny, nz = self.get_Shape(curVoxelModel)
        endspec = '\n'
        header = ''
        end = '##}' + endspec
        mid = ''
        header = header + '#INRIMAGE-4#{' + endspec
        header = header + 'XDIM=' + repr(nx) + endspec
        header = header + 'YDIM=' + repr(ny) + endspec
        header = header + 'ZDIM=' + repr(nz) + endspec
        header = header + 'VDIM=1' + endspec
        header = header + 'VX=' + repr(dimList[0]) + endspec
        header = header + 'VY=' + repr(dimList[1]) + endspec
        header = header + 'VZ=' + repr(dimList[2]) + endspec
        header = header + 'TYPE=' + vType + endspec
        header = header + 'PIXSIZE=' + bits + ' bits' + endspec
        header = header + 'SCALE=2**0' + endspec
        header = header + 'CPU=pc' + endspec
        for i in range(256 - len(header) - len(end)):
            mid += endspec

        header = header + mid + end
        fos.write('%s' % header)
        time2 = time.clock()
        sum = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for k in range(nz):
            for j in range(ny):
                progress = float(sum) / float(nx * ny * nz) * 10.0
                for i in range(nx):
                    sum = sum + 1
                    fos.write(struct.pack(format, curVoxelModel[k, j, i]))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        fos.close()
        time3 = time.clock()
        stdout.write('     -> write to file in    : %8.2f sec\n' % (time3 - time2))
        stdout.flush()
        return

    def writeVista(self, outFileName, curVoxelModel, thresList, dimList):
        """
        Vista file writer for VGRID.
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  thresList: Threshold list
            - TYPE: list[treshID] = value
            - int thresID ... threshold ID, 0,1,2 .....
            - int value   ... threshold value for this ID, 0...255
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        @return:
          no return value
        """
        stdout.write(' ... write vista data file\n')
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'B')
        time1 = time.clock()
        sum = 0
        os = open(outFileName, 'w')
        if thresList == None:
            stdout.write('\n **ERROR** writeVista(): Threshold not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dimList == None:
            stdout.write('\n **ERROR** writeVista(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx, ny, nz = self.get_Shape(curVoxelModel)
        os.write('V-data 2 {\n')
        os.write('       image: image {\n')
        os.write('        data: 0\n')
        os.write('        length: %15i\n' % (nx * ny * nz))
        os.write('        nbands: %15i\n' % nz)
        os.write('        nframes: %15i\n' % nz)
        os.write('        nrows: %15i\n' % nx)
        os.write('        ncolumns: %15i\n' % ny)
        os.write('        repn: ubyte\n')
        os.write('        voxel: "%15.9f %15.9f %15.9f"\n' % (dimList[0], dimList[1], dimList[2]))
        os.write('        }\n')
        os.write('}\n')
        os.write('\x0c\n')
        dpUtils.progressStart('     -> Processed Data      : ')
        for k in range(nz):
            for j in range(ny):
                progress = float(sum) / float(nx * ny * nz) * 10.0
                for i in range(nx):
                    sum = sum + 1
                    value = int(curVoxelModel[k, j, i])
                    if value < thresList[0]:
                        cv = '0'
                    else:
                        cv = '1'
                    os.write(struct.pack('c', cv))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        time2 = time.clock()
        stdout.write('     -> Processed Data      : %10i      \n' % sum)
        stdout.flush()
        stdout.write('     -> write finished in   :   %8.1f sec\n' % (time2 - time1))
        stdout.flush()
        os.close
        return

    def writeAbaqusGeneral(self, outFileName, curVoxelModel, dimList, templateFile, smooth):
        """
        General Abaqus *.inp file writer. For these materials a default material will be 
        applied. Supported commands: 
          *USER NODE
          *USER ELEMENT 
          *USER NSET, type=point, location=arbitrary
            generate NSET: ARB_NODE_S, ARB_NODE_N, ARB_NODE_E, ARB_NODE_W, ARB_NODE_T, ARB_NODE_B 
          *USER NSET, type=point, location=addcorner
            generate NSET: ACOR_NODE_SWB, ACOR_NODE_SEB, ACOR_NODE_NEB, ACOR_NODE_NWB,
                           ACOR_NODE_SWT, ACOR_NODE_SET, ACOR_NODE_NET, ACOR_NODE_NWT
          *USER NSET, type=face 
            generate NSET: ALL_NODE_S, ALL_NODE_N, ALL_NODE_E, ALL_NODE_W, ALL_NODE_T, ALL_NODE_B 
          *USER ELSET, type=face
            generate ELSET: ALL_S, ALL_ELEM_N, ALL_ELEM_E, ALL_ELEM_W, ALL_ELEM_T, ALL_ELEM_B      
          *USER PROPERTY, file=property_temp.inp, range=5:367
            generate multiple material cards, internal variables are "SetName, CardName, GrayValue"
            for the given example: GrayValues > 5 and GrayValues <= 367 are written 
            This card can be used multiple times
            If range=... is not given, material cards for all GrayValues are written
        Elements are only written for the given ranges in *USER PROPERTY
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param  templateFile: name of the template file           
        @param  smoothParam: taubin voxel model smoothing parameter 
            - TYPE: list[0] = iter, list[1] = lambda, list[2] = kPB, 
                    list[3] = nearIntf, list[4] = bcid, list[5] = shrink 
            - int iter, float lambda, float kPB, int nearIntf      
        
        @return:
          no return value
        """
        stdout.write(' ... setup ABAQUS *.inp file from template\n')
        stdout.write("     -> recast model from '%s' to 'i'\n" % curVoxelModel.dtype.char)
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'i')
        time1 = time.clock()
        if dimList.all() == None:                                                                         # 12.01.01 change: if dimList == None:     to if dimList.all() == None
            print '\n **ERROR** writeAbaqusGeneral(): Voxel size not optional for this function!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        xvox = dimList[0]
        yvox = dimList[1]
        zvox = dimList[2]
        nx, ny, nz = self.get_Shape(curVoxelModel)
        minVox, maxVox = self.computeNumpyMinMax(curVoxelModel, 0)
        minVox = int(minVox + 0.5)
        maxVox = int(maxVox + 0.5)
        noid = 0
        activeNodes = {}
        nodeSets = {}
        nodeSets['ALL_NODE_S'] = []
        nodeSets['ALL_NODE_N'] = []
        nodeSets['ALL_NODE_E'] = []
        nodeSets['ALL_NODE_W'] = []
        nodeSets['ALL_NODE_T'] = []
        nodeSets['ALL_NODE_B'] = []
        elemSets = {}
        elemSets['ALL_ELEM_S'] = []
        elemSets['ALL_ELEM_N'] = []
        elemSets['ALL_ELEM_E'] = []
        elemSets['ALL_ELEM_W'] = []
        elemSets['ALL_ELEM_T'] = []
        elemSets['ALL_ELEM_B'] = []
        tempflag = False
        if templateFile == None:
            tempflag = True
            OS = open('temp.inp', 'w')
            OS.write('*USER NODE\n*USER ELEMENT\n*USER PROPERTY, file=prop.inp, range=1:255\n')
            templateFile = 'temp.inp'
            OS.close()
            OS = open('prop.inp', 'w')
            OS.write('*SOLID SECTION, ELSET=SetName, MATERIAL=CardName\n1.\n')
            OS.write('*MATERIAL,NAME=CardName\n')
            OS.write('*ELASTIC\n')
            OS.write('20000., 0.3\n')
            OS.close()
            templateFile = 'temp.inp'
        OS = open(outFileName, 'w')
        try:
            osTempFile = open(templateFile, 'r')
        except IOError:
            stdout.write("\n **ERROR** mic.writeAbaqusGeneral(): Abaqus Template file '%s' not found!\n\n" % templateFile)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        lines = osTempFile.readlines()
        elsetNodes = {}
        ranges = {}
        thresList = []
        rangeMin = 0
        rangeMax = 255
        outFlag = False
        overlap = numpy.zeros(rangeMax + 1, numpy.int)
        for line in lines:
            line = line.replace('\n', '')
            if line.upper().find('*USER PROPERTY') == 0:
                line = line.replace(' ', '')
                args = line.split(',')
                matTemplateFilename = ''
                for arg in args:
                    if arg.upper().find('RANGE') == 0:
                        dummy, rangeStr = arg.split('=')
                        rangeMin, rangeMax = dpUtils.userSplit(rangeStr)
                        rangeMin = int(rangeMin)
                        rangeMax = int(rangeMax)
                        if rangeMin < 1:
                            stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): Minimum Range < 1!\n\n')
                            stdout.write('\n E N D E D  with ERRORS \n\n')
                            stdout.flush()
                            exit(1)
                        if rangeMax > maxVox:
                            outFlag = True
                        for ii in range(rangeMax - rangeMin + 1):
                            overlap[rangeMin + ii] += 1

                    if arg.upper().find('FILE') == 0:
                        dummy, matTemplateFilename = arg.split('=')

                ranges[matTemplateFilename] = (
                 rangeMin, rangeMax)

        if len(ranges) == 0:
            stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): *USER PROPERTY: keyword missing!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if rangeMax > maxVox:
            stdout.write('\n **WARNING** mic.writeAbaqusGeneral(): *USER PROPERTY: Max GV Range (%i) > Max Image GV (%i)!\n\n' % (rangeMax, maxVox))
        if numpy.sum(numpy.greater(overlap, 1)) > 0:
            stdout.write('\n **ERROR** mic.writeAbaqusGeneral(): *USER PROPERTY: Ranges in property template overlap!\n\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        for crange in ranges:
            for matId in range(ranges[crange][0], ranges[crange][1] + 1):
                elsetNodes[repr(matId)] = []
                thresList.append(matId)

        elid = 0
        nx1 = nx + 1
        nxy1 = (ny + 1) * (nx + 1)
        sum = 0
        dpUtils.progressStart('     -> setup Element Data  : ')
        for k in range(nz):
            sum += 1
            progress = float(sum) / float(nz) * 10.0
            for j in range(ny):
                for i in range(nx):
                    grayValue = curVoxelModel[k, j, i]
                    if elsetNodes.has_key(repr(grayValue)):
                        elid = elid + 1
                        elnds = [nxy1 * k + nx1 * j + (i + 1),
                         nxy1 * k + nx1 * j + (i + 2),
                         nxy1 * k + nx1 * (j + 1) + (i + 2),
                         nxy1 * k + nx1 * (j + 1) + (i + 1),
                         nxy1 * (k + 1) + nx1 * j + (i + 1),
                         nxy1 * (k + 1) + nx1 * j + (i + 2),
                         nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2),
                         nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)]
                        elsetNodes[repr(grayValue)].append((elid, elnds))
                        if k == 0:
                            elemSets['ALL_ELEM_B'].append(elid)
                        if k == nz - 1:
                            elemSets['ALL_ELEM_T'].append(elid)
                        if j == 0:
                            elemSets['ALL_ELEM_S'].append(elid)
                        if j == ny - 1:
                            elemSets['ALL_ELEM_N'].append(elid)
                        if i == 0:
                            elemSets['ALL_ELEM_W'].append(elid)
                        if i == nx - 1:
                            elemSets['ALL_ELEM_E'].append(elid)

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        stdout.write('     -> setup Node Data     :')
        for matid in thresList:
            if len(elsetNodes[repr(matid)]) > 0:
                matidStr = 'SET' + repr(matid)
                for elnds in elsetNodes[repr(matid)]:
                    elid = elnds[0]
                    for elnd in elnds[1]:
                        activeNodes[elnd] = 1

        noid = 0
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    noid = noid + 1
                    if activeNodes.has_key(noid):
                        if k == 0:
                            nodeSets['ALL_NODE_B'].append(noid)
                        if k == nz:
                            nodeSets['ALL_NODE_T'].append(noid)
                        if j == 0:
                            nodeSets['ALL_NODE_S'].append(noid)
                        if j == ny:
                            nodeSets['ALL_NODE_N'].append(noid)
                        if i == 0:
                            nodeSets['ALL_NODE_W'].append(noid)
                        if i == nx:
                            nodeSets['ALL_NODE_E'].append(noid)

        stdout.write(' Done\n')
        stdout.flush()
        nodeCoord = {}
        nodeCoordOrig = {}
        if smooth != None:
            activeNodes2, nElem, nNode, nIntElem, nIntFaces, nIntNode = dpMesherf77.check_voxmesh2d(curVoxelModel, smooth[4], 2)
            nodeCoordF77, nodeCoordInt, noidF77PY, noidIntVoxF77 = dpMesherf77.get_voxmesh2d_nodes(curVoxelModel, dimList, smooth, activeNodes2, nNode, nIntElem, nIntNode, 2)
            for noidF77 in range(len(noidF77PY)):
                noid = noidF77PY[noidF77]
                nodeCoord[noid] = (nodeCoordF77[noidF77][0], nodeCoordF77[noidF77][1], nodeCoordF77[noidF77][2])

        else:
            noid = 0
            for k in range(nz + 1):
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        noid = noid + 1
                        if activeNodes.has_key(noid):
                            nodeCoord[noid] = (
                             float(xvox * i), float(yvox * j), float(zvox * k))

        curPathFilename, ext = self.getFilenameAndExtension(outFileName)
        curFilename, ext = self.getShortFilenameAndExtension(outFileName)
        stdout.write(' ... write ABAQUS *.inp file from template\n')
        for line in lines:
            line = line.replace('\n', '')
            line = line.replace('$filename', curFilename)
            line = line.replace('$pathfilename', curPathFilename)
            if line.upper().find('*USER NODE') > -1:
                OS.write('*NODE\n')
                noid2 = 0
                noid = 0
                dpUtils.progressStart('     -> process Node IDs    : ')
                for k in range(nz + 1):
                    progress = float(k + 1) / float(nz + 1) * 10.0
                    for j in range(ny + 1):
                        for i in range(nx + 1):
                            noid = noid + 1
                            if activeNodes.has_key(noid):
                                noid2 = noid2 + 1
                                OS.write('%12i,%13.6g,%13.6g,%13.6g\n' % (noid, nodeCoord[noid][0], nodeCoord[noid][1], nodeCoord[noid][2]))

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
                stdout.write('     -> write Nodes         : %10i \n' % noid2)
                stdout.flush()
            elif line.upper().find('*USER ELEMENT') > -1:
                count = 0
                dpUtils.progressStart('     -> process Elements    : ')
                for matid in thresList:
                    count += 1
                    progress = count / float(len(thresList)) * 10.0
                    if len(elsetNodes[repr(matid)]) > 0:
                        matidStr = 'SET' + repr(matid)
                        OS.write('*ELEMENT, TYPE=C3D8, ELSET=%s\n' % matidStr)
                        for elnds in elsetNodes[repr(matid)]:
                            elid = elnds[0]
                            OS.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (elid, elnds[1][0], elnds[1][1], elnds[1][2], elnds[1][3], elnds[1][4], elnds[1][5], elnds[1][6], elnds[1][7]))
                            for elnd in elnds[1]:
                                activeNodes[elnd] = 1

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
                stdout.write('     -> write Elements      : %10i             \n' % elid)
                stdout.flush()
            elif line.upper().find('*USER NSET') > -1:
                if line.upper().find('TYPE=FACE') > -1:
                    stdout.write('     -> write BCs Node Sets     \n')
                    stdout.flush()
                    for nsetName in nodeSets:
                        OS.write('*NSET, NSET=%s\n' % nsetName)
                        entry = 0
                        for noid in nodeSets[nsetName]:
                            entry = entry + 1
                            if entry == 16:
                                OS.write('%s' % repr(noid))
                                entry = 0
                                OS.write('\n')
                            else:
                                OS.write('%s,' % repr(noid))

                        OS.write('\n')

                if line.upper().find('TYPE=POINT') > -1:
                    if line.upper().find('LOCATION=ARBITRARY') > -1:
                        for nsetName in nodeSets:
                            if len(nodeSets[nsetName]) > 0:
                                nid = nodeSets[nsetName][0]
                                name = nsetName.replace('ALL_NODE_', 'ARB_NODE_')
                                OS.write('*NSET, NSET=%s\n' % name)
                                OS.write('%s\n' % repr(nid))

                    if line.upper().find('LOCATION=ADDCORNER') > -1:
                        nid = (nx + 1) * (ny + 1) * (nz + 1)
                        OS.write('*NODE, NSET=ACOR_NODE_SWB\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 1, 0.0, 0.0, 0.0))
                        OS.write('*NODE, NSET=ACOR_NODE_SEB\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 2, nx * xvox, 0.0, 0.0))
                        OS.write('*NODE, NSET=ACOR_NODE_NEB\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 3, nx * xvox, ny * yvox, 0.0))
                        OS.write('*NODE, NSET=ACOR_NODE_NWB\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 4, 0.0, ny * yvox, 0.0))
                        OS.write('*NODE, NSET=ACOR_NODE_SWT\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 5, 0.0, 0.0, nz * zvox))
                        OS.write('*NODE, NSET=ACOR_NODE_SET\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 6, nx * xvox, 0.0, nz * zvox))
                        OS.write('*NODE, NSET=ACOR_NODE_NET\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 7, nx * xvox, ny * yvox, nz * zvox))
                        OS.write('*NODE, NSET=ACOR_NODE_NWT\n')
                        OS.write('%i, %13.6g,  %13.6g, %13.6g\n' % (nid + 8, 0.0, ny * yvox, nz * zvox))
            elif line.upper().find('*USER ELSET') > -1:
                if line.upper().find('TYPE=FACE') > -1:
                    stdout.write('     -> Write BCs Elem Sets          \n')
                    stdout.flush()
                    for elsetName in elemSets:
                        OS.write('*ELSET, ELSET=%s\n' % elsetName)
                        entry = 0
                        for elid in elemSets[elsetName]:
                            entry = entry + 1
                            if entry == 16:
                                OS.write('%s' % repr(elid))
                                entry = 0
                                OS.write('\n')
                            else:
                                OS.write('%s,' % repr(elid))

                        OS.write('\n')

            elif line.upper().find('*USER PROPERTY') > -1:
                line = line.replace(' ', '')
                args = line.split(',')
                rangeMin = minVox
                rangeMax = maxVox
                matTemplateFilename = ''
                for arg in args:
                    if arg.upper().find('RANGE') == 0:
                        dummy, rangeStr = arg.split('=')
                        rangeMin, rangeMax = dpUtils.userSplit(rangeStr)
                        rangeMin = int(rangeMin)
                        rangeMax = int(rangeMax)
                    if arg.upper().find('FILE') == 0:
                        dummy, matTemplateFilename = arg.split('=')

                stdout.write('     -> Write Property      : %s \n' % matTemplateFilename)
                try:
                    osMatCard = open(matTemplateFilename, 'r')
                except IOError:
                    stdout.write("\n **ERROR** writeAbaqusGeneral(): Material template file '%s' not found!\n\n" % matTemplateFilename)
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

                lines = osMatCard.readlines()
                for matid in thresList:
                    GrayValue = matid
                    if len(elsetNodes[repr(matid)]) > 0:
                        if matid >= rangeMin and matid <= rangeMax:
                            matidStr = 'SET' + repr(matid)
                            GrayValue = matid
                            for line in lines:
                                line = line.replace('\n', '')
                                if line.find('SetName') > -1:
                                    line = line.replace('SetName', matidStr)
                                if line.find('CardName') > -1:
                                    line = line.replace('CardName', 'MAT' + matidStr)
                                if line.find('GrayValue') > -1:
                                    exprList = line.split(',')
                                    count = 0
                                    for expr in exprList:
                                        if expr.find('GrayValue') > -1:
                                            compValue = eval(expr)
                                            OS.write('%s' % repr(compValue))
                                        else:
                                            OS.write('%s' % expr)
                                        if count < len(exprList) - 1:
                                            OS.write(',')
                                        count += 1

                                    OS.write('\n')
                                else:
                                    OS.write('%s\n' % line)

                osMatCard.close()
            else:
                OS.write('%s\n' % line)

        osTempFile.close()
        if tempflag:
            os.remove('temp.inp')
            os.remove('prop.inp')
        time2 = time.clock()
        stdout.write('     -> Write finished in   :   %8.1f sec  \n' % (time2 - time1))
        stdout.flush()
        OS.close
        return

    def writeParfe(self, outFileName, curVoxelModel, dimList, templateFile, smooth):
        """
        General Parfe *.mesh file writer. For these materials a default material will be 
        applied. Supported commands: 
          *USER HEADER
            Generates header with: 
              # l1: elem_type, nElems, nNodes, nElemIp, nNodeDof, nElemNodes, nDimSE, nDim
              # l2: lx, ly, lz (isotropic element dimenions in mm)
          *USER PROPERTY
             generate multiple material cards, internal variables is GrayValue
             Example: 
              *USER PROPERTY, range=70:80, param=5400.*(GrayValue/250.)**2.5:0.3
             means 
               -> GrayValues > 70 and GrayValues <= 80 are written
               -> 
             This card can be used multiple times
             If range=... is not given, material cards for all GrayValues are written
             Elements are only written for the given ranges in *USER PROPERTY
             
          *USER NODE
          *USER ELEMENT
          *USER BC1
          *USER BC2
          *USER BC3
             generate BCs, internal sets are  ALL_NODE_S, ALL_NODE_N, ALL_NODE_E, 
             ALL_NODE_W, ALL_NODE_T, ALL_NODE_B, ID_NODE_???. '???' is the node ID
             which follows from : 
                ID_Node_1 = nxy1*(k)   + nx1*(j)   + (i+1)
                ID_Node_2 = nxy1*(k)   + nx1*(j)   + (i+2)
                ID_Node_3 =  nxy1*(k)  + nx1*(j+1) + (i+2)
                ID_Node_4 = nxy1*(k)   + nx1*(j+1) + (i+1)
                ID_Node_5 = nxy1*(k+1) + nx1*(j)   + (i+1)
                ID_Node_6 = nxy1*(k+1) + nx1*(j)   + (i+2)
                ID_Node_7 = nxy1*(k+1) + nx1*(j+1) + (i+2)
                ID_Node_8 = nxy1*(k+1) + nx1*(j+1) + (i+1)
              with nxy1=(nx+1)*(ny+1), nx1=(nx+1), i=0..(nx-1), j=0..(ny-1), k=0..(nz-1)
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param  templateFile: name of the template file           
        @param  smoothParam: taubin voxel model smoothing parameter 
            - TYPE: list[0] = iter, list[1] = lambda, list[2] = kPB, 
                    list[3] = nearIntf, list[4] = bcid, list[5] = shrink 
            - int iter, float lambda, float kPB, int nearIntf      
        
        @return:
          no return value
        """
        stdout.write(' ... setup ParFE ASCII file from template\n')
        stdout.write("     -> recast model from '%s' to 'i'\n" % curVoxelModel.dtype.char)
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'i')
        time1 = time.clock()
        if dimList == None:
            print '\n **ERROR** writeParfe(): Voxel size not optional for this function!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        xvox = dimList[0]
        yvox = dimList[1]
        zvox = dimList[2]
        nx, ny, nz = self.get_Shape(curVoxelModel)
        minVox, maxVox = self.computeNumpyMinMax(curVoxelModel, 0)
        minVox = int(minVox + 0.5)
        maxVox = int(maxVox + 0.5)
        noid = 0
        activeNodes = {}
        nodeSets = {}
        nodeSets['ALL_NODE_S'] = []
        nodeSets['ALL_NODE_N'] = []
        nodeSets['ALL_NODE_E'] = []
        nodeSets['ALL_NODE_W'] = []
        nodeSets['ALL_NODE_T'] = []
        nodeSets['ALL_NODE_B'] = []
        elemSets = {}
        elemSets['ALL_ELEM_S'] = []
        elemSets['ALL_ELEM_N'] = []
        elemSets['ALL_ELEM_E'] = []
        elemSets['ALL_ELEM_W'] = []
        elemSets['ALL_ELEM_T'] = []
        elemSets['ALL_ELEM_B'] = []
        if templateFile == None:
            stdout.write('\n **ERROR** mic.writeParfe(): parFe Template file not given!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        OS = open(outFileName, 'w')
        try:
            osTempFile = open(templateFile, 'r')
        except IOError:
            stdout.write("\n **ERROR** mic.writeParfe(): parFe Template file '%s' not found!\n\n" % matCardFileName)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        lines = osTempFile.readlines()
        elsetNodes = {}
        ranges = {}
        thresList = []
        rangeMin = minVox
        rangeMax = maxVox
        for line in lines:
            line = line.replace('\n', '')
            if line.upper().find('*USER PROPERTY') > -1:
                line = line.replace(' ', '')
                args = line.split(',')
                matTemplateFilename = ''
                for arg in args:
                    if arg.upper().find('RANGE') == 0:
                        dummy, rangeStr = arg.split('=')
                        rangeMin, rangeMax = dpUtils.userSplit(rangeStr)
                        rangeMin = int(rangeMin)
                        rangeMax = int(rangeMax)
                    if arg.upper().find('PARAM') == 0:
                        dummy, paramStr = arg.split('=')
                        propertyStr = paramStr.replace(':', ' ')

                ranges[matTemplateFilename] = (
                 rangeMin, rangeMax)
            if line.upper().find('*USER DIRECT PROPERTY') > -1:
                if smooth == None:
                    nodeCoordF77 = numpy.zeros((1, 3), numpy.float)
                    error = micf77.writeparfeascii(outFileName, curVoxelModel, dimList, templateFile, 0, nodeCoordF77, 2)
                else:
                    activeNodes2, nElem, nNode, nIntElem, nIntFaces, nIntNode = dpMesherf77.check_voxmesh2d(curVoxelModel, smooth[4], 2)
                    nodeCoordF77, nodeCoordInt, noidF77PY, noidIntVoxF77 = dpMesherf77.get_voxmesh2d_nodes(curVoxelModel, dimList, smooth, activeNodes2, nNode, nIntElem, nIntNode, 2)
                    error = micf77.writeparfeascii(outFileName, curVoxelModel, dimList, templateFile, 1, nodeCoordF77, 2)
                if error:
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                return

        for crange in ranges:
            for matId in range(ranges[crange][0], ranges[crange][1] + 1):
                elsetNodes[repr(matId)] = []
                thresList.append(matId)

        sum = 0
        elid = 0
        nx1 = nx + 1
        nxy1 = (ny + 1) * (nx + 1)
        nElems = 0
        dpUtils.progressStart('     -> setup Element Data  : ')
        for k in range(nz):
            progress = float(k + 1) / float(nz) * 10.0
            for j in range(ny):
                for i in range(nx):
                    grayValue = curVoxelModel[k, j, i]
                    if grayValue != 0 and grayValue <= rangeMax and grayValue > rangeMin:
                        elid = elid + 1
                        elnds = [nxy1 * k + nx1 * j + (i + 1),
                         nxy1 * k + nx1 * j + (i + 2),
                         nxy1 * k + nx1 * (j + 1) + (i + 2),
                         nxy1 * k + nx1 * (j + 1) + (i + 1),
                         nxy1 * (k + 1) + nx1 * j + (i + 1),
                         nxy1 * (k + 1) + nx1 * j + (i + 2),
                         nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2),
                         nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)]
                        elsetNodes[repr(grayValue)].append((elid, elnds))
                        if k == 0:
                            elemSets['ALL_ELEM_B'].append(elid)
                        if k == nz - 1:
                            elemSets['ALL_ELEM_T'].append(elid)
                        if j == 0:
                            elemSets['ALL_ELEM_S'].append(elid)
                        if j == ny - 1:
                            elemSets['ALL_ELEM_N'].append(elid)
                        if i == 0:
                            elemSets['ALL_ELEM_W'].append(elid)
                        if i == nx - 1:
                            elemSets['ALL_ELEM_E'].append(elid)

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        nElems = elid
        stdout.write('     -> Setup Node Data     :')
        for matid in thresList:
            if len(elsetNodes[repr(matid)]) > 0:
                matidStr = 'SET' + repr(matid)
                for elnds in elsetNodes[repr(matid)]:
                    elid = elnds[0]
                    for elnd in elnds[1]:
                        activeNodes[elnd] = 1

        noid = 0
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    noid = noid + 1
                    if activeNodes.has_key(noid):
                        if k == 0:
                            nodeSets['ALL_NODE_B'].append(noid)
                        if k == nz:
                            nodeSets['ALL_NODE_T'].append(noid)
                        if j == 0:
                            nodeSets['ALL_NODE_S'].append(noid)
                        if j == ny:
                            nodeSets['ALL_NODE_N'].append(noid)
                        if i == 0:
                            nodeSets['ALL_NODE_W'].append(noid)
                        if i == nx:
                            nodeSets['ALL_NODE_E'].append(noid)

        stdout.write(' Done\n')
        stdout.flush()
        nodeCoord = {}
        nodeCoordOrig = {}
        if smooth != None:
            noid = 0
            for k in range(nz + 1):
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        noid = noid + 1
                        if activeNodes.has_key(noid):
                            nodeCoordOrig[noid] = (
                             float(xvox * i), float(yvox * j), float(zvox * k))

            activeNodes2, nElem, nNode, nIntElem, nIntNode = micf77.check_voxelsmoother(curVoxelModel, smooth[4], 2)
            nodeCoordF77, noidF77PY = micf77.computesmoothvox(curVoxelModel, dimList, smooth, activeNodes2, nElem, nNode, nIntElem, nIntNode, 2)
            for noidF77 in range(len(noidF77PY)):
                noid = noidF77PY[noidF77]
                nodeCoord[noid] = (nodeCoordF77[noidF77][0], nodeCoordF77[noidF77][1], nodeCoordF77[noidF77][2])

        else:
            noid = 0
            for k in range(nz + 1):
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        noid = noid + 1
                        if activeNodes.has_key(noid):
                            nodeCoord[noid] = (
                             float(xvox * i), float(yvox * j), float(zvox * k))

            lineId = 0
            linesToRemove = {}
            bcDict = {}
            while lineId < len(lines):
                line = lines[lineId]
                line = line.replace('\n', '')
                if line.upper().find('*USER BC') == 0:
                    bcIdName = line
                    bcDict[bcIdName] = []
                    lineId += 1
                    nBC = int(lines[lineId])
                    linesToRemove[lineId] = True
                    for iBC in range(nBC):
                        lineId += 1
                        line = lines[lineId]
                        linesToRemove[lineId] = True
                        lineList = line.split()
                        appendStr = ''
                        for eId in range(1, len(lineList)):
                            appendStr = appendStr + lineList[eId] + ' '

                        nsetName = lineList[0].replace(' ', '')
                        if nsetName.find('ALL_NODE') > -1:
                            for noid in nodeSets[nsetName]:
                                bcDict[bcIdName].append((noid, appendStr))

                        if nsetName.find('ID_NODE') > -1:
                            idList = nsetName.split('_')
                            noid = int(idList[2])
                            bcDict[bcIdName].append((noid, appendStr))

                    bcDict[bcIdName]
                lineId += 1

        stdout.write(' ... write ParFE ASCII file from template\n')
        nodeIdMap = {}
        nMat = 0
        lineId = 0
        for line in lines:
            if linesToRemove.has_key(lineId):
                lineId += 1
                continue
            lineId += 1
            line = line.replace('\n', '')
            if line.upper().find('*USER HEADER') == 0:
                OS.write("'hexahedron' %i %i 8 3 8 6 3\n" % (nElems, len(nodeCoord)))
                OS.write('%13.6g %13.6g %13.6g\n' % (xvox, yvox, zvox))
                continue
            elif line.upper().find('*USER PROPERTY') == 0:
                for matid in elsetNodes:
                    if len(elsetNodes[matid]) > 0:
                        nMat += 1

                OS.write('2 %i\n' % nMat)
                for matid in elsetNodes:
                    if len(elsetNodes[matid]) > 0:
                        GrayValue = int(matid)
                        if nMat > 1:
                            OS.write('%i ' % int(matid))
                        if propertyStr.find('GrayValue') > -1:
                            exprList = propertyStr.split()
                            count = 0
                            for expr in exprList:
                                if expr.find('GrayValue') > -1:
                                    compValue = eval(expr)
                                    OS.write('%s' % repr(compValue))
                                else:
                                    OS.write('%s' % expr)
                                if count < len(exprList) - 1:
                                    OS.write(' ')
                                count += 1

                            OS.write('\n')
                        else:
                            OS.write('%s\n' % propertyStr)

            elif line.upper().find('*USER NODE') == 0:
                noid2 = 0
                noid = 0
                dpUtils.progressStart('     -> Write FEM Node ID   : ')
                for k in range(nz + 1):
                    progress = float(k + 1) / float(nz + 1) * 10.0
                    for j in range(ny + 1):
                        for i in range(nx + 1):
                            noid = noid + 1
                            if activeNodes.has_key(noid):
                                noid2 = noid2 + 1
                                nodeIdMap[noid] = noid2
                                OS.write('%13.6g %13.6g %13.6g\n' % (nodeCoord[noid][0], nodeCoord[noid][1], nodeCoord[noid][2]))

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
                stdout.write('     -> write Nodes         : %10i \n' % noid2)
                stdout.flush()
            elif line.upper().find('*USER ELEMENT') == 0:
                count = 0
                dpUtils.progressStart('     -> process Elements    : ')
                for matid in thresList:
                    count += 1
                    progress = count / float(len(thresList)) * 10.0
                    if len(elsetNodes[repr(matid)]) > 0:
                        for elnds in elsetNodes[repr(matid)]:
                            OS.write('%s %s %s %s %s %s %s %s\n' % (nodeIdMap[elnds[1][0]],
                             nodeIdMap[elnds[1][4]],
                             nodeIdMap[elnds[1][5]],
                             nodeIdMap[elnds[1][1]],
                             nodeIdMap[elnds[1][3]],
                             nodeIdMap[elnds[1][7]],
                             nodeIdMap[elnds[1][6]],
                             nodeIdMap[elnds[1][2]]))

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
                stdout.write('     -> write Elements      : %10i             \n' % elid)
                stdout.flush()
                if nMat > 1:
                    count = 0
                    dpUtils.progressStart('     -> write Material IDs  : ')
                    for matid in thresList:
                        count += 1
                        progress = count / float(len(thresList)) * 10.0
                        if len(elsetNodes[repr(matid)]) > 0:
                            for elnds in elsetNodes[repr(matid)]:
                                OS.write('%i\n' % matid)

                        stdout.write('     -> write Material IDs  : %10i             \n' % nElems)
                        dpUtils.progressNext(progress)

                    dpUtils.progressEnd()
            elif line.upper().find('*USER BC') == 0:
                OS.write('%i\n' % len(bcDict[line]))
                for entry in bcDict[line]:
                    if nodeIdMap.has_key(entry[0]):
                        OS.write('%i %s\n' % (nodeIdMap[entry[0]], entry[1]))
                    else:
                        stdout.write('\n **ERROR** mic.writeParfe(): %s: Node ID %i not found!\n\n' % (line, entry[0]))
                        stdout.flush()
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)

            else:
                OS.write('%s\n' % line)

        osTempFile.close()
        time2 = time.clock()
        stdout.write('     -> Write finished in   :   %8.1f sec  \n' % (time2 - time1))
        stdout.flush()
        OS.close
        return

    def writeOOFEM(self, outFileName, curVoxelModel, dimList, templateFile):
        """
        General Abaqus *.inp file writer. For these materials a default material will be 
        applied. Supported commands: 
          *USER FILENAME
          *USER SIZE RECORD
          *USER NODE
          *USER ELEMENT 
          *USER PROPERTY, file=property_temp.in, range=5:367
            generate multiple material cards, internal variables are "SetName, CardName, GrayValue"
            for the given example: GrayValues > 5 and GrayValues <= 367 are written 
            This card can be used multiple times
            If range=... is not given, material cards for all GrayValues are written
        Elements are only written for the given ranges in *USER PROPERTY
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param  templateFile: name of the template file           
        
        @return:
          no return value
        """
        stdout.write(' ... setup OOFEM *.in file from template\n')
        stdout.write("     -> recast model from '%s' to 'i'\n" % curVoxelModel.dtype.char)
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'i')
        time1 = time.clock()
        if dimList == None:
            print '\n **ERROR** writeOOFEM(): Voxel size not optional for this function!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        xvox = dimList[0]
        yvox = dimList[1]
        zvox = dimList[2]
        nx, ny, nz = self.get_Shape(curVoxelModel)
        minVox, maxVox = self.computeNumpyMinMax(curVoxelModel, 0)
        minVox = int(minVox + 0.5)
        maxVox = int(maxVox + 0.5)
        noid = 0
        activeNodes = {}
        nodeIdName = {}
        allNodeSetsName = {}
        nodeIdName['T'] = 1
        nodeIdName['B'] = 2
        nodeIdName['N'] = 3
        nodeIdName['S'] = 4
        nodeIdName['E'] = 5
        nodeIdName['W'] = 6
        nodeIdName['NT'] = 7
        nodeIdName['ST'] = 8
        nodeIdName['ET'] = 9
        nodeIdName['WT'] = 10
        nodeIdName['NB'] = 11
        nodeIdName['SB'] = 12
        nodeIdName['EB'] = 13
        nodeIdName['WB'] = 14
        nodeIdName['NE'] = 15
        nodeIdName['NW'] = 16
        nodeIdName['SE'] = 17
        nodeIdName['SW'] = 18
        nodeIdName['NET'] = 19
        nodeIdName['NWT'] = 20
        nodeIdName['SET'] = 21
        nodeIdName['SWT'] = 22
        nodeIdName['NEB'] = 23
        nodeIdName['NWB'] = 24
        nodeIdName['SEB'] = 25
        nodeIdName['SWB'] = 26
        tempflag = False
        if templateFile == None:
            tempflag = True
            OS = open('temp.in', 'w')
            OS.write('*USER NODE\n*USER ELEMENT\n*USER PROPERTY, file=prop.in\n')
            templateFile = 'temp.inp'
            OS.close()
            OS = open('prop.inp', 'w')
            OS.write('IsoLE GrayValue d 1. E 5400.*(GrayValue/250.)**2.5 n 0.2\n')
            OS.close()
            templateFile = 'temp.in'
        OS = open(outFileName, 'w')
        try:
            osTempFile = open(templateFile, 'r')
        except IOError:
            stdout.write("\n **ERROR** mic.writeOOFEM(): Abaqus Template file '%s' not found!\n\n" % templateFile)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        lines = osTempFile.readlines()
        elsetNodes = {}
        ranges = {}
        thresList = []
        rangeMin = minVox
        rangeMax = maxVox
        for line in lines:
            line = line.replace('\n', '')
            if line.upper().find('*USER PROPERTY') == 0:
                line = line.replace(' ', '')
                args = line.split(',')
                matTemplateFilename = ''
                for arg in args:
                    if arg.upper().find('RANGE') == 0:
                        dummy, rangeStr = arg.split('=')
                        rangeMin, rangeMax = dpUtils.userSplit(rangeStr)
                        rangeMin = int(rangeMin)
                        rangeMax = int(rangeMax)
                    if arg.upper().find('FILE') == 0:
                        dummy, matTemplateFilename = arg.split('=')

                ranges[matTemplateFilename] = (
                 rangeMin, rangeMax)

        for crange in ranges:
            for matId in range(ranges[crange][0] + 1, ranges[crange][1] + 1):
                elsetNodes[repr(matId)] = []
                thresList.append(matId)

        elid = 0
        nx1 = nx + 1
        nxy1 = (ny + 1) * (nx + 1)
        nElems = 0
        dpUtils.progressStart('     -> setup Element Data  : ')
        for k in range(nz):
            progress = float(k + 1) / float(nz) * 10.0
            for j in range(ny):
                for i in range(nx):
                    grayValue = curVoxelModel[k, j, i]
                    if elsetNodes.has_key(repr(grayValue)):
                        elid += 1
                        nElems += 1
                        elnds = [nxy1 * k + nx1 * j + (i + 1),
                         nxy1 * k + nx1 * j + (i + 2),
                         nxy1 * k + nx1 * (j + 1) + (i + 2),
                         nxy1 * k + nx1 * (j + 1) + (i + 1),
                         nxy1 * (k + 1) + nx1 * j + (i + 1),
                         nxy1 * (k + 1) + nx1 * j + (i + 2),
                         nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2),
                         nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)]
                        elsetNodes[repr(grayValue)].append((elid, elnds))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        nMatCards = 0
        dictGrayMat = {}
        for elsets in elsetNodes:
            if len(elsetNodes[elsets]) > 0:
                nMatCards += 1
                dictGrayMat[elsets] = nMatCards

        stdout.write('     -> Setup Element Data  : Done                    \n')
        stdout.write('     -> Setup Node Data     :')
        for matid in thresList:
            if len(elsetNodes[repr(matid)]) > 0:
                for elnds in elsetNodes[repr(matid)]:
                    elid = elnds[0]
                    for elnd in elnds[1]:
                        activeNodes[elnd] = 1

        noid = 0
        noid2 = 0
        noidMicOOFEM = {}
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    noid = noid + 1
                    if activeNodes.has_key(noid):
                        noid2 = noid2 + 1
                        noidMicOOFEM[noid] = noid2
                        foundCorner = False
                        if i == 0 and j == 0 and k == 0:
                            foundCorner = True
                            allNodeSetsName[noid] = 'SWB'
                        if i == nx and j == 0 and k == 0:
                            foundCorner = True
                            allNodeSetsName[noid] = 'SEB'
                        if i == nx and j == ny and k == 0:
                            foundCorner = True
                            allNodeSetsName[noid] = 'NEB'
                        if i == 0 and j == ny and k == 0:
                            foundCorner = True
                            allNodeSetsName[noid] = 'NWB'
                        if i == 0 and j == 0 and k == nz:
                            foundCorner = True
                            allNodeSetsName[noid] = 'SWT'
                        if i == nx and j == 0 and k == nz:
                            foundCorner = True
                            allNodeSetsName[noid] = 'SET'
                        if i == nx and j == ny and k == nz:
                            foundCorner = True
                            allNodeSetsName[noid] = 'NET'
                        if i == 0 and j == ny and k == nz:
                            foundCorner = True
                            allNodeSetsName[noid] = 'NWT'
                        foundEdge = False
                        if foundCorner == False:
                            if i == 0 and j == 0:
                                foundEdge = True
                                allNodeSetsName[noid] = 'SW'
                            if i == nx and j == 0:
                                foundEdge = True
                                allNodeSetsName[noid] = 'SE'
                            if i == nx and j == ny:
                                foundEdge = True
                                allNodeSetsName[noid] = 'NE'
                            if i == 0 and j == ny:
                                foundEdge = True
                                allNodeSetsName[noid] = 'NW'
                            if i == 0 and k == 0:
                                foundEdge = True
                                allNodeSetsName[noid] = 'WB'
                            if i == nx and k == 0:
                                foundEdge = True
                                allNodeSetsName[noid] = 'EB'
                            if j == 0 and k == 0:
                                foundEdge = True
                                allNodeSetsName[noid] = 'SB'
                            if j == ny and k == 0:
                                foundEdge = True
                                allNodeSetsName[noid] = 'NB'
                            if i == 0 and k == nz:
                                foundEdge = True
                                allNodeSetsName[noid] = 'WT'
                            if i == nx and k == nz:
                                foundEdge = True
                                allNodeSetsName[noid] = 'ET'
                            if j == 0 and k == nz:
                                foundEdge = True
                                allNodeSetsName[noid] = 'ST'
                            if j == ny and k == nz:
                                foundEdge = True
                                allNodeSetsName[noid] = 'NT'
                        if foundCorner == False and foundEdge == False:
                            if k == 0:
                                allNodeSetsName[noid] = 'B'
                            if k == nz:
                                allNodeSetsName[noid] = 'T'
                            if j == 0:
                                allNodeSetsName[noid] = 'S'
                            if j == ny:
                                allNodeSetsName[noid] = 'N'
                            if i == 0:
                                allNodeSetsName[noid] = 'W'
                            if i == nx:
                                allNodeSetsName[noid] = 'E'

        stdout.write(' Done\n')
        stdout.flush()
        nodeCoord = {}
        noid = 0
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    noid = noid + 1
                    if activeNodes.has_key(noid):
                        nodeCoord[noid] = (
                         float(xvox * i), float(yvox * j), float(zvox * k))

        nbc = 0
        for line in lines:
            line = line.replace('\n', '')
            if line.find('BoundaryCondition ') == 0:
                nbc += 1

        stdout.write(' ... write OOFEM *.in file from template\n')
        noShift = 26
        for line in lines:
            line = line.replace('\n', '')
            if line.upper().find('*USER OUTFILE') > -1:
                outList = outFileName.split('.')
                oofemoutfile = outFileName.replace(outList.pop(), 'out')
                shortfile, ext = self.getShortFilenameAndExtension(oofemoutfile)
                OS.write('%s\n' % (shortfile + '.' + ext))
            elif line.upper().find('*USER SIZE RECORD') > -1:
                OS.write('## USER SIZE RECORD\n')
                OS.write('ndofman %i nelem %i ncrosssect 1 nmat %i nbc %i nic 0 nltf 1\n' % (len(activeNodes) + noShift, nElems, nMatCards, nbc))
            elif line.upper().find('*USER NODE') > -1:
                OS.write('## USER NODES\n')
                noid2 = 0
                noid = 0
                dpUtils.progressStart('     -> Write FEM Node ID   : ')
                for k in range(nz + 1):
                    progress = (k + 1) / float(nz + 1) * 10.0
                    for j in range(ny + 1):
                        for i in range(nx + 1):
                            noid = noid + 1
                            if activeNodes.has_key(noid):
                                noid2 = noid2 + 1
                                OS.write('Node %i coords 3 %g %g %g' % (noidMicOOFEM[noid] + noShift, nodeCoord[noid][0], nodeCoord[noid][1], nodeCoord[noid][2]))
                                noidMicOOFEM
                                if allNodeSetsName.has_key(noid):
                                    nsetId = nodeIdName[allNodeSetsName[noid]]
                                    OS.write(' dofType 3 0 0 0 masterMask 3 %i %i %i\n' % (nsetId, nsetId, nsetId))
                                else:
                                    OS.write('\n')

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
                stdout.write('     -> Write FEM Nodes     : %10i     \n' % noid2)
                stdout.flush()
            elif line.upper().find('*USER ELEMENT') > -1:
                elid = 0
                OS.write('## USER ELEMENT\n')
                count = 0
                dpUtils.progressStart('     -> write FEM Elements  : ')
                for matid in thresList:
                    count += 1
                    progress = count / float(len(thresList)) * 10.0
                    if len(elsetNodes[repr(matid)]) > 0:
                        for elnds in elsetNodes[repr(matid)]:
                            elid += 1
                            OS.write('lspace %i nodes 8 %i %i %i %i %i %i %i %i mat %i crossSect 1\n' % (
                             elid,
                             noidMicOOFEM[elnds[1][7]] + noShift,
                             noidMicOOFEM[elnds[1][6]] + noShift,
                             noidMicOOFEM[elnds[1][5]] + noShift,
                             noidMicOOFEM[elnds[1][4]] + noShift,
                             noidMicOOFEM[elnds[1][3]] + noShift,
                             noidMicOOFEM[elnds[1][2]] + noShift,
                             noidMicOOFEM[elnds[1][1]] + noShift,
                             noidMicOOFEM[elnds[1][0]] + noShift,
                             dictGrayMat[repr(matid)]))

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
                stdout.write('     -> write FEM Elements  : %10i             \n' % elid)
                stdout.flush()
            elif line.upper().find('*USER PROPERTY') > -1:
                OS.write('## USER PROPERTY\n')
                line = line.replace(' ', '')
                args = line.split(',')
                rangeMin = minVox
                rangeMax = maxVox
                matTemplateFilename = ''
                for arg in args:
                    if arg.upper().find('RANGE') == 0:
                        dummy, rangeStr = arg.split('=')
                        rangeMin, rangeMax = dpUtils.userSplit(rangeStr)
                        rangeMin = int(rangeMin)
                        rangeMax = int(rangeMax)
                    if arg.upper().find('FILE') == 0:
                        dummy, matTemplateFilename = arg.split('=')

                stdout.write('     -> Write Property      : %s \n' % matTemplateFilename)
                try:
                    osMatCard = open(matTemplateFilename, 'r')
                except IOError:
                    stdout.write("\n **ERROR** writeOOFEM(): Material template file '%s' not found!\n\n" % matTemplateFilename)
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

                lines = osMatCard.readlines()
                for matid in thresList:
                    if len(elsetNodes[repr(matid)]) > 0:
                        if matid >= rangeMin and matid <= rangeMax:
                            GrayValue = matid
                            for line in lines:
                                line = line.replace('\n', '')
                                if line.find('GrayValue') > -1:
                                    exprList = line.split(' ')
                                    count = 0
                                    for expr in exprList:
                                        if expr.find('CardName') > -1:
                                            expr = expr.replace('CardName', repr(dictGrayMat[repr(matid)]))
                                        if expr.find('GrayValue') > -1:
                                            compValue = eval(expr)
                                            OS.write('%s' % repr(compValue))
                                        else:
                                            OS.write('%s' % expr)
                                        if count < len(exprList) - 1:
                                            OS.write(' ')
                                        count += 1

                                    OS.write('\n')
                                else:
                                    OS.write('%s\n' % line)

                osMatCard.close()
            else:
                OS.write('%s\n' % line)

        osTempFile.close()
        if tempflag:
            os.remove('rm temp.in')
            os.remove('rm prop.in')
        time2 = time.clock()
        stdout.write('     -> Write finished in   :   %8.1f sec  \n' % (time2 - time1))
        stdout.flush()
        OS.close
        return

    def writeOptistruct(self, outFileName, curVoxelModel, thresList, dimList):
        """
        Altair Optistruct file writer. Writes HyperMesh Optistruct *.fem and *.sh file.
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  thresList: Threshold list
            - TYPE: list[treshID] = value
            - int thresID ... threshold ID, 0,1,2 .....
            - int value   ... threshold value for this ID, 0...255
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        @return:
          no return value
        """
        time1 = time.clock()
        stdout.write(' ... write HyperMesh Optistruct *.fem *.sh file\n')
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'B')
        if thresList == None:
            stdout.write('\n **ERROR** writeOptistruct(): Threshold not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dimList == None:
            stdout.write('\n **ERROR** writeOptistruct(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        fileName, ext = self.getFilenameAndExtension(outFileName)
        os = open(fileName + '.fem', 'w')
        ossh = open(fileName + '.sh', 'w')
        os.write('$$\n')
        os.write('$$ Optistruct Input Deck Generated by : mic.py for HyperMesh Version 7.0\n')
        os.write('$$ Dr D. H. Pahr, ILSB, TU-Wien, 2005\n')
        os.write('FORMAT HM\n')
        os.write('FORMAT H3D\n')
        os.write('SCREEN OUT\n')
        os.write('$$------------------------------------------------------------------------------$\n')
        os.write('$$                      Case Control Cards                                      $\n')
        os.write('$$------------------------------------------------------------------------------$\n')
        os.write('$$\n')
        os.write('BEGIN BULK\n')
        os.write('$\n')
        os.write('$  GRID Data\n')
        os.write('$\n')
        noid = 0
        elid = 0
        sumid = 0
        coord = ' '
        nx, ny, nz = self.get_Shape(curVoxelModel)
        xvox = dimList[0]
        yvox = dimList[1]
        zvox = dimList[2]
        dpUtils.progressStart('     -> write FEM Nodes     : ')
        for k in range(nz + 1):
            progress = float(noid) / float((nx + 1) * (ny + 1) * (nz + 1)) * 10.0
            for j in range(ny + 1):
                for i in range(nx + 1):
                    noid = noid + 1
                    os.write('%8s%8i%8s%8.1f%8.1f%8.1f\n' % ('GRID    ', noid, coord, float(xvox * i), float(yvox * j), float(zvox * k)))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        stdout.write('     -> write FEM Nodes     : %10i     \n' % noid)
        stdout.flush()
        os.write('$\n')
        os.write('$  CHEXA Elements: First Order\n')
        os.write('$\n')
        iterations = 1
        ossh.write(' 7.0\n')
        ossh.write('%10i%13i\n' % (nx * ny * nz, iterations))
        artifical_mat = 1
        sumid = 0
        dpUtils.progressStart('     -> write FEM Elements  : ')
        for k in range(nz):
            progress = float(sumid) / float(nx * ny * nz) * 10.0
            for j in range(ny):
                for i in range(nx):
                    density = curVoxelModel[k, j, i] / 255.0
                    sumid = sumid + 1
                    node1 = (ny + 1) * (nx + 1) * k + (nx + 1) * j + (i + 1)
                    node2 = (ny + 1) * (nx + 1) * k + (nx + 1) * j + (i + 2)
                    node3 = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + (i + 2)
                    node4 = (ny + 1) * (nx + 1) * k + (nx + 1) * (j + 1) + (i + 1)
                    node5 = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + (i + 1)
                    node6 = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * j + (i + 2)
                    node7 = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + (i + 2)
                    node8 = (ny + 1) * (nx + 1) * (k + 1) + (nx + 1) * (j + 1) + (i + 1)
                    os.write('%8s%8i%8i%8i%8i%8i%8i%8i%8i\n' % ('CHEXA   ', sumid, artifical_mat, node1, node2, node3, node4, node5, node6))
                    os.write('%8s%8i%8i\n' % ('+       ', node7, node8))
                    ossh.write('%10i%13.6e\n' % (sumid, density))

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        stdout.write('     -> write FEM Elements  : %10i     \n' % elid)
        stdout.flush()
        os.write('ENDDATA \n')
        time2 = time.clock()
        stdout.write('     -> write finished in   :   %8.1f sec  \n' % (time2 - time1))
        stdout.flush()
        os.close
        ossh.close
        return

    def writeXml(self, outfilename, curVoxelModel, dimList, rawFormat, history):
        """
        Julius Image file writer. Two files are written. One xml file and one
        data file. 
        
        @param outfilename: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z   
        @param rawFormat: Binary format specifier.
            - TYPE: char
        @param  history: history of the performed modifications on the model
            - TYPE: dict[key] = values
            - int key    ... identifier of the applied option (e.g. -in, -thres)
            - int values ... values for this identifier
          
        @return:
          no return value
        """
        stdout.write(' ... write xml files\n')
        stdout.flush()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        time1 = time.clock()
        if outfilename.upper().find('.RAW.XML') == -1:
            stdout.write('\n **ERROR**: mic.writeXml() a filename of the form "name.raw.xml" is required (1)! \n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        xmlpathRawFile = os.path.realpath(outfilename)
        path, RawXmlFile = os.path.split(xmlpathRawFile)
        filename, ext, xmlext = RawXmlFile.split('.')
        rawfilename = path + '/' + filename + '.raw'
        upext = ext.upper()
        if upext != 'RAW':
            stdout.write('\n **ERROR**: mic.writeXml() a filename of the form "name.raw.xml" is required (2)! \n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if filename.find('<') > -1 or filename.find('>') > -1 or filename.find('&') > -1 or filename.find('"') > -1 or filename.find("'") > -1:
            stdout.write('\n **ERROR**: mic.writeXml() Following characters cannot be used in normal XML strings')
            stdout.write('\n            <,>,&,",\' use &lt,&gt,&quot,&#39 instead!\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if rawFormat == None:
            rawFormat = curVoxelModel.dtype.char
        if rawFormat == 'B':
            format = 'uchar'
            bits = '8'
        elif rawFormat == 'h':
            format = 'short'
            bits = '16'
        elif rawFormat == 'i':
            format = 'integer'
            bits = '32'
        elif rawFormat == 'f':
            format = 'float'
            bits = '32'
        else:
            stdout.write("\n **ERROR** mic.writeXml(): '-raw' ... Format '%s' not supported!\n\n" % rawFormat)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dimList == None:
            stdout.write("\n **ERROR**: mic.writeXml() '-ldim' option not optional for writing RAW/XML data! \n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lx = float(dimList[0])
        ly = float(dimList[1])
        lz = float(dimList[2])
        if history != None:
            shisto = "\n  <modifications date='" + history['Session'] + "'>"
            for item in history:
                nitem = item.replace('-', '')
                if item != 'Session':
                    shisto += '\n   <' + nitem + '>' + history[item] + '</' + nitem + '>'

            shisto += '\n  </modifications>'
        else:
            shisto = ''
        string = '<jxml>\n <patient showroottypenameseparately="true" >\n  <initials>S.B.S.</initials>\n  <dateofbirth>\n   <day>1</day>\n   <month>1</month>\n   <year>1998</year>\n  </dateofbirth>\n  <comments>None</comments>\n  <modality>NA</modality>\n </patient>\n <volume showroottypenameseparately="true" >\n  <name>Schwammkopf</name>\n  <description>data from TU Wien</description>\n  <dimensions textchildnodesononeline="true" >\n   <x>' + repr(nx) + '</x>\n   <y>' + repr(ny) + '</y>\n   <z>' + repr(nz) + '</z>\n  </dimensions>\n  <bigendian>no</bigendian>\n  <numscalar>1</numscalar>\n  <volumeasslicefiles>no</volumeasslicefiles>\n  <spacing textchildnodesononeline="true" >\n   <x>' + repr(lx) + '</x>\n   <y>' + repr(ly) + '</y>\n   <z>' + repr(lz) + '</z>\n  </spacing>\n  <scalartype>' + format + '</scalartype>\n  <bitsused>' + bits + '</bitsused>\n  <headersize>0</headersize>\n  <volumeinonefile_info>\n   <rawfilepath>' + '.' + '</rawfilepath>\n   <rawfilename>' + filename + '.raw' + '</rawfilename>\n  </volumeinonefile_info>\n  <volumeinslicefiles_info>\n   <slicefilepath/>\n   <slicefilepattern/>\n   <numberofslicefiles/>\n  </volumeinslicefiles_info>' + shisto + '\n </volume>\n</jxml>'
        xmlDoc = xml.dom.minidom.parseString(string)
        outs = open(rawfilename + '.xml', 'w')
        outs.write(xmlDoc.toxml())
        outs.close()
        outs = open(rawfilename + '.xml', 'r')
        lines = outs.readlines()
        outs.close()
        outs = open(rawfilename + '.xml', 'w')
        i = 0
        for line in lines:
            if i == 0:
                outs.write('<!DOCTYPE JuliusXML>\n')
            else:
                outs.write('%s' % line)
            i = i + 1

        outs.close()
        self.writeGeneralRaw3(rawfilename, curVoxelModel, rawFormat=rawFormat, bigEndian=False)
        return

    def writeTxt(self, outfilename, curVoxelModel, dimList, echo=True):
        if echo:
            stdout.write(' ... write txt files\n')
            stdout.flush()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        if dimList == None:
            dpUtils.throwError("mic.writeTxt(): '-ldim' option not optional for writing text data!")
        lx = float(dimList[0])
        ly = float(dimList[1])
        lz = float(dimList[2])
        if not (numpy.allclose(lx, ly) and numpy.allclose(lx, lz) and numpy.allclose(ly, lz)):
            dpUtils.throwError('mic.writeTxt(): Only isotropic data are allowed !')
        outs = open(outfilename, 'w')
        outs.write('%i\n' % (nx * ny * nz))
        outs.write('%g %g %g\n' % (nx, ny, nz))
        outs.write('%g\n' % lx)
        outs.write('%s\n' % outfilename)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    outs.write('%g\n' % curVoxelModel[k, j, i])

        outs.close()
        return

    def writeMhd(self, outfilename, curVoxelModel, dimList, rawFormat, history, geomList=None, echo=True):
        """
              Meta Image file writer. Two files are written. One mhd file and one
              data file. 
        
              @param outfilename: name of the output file
              @param curVoxelModel: voxel model of the RVE
                  - TYPE: numpy.array[iZ, jY, kX] = grayValue
                  - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
                  - int grayValue  ... value of voxel
              @param  dimList: list of voxel dimension's
                  - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
                  - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z   
              @param rawFormat: Binary format specifier.
                  - TYPE: char
              @param  history: history of the performed modifications on the model
                  - TYPE: dict[key] = values
                  - int key    ... identifier of the applied option (e.g. -in, -thres)
                  - int values ... values for this identifier
              outs = open(rawfilename+'.xml','r')
              lines = outs.readlines()
              outs.close()
              outs = open(rawfilename+'.xml','w')
              i = 0 
              for line in lines :
                if i == 0 :
                  outs.write('<!DOCTYPE JuliusXML>
        ')
                else:
                  outs.write('%s' % line)
                i=i+1
              outs.close()
        
              @return:
                no return value
              """
        if echo:
            stdout.write(' ... write mhd files\n')
            stdout.flush()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        time1 = time.clock()
        if outfilename.upper().find('.MHD') == -1:
            stdout.write('\n **ERROR**: mic.writeMhd() a filename with *.mhd extension is required (1)! \n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        mhdPathRawFile = os.path.realpath(outfilename)
        path, RawMhaFile = os.path.split(mhdPathRawFile)
        filename, ext = RawMhaFile.split('.')
        rawfilename = path + '/' + filename
        if rawFormat == None:
            rawFormat = curVoxelModel.dtype.char
        if rawFormat == 'B':
            format = 'MET_UCHAR'
        elif rawFormat == 'h':
            format = 'MET_SHORT'
        elif rawFormat == 'H':
            format = 'MET_USHORT'
        elif rawFormat == 'i':
            format = 'MET_INT'
        elif rawFormat == 'f':
            format = 'MET_FLOAT'
        else:
            stdout.write("\n **ERROR** mic.writeMhd(): '-raw' ... Format '%s' not supported!\n\n" % rawFormat)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dimList == None:
            stdout.write("\n **ERROR**: mic.writeMhd() '-ldim' option not optional for writing Mha data! \n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lx = float(dimList[0])
        ly = float(dimList[1])
        lz = float(dimList[2])
        if geomList != None:
            TM = geomList['TransformMatrix']
            TransformMatrix = repr(TM[0]) + ' ' + repr(TM[1]) + ' ' + repr(TM[2]) + ' ' + repr(TM[3]) + ' ' + repr(TM[4]) + ' ' + repr(TM[5]) + ' ' + repr(TM[6]) + ' ' + repr(TM[7]) + ' ' + repr(TM[8])
            OF = geomList['Offset']
            Offset = repr(OF[0]) + ' ' + repr(OF[1]) + ' ' + repr(OF[2])
            CR = geomList['CenterOfRotation']
            CenterOfRotation = repr(CR[0]) + ' ' + repr(CR[1]) + ' ' + repr(CR[2])
            AnatomicalOrientation = geomList['AnatomicalOrientation']
        else:
            TransformMatrix = '1 0 0 0 1 0 0 0 1'
            Offset = '0 0 0'
            CenterOfRotation = '0 0 0'
            AnatomicalOrientation = 'LPS'
        outs = open(rawfilename + '.mhd', 'w')
        outs.write('ObjectType = Image\n')
        outs.write('NDims = 3\n')
        outs.write('BinaryData = True\n')
        outs.write('BinaryDataByteOrderMSB = False\n')
        outs.write('CompressedData = False\n')
        outs.write('TransformMatrix = %s \n' % TransformMatrix)
        outs.write('Offset = %s \n' % Offset)
        outs.write('CenterOfRotation = %s \n' % CenterOfRotation)
        outs.write('AnatomicalOrientation = %s \n' % AnatomicalOrientation)
        outs.write('ElementSpacing = %g %g %g\n' % (lx, ly, lz))
        outs.write('DimSize = %i %i %i\n' % (nx, ny, nz))
        outs.write('ElementType = %s\n' % format)
        outs.write('ElementDataFile = %s\n' % (filename + '.raw'))
        outs.close()
        self.writeGeneralRaw3(rawfilename + '.raw', curVoxelModel, rawFormat=rawFormat, bigEndian=False, echo=echo)
        return

    def writeNhdr(self, outfilename, curVoxelModel, dimList, rawFormat, history, geomList=None, echo=True):
        """
              NRRD Image file writer. Two files are written. One header file and one
              data file. 
        
              @param outfilename: name of the output file
              @param curVoxelModel: voxel model of the RVE
                  - TYPE: numpy.array[iZ, jY, kX] = grayValue
                  - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
                  - int grayValue  ... value of voxel
              @param  dimList: list of voxel dimension's
                  - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
                  - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z   
              @param rawFormat: Binary format specifier.
                  - TYPE: char
              @param  history: history of the performed modifications on the model
                  - TYPE: dict[key] = values
                  - int key    ... identifier of the applied option (e.g. -in, -thres)
                  - int values ... values for this identifier
              outs = open(rawfilename+'.xml','r')
              lines = outs.readlines()
              outs.close()
              outs = open(rawfilename+'.xml','w')
              i = 0 
              for line in lines :
                if i == 0 :
                  outs.write('<!DOCTYPE JuliusXML>
        ')
                else:
                  outs.write('%s' % line)
                i=i+1
              outs.close()
        
              @return:
                no return value
              """
        if echo:
            stdout.write(' ... write nhdr files\n')
            stdout.flush()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        time1 = time.clock()
        if outfilename.upper().find('.NHDR') == -1:
            stdout.write('\n **ERROR**: mic.writeNhdr() a filename with *.nhdr extension is required (1)! \n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nhdrPathRawFile = os.path.realpath(outfilename)
        path, RawNhdrFile = os.path.split(nhdrPathRawFile)
        filename, ext = RawNhdrFile.split('.')
        rawfilename = path + '/' + filename
        if rawFormat == None:
            rawFormat = curVoxelModel.dtype.char
        if rawFormat == 'B':
            format = 'unsigned char'
        elif rawFormat == 'h':
            format = 'short'
        elif rawFormat == 'H':
            format = 'unsigned short'
        elif rawFormat == 'i':
            format = 'int'
        elif rawFormat == 'f':
            format = 'float'
        else:
            stdout.write("\n **ERROR** mic.writeNhdr(): '-raw' ... Format '%s' not supported!\n\n" % rawFormat)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dimList == None:
            stdout.write("\n **ERROR**: mic.writeNhdr(): '-ldim' option not optional for writing Nhdr data! \n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lx = float(dimList[0])
        ly = float(dimList[1])
        lz = float(dimList[2])
        if geomList != None:
            TM = geomList['TransformMatrix']
            TransformMatrix = '(' + repr(TM[0] * lx) + ',' + repr(TM[1] * lx) + ',' + repr(TM[2] * lx) + ') (' + repr(TM[3] * ly) + ',' + repr(TM[4] * ly) + ',' + repr(TM[5] * ly) + ') (' + repr(TM[6] * lz) + ',' + repr(TM[7] * lz) + ',' + repr(TM[8] * lz) + ')'
            OF = geomList['Offset']
            Offset = '(' + repr(OF[0]) + ',' + repr(OF[1]) + ',' + repr(OF[2]) + ')'
            AO = geomList['AnatomicalOrientation'].replace(' ', '')
            CR = geomList['CenterOfRotation']
            if not numpy.allclose(CR, [0.0, 0.0, 0.0], rtol=1e-05, atol=1e-08):
                stdout.write('\n **ERROR** mic.writeNhdr(): center of rotation have to be at [0 0 0]!\n\n')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if AO == 'LPS' or AO == 'LAS' or AO == 'RAS':
                AnatomicalOrientation = AO
            else:
                stdout.write("\n **ERROR** mic.writeNhdr(): space '%s' not supported!\n\n" % AO)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            TransformMatrix = '(' + repr(lx) + ', 0, 0)  (0,' + repr(ly) + ', 0) (0, 0,' + repr(lz) + ')'
            Offset = '(0, 0, 0)'
            AnatomicalOrientation = 'LPS'
        outs = open(rawfilename + '.nhdr', 'w')
        outs.write('NRRD0005\n')
        outs.write('# Complete NRRD file format specification at:\n')
        outs.write('# http://teem.sourceforge.net/nrrd/format.html\n')
        outs.write('# generated by medtool: www.dr-pahr.at\n')
        outs.write('type: %s\n' % format)
        outs.write('dimension: 3\n')
        outs.write('space: %s\n' % AnatomicalOrientation)
        outs.write('sizes: %i %i %i\n' % (nx, ny, nz))
        outs.write('space directions: %s \n' % TransformMatrix)
        outs.write('kinds: domain domain domain\n')
        outs.write('endian: little\n')
        outs.write('encoding: raw\n')
        outs.write('space origin: %s \n' % Offset)
        outs.write('data file: %s\n' % (filename + '.raw'))
        outs.close()
        self.writeGeneralRaw3(rawfilename + '.raw', curVoxelModel, rawFormat=rawFormat, bigEndian=False, echo=echo)
        return

    def writeSvx(self, outfilename, curVoxelModel, dimList, rawFormat, history, geomList=None, echo=True):
        """
        Shapeways Svx image file writer. 
        """
        if echo:
            stdout.write(' ... write svx files\n')
            stdout.flush()
        time1 = time.clock()
        nx, ny, nz = self.get_Shape(curVoxelModel)
        lx, ly, lz = float(dimList[0]), float(dimList[1]), float(dimList[2])
        if geomList == None:
            ox, oy, oz = (0.0, 0.0, 0.0)
        elif 'Offset' in geomList:
            ox, oy, oz = float(geomList['Offset'][0]), float(geomList['Offset'][1]), float(geomList['Offset'][2])
        else:
            ox, oy, oz = (0.0, 0.0, 0.0)
        if outfilename.upper().find('.SVX') == -1:
            stdout.write('\n **ERROR**: mic.writeSvx() a filename with *.svx extension is required (1)! \n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        path, filenameext = os.path.split(outfilename)
        filename, ext = os.path.splitext(filenameext)
        if path == '':
            path = os.getcwd()
        pathToInstaller = path + os.sep + 'temp'
        os.mkdir(pathToInstaller)
        os.mkdir(pathToInstaller + os.sep + 'density')
        curVoxelModel = self.castType(curVoxelModel, 'B')
        curVoxelModel = self.scaleModel(curVoxelModel, 0, 255, echo=False)
        if len(numpy.unique(curVoxelModel)) > 2:
            dpUtils.throwError('mic.writeSvx(): Only one material (grayvalue) is allowed!')
        self.writeImages(pathToInstaller + os.sep + 'density' + os.sep + 'slice.png', curVoxelModel, echo=echo, zero=True)
        outs = open(pathToInstaller + os.sep + 'manifest.xml', 'w')
        outs.write('<?xml version="1.0"?>\n')
        outs.write('<grid ')
        outs.write('gridSizeX="%i" gridSizeY="%i" gridSizeZ="%i" ' % (nx, ny, nz))
        outs.write('originX="%g" originY="%g" originZ="%g" ' % (ox / 1000.0, oy / 1000.0, oz / 1000.0))
        if not (numpy.allclose(lx, ly) and numpy.allclose(lx, lz) and numpy.allclose(ly, lz)):
            dpUtils.throwError('mic.writeSvx(): Only isotropic data are allowed !')
        outs.write('voxelSize="%g" ' % (lx / 1000.0))
        outs.write('subvoxelBits="8" ')
        outs.write('slicesOrientation="Z" ')
        outs.write('>\n')
        outs.write('\n')
        outs.write('    <channels>\n')
        outs.write('        <channel type="DENSITY" bits="8" slices="density/slice%04d.png" />\n')
        outs.write('    </channels>\n')
        outs.write('\n')
        outs.write('    <metadata>\n')
        outs.write('        <entry key="author" value="medool" />\n')
        outs.write('        <entry key="creationDate" value="%s" />\n' % time.strftime('%Y/%m/%d'))
        outs.write('    </metadata>\n')
        outs.write('\n')
        outs.write('</grid> ')
        outs.close()
        shutil.make_archive(path + os.sep + filename, format='zip', root_dir=pathToInstaller)
        shutil.rmtree(pathToInstaller)
        shutil.move(path + os.sep + filename + '.zip', path + os.sep + filenameext)
        return

    def writeMesh2d(self, outFileName, curVoxelModel, dimList, smooth, eltype):
        """
        Write 2D surface mesh file. 
        
        @param outFileName: name of the output file
        @param curVoxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param  dimList: list of voxel dimension's
            - TYPE: list[0] = lenX, list[0] = lenY, list[0] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        @param  smoothParam: taubin voxel model smoothing parameter 
            - TYPE: list[0] = iter, list[1] = lambda, list[2] = kPB, 
                    list[3] = nearIntf, list[4] = bcid, list[5] = shrink
            - int iter, float lambda, float kPB, int nearIntf      
        @param  eltype: element type  
            - TYPE: string
        @return:
          no return value
        """
        stdout.write(' ... write 2D mesh file \n')
        stdout.write("     -> recast model from '%s' to 'i'\n" % curVoxelModel.dtype.char)
        stdout.flush()
        curVoxelModel = self.castType(curVoxelModel, 'i')
        time1 = time.clock()
        if dimList == None:
            print '\n **ERROR** writeMesh2d(): Voxel size not optional for this function!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if not (eltype == 'tria3' or eltype == 'quad4'):
            print '\n **ERROR** writeMesh2d(): element type ', eltype, ' not implemented!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if numpy.min(curVoxelModel) < 0:
            print '\n **ERROR** writeMesh2d(): voxel values have to be bigger than 0!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        xvox = dimList[0]
        yvox = dimList[1]
        zvox = dimList[2]
        nx, ny, nz = self.get_Shape(curVoxelModel)
        activeNodes = {}
        nx1 = nx + 1
        nxy1 = (ny + 1) * (nx + 1)
        dpUtils.progressStart('     -> Setup Node Data     : ')
        for k in range(nz):
            progress = float(k + 1) / float(nz) * 10.0
            for j in range(ny):
                for i in range(nx):
                    grayValue = curVoxelModel[k, j, i]
                    if grayValue > 0:
                        activeNodes[nxy1 * k + nx1 * j + (i + 1)] = 1
                        activeNodes[nxy1 * k + nx1 * j + (i + 2)] = 1
                        activeNodes[nxy1 * k + nx1 * (j + 1) + (i + 2)] = 1
                        activeNodes[nxy1 * k + nx1 * (j + 1) + (i + 1)] = 1
                        activeNodes[nxy1 * (k + 1) + nx1 * j + (i + 1)] = 1
                        activeNodes[nxy1 * (k + 1) + nx1 * j + (i + 2)] = 1
                        activeNodes[nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2)] = 1
                        activeNodes[nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)] = 1

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        nodeCoord = {}
        if smooth != None:
            activeNodes2, nElem, nNode, nIntElem, nIntFaces, nIntNode = dpMesherf77.check_voxmesh2d(curVoxelModel, smooth[4], 2)
            nodeCoordF77, nodeCoordInt, noidF77PY, noidIntVoxF77 = dpMesherf77.get_voxmesh2d_nodes(curVoxelModel, dimList, smooth, activeNodes2, nNode, nIntElem, nIntNode, 2)
            for noidF77 in range(len(noidF77PY)):
                noid = noidF77PY[noidF77]
                nodeCoord[noid] = (nodeCoordF77[noidF77][0], nodeCoordF77[noidF77][1], nodeCoordF77[noidF77][2])

        else:
            noid = 0
            for k in range(nz + 1):
                for j in range(ny + 1):
                    for i in range(nx + 1):
                        noid = noid + 1
                        if activeNodes.has_key(noid):
                            nodeCoord[noid] = (
                             float(xvox * i), float(yvox * j), float(zvox * k))

        fecModel = fec.fec()
        Nodes = {}
        Elems = {}
        ElemsSets = []
        EResults = {}
        sum = 0
        surfElems = {}
        elid = 0
        dpUtils.progressStart('     -> Setup Surface Elems : ')
        for k in range(-1, nz):
            sum += 1
            progress = float(sum) / float(nz + 1) * 10.0
            for j in range(-1, ny):
                for i in range(-1, nx):
                    if k > -1 and j > -1:
                        bset = ''
                        if i == -1:
                            grayValue = 0
                            bset = 'B'
                        else:
                            grayValue = curVoxelModel[k, j, i]
                        if i == nx - 1:
                            grayValue2 = 0
                            bset = 'B'
                        else:
                            grayValue2 = curVoxelModel[k, j, i + 1]
                        if grayValue2 != grayValue:
                            if grayValue < grayValue2:
                                set = bset + 'SET' + repr(grayValue) + '_' + repr(grayValue2)
                            else:
                                set = bset + 'SET' + repr(grayValue2) + '_' + repr(grayValue)
                            if not surfElems.has_key(set):
                                surfElems[set] = []
                                ElemsSets.append(set)
                            nid1 = nxy1 * k + nx1 * j + (i + 2)
                            nid2 = nxy1 * k + nx1 * (j + 1) + (i + 2)
                            nid3 = nxy1 * (k + 1) + nx1 * j + (i + 2)
                            nid4 = nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2)
                            if not Nodes.has_key(nid1):
                                Nodes[nid1] = dpFem.node(nid1, nodeCoord[nid1][0], nodeCoord[nid1][1], nodeCoord[nid1][2])
                            if not Nodes.has_key(nid2):
                                Nodes[nid2] = dpFem.node(nid2, nodeCoord[nid2][0], nodeCoord[nid2][1], nodeCoord[nid2][2])
                            if not Nodes.has_key(nid3):
                                Nodes[nid3] = dpFem.node(nid3, nodeCoord[nid3][0], nodeCoord[nid3][1], nodeCoord[nid3][2])
                            if not Nodes.has_key(nid4):
                                Nodes[nid4] = dpFem.node(nid4, nodeCoord[nid4][0], nodeCoord[nid4][1], nodeCoord[nid4][2])
                            if eltype == 'tria3':
                                dist14 = numpy.sqrt((nodeCoord[nid1][0] - nodeCoord[nid4][0]) ** 2 + (nodeCoord[nid1][1] - nodeCoord[nid4][1]) ** 2 + (nodeCoord[nid1][2] - nodeCoord[nid4][2]) ** 2)
                                dist23 = numpy.sqrt((nodeCoord[nid2][0] - nodeCoord[nid3][0]) ** 2 + (nodeCoord[nid2][1] - nodeCoord[nid3][1]) ** 2 + (nodeCoord[nid2][2] - nodeCoord[nid3][2]) ** 2)
                                if dist14 < dist23:
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid2, nid4], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid4, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                else:
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid2, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid2, nid4, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                            else:
                                elid += 1
                                Elems[elid] = dpFem.element(elid, [nid1, nid2, nid4, nid3], 'quad4')
                                EResults[elid] = ElemsSets.index(set)
                                surfElems[set].append(elid)
                    if k > -1 and i > -1:
                        bset = ''
                        if j == -1:
                            grayValue = 0
                            bset = 'B'
                        else:
                            grayValue = curVoxelModel[k, j, i]
                        if j == ny - 1:
                            grayValue2 = 0
                            bset = 'B'
                        else:
                            grayValue2 = curVoxelModel[k, j + 1, i]
                        if grayValue2 != grayValue:
                            if grayValue < grayValue2:
                                set = bset + 'SET' + repr(grayValue) + '_' + repr(grayValue2)
                            else:
                                set = bset + 'SET' + repr(grayValue2) + '_' + repr(grayValue)
                            if not surfElems.has_key(set):
                                surfElems[set] = []
                                ElemsSets.append(set)
                            nid1 = nxy1 * k + nx1 * (j + 1) + (i + 2)
                            nid2 = nxy1 * k + nx1 * (j + 1) + (i + 1)
                            nid3 = nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2)
                            nid4 = nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)
                            if not Nodes.has_key(nid1):
                                Nodes[nid1] = dpFem.node(nid1, nodeCoord[nid1][0], nodeCoord[nid1][1], nodeCoord[nid1][2])
                            if not Nodes.has_key(nid2):
                                Nodes[nid2] = dpFem.node(nid2, nodeCoord[nid2][0], nodeCoord[nid2][1], nodeCoord[nid2][2])
                            if not Nodes.has_key(nid3):
                                Nodes[nid3] = dpFem.node(nid3, nodeCoord[nid3][0], nodeCoord[nid3][1], nodeCoord[nid3][2])
                            if not Nodes.has_key(nid4):
                                Nodes[nid4] = dpFem.node(nid4, nodeCoord[nid4][0], nodeCoord[nid4][1], nodeCoord[nid4][2])
                            if eltype == 'tria3':
                                dist14 = numpy.sqrt((nodeCoord[nid1][0] - nodeCoord[nid4][0]) ** 2 + (nodeCoord[nid1][1] - nodeCoord[nid4][1]) ** 2 + (nodeCoord[nid1][2] - nodeCoord[nid4][2]) ** 2)
                                dist23 = numpy.sqrt((nodeCoord[nid2][0] - nodeCoord[nid3][0]) ** 2 + (nodeCoord[nid2][1] - nodeCoord[nid3][1]) ** 2 + (nodeCoord[nid2][2] - nodeCoord[nid3][2]) ** 2)
                                if dist14 < dist23:
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid2, nid4], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid4, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                else:
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid2, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid2, nid4, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                            else:
                                elid += 1
                                Elems[elid] = dpFem.element(elid, [nid1, nid2, nid4, nid3], 'quad4')
                                EResults[elid] = ElemsSets.index(set)
                                surfElems[set].append(elid)
                    if i > -1 and j > -1:
                        bset = ''
                        if k == -1:
                            grayValue = 0
                            bset = 'B'
                        else:
                            grayValue = curVoxelModel[k, j, i]
                        if k == nz - 1:
                            grayValue2 = 0
                            bset = 'B'
                        else:
                            grayValue2 = curVoxelModel[k + 1, j, i]
                        if grayValue2 != grayValue:
                            if grayValue < grayValue2:
                                set = bset + 'SET' + repr(grayValue) + '_' + repr(grayValue2)
                            else:
                                set = bset + 'SET' + repr(grayValue2) + '_' + repr(grayValue)
                            if not surfElems.has_key(set):
                                surfElems[set] = []
                                ElemsSets.append(set)
                            nid1 = nxy1 * (k + 1) + nx1 * j + (i + 1)
                            nid2 = nxy1 * (k + 1) + nx1 * j + (i + 2)
                            nid4 = nxy1 * (k + 1) + nx1 * (j + 1) + (i + 2)
                            nid3 = nxy1 * (k + 1) + nx1 * (j + 1) + (i + 1)
                            if not Nodes.has_key(nid1):
                                Nodes[nid1] = dpFem.node(nid1, nodeCoord[nid1][0], nodeCoord[nid1][1], nodeCoord[nid1][2])
                            if not Nodes.has_key(nid2):
                                Nodes[nid2] = dpFem.node(nid2, nodeCoord[nid2][0], nodeCoord[nid2][1], nodeCoord[nid2][2])
                            if not Nodes.has_key(nid3):
                                Nodes[nid3] = dpFem.node(nid3, nodeCoord[nid3][0], nodeCoord[nid3][1], nodeCoord[nid3][2])
                            if not Nodes.has_key(nid4):
                                Nodes[nid4] = dpFem.node(nid4, nodeCoord[nid4][0], nodeCoord[nid4][1], nodeCoord[nid4][2])
                            if eltype == 'tria3':
                                dist14 = numpy.sqrt((nodeCoord[nid1][0] - nodeCoord[nid4][0]) ** 2 + (nodeCoord[nid1][1] - nodeCoord[nid4][1]) ** 2 + (nodeCoord[nid1][2] - nodeCoord[nid4][2]) ** 2)
                                dist23 = numpy.sqrt((nodeCoord[nid2][0] - nodeCoord[nid3][0]) ** 2 + (nodeCoord[nid2][1] - nodeCoord[nid3][1]) ** 2 + (nodeCoord[nid2][2] - nodeCoord[nid3][2]) ** 2)
                                if dist14 < dist23:
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid2, nid4], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid4, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                else:
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid1, nid2, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                                    elid += 1
                                    Elems[elid] = dpFem.element(elid, [nid2, nid4, nid3], 'tria3')
                                    EResults[elid] = ElemsSets.index(set)
                                    surfElems[set].append(elid)
                            else:
                                elid += 1
                                Elems[elid] = dpFem.element(elid, [nid1, nid2, nid4, nid3], 'quad4')
                                EResults[elid] = ElemsSets.index(set)
                                surfElems[set].append(elid)

            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        fecModel.write(outFileName, None, Nodes, None, Elems, surfElems, EscaResults=[EResults])
        return

    def writeMeshStlBin(self, outFileName, voxelModel, dimList, smooth, fast=True, additData=None):
        """
        Write 2D surface mesh file in stl binary format
        """
        stdout.write(' ... write stl mesh file \n')
        stdout.write("     -> recast model from '%s' to 'i'\n" % voxelModel.dtype.char)
        stdout.flush()
        curVoxelModel = self.castType(voxelModel, 'i')
        time1 = time.clock()
        if dimList == None:
            print '\n **ERROR** writeMeshStlBin(): Voxel size not optional for this function!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if numpy.min(voxelModel) < 0:
            print '\n **ERROR** writeMeshStlBin(): voxel values have to be bigger than 0!\n'
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if smooth == None:
            smooth = [0, 0.6307, 0.1, 0, 1, 0.05, 0.6]
            fast = True
        else:
            fast = False
        nodeCoord, faces, facesSets = dpMesher.computeImageToSmoothMesh2D(voxelModel, dimList, smooth=smooth, fast=fast, echo=1)
        offset = None
        if additData != None:
            if 'Offset' in additData:
                offset = numpy.array(additData['Offset'])
                if numpy.allclose(offset, [0.0, 0.0, 0.0]):
                    offset = None
        outFile = open(outFileName, 'wb')
        for i in range(80):
            outFile.write(struct.pack('B', 0))

        outFile.write(struct.pack('I', len(faces)))
        sumId = 0
        maxId = len(faces)
        dpUtils.progressStart('     -> writing binary : ')
        for face in faces:
            sumId += 1
            progress = float(sumId) / float(maxId) * 10.0
            if offset == None:
                n1 = numpy.array(nodeCoord[face[0]])
                n2 = numpy.array(nodeCoord[face[1]])
                n3 = numpy.array(nodeCoord[face[2]])
            else:
                n1 = numpy.array(nodeCoord[face[0]]) + offset
                n2 = numpy.array(nodeCoord[face[1]]) + offset
                n3 = numpy.array(nodeCoord[face[2]]) + offset
            a = n2 - n1
            b = n3 - n1
            nVec = numpy.array([a[1] * b[2] - a[2] * b[1],
             a[2] * b[0] - a[0] * b[2],
             a[0] * b[1] - a[1] * b[0]])
            outFile.write(struct.pack('12fH', nVec[0], nVec[1], nVec[2], n1[0], n1[1], n1[2], n2[0], n2[1], n2[2], n3[0], n3[1], n3[2], 0))
            dpUtils.progressNext(progress)

        dpUtils.progressEnd()
        outFile.close()
        return

    def writeMidplaneImages(self, outFileName, voxelModel):
        """
        Write three images of the midplane.
        
        @param outFileName: name of the output file
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        
        @return:
          no return value
        """
        stdout.write(' ... write midplane images                       \n')
        stdout.flush()
        minValue, maxValue = self.computeNumpyMinMax(voxelModel, 0)
        delta = float(abs(maxValue - minValue))
        if numpy.allclose(delta, 0.0, rtol=1e-05, atol=1e-08):
            if minValue >= 0 and maxValue <= 255:
                scale = 1.0
            else:
                scale = 0.0
        else:
            scale = 255 / delta
        nx, ny, nz = self.get_Shape(voxelModel)
        filename, ext = self.getFilenameAndExtension(outFileName)
        upext = ext.upper()
        if upext == 'JPG':
            upext = 'JPEG'
        name = filename + '-XM' + '.' + ext
        misc.toimage((voxelModel[:, :, int(nx / 2)] - minValue) * scale, cmin=0, cmax=255).save(name)
        name = filename + '-YM' + '.' + ext
        misc.toimage((voxelModel[:, int(ny / 2), :] - minValue) * scale, cmin=0, cmax=255).save(name)
        name = filename + '-ZM' + '.' + ext
        misc.toimage((voxelModel[int(nz / 2), :, :] - minValue) * scale, cmin=0, cmax=255).save(name)

    def writeSingleImages(self, voxelModel, direct, layer, outFileName):
        """
        Write a single image
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param outFileName: name of the output file
        @param direct: direction 1+,2,+,3+,1-,2-,3-
        @param layer: slice number to be written, started with 1
        @param outFileName: name of the output file
        
        @return:
          no return value
        """
        stdout.write(' ... write single images                       \n')
        stdout.flush()
        minValue, maxValue = self.computeNumpyMinMax(voxelModel, 0)
        delta = float(abs(maxValue - minValue))
        scale = 255 / delta
        nx, ny, nz = self.get_Shape(voxelModel)
        filename, ext = self.getFilenameAndExtension(outFileName)
        upext = ext.upper()
        if upext == 'JPG':
            upext = 'JPEG'
        layer0 = layer - 1
        sliceModel = None
        if direct == '1+' or direct == '1-':
            if layer0 >= 0 and layer0 < nx:
                if direct == '1+':
                    sliceModel = (voxelModel[:, :, layer0] - minValue) * scale
                else:
                    sliceModel = numpy.flipud((voxelModel[:, :, layer0] - minValue) * scale)
            else:
                stdout.write('\n **ERROR** : writeSingleImages(): sliceNo=%s outside voxel model!\n' % repr(layer))
                stdout.write('             Valid range:  1 <= sliceNo <= %i!\n' % nx)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        elif direct == '2+' or direct == '2-':
            if layer0 >= 0 and layer0 < ny:
                if direct == '2+':
                    sliceModel = (voxelModel[:, layer0, :] - minValue) * scale
                else:
                    sliceModel = numpy.flipud((voxelModel[:, layer0, :] - minValue) * scale)
            else:
                stdout.write('\n **ERROR** : writeSingleImages(): sliceNo=%s outside voxel model!\n' % repr(layer))
                stdout.write('             Valid range:  1 <= sliceNo <= %i!\n' % ny)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        elif direct == '3+' or direct == '3-':
            if layer0 >= 0 and layer0 < nz:
                if direct == '3+':
                    sliceModel = (voxelModel[layer0, :, :] - minValue) * scale
                else:
                    sliceModel = numpy.flipud((voxelModel[layer0, :, :] - minValue) * scale)
            else:
                stdout.write('\n **ERROR** : writeSingleImages(): sliceNo=%s outside voxel model!\n' % repr(layer))
                stdout.write('             Valid range:  1 <= sliceNo <= %i!\n' % nz)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        else:
            stdout.write('\n **ERROR** : writeSingleImages(): direction=%s not supported!\n' % repr(direct))
            stdout.write('             Possible directions: 1+, 1-, 2+, 2-, 3+, 3-!\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        name = filename + '_' + direct + '_' + repr(layer) + '.' + ext
        misc.toimage(sliceModel, cmin=0, cmax=255).save(name)
        return

    def writeGeomObject(self, parList):
        """
        Function writes a geometric object to the voxel model
        """
        stdout.write(' ... write geometric object \n')
        stdout.flush()
        gType = parList[0]
        grayValue = float(parList[1])
        phyLen = float(parList[2])
        filename = parList[3]
        geomList = None
        dimList = [
         phyLen, phyLen, phyLen]
        if gType.lower() == 'sph':
            diam = int(parList[4])
            voxelModel = self.make_sphere(diam, grayValue)
        elif gType.lower() == 'cyl':
            diam = int(parList[4])
            height = int(parList[5])
            voxelModel = self.make_zylinder(diam, height, grayValue)
        else:
            stdout.write("\n **ERROR** writeGeomObject(): Object type '%s' not known!\n\n" % gType)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        self.write(filename, voxelModel, dimList=dimList, format='f', geomList=geomList)
        return

    def writeHistogram(self, outFileName, voxelModel, normalize='n', noInterval=None):
        """
        Writes histogram data into an output file. 
        
        @param outFileName: name of the output file
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... value of voxel
        @param normalize : if 'y' or 'Y' it will normalize the histogram to 
                           0-100%
        @param noInterval: number of intervals 
            - int noIntervall
            
        @return:
          no return value
        """
        print ' ... write Histogram '
        norm = False
        if normalize == 'y' or normalize == 'Y':
            norm = True
        nz, ny, nx = voxelModel.shape
        noVox = nz * ny * nx
        histogram = {}
        minVox, maxVox = self.computeNumpyMinMax(voxelModel, 0)
        if noInterval == None:
            noInterval = int(maxVox + 1 - minVox)
        else:
            noInterval = int(noInterval)
        dInt = (maxVox + 1 - minVox) / float(noInterval)
        for i in range(noInterval):
            histogram[i] = 0

        for voxel in voxelModel.flat:
            i = int((voxel - minVox) / dInt)
            histogram[i] += 1

        os = open(outFileName, 'w')
        sum = 0.0
        for i in range(noInterval):
            if norm:
                sum += histogram[i] / float(noVox) * 100.0
                os.write('%12i %13.6g  %13.6g\n' % (i * dInt + minVox,
                 histogram[i] / float(noVox) * 100.0, sum))
            else:
                sum += histogram[i]
                os.write('%12i %12i %13i\n' % (i * dInt + minVox, histogram[i], sum))

        os.close()
        return

    def changeResolution2(self, voxelModel, resolut, ldim=None):
        """
        Function changes the resultion of a given voxel model. IF "ldim" 
        is given the dimension will be updated by reference.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param resolut: New image resolution
            - TYPE: list[0] = resX, list[1] = resY, list[2] = resZ
                list[3] = enhance, list[4] = repeat
            - int resX, resY, resZ ... new number of voxels in x,y,z
        @param ldim: voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        """
        stdout.write(' ... Change Resolution 2\n')
        stdout.flush()
        resX = int(resolut[0])
        resY = int(resolut[1])
        resZ = int(resolut[2])
        nx, ny, nz = self.get_Shape(voxelModel)
        resConnectDictX = {}
        resConnectDictY = {}
        resConnectDictZ = {}
        scaleX = float(nx) / float(resX)
        scaleY = float(ny) / float(resY)
        scaleZ = float(nz) / float(resZ)
        for nneu in range(resX):
            sliceValue = []
            d1n = nneu * scaleX
            d2n = (nneu + 1) * scaleX
            for n in range(nx):
                d1a = n
                d2a = n + 1
                if d1a >= d1n and d1a < d2n and d2a <= d2n:
                    sliceValue.append((d1a, 1.0))
                elif d1a >= d1n and d1a < d2n and d2a > d2n:
                    ratio = float(d2n - d1a)
                    sliceValue.append((d1a, ratio))
                elif d1a < d1n and d2a > d1n:
                    ratio = float(d2a - d1n)
                    sliceValue.append((d1a, ratio))

            resConnectDictX[nneu] = sliceValue

        for nneu in range(resY):
            sliceValue = []
            d1n = nneu * scaleY
            d2n = (nneu + 1) * scaleY
            for n in range(ny):
                d1a = n
                d2a = n + 1
                if d1a >= d1n and d1a < d2n and d2a <= d2n:
                    sliceValue.append((d1a, 1.0))
                elif d1a >= d1n and d1a < d2n and d2a > d2n:
                    ratio = float(d2n - d1a)
                    sliceValue.append((d1a, ratio))
                elif d1a < d1n and d2a > d1n:
                    ratio = float(d2a - d1n)
                    sliceValue.append((d1a, ratio))

            resConnectDictY[nneu] = sliceValue

        for nneu in range(resZ):
            sliceValue = []
            d1n = nneu * scaleZ
            d2n = (nneu + 1) * scaleZ
            for n in range(nz):
                d1a = n
                d2a = n + 1
                if d1a >= d1n and d1a < d2n and d2a <= d2n:
                    sliceValue.append((d1a, 1.0))
                elif d1a >= d1n and d1a < d2n and d2a > d2n:
                    ratio = float(d2n - d1a)
                    sliceValue.append((d1a, ratio))
                elif d1a < d1n and d2a > d1n:
                    ratio = float(d2a - d1n)
                    sliceValue.append((d1a, ratio))

            resConnectDictZ[nneu] = sliceValue

        helpVoxelModel = numpy.zeros((resZ, resY, resX), numpy.float32)
        for nnZ in range(resZ):
            for nnY in range(resY):
                for nnX in range(resX):
                    listX = resConnectDictX[nnX]
                    listY = resConnectDictY[nnY]
                    listZ = resConnectDictZ[nnZ]
                    sum = 0.0
                    graySum = 0.0
                    for oldZ in listZ:
                        for oldY in listY:
                            for oldX in listX:
                                sum += oldX[1] * oldY[1] * oldZ[1]
                                graySum += voxelModel[oldZ[0], oldY[0], oldX[0]] * oldX[1] * oldY[1] * oldZ[1]

                    helpVoxelModel[nnZ, nnY, nnX] = graySum / sum

        if not ldim == None:
            ldim[0] = ldim[0] * float(nx) / float(resolut[0])
            ldim[1] = ldim[1] * float(ny) / float(resolut[1])
            ldim[2] = ldim[2] * float(nz) / float(resolut[2])
        return helpVoxelModel

    def changeResolution3(self, voxelModel, resolut, ldim):
        """
        Function changes the resultion of a given voxel model. "ldim" 
        will be updated by reference.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param resolut: approximate size of new image voxels 
            - TYPE: list[0] = sizeX, list[1] = sizeY, list[2] = sizeZ
                list[3] = enhance, list[4] = repeat
            - float sizeX, sizeY, sizeZ ... new number of voxels in x,y,z
            - float enhance: enhancement factor 0...1.0 -> blur,
                1.0...2.0 -> sharpen.
            - int repeat: factor how often the enhancement should be repeated
        @param ldim: voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
            
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
          ldim: new voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z          
        """
        stdout.write(' ... Change Resolution 3\n')
        stdout.flush()
        if ldim == None:
            stdout.write('\n **ERROR** changeResolution3(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lXn = resolut[0]
        lYn = resolut[1]
        lZn = resolut[2]
        lX0 = ldim[0]
        lY0 = ldim[1]
        lZ0 = ldim[2]
        nx, ny, nz = self.get_Shape(voxelModel)
        resX = int(nx * lX0 / lXn)
        resY = int(ny * lY0 / lYn)
        resZ = int(nz * lZ0 / lZn)
        helpVoxelModel = self.changeResolution2(voxelModel, [resX, resY, resZ], ldim)
        return helpVoxelModel

    def changeResolution3(self, voxelModel, resolut, ldim):
        """
        Function changes the resultion of a given voxel model. "ldim" 
        will be updated by reference.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param resolut: approximate size of new image voxels 
            - TYPE: list[0] = sizeX, list[1] = sizeY, list[2] = sizeZ
                list[3] = enhance, list[4] = repeat
            - float sizeX, sizeY, sizeZ ... new number of voxels in x,y,z
            - float enhance: enhancement factor 0...1.0 -> blur,
                1.0...2.0 -> sharpen.
            - int repeat: factor how often the enhancement should be repeated
        @param ldim: voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
            
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
          ldim: new voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z          
        """
        stdout.write(' ... Change Resolution 3\n')
        stdout.flush()
        if ldim == None:
            stdout.write('\n **ERROR** changeResolution3(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        lXn = resolut[0]
        lYn = resolut[1]
        lZn = resolut[2]
        lX0 = ldim[0]
        lY0 = ldim[1]
        lZ0 = ldim[2]
        nx, ny, nz = self.get_Shape(voxelModel)
        resX = int(nx * lX0 / lXn)
        resY = int(ny * lY0 / lYn)
        resZ = int(nz * lZ0 / lZn)
        helpVoxelModel = self.changeResolution2(voxelModel, [resX, resY, resZ], ldim)
        return helpVoxelModel

    def changeResolutionF(self, voxelModel, factor, ldim):
        """
        Function changes the resolution of a given voxel model. "ldim" 
        will be updated by reference.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param resolut: recoarsening factor  
            - TYPE: int
        @param ldim: voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
            
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
          ldim: new voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z          
        """

        if ldim.all() == None:                                                                         # 12.01.01 change: if ldim == None:     to if ldim.all() == None
            stdout.write('\n **ERROR** changeResolutionF(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx, ny, nz = self.get_Shape(voxelModel)
        if factor < 1:
            stdout.write(' ... refine resolution with factor %i \n' % (-1 * factor))
            stdout.flush()
            nnx = -nx * factor
            nny = -ny * factor
            nnz = -nz * factor
            ldim[0] = -ldim[0] / float(factor)
            ldim[1] = -ldim[1] / float(factor)
            ldim[2] = -ldim[2] / float(factor)
        elif factor > 1:
            stdout.write(' ... recoarse resolution with factor %i \n' % factor)
            stdout.flush()
            nnx = int(nx / factor)
            nny = int(ny / factor)
            nnz = int(nz / factor)
            ldim[0] = ldim[0] * factor
            ldim[1] = ldim[1] * factor
            ldim[2] = ldim[2] * factor
        else:
            stdout.write('\n **ERROR** changeResolutionF(): Factor needs to be lower than -1 or bigger than 1!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        modVoxelModel = self.createVoxelModel(nnx, nny, nnz, 'f')
        if factor > 1:
            voxelModel, modVoxelModel = micf77.changeresolution(voxelModel, modVoxelModel, 0)
        else:
            voxelModel, modVoxelModel = micf77.changeresolution2(voxelModel, modVoxelModel, 0)
        del voxelModel
        return modVoxelModel

    def changeResolutionDF(self, voxelModel, direct, factor, ldim):
        """
        Function changes the resultion of a given voxel model. "ldim" 
        will be updated by reference.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param resolut: recoarsening factor  
            - TYPE: int
        @param ldim: voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
            
        @return:
          voxelModel: voxel model which was read from file
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
          ldim: new voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z          
        """
        stdout.write(' ... change resolution with direct=%i factor=%i \n' % (direct, factor))
        stdout.flush()
        if ldim == None:
            stdout.write('\n **ERROR** changeResolutionDF(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx, ny, nz = self.get_Shape(voxelModel)
        nnx = nx
        nny = ny
        nnz = nz
        if direct == 1:
            nnx = int(nx * factor)
        elif direct == 2:
            nny = int(ny * factor)
        elif direct == 3:
            nnz = int(nz * factor)
        else:
            stdout.write('\n **ERROR** changeResolutionDF(): dir=%i not implemented! \n\n' % direct)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        modVoxelModel = self.createVoxelModel(nnx, nny, nnz, 'f')
        if direct == 1:
            for i in range(nx):
                for k in range(nz):
                    for j in range(ny):
                        for ii in range(factor):
                            modVoxelModel[k, j, i * factor + ii] = voxelModel[k, j, i]

        if direct == 2:
            for j in range(ny):
                for k in range(nz):
                    for i in range(nx):
                        for jj in range(factor):
                            modVoxelModel[k, j * factor + jj, i] = voxelModel[k, j, i]

        if direct == 3:
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        for kk in range(factor):
                            modVoxelModel[k * factor + kk, j, i] = voxelModel[k, j, i]

        del voxelModel
        if direct == 1:
            ldim[0] = ldim[0] / float(factor)
        elif direct == 2:
            ldim[1] = ldim[1] / float(factor)
        elif direct == 3:
            ldim[2] = ldim[2] / float(factor)
        return modVoxelModel

    def computeCgalCenter(self, voxelModel, threshold, ldim, file=None, mode=None):
        """
        Function cuts a block/RVE from a big voxel models and returns it.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param  threshold: Threshold value
            - TYPE: float
        @param ldim: voxel dimensions in x,y,z
            - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
            - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
                      
        @return:
          cgalCenterRadius: information about center and radius of  bounding 
            sphere.  
            - TYPE:list[0] = cenX, list[1] = cenY, list[2] = cenZ, 
                    list[3] = radius
        """
        stdout.write(' ... compute Cgal center/radius with threshold %s   \n' % repr(threshold))
        stdout.flush()
        if ldim == None:
            stdout.write('\n **ERROR** computeCgalCenter(): Voxel size not optional for this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        nx, ny, nz = self.get_Shape(voxelModel)
        isOutside = True
        i = int(nx / 2)
        j = int(ny / 2)
        k = int(nz / 2)
        n = 0
        ni = -1
        nj = -1
        nk = -1
        cenX = -1.0
        cenY = -1.0
        cenZ = -1.0
        radius = -1.0
        while isOutside:
            if i - n < 0 or i + n >= nx or j - n < 0 or j + n >= ny or k - n < 0 or k + n >= nz:
                isOutside = False
                continue
            else:
                if voxelModel[k, j, i + n] >= threshold:
                    isOutside = False
                    cenX = (i + n + 0.5) * ldim[0]
                    cenY = (j + 0.5) * ldim[1]
                    cenZ = (k + 0.5) * ldim[2]
                    radius = numpy.sqrt(cenX * cenX + cenY * cenY + cenZ * cenZ)
                    ni = i + n
                    nj = j
                    nk = k
                    continue
                if voxelModel[k, j + n, i] >= threshold:
                    isOutside = False
                    cenX = (i + 0.5) * ldim[0]
                    cenY = (j + n + 0.5) * ldim[1]
                    cenZ = (k + 0.5) * ldim[2]
                    radius = numpy.sqrt(cenX * cenX + cenY * cenY + cenZ * cenZ)
                    ni = i
                    nj = j + n
                    nk = k
                    continue
                if voxelModel[k + n, j, i] >= threshold:
                    isOutside = False
                    cenX = (i + 0.5) * ldim[0]
                    cenY = (j + 0.5) * ldim[1]
                    cenZ = (k + n + 0.5) * ldim[2]
                    radius = numpy.sqrt(cenX * cenX + cenY * cenY + cenZ * cenZ)
                    ni = i
                    nj = j
                    nk = k + n
                    continue
                if voxelModel[k, j, i - n] >= threshold:
                    isOutside = False
                    cenX = (i - n + 0.5) * ldim[0]
                    cenY = (j + 0.5) * ldim[1]
                    cenZ = (k + 0.5) * ldim[2]
                    dist = nx * ldim[0] - cenX
                    radius = numpy.sqrt(dist * dist + cenY * cenY + cenZ * cenZ)
                    ni = i - n
                    nj = j
                    nk = k
                    continue
                if voxelModel[k, j - n, i] >= threshold:
                    isOutside = False
                    cenX = (i + 0.5) * ldim[0]
                    cenY = (j - n + 0.5) * ldim[1]
                    cenZ = (k + 0.5) * ldim[2]
                    dist = ny * ldim[1] - cenY
                    radius = numpy.sqrt(cenX * cenX + dist * dist + cenZ * cenZ)
                    ni = i
                    nj = j - n
                    nk = k
                    continue
                if voxelModel[k - n, j, i] >= threshold:
                    isOutside = False
                    cenX = (i + 0.5) * ldim[0]
                    cenY = (j + 0.5) * ldim[1]
                    cenZ = (k - n + 0.5) * ldim[2]
                    dist = nz * ldim[2] - cenZ
                    radius = numpy.sqrt(cenX * cenX + cenY * cenY + dist * dist)
                    ni = i
                    nj = j
                    nk = k - n
                    continue
            n += 1

        if radius < 0.0:
            stdout.write('     -> no valid voxel found !              \n')
            stdout.flush()
        else:
            stdout.write('     -> voxX:voxY:voxZ  thres : %i:%i:%i  %f\n' % (ni, nj, nk, voxelModel[nk, nj, ni]))
            stdout.write('     -> cenX:cenY:cenZ        : %f:%f:%f \n' % (cenX, cenY, cenZ))
            stdout.write('     -> radius                : %f       \n' % radius)
            stdout.flush()
        if file != None:
            if mode == 'a':
                osfil = open(file, 'a')
            elif mode == 'w':
                osfil = open(file, 'w')
                osfil.write('#%s;%s;%s;%s;%s;%s;%s;%s\n' % ('voxX', 'voxY', 'voxZ',
                                                            'Threshold', 'centerX',
                                                            'centerY', 'centerZ',
                                                            'radius'))
            else:
                stdout.write("\n **ERROR** computeCgalCenter(): Output mode '%s' not known!\n\n" % mode)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            osfil.write('%i;%i;%i;%g;%g;%g;%g;%g\n' % (ni, nj, nk, voxelModel[nk, nj, ni], cenX, cenY, cenZ, radius))
            osfil.close()
        return [
         cenX, cenY, cenZ, radius]

    def cutBlockROI(self, voxelModel, cutList, additData=None, echo=True):
        """
        Function cuts a block/RVE from a big voxel models and returns it. 
        It uses a ROI given in physical dimensions- 
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param  cutList: Cutting/cropping parameters
            - TYPE: list[0] = n0X, list[1] = n0Y, list[2] = n0Z
                 list[3] = dnX, list[4] = dnY, list[5] = dnZ
            - int x0, y0, z0 ... x,y,z coordinate of ROI start 
            - int x1, y1, z1 ... x,y,z coordinate of ROI end 
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
        
        @return:
          newVoxelModelvoxel: cut/cropped model
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        """
        if echo == True:
            stdout.write(' ... cut block by physical dimensions\n')
            stdout.flush()
        x0c = float(cutList[0])
        y0c = float(cutList[1])
        z0c = float(cutList[2])
        x1c = float(cutList[3])
        y1c = float(cutList[4])
        z1c = float(cutList[5])
        vOFF = [
         0.0, 0.0, 0.0]
        if additData.has_key('Offset'):
            vOFF = additData['Offset']
        vnx, vny, vnz = self.get_Shape(voxelModel)
        vlx, vly, vlz = additData['-ldim']
        cOFF = [
         x0c, y0c, z0c]
        nx0 = int((cOFF[0] - vOFF[0]) / float(vlx))
        ny0 = int((cOFF[1] - vOFF[1]) / float(vly))
        nz0 = int((cOFF[2] - vOFF[2]) / float(vlz))
        prec = 1e-07
        dnx = int((x1c - x0c + prec) / vlx)
        dny = int((y1c - y0c + prec) / vly)
        dnz = int((z1c - z0c + prec) / vlz)
        return self.cutBlock(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], additData=additData, echo=False)

    def cutBlock(self, voxelModel, cutList, additData=None, echo=True):
        """
        Function cuts a block/RVE from a big voxel models and returns it.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param  cutList: Cutting/cropping parameters
            - TYPE: list[0] = n0X, list[1] = n0Y, list[2] = n0Z
                 list[3] = dnX, list[4] = dnY, list[5] = dnZ
            - int n0X, n0Y, n0Z ... start voxel in x,y,z
                first voxel in old voxel model =(0,0,0) and in new (n0X,n0Y,n0Z)
            - int dnX, dnY, dnZ ... number of voxels of the new voxel model
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
            
        @return:
          newVoxelModelvoxel: cut/cropped model
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        """
        if echo == True:
            stdout.write(' ... cut block\n')
            stdout.flush()
        nx0 = cutList[0]
        ny0 = cutList[1]
        nz0 = cutList[2]
        dnx = cutList[3]
        dny = cutList[4]
        dnz = cutList[5]
        nx, ny, nz = self.get_Shape(voxelModel)
        if nx0 + dnx > nx:
            stdout.write('\n **ERROR** cutBlock(): Out of model!\n   nx0 + dnx > nx  <->  %s + %s > %s\n' % (nx0, dnx, nx))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if ny0 + dny > ny:
            stdout.write('\n **ERROR** cutBlock(): Out of model!\n   ny0 + dny > ny  <->  %s + %s > %s\n' % (ny0, dny, ny))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if nz0 + dnz > nz:
            stdout.write('\n **ERROR** cutBlock():  Out of model!\n   nz0 + dnz > nz  <->  %s + %s > %s\n' % (nz0, dnz, nz))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if nx0 < 0:
            stdout.write('\n **ERROR** cutBlock():  Out of model!\n   nx0 = %s < 0 \n' % nx0)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if ny0 < 0:
            stdout.write('\n **ERROR** cutBlock():  Out of model!\n   ny0 = %s < 0 \n' % ny0)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if nz0 < 0:
            stdout.write('\n **ERROR** cutBlock():  Out of model!\n   nz0 = %s < 0 \n' % nz0)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        sum = 0
        time1 = time.clock()
        allElems = dnx * dny * dnz
        newVoxelModel = self.createVoxelModel(dnx, dny, dnz, 'f')
        newVoxelModel = voxelModel[nz0:nz0 + dnz, ny0:ny0 + dny, nx0:nx0 + dnx]
        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox + nx0 * lx, oy + ny0 * ly, oz + nz0 * lz]
        return newVoxelModel

    def mcutBlock(self, voxelModel, cutList, additData=None, echo=True):
        """
        Function cuts a block/RVE from the middle of a voxel models and returns it.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param  cutList: Cutting/cropping parameters
            - TYPE: list[0] = dnX, list[1] = dnY, list[2] = dnZ
            - int dnX, dnY, dnZ ... number of voxels of the new voxel model
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
            
        @return:
          newVoxelModelvoxel: cut/cropped model
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        """
        if echo == True:
            stdout.write(' ... Cut Block from Middle\n')
            stdout.flush()
        dnx = cutList[0]
        dny = cutList[1]
        dnz = cutList[2]
        nx, ny, nz = self.get_Shape(voxelModel)
        if dnx > nx:
            stdout.write('\n **ERROR** mcutBlock(): Out of model!\n   dnx > nx  <->  %s > %s\n' % (dnx, nx))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dny > ny:
            stdout.write('\n **ERROR** mcutBlock(): Out of model!\n   dny > ny  <->  %s > %s\n' % (dny, ny))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if dnz > nz:
            stdout.write('\n **ERROR** mcutBlock():  Out of model!\n   dnz > nz  <->  %s > %s\n' % (dnz, nz))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        sum = 0
        time1 = time.clock()
        allElems = dnx * dny * dnz
        newVoxelModel = self.createVoxelModel(dnx, dny, dnz, 'f')
        nx0 = int((nx - dnx) / 2)
        ny0 = int((ny - dny) / 2)
        nz0 = int((nz - dnz) / 2)
        newVoxelModel = voxelModel[nz0:nz0 + dnz, ny0:ny0 + dny, nx0:nx0 + dnx]
        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox + nx0 * lx, oy + ny0 * ly, oz + nz0 * lz]
        return newVoxelModel

    def threshold(self, voxelModel, thres, echo=True):
        """
        Function segment the given voxel model by using a given threshold value list.
        The smallest data value is the 'external' region.
        
        @param voxelModel: voxel model of the RVE
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        @param  thres: Threshold list
            - TYPE: list[treshID] = value
            - int thresID ... threshold ID, 0,1,2 .....
            - int value   ... threshold value
        @param  echo: Flag for printing the  " ... Threshold Data" on stdout
            - TYPE: bool
        
        @return:
          voxelModel: thresholded voxel model
            - TYPE: numpy.array[iZ, jY, kX] = grayValue
            - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
            - int grayValue  ... gray value of voxel,
        """
        if echo == True:
            stdout.write(' ... threshold data\n')
            stdout.flush()
        minValue, maxThres = self.computeNumpyMinMax(voxelModel, 2)
        if minValue > thres[0]:
            stdout.write('\n **ERROR** threshold(): Lower threshold values (=%s) is smaller then ' % repr(thres[0]))
            stdout.write('\n                        lowest grayvalue (=%s)!\n\n' % repr(minValue))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        sum = 0
        time1 = time.clock()
        ctime1 = time.time()
        nx, ny, nz = self.get_Shape(voxelModel)
        if len(thres) == 1:
            voxelModel = micf77.computethreshold(voxelModel, float(thres[0]), 2)
            sum = nx * ny * nz
        elif len(thres) == 2:
            voxelModel = micf77.computethreshold2(voxelModel, float(thres[0]), float(thres[1]), 2)
            sum = nx * ny * nz
        else:
            rangeList = []
            if minValue != thres[0]:
                rangeList.append(minValue)
            for matid in thres:
                rangeList.append(matid)

            rangeList.sort()
            flatArray = numpy.reshape(voxelModel, (nz * ny * nx,))
            for k in xrange(nz * ny * nx):
                sum += 1
                value = flatArray[k]
                for pos in range(len(rangeList)):
                    if pos + 1 != len(rangeList):
                        if value >= rangeList[pos] and value < rangeList[pos + 1]:
                            flatArray[k] = rangeList[pos]
                            break
                    else:
                        flatArray[k] = rangeList[pos]

            voxelModel = numpy.reshape(flatArray, (nz, ny, nx))
            boneStatistics, sumVox = myModel.boneStatistics(voxelModel, thres)
        time2 = time.clock()
        ctime2 = time.time()
        if echo == True:
            stdout.write('     -> Processed Data      : %10i      \n' % sum)
            stdout.flush()
            stdout.write('     -> Threshold Data in   :   %8.1f / %8.1fsec    \n' % (time2 - time1, ctime2 - ctime1))
            stdout.flush()
        return voxelModel

    def showHistogramFit(self, voxelModel, fitNo, infile=None, outfile=None, mode='w'):
        """
        Function shows the histogram and 1,2,3 fit functions
        """
        nz, ny, nx = voxelModel.shape
        noVox = nz * ny * nx
        stdout.write(' ... show histogram and fit \n')
        fit = int(fitNo)
        interval = []
        noInterval = 256
        histogram, minVox, maxVox = micf77.computehistogram2(voxelModel, noInterval - 1, 2)
        dInt = (maxVox + 1 - minVox) / float(noInterval)
        for i in range(noInterval):
            interval.append(dInt * i + minVox)

        maxHisto = -100000000
        maxI = 0
        for i in range(noInterval):
            if histogram[i] > maxHisto:
                maxHisto = histogram[i]
                maxI = i

        stdout.write('     -> compute fit \n')
        stdout.flush()
        x = interval
        y = histogram
        if fit == 1:
            fp = lambda v, x: v[0] * numpy.exp(-((x - v[1]) / v[2]) ** 2)
            e = lambda v, x, y: fp(v, x) - y
            v0 = [
             10000.0, maxI, 10.0]
            v, success = leastsq(e, v0, args=(x, y), maxfev=100000)
            print '    -> Estimater parameters: ', v
            X = numpy.linspace(minVox, maxVox, 200)
            p1, = plt.plot(interval, histogram, 'ro', label='Histogram')
            p2, = plt.plot([v[1] + 2 * v[2], v[1] + 2 * v[2]], [0, v[0]], 'y-', label='2*sigma')
            p3, = plt.plot(X, fp(v, X), label='SumFit')
            plt.legend()
            plt.show()
            if outfile != None:
                if mode == 'a':
                    osfil = open(outfile, 'a')
                elif mode == 'w':
                    osfil = open(outfile, 'w')
                    osfil.write('#%s;%s;%s;%s\n' % ('inFile', 'A1', 'B1', 'C1'))
                else:
                    stdout.write("\n **ERROR** showHistogramFit(): Output mode '%s' not known!\n\n" % mode)
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                osfil.write('%s;%g;%g;%g\n' % (repr(infile), v[0], v[1], v[2]))
        elif fit == 2:
            fp = lambda v, x: v[0] * numpy.exp(-((x - v[1]) / v[2]) ** 2) + v[3] * numpy.exp(-((x - v[4]) / v[5]) ** 2)
            fp1 = lambda v, x: v[0] * numpy.exp(-((x - v[1]) / v[2]) ** 2)
            fp2 = lambda v, x: v[3] * numpy.exp(-((x - v[4]) / v[5]) ** 2)
            e = lambda v, x, y: fp(v, x) - y
            v0 = [
             10000.0, maxI, 10.0, 1000.0, maxI * 1.5, 10.0]
            v, success = leastsq(e, v0, args=(x, y), maxfev=100000)
            print '    -> Estimater parameters: ', v
            X = numpy.linspace(minVox, maxVox, 200)
            p1, = plt.plot(interval, histogram, 'ro', label='Histogram')
            p2, = plt.plot(X, fp(v, X), label='SumFit')
            p3, = plt.plot(X, fp1(v, X), label='Mat1')
            p4, = plt.plot(X, fp2(v, X), label='Mat2')
            plt.legend()
            plt.show()
            if outfile != None:
                if mode == 'a':
                    osfil = open(outfile, 'a')
                elif mode == 'w':
                    osfil = open(outfile, 'w')
                    osfil.write('#%s;%s;%s;%s;%s;%s;%s\n' % ('inFile', 'A1', 'B1',
                                                             'C1', 'A2', 'B2', 'C2'))
                else:
                    stdout.write("\n **ERROR** showHistogramFit(): Output mode '%s' not known!\n\n" % mode)
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                osfil.write('%s;%g;%g;%g;%g;%g;%g\n' % (repr(infile), v[0], v[1], v[2], v[3], v[4], v[5]))
        elif fit == 3:
            fp = lambda v, x: v[0] * numpy.exp(-((x - v[1]) / v[2]) ** 2) + v[3] * numpy.exp(-((x - v[4]) / v[5]) ** 2) + v[6] * numpy.exp(-((x - v[7]) / v[8]) ** 2)
            fp1 = lambda v, x: v[0] * numpy.exp(-((x - v[1]) / v[2]) ** 2)
            fp2 = lambda v, x: v[3] * numpy.exp(-((x - v[4]) / v[5]) ** 2)
            fp3 = lambda v, x: v[6] * numpy.exp(-((x - v[7]) / v[8]) ** 2)
            e = lambda v, x, y: fp(v, x) - y
            v0 = [
             10000.0, (maxI + minVox) / 2.0, 10.0, 1000.0, maxI, 100.0, 1000.0, (maxI + maxVox) / 2, 10.0]
            v, success = leastsq(e, v0, args=(x, y), maxfev=100000)
            print '    -> Estimater parameters: ', v
            X = numpy.linspace(minVox, maxVox, 200)
            p1, = plt.plot(interval, histogram, 'ro', label='Histogram')
            p2, = plt.plot(X, fp(v, X), label='SumFit')
            p3, = plt.plot(X, fp1(v, X), label='Mat1')
            p4, = plt.plot(X, fp2(v, X), label='Mat2')
            p5, = plt.plot(X, fp3(v, X), label='Mat3')
            plt.legend()
            plt.show()
            if outfile != None:
                if mode == 'a':
                    osfil = open(outfile, 'a')
                elif mode == 'w':
                    osfil = open(outfile, 'w')
                    osfil.write('#%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n' % ('inFile', 'A1',
                                                                      'B1', 'C1',
                                                                      'A2', 'B2',
                                                                      'C2', 'A3',
                                                                      'B3', 'C3'))
                else:
                    stdout.write("\n **ERROR** showHistogramFit(): Output mode '%s' not known!\n\n" % mode)
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                osfil.write('%s;%g;%g;%g;%g;%g;%g;%g;%g;%g\n' % (repr(infile), v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]))
        else:
            stdout.write('\n **ERROR**  showHistogramFit(): fit has to be 1,2 or 3 not %i!\n\n' % fitNo)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            exit(1)
        return

    def computeFixThreshold(self, voxelModel, thresRatio, minVal, maxVal, echo=False):
        """
        Function segment the given voxel model by using a fixed threshold 
        using the Laplace Hamming filtering
        """
        if echo == True:
            stdout.write(' ... compute fix threshold \n')
            stdout.flush()
        time1 = time.clock()
        maxValue = 0.0
        minValue = 0.0
        minVox, maxVox = self.computeNumpyMinMax(voxelModel, 0)
        if minVal.lower() == 'minval':
            minValue = minVox
        else:
            minValue = float(minVal)
        if maxVal.lower() == 'maxval':
            maxValue = maxVox
        else:
            maxValue = float(maxVal)
        newThreshold = (maxValue - minValue) * float(thresRatio)
        voxelModel = self.scaleModel(voxelModel, minValue, maxValue, echo=False)
        voxelModel = micf77.computethreshold(voxelModel, newThreshold, 2)
        voxelModel = self.castType(voxelModel, 'f')
        time2 = time.clock()
        stdout.write('     -> finished in         :   %8.1f sec         \n' % (time2 - time1))
        stdout.flush()
        return (
         voxelModel, newThreshold)

    def computeLaplaceHamming(self, voxelModel, laplace_eps, lp_cutoff_freq, hamming_amp, ldim, echo=False):
        """
        Function computes the laplace-hamming filter (Scanco implementation)
        """
        time1 = time.clock()
        if echo == True:
            stdout.write(' ... compute Laplace-Hamming (Scanco)\n')
        nx, ny, nz = self.get_Shape(voxelModel)
        lx, ly, lz = ldim[0], ldim[1], ldim[2]
        px = 0
        while 2 ** px < nx:
            px += 1

        nfftx = 2 ** px
        py = 0
        while 2 ** py < ny:
            py += 1

        nffty = 2 ** py
        pz = 0
        while 2 ** pz < nz:
            pz += 1

        nfftz = 2 ** pz
        if echo == True:
            stdout.write('     -> setup model for FFT  \n')
            stdout.flush()
        max_phys_freq = 1.0 / numpy.min([lx, ly, lz])
        lp_phys_freq2 = (max_phys_freq * lp_cutoff_freq) ** 2
        nyquist_phys = 0.5 * max_phys_freq
        if echo == True:
            stdout.write('        dims         = %6i %6i %6i\n' % (nx, ny, nz))
            stdout.write('        dims_fft     = %6i %6i %6i\n' % (nfftx, nffty, nfftz))
            stdout.write('        phys_len     = %6.3f %6.3f %6.3f\n' % (nx * lx, ny * ly, nz * lz))
            stdout.write('        phys_len_fft = %6.3f %6.3f %6.3f\n' % (nfftx * lx, nffty * ly, nfftz * lz))
            stdout.write('        nyquist phys = %9.6f [lp/mm] \n' % nyquist_phys)
            stdout.write('        lp_phys_freq = %9.6f [lp/mm] \n' % numpy.sqrt(lp_phys_freq2))
        mirrVoxelModel = self.createVoxelModel(nfftx, nffty, nfftz, 'f')
        shift = [int((nfftx - nx) / 2), int((nffty - ny) / 2), int((nfftz - nz) / 2)]
        voxelModel, mirrVoxelModel = micf77.computemirror2(voxelModel, mirrVoxelModel, shift, 2)
        del voxelModel
        voxelModel = mirrVoxelModel
        if echo == True:
            stdout.write('     -> FFT transform        \n')
            stdout.flush()
        freqMod = numpy.fft.fftn(voxelModel)
        lx2 = nfftx * nfftx * lx * lx
        ly2 = nffty * nffty * ly * ly
        lz2 = nfftz * nfftz * lz * lz
        phys_freq2 = self.createVoxelModel(nfftx, nffty, nfftz, 'f')
        hamming = self.createVoxelModel(nfftx, nffty, nfftz, 'f')
        for k in range(nfftz / 2):
            for j in range(nffty / 2):
                for i in range(nfftx / 2):
                    i2 = (i + 0.5) * (i + 0.5)
                    j2 = (j + 0.5) * (j + 0.5)
                    k2 = (k + 0.5) * (k + 0.5)
                    sphys_freq2 = i2 / lx2 + j2 / ly2 + k2 / lz2
                    phys_freq2[k, j, i] = sphys_freq2
                    phys_freq2[k, j, nfftx - 1 - i] = sphys_freq2
                    phys_freq2[k, nffty - 1 - j, i] = sphys_freq2
                    phys_freq2[k, nffty - 1 - j, nfftx - 1 - i] = sphys_freq2
                    phys_freq2[nfftz - 1 - k, j, i] = sphys_freq2
                    phys_freq2[nfftz - 1 - k, j, nfftx - 1 - i] = sphys_freq2
                    phys_freq2[nfftz - 1 - k, nffty - 1 - j, i] = sphys_freq2
                    phys_freq2[nfftz - 1 - k, nffty - 1 - j, nfftx - 1 - i] = sphys_freq2
                    if sphys_freq2 > lp_phys_freq2:
                        hamm = 0.0
                    else:
                        hamm = 1.0 + hamming_amp / 2.0 * (numpy.cos(numpy.pi * numpy.sqrt(sphys_freq2 / lp_phys_freq2)) - 1.0)
                    hamming[k, j, i] = hamm
                    hamming[k, j, nfftx - 1 - i] = hamm
                    hamming[k, nffty - 1 - j, i] = hamm
                    hamming[k, nffty - 1 - j, nfftx - 1 - i] = hamm
                    hamming[nfftz - 1 - k, j, i] = hamm
                    hamming[nfftz - 1 - k, j, nfftx - 1 - i] = hamm
                    hamming[nfftz - 1 - k, nffty - 1 - j, i] = hamm
                    hamming[nfftz - 1 - k, nffty - 1 - j, nfftx - 1 - i] = hamm

        if echo == True:
            stdout.write('     -> FFT back transform \n')
            stdout.flush()
        n_factor = 4.0 * numpy.pi ** 2
        newMod = numpy.fft.ifftn(freqMod * n_factor * (1.0 + laplace_eps * (phys_freq2 - 1.0)) * hamming)
        if echo == True:
            stdout.write('     -> cut model          \n')
            stdout.flush()
        voxelModel = self.cutBlock(newMod.real, [shift[0], shift[1], shift[2], nx, ny, nz], echo=False)
        voxelModel = self.castType(voxelModel, 'f')
        del newMod
        time2 = time.clock()
        stdout.write('     -> finished in         :   %8.1f sec         \n' % (time2 - time1))
        stdout.flush()
        return voxelModel

    def computeLaplaceHammingTesting(self, voxelModel, weight, cutoff, ampl, ldim, echo=False):
        """
        Function computs the laplace-hamming filter
        """
        time1 = time.clock()
        if echo == True:
            stdout.write(' ... compute Laplace-Hamming \n')
        if weight < 0 or weight > 1:
            stdout.write('\n **ERROR**  computeLaplaceHamming(): weight=%g (0 < weight < 1)!\n\n' % weight)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            exit(1)
        if cutoff < 0 or cutoff > 1:
            stdout.write('\n **ERROR**  computeLaplaceHamming(): cutoff=%g (0 < cutoff < 1)!\n\n' % cutoff)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            exit(1)
        if ampl < 0 or ampl > 1:
            stdout.write('\n **ERROR**  computeLaplaceHamming(): amplitude=%g (0 < amplitude < 1)!\n\n' % ampl)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            exit(1)
        stdout.flush()
        nx, ny, nz = self.get_Shape(voxelModel)
        lx, ly, lz = ldim[0], ldim[1], ldim[2]
        nmax = numpy.max([nx, ny, nz])
        p = 0
        while 2 ** p < nmax:
            p += 1

        nfft = 2 ** p
        if echo == True:
            stdout.write('     -> make array %ix%ix%s  \n' % (nfft, nfft, nfft))
            stdout.flush()
        if nfft != nx and nfft != ny and nfft != nz:
            mirrVoxelModel = self.createVoxelModel(nfft, nfft, nfft, 'f')
            voxelModel, mirrVoxelModel = micf77.computemirror2(voxelModel, mirrVoxelModel, 1)
            del voxelModel
            voxelModel = mirrVoxelModel
        if echo == True:
            stdout.write('     -> FFT transform        \n')
            stdout.flush()
        freqMod = numpy.fft.fftn(voxelModel)
        hamfft = numpy.zeros(nfft, numpy.float32) + (1.0 - ampl)
        fcutoff = int(nfft / 2.0 * cutoff)
        for k in range(fcutoff):
            hamfft[k] = 1.0 - ampl / 2.0 + ampl / 2.0 * numpy.cos(2.0 * numpy.pi * k / (nfft * cutoff - 1.0))
            hamfft[nfft - 1 - k] = 1.0 - ampl / 2.0 + ampl / 2.0 * numpy.cos(2.0 * numpy.pi * k / (nfft * cutoff - 1.0))

        pi2N = numpy.pi / float(nfft)
        lx2 = lx * lx
        ly2 = ly * ly
        lz2 = lz * lz
        hamming = numpy.zeros((nfft, nfft, nfft), numpy.float32)
        fft_d = numpy.zeros((nfft, nfft, nfft), numpy.float32)
        for k in range(nfft):
            for j in range(nfft):
                for i in range(nfft):
                    hamming[k, j, i] = hamfft[k] * hamfft[j] * hamfft[i]
                    sin_i = numpy.sin(pi2N * i)
                    sin_j = numpy.sin(pi2N * j)
                    sin_k = numpy.sin(pi2N * k)
                    fft_d[k, j, i] = -4.0 * (sin_i * sin_i / lx2 + sin_j * sin_j / ly2 + sin_k * sin_k / lz2)

        if echo == True:
            stdout.write('     -> FFT back transform \n')
            stdout.flush()
        newMod = numpy.fft.ifftn((-freqMod * fft_d * weight + 39 * freqMod * (1.0 - weight)) * hamming)
        import pylab
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Signal along x-axis', fontsize=18)
        ax.set_xlabel('Voxel k', fontsize=16)
        ax.set_ylabel('Amplitude w(k)', fontsize=16)
        pylab.plot(numpy.fft.fftshift(numpy.hamming(nfft)), color='blue', label='Hamming Window')
        pylab.plot(numpy.fft.fftshift(numpy.hanning(nfft)), color='green', label='Hanning Window')
        pylab.plot(hamfft, color='red', label='Used Scaled/Cut Window')
        signal = freqMod
        ax.set_xlim(0, nfft)
        pylab.legend(loc=(0.25, 0.75))
        pylab.show()
        if echo == True:
            stdout.write('     -> cut model          \n')
            stdout.flush()
        voxelModel = self.cutBlock(newMod, [0, 0, 0, nx, ny, nz], echo=False)
        voxelModel = self.castType(voxelModel, 'f')
        time2 = time.clock()
        stdout.write('     -> finished in         :   %8.1f sec         \n' % (time2 - time1))
        stdout.flush()
        return voxelModel

    def computeCubeFill(self, voxelModel, nlevel):
        """
        Functions fills the voxel model with cubes of increasing size (octree)
        """
        stdout.write(' ... compute cube fill \n')
        stdout.flush()
        voxelModel = numpy.greater_equal(voxelModel, 1) * (nlevel + 1)
        nx, ny, nz = self.get_Shape(voxelModel)
        auxVoxelModel = self.createVoxelModel(nx, ny, nz, 'f')
        for i in range(1, nlevel + 1):
            auxVoxelModel = micf77.fillcubes(voxelModel, auxVoxelModel, i * 2, 2)
            voxelModel = voxelModel - auxVoxelModel

        del auxVoxelModel
        return self.castType(voxelModel, 'f')

    def clean(self, curVoxelModel, thresList, cleanList, echo=False):
        """
        Functions cleans the voxel model - removes islands and not proper connected regions.
        
        @param curVoxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        @param  thresList:  list of thresholds - only one threshold implemented!
              - TYPE: list[thresId] = thres
              - int thresID ... threshold id, sorted ascending
              - int thres   ... threshold value 0..25
        @param  cleanList: list of clean parameters
              - TYPE: list[0] = nodeShared, list[1] = islandSize
              - int nodeShared ... number of shared nodes of a bone to bone connection
              - int islandSize ... number of Voxel isolated bone/marrow voxel region
                  which should be removed.
        @param echo: Flag if extended echo should be written on stdout.
              - TYPE: bool
        
        @return:
           cleanVoxelModel: cleaned voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        """
        print ' ... clean model'
        nodeShared = cleanList[0]
        islandSize = cleanList[1]
        time1 = time.clock()
        ok = False
        nx, ny, nz = self.get_Shape(curVoxelModel)
        cleanVoxelModel = self.createVoxelModel(nx, ny, nz, 'f')
        cleanVoxelModel = curVoxelModel
        step = 0
        while not ok:
            step = step + 1
            stdout.write('     Analyses Step          : %10i       \n' % step)
            stdout.flush()
            sumid = 0
            if thresList == None or len(thresList) != 1:
                stdout.write('\n **ERROR** find_parts(): Exactly one threshold needed for this function!\n\n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            partMatNo = 0
            checkList = []
            partIdElemDict = {}
            partIdsVolume = numpy.zeros((nz, ny, nx))
            nbPartIdList = []
            matid = thresList[0]
            dpUtils.progressStart('     -> look for neighbours : ')
            for k in range(nz):
                progress = float(sumid) / float(nx * ny * nz) * 10.0
                for j in range(ny):
                    for i in range(nx):
                        sumid = sumid + 1
                        if i - 1 >= 0:
                            checkList.append((k, j, i - 1))
                        if j - 1 >= 0:
                            checkList.append((k, j - 1, i))
                        if k - 1 >= 0:
                            checkList.append((k - 1, j, i))
                        if nodeShared <= 2:
                            if i - 1 >= 0 and j - 1 >= 0:
                                checkList.append((k, j - 1, i - 1))
                            if i - 1 >= 0 and k - 1 >= 0:
                                checkList.append((k - 1, j, i - 1))
                            if k - 1 >= 0 and j - 1 >= 0:
                                checkList.append((k - 1, j - 1, i))
                        if nodeShared <= 1:
                            if i - 1 >= 0 and j - 1 >= 0 and k - 1 >= 0:
                                checkList.append((k - 1, j - 1, i - 1))
                        for checkElem in checkList:
                            iN = checkElem[2]
                            jN = checkElem[1]
                            kN = checkElem[0]
                            if cleanVoxelModel[k, j, i] == cleanVoxelModel[kN, jN, iN]:
                                partId = partIdsVolume[kN, jN, iN]
                                if partId not in nbPartIdList:
                                    nbPartIdList.append(partId)

                        partIdNo = len(nbPartIdList)
                        if partIdNo == 0:
                            partMatNo = partMatNo + 1
                            partIdsVolume[k, j, i] = partMatNo
                            partIdElemDict[partMatNo] = [(k, j, i)]
                        elif partIdNo == 1:
                            partIdsVolume[k, j, i] = nbPartIdList[0]
                            partIdElemDict[nbPartIdList[0]].append((k, j, i))
                        elif partIdNo > 1:
                            maxElem = 0
                            newPartId = 0
                            for partId in nbPartIdList:
                                if len(partIdElemDict[partId]) > maxElem:
                                    maxElem = len(partIdElemDict[partId])
                                    newPartId = partId

                            nbPartIdList.remove(newPartId)
                            for partId in nbPartIdList:
                                for element in partIdElemDict[partId]:
                                    partIdElemDict[newPartId].append(element)

                                for elemTuple in partIdElemDict[partId]:
                                    partIdsVolume[elemTuple[0], elemTuple[1], elemTuple[2]] = newPartId

                                del partIdElemDict[partId]

                            partIdsVolume[k, j, i] = newPartId
                            partIdElemDict[newPartId].append((k, j, i))
                        del nbPartIdList[0:len(nbPartIdList)]
                        del checkList[0:len(checkList)]

                dpUtils.progressNext(progress)

            dpUtils.progressEnd()
            ok = self.updateModel(partIdElemDict, cleanVoxelModel, thresList, islandSize)
            if echo == True or ok == True:
                for partId in partIdElemDict.keys():
                    ck = partIdElemDict[partId][0][0]
                    cj = partIdElemDict[partId][0][1]
                    ci = partIdElemDict[partId][0][2]
                    stdout.write('     -> Part %4i;  Elements in Part = %8i;  Material = %3i\n' % (partId, len(partIdElemDict[partId]), cleanVoxelModel[ck, cj, ci]))
                    stdout.flush()

        time2 = time.clock()
        stdout.write('     -> clean finished in   :   %8.1f sec         \n' % (time2 - time1))
        stdout.flush()
        return cleanVoxelModel

    def clean_numpy(self, voxelModel, N=1, min_size=0, invert=False, labeled=False, kernel='6'):
        """
        @author: lsteiner/dpahr
        
        Image cleaner that leaves the N largest structures.
        Uses numpy's label and bincount functions.
        
        Parameters:
        -----------
        
        voxelModel: int, float
            Input binary image array, 2- or 3-dimensional with exactly 2 grayvalues (BG/FG).
        
        N: int
            Number of contecting structures in the output image. Default: 1.
        
        min_size: int
            Minimal size of structures in voxel in the output image. 
            Structures are deleted in any case if they are smaller than min_size. 
            If min_size = 0 all structures are considered ie. noting is deleted.
            Default: 0. 
        
        invert: Bool
            If this is True, all operations are executed on the background value
        
        labeled: Bool
            If set True, the found structures are labeled consecutively, with 1 being the largest structure
        
        kernel_str: string, array
            A provided string defines the labeling kernel using kernel3d. Default: '6'. For further info read the corresponding help.
            A provided kernel array is used directly.
        
        Returns:
        --------
        
        array of same shape as input array with N largest structures optionally labeled with size-ranking
        """
        stdout.write('     -> check image\n')
        stdout.flush()
        values = numpy.unique(voxelModel)
        if len(values) != 2:
            dpUtils.throwError('clean_numpy(): Binary image not OK, exactly two values required!')
        if invert:
            bg = values.max()
        else:
            bg = values.min()
        if type(kernel) == str:
            kernel = self.kernel3d(kernel)
        if len(voxelModel.shape) == 2:
            kernel = kernel[(len(kernel) - 1) / 2]
        stdout.write('     -> compute labels\n')
        stdout.flush()
        import scipy.ndimage as ndi
        label, nlabels = ndi.label(voxelModel != bg, structure=kernel)
        sizes = numpy.bincount(label.ravel())
        stdout.write('     -> find largest region\n')
        stdout.flush()
        if min_size == 0:
            Nlargest_labels = numpy.argsort(sizes[1:])[-N:] + 1
            mask = numpy.in1d(label, Nlargest_labels).reshape(label.shape)
        elif min_size != 0:
            Nlargest_labels = numpy.argsort(sizes[1:])[-N:] + 1
            Nlargest_sizes = numpy.sort(sizes[1:])[-N:] + 1
            Nlargest_minsized_loc = numpy.where(Nlargest_sizes > min_size)[0]
            Nlargest_labels = Nlargest_labels[Nlargest_minsized_loc]
            mask = numpy.in1d(label, Nlargest_labels).reshape(label.shape)
        stdout.write('     -> update voxel model\n')
        stdout.flush()
        if not labeled:
            cleanVoxelModel = numpy.copy(voxelModel)
            cleanVoxelModel[mask == False] = bg
        elif labeled:
            cleanVoxelModel = numpy.copy(label)
            features = numpy.arange(N, dtype=label.dtype)[::-1] + bg + 1
            pairs = numpy.voxelModel((Nlargest_labels, features)).T
            for i, j in pairs:
                numpy.place(cleanVoxelModel, cleanVoxelModel == i, j)

            cleanVoxelModel[mask == False] = bg
        return cleanVoxelModel

    def boneStatistics(self, voxelModel, thresList, echo=True):
        """
        Computes bone statistics of a given voxel model (number of voxels for
        each material phase/threshold).
        
        @param  voxelModel: segmented voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel,
        @param  thresList:  list of thresholds
          - TYPE: list[thresId] = thres
          - int thresID ... threshold id, sorted ascending
          - int thres   ... threshold value
        @param echo: Flag if extended echo should be written on stdout.
          - TYPE: bool
            
        @return:
         BVlist: No of voxel within each threshold range e.g th1 <= X < th2
          - TYPE: dict[thres] = noVox
          - int thres ... threshold value  , first threshold value is "0"
          - int noVox ... number of voxels for this value
        """
        if echo == True:
            stdout.write(' ... compute BV/TV\n')
            stdout.flush()
        nx, ny, nz = self.get_Shape(voxelModel)
        BVTVlist = {}
        BVTVlist[0] = 0
        rangeList = []
        rangeList.append(0)
        for ID in thresList:
            rangeList.append(ID)
            BVTVlist[ID] = 0

        if len(thresList) == 1:
            ncount = micf77.computegreater(voxelModel, thresList[0] - 1.0, 0)
            BVTVlist[thresList[0]] = ncount
            BVTVlist[0] = nz * ny * nx - BVTVlist[thresList[0]]
        else:
            minThres, maxThres = self.computeNumpyMinMax(voxelModel, 0)
            if minThres < 0:
                stdout.write('\n **ERROR** : boneStatistics(): Multiple thresholds and negative values not allowed!\n\n')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            flatArray = numpy.reshape(voxelModel, (nz * ny * nx,))
            for k in range(nz * ny * nx):
                value = flatArray[k]
                for pos in range(len(rangeList)):
                    if pos + 1 != len(rangeList):
                        if value >= rangeList[pos] and value < rangeList[pos + 1]:
                            BVTVlist[rangeList[pos]] += 1
                            break
                    else:
                        BVTVlist[rangeList[pos]] += 1

            del flatArray
        sumVox = 0
        for region in BVTVlist:
            sumVox += BVTVlist[region]

        if sumVox != nx * ny * nz:
            stdout.write('\n **ERROR** : boneStatistics(): Problems with the Voxel sum!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if echo == True:
            for region in BVTVlist:
                if int(region) > 0:
                    BVTV = float(BVTVlist[region] / float(sumVox) * 100.0)
                    stdout.write('     -> threshold / BVTV    : %13i /%13.3f %s \n' % (int(region), BVTV, '%'))
                    stdout.flush()
                    stdout.write('     -> voxels    / sumVox  : %13i /%13i   \n' % (BVTVlist[region], sumVox))
                    stdout.flush()

        return (
         BVTVlist, sumVox)

    def findThreshold(self, BVTVErrorList, voxelModel):
        """
        Function finds a threshold for a given BVTV value within a given error.
        
        @param BVTVErrorList: Function parameters
           - TYPE: list[0] = BVTV, list[1]=error, list[2]=thresholdEstimate
           - float BVTV  ... BVTV/density which should be reached. Absolute value not in %
           - float error ... maximum error allowed. Absolute value not in %
           - int thresholdEstimate ... Estimation of threshold.
        @param  voxelModel: voxel model
           - TYPE: numpy.array[kZ, jY, iX] = grayValue
           - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
           - int grayValue  ... gray value of voxel,
        
        @return:
         estThreshold: Estimation of new threshold
           - TYPE: int
        """
        stdout.write(' ... find Threshold\n')
        stdout.flush()
        newThresholdList = []
        BVTV = BVTVErrorList[0]
        error = BVTVErrorList[1]
        thresestimate = BVTVErrorList[2]
        curError = 100.0
        deltaThres = abs(thresestimate) * 0.1
        curThresList1 = [thresestimate - deltaThres]
        curThresList2 = [thresestimate + deltaThres]
        estThreshold = 0
        curBVTV = 0.0
        step = 0
        minThres, maxThres = self.computeNumpyMinMax(voxelModel, 2)
        stdout.write('     -> start search: min=%s, max=%s ... \n' % (repr(minThres), repr(maxThres)))
        stdout.flush()
        while abs(curError) > error:
            step += 1
            thresVoxelModel = micf77.computethreshold(voxelModel.copy(), curThresList1[0], 0)
            BVTVlist, sumVox = self.boneStatistics(thresVoxelModel, curThresList1, echo=False)
            curBVTV1 = float(BVTVlist[curThresList1[0]]) / float(sumVox)
            thresVoxelModel = micf77.computethreshold(voxelModel.copy(), curThresList2[0], 0)
            BVTVlist, sumVox = self.boneStatistics(thresVoxelModel, curThresList2, echo=False)
            curBVTV2 = float(BVTVlist[curThresList2[0]]) / float(sumVox)
            y0 = curThresList1[0]
            y1 = curThresList2[0]
            x0 = curBVTV1
            x1 = curBVTV2
            x = BVTV
            if x1 == x0:
                print '\n **WARNING** findThreshold(): BVTV does not change during iteration!'
                print '             Error tolerance can not be reached.\n'
                break
            else:
                estThreshold = y0 + (y1 - y0) / (x1 - x0) * (x - x0)
            if estThreshold < minThres:
                print '\n\n **ERROR** findThreshold(): Given BVTV to high! No Threshold found.\n'
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if estThreshold > maxThres:
                print '\n\n **ERROR** findThreshold(): Given BVTV to low! No Threshold found.\n'
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            curThresList1[0] = curThresList1[0] + (estThreshold - curThresList1[0]) / 2.0
            curThresList2[0] = curThresList2[0] - (curThresList2[0] - estThreshold) / 2.0
            estThresholdList = [
             estThreshold]
            thresVoxelModel = micf77.computethreshold(voxelModel, estThresholdList, 0)
            BVTVlist, sumVox = self.boneStatistics(thresVoxelModel, estThresholdList, echo=False)
            curBVTV = float(BVTVlist[estThresholdList[0]]) / float(sumVox)
            curError = abs(BVTV - curBVTV) / BVTV
            stdout.write('     -> Step/Threshold/Error:   %3i/%8.2f/%8.4f \n' % (step, estThreshold, curError))
            stdout.flush()

        stdout.write('     -> current BVTV        :   %8.4f \n' % curBVTV)
        stdout.flush()
        return int(estThreshold + 0.5)

    def scaleModel(self, voxelModel, newMinValue, newMaxValue, oldMinValue=None, oldMaxValue=None, format=None, echo=True):
        """
        Function changes the min/max values of a voxel model. If old values are not
        given the min and max values of the voxel model will taken as these values.
        
        @param  voxelModel: voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel,
        @param  newMinValue: new mininum value
          - TYPE: int
        @param  newMaxValue: new maximum value
          - TYPE: int
        @param  oldMinValue: old mininum value. If not given the minimum data
          value of the voxel model will be used.
            - TYPE: int
        @param  oldMaxValue: old maximum value. If not given the minimum data
          value of the voxel model will be used.
            - TYPE: int
        @param  format: data format of the resulting array (B,h,f)
            - TYPE: string
        @param echo: Flag if extended echo should be written on stdout.
            - TYPE: bool
            
        @return:
         voxelModel: scaled voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel, 0..25
        """
        if echo == True:
            stdout.write(' ... scale model   \n')
            stdout.flush()
        minVox, maxVox = self.computeNumpyMinMax(voxelModel, 2)
        if oldMinValue == None or oldMinValue.lower() == 'oldmin':
            minValue = minVox
        else:
            minValue = float(oldMinValue)
        if oldMaxValue == None or oldMaxValue.lower() == 'oldmax':
            maxValue = maxVox
        else:
            maxValue = float(oldMaxValue)
        if repr(newMinValue).lower() == "'oldmin'":
            newMinValue = minVox
        else:
            newMinValue = float(newMinValue)
        if repr(newMaxValue).lower() == "'oldmax'":
            newMaxValue = maxVox
        else:
            newMaxValue = float(newMaxValue)
        if minValue > minVox:
            stdout.write('\n **ERROR** : scaleModel(): Lower scale range (=%s) is bigger than \n' % repr(minValue))
            stdout.flush()
            stdout.write('\n             minimum value in voxel model (=%s)!\n\n' % repr(minVox))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if maxValue < maxVox:
            stdout.write('\n **ERROR** : scaleModel(): Upper scale range (=%s) is smaller than \n' % repr(maxValue))
            stdout.flush()
            stdout.write('\n             maximum value in voxel model (=%s)!\n\n' % repr(maxVox))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        delta = float(abs(maxValue - minValue))
        newDelta = float(abs(newMaxValue - newMinValue))
        scale = newDelta / delta
        voxelModel = newMinValue + (voxelModel - minValue) * scale
        if echo == True:
            newMinVox, newMaxVox = self.computeNumpyMinMax(voxelModel, 0)
            stdout.write('     -> from  min / max     : %12.5g / %12.5g\n' % (minValue, maxValue))
            stdout.flush()
            stdout.write('     -> to    min / max     : %12.5g / %12.5g\n' % (newMinVox, newMaxVox))
            stdout.flush()
        return voxelModel

    def computeEquation(self, voxelModel, arithList, echo=1):
        """
        Algebraic manipulation of the voxelModel. Currently implemented::
        +A   ... add value A
        -A   ... subtract value A
        *A   ... multiply with A
        /A   ... divide by A
        <B   ... store current value in B
        >B   ... load B in current value
        ^I   ... power I (only integers!) of current value
        >>C  ... load C from file
        ?    ... logical condition
        &sin ... compute sinus of array in memory
        &cos ... compute cosinus of array in memory
        &arcsin ... compute arcus sinus of array in memory
        &arccos ... compute arcus cosinus of array in memory
        &tan    ... compute tangens of array in memory
        &arctan ... compute arcus tangens of array in memory
        
        Note all operations are bitwise operations!
        
        Example:: ?<0=0:?>5000=5000:/8000:*1500:-200:<BMD:+225:<BB:>BMD:+100:^2:/1047:/BB:?==75=80
        
        @param  voxelModel: voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param  arithList: list of arithmetic operations
             - TYPE: list[0] = Operation 1, list[1]= Operation 2, ....
        
        @return:
            voxelModel: modified voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel, 0..25
        """
        nx, ny, nz = self.get_Shape(voxelModel)
        helpVox = {}
        if echo > 0:
            stdout.write(' ... apply arithmetic operations \n')
            stdout.flush()
        for oper in arithList:
            if oper[0] == '*':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                if helpVox.has_key(oper[1:]):
                    voxelModel = voxelModel * helpVox[oper[1:]]
                else:
                    voxelModel = voxelModel * float(oper[1:])
            elif oper[0] == '/':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                if helpVox.has_key(oper[1:]):
                    voxelModel = voxelModel / helpVox[oper[1:]]
                else:
                    voxelModel = voxelModel / float(oper[1:])
            elif oper[0] == '+':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                if helpVox.has_key(oper[1:]):
                    voxelModel = voxelModel + helpVox[oper[1:]]
                else:
                    voxelModel = voxelModel + float(oper[1:])
            elif oper[0] == '-':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                if helpVox.has_key(oper[1:]):
                    voxelModel = voxelModel - helpVox[oper[1:]]
                else:
                    voxelModel = voxelModel - float(oper[1:])
            elif oper[0] == '<':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                key = oper[1:]
                helpVox[key] = numpy.zeros((nz, ny, nx), numpy.float32) + 1.0
                helpVox[key] = helpVox[key] * voxelModel
            elif oper[0] == '>':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                if oper[1] == '>':
                    filename = oper[2:]
                    if filename.find('.xml') > -1:
                        voxelModel, additData = self.readXmlFile(filename)
                    elif filename.find('.mhd') > -1:
                        voxelModel, additData = self.readMhdFile(filename)
                    elif filename.find('.nhdr') > -1:
                        voxelModel, additData = self.readNhdrFile(filename)
                    elif filename.find('.isq') > -1:
                        voxelModel, additData = self.readIsqFile(filename)
                    else:
                        stdout.write('\n **ERROR** : computeEquation(): Filetype not supported!\n')
                        stdout.flush()
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                else:
                    key = oper[1:]
                    voxelModel = helpVox[key]
            elif oper[0] == '^':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                voxelModel = numpy.power(voxelModel, float(oper[1:]))
            elif oper[0] == '&':
                function = oper[1:]
                if function.lower().find('sin') == 0:
                    voxelModel = numpy.sin(voxelModel)
                if function.lower().find('cos') == 0:
                    voxelModel = numpy.cos(voxelModel)
                if function.lower().find('tan') == 0:
                    voxelModel = numpy.tan(voxelModel)
                if function.lower().find('arcsin') == 0:
                    voxelModel = numpy.arcsin(voxelModel)
                if function.lower().find('arccos') == 0:
                    voxelModel = numpy.arccos(voxelModel)
                if function.lower().find('arctan') == 0:
                    voxelModel = numpy.arctan(voxelModel)
            elif oper[0] == '?':
                if echo > 0:
                    stdout.write('     -> %s\n' % oper)
                    stdout.flush()
                oper = oper.replace('==', '~')
                cond, newValue = oper[1:].split('=')
                checkValue = float(cond[1:])
                newValue = float(newValue)
                if cond[0] == '<':
                    voxelModel = numpy.greater_equal(voxelModel, checkValue) * voxelModel + numpy.less(voxelModel, checkValue) * float(newValue)
                elif cond[0] == '>':
                    voxelModel = numpy.less_equal(voxelModel, checkValue) * voxelModel + numpy.greater(voxelModel, checkValue) * float(newValue)
                elif cond[0] == '~':
                    prec = 1e-05
                    voxelModel = voxelModel + numpy.less_equal(voxelModel, checkValue + prec) * numpy.greater_equal(voxelModel, checkValue - prec) * float(newValue - checkValue)
                else:
                    stdout.write('\n **ERROR** : computeEquation(): Operator %s in logical equation not known!\n' % cond[0])
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
            else:
                stdout.write('\n **ERROR** : computeEquation(): Operator %s not known!\n' % oper[0])
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            voxelModel = self.castType(voxelModel, 'f')

        return voxelModel

    def rotateModel(self, voxelModel, rotList, additData, echo=True):
        """
        Rotate the voxel model, currently only isotropic resolutions
        
        @param  voxelModel: voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ...u voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param  rotList: list of arithmetic operations
             - TYPE: list[0] = type, list[1]= Option1, list[1]= Option2
        @param ldim: voxel dimensions in x,y,z
               - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
               - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        @return:
            voxelModel: modified voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel, 0..25
        """
        nx = voxelModel.shape[2]
        ny = voxelModel.shape[1]
        nz = voxelModel.shape[0]
        stdout.write(' ... rotate model \n')
        stdout.flush()
        ldim = additData['-ldim']
        if int(ldim[0] * 1000000) != int(ldim[1] * 1000000) or int(ldim[0] * 1000000) != int(ldim[2] * 1000000):
            stdout.write('\n **ERROR** : rotateModel(): Only isotropic resolutions implemented!\n')
            stdout.flush()
            stdout.write('             -ldim = %s \n' % repr(ldim))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        intpol = None
        debug = 0
        if rotList[0] == 'angle':
            ptN = numpy.zeros(3, numpy.float32)
            ptO = numpy.zeros(3, numpy.float32)
            ROf = numpy.array([int(rotList[3]), int(rotList[2]), int(rotList[1])])
            RAf = numpy.array([-float(rotList[4]), -float(rotList[5]), -float(rotList[6])])
            a1 = RAf[0] * numpy.pi / 180.0
            a2 = RAf[1] * numpy.pi / 180.0
            a3 = RAf[2] * numpy.pi / 180.0
            sa1 = numpy.sin(a1)
            sa2 = numpy.sin(a2)
            sa3 = numpy.sin(a3)
            ca1 = numpy.cos(a1)
            ca2 = numpy.cos(a2)
            ca3 = numpy.cos(a3)
            R = numpy.zeros((3, 3, 3), numpy.float32)
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
            Rot = numpy.dot(R[ROf[2] - 1], numpy.dot(R[ROf[1] - 1], R[ROf[0] - 1]))
            intpol = rotList[7].upper()
        elif rotList[0] == 'matrix':
            Rot = numpy.matrix([[float(rotList[1]), float(rotList[2]), float(rotList[3])],
             [
              float(rotList[4]), float(rotList[5]), float(rotList[6])],
             [
              float(rotList[7]), float(rotList[8]), float(rotList[9])]])
            intpol = rotList[10].upper()
        elif rotList[0] == 'file':
            exponent = 1
            tFilename = rotList[1]
            myTransformer = dpTransform.Transformer()
            myTransformer.readTransformMatrix(tFilename)
            myTransformer.setTransformDirection(exponent)
            M = myTransformer.getTransformMatrix4()
            Rot = numpy.matrix([[M[(0, 0)], M[(0, 1)], M[(0, 2)]],
             [
              M[(1, 0)], M[(1, 1)], M[(1, 2)]],
             [
              M[(2, 0)], M[(2, 1)], M[(2, 2)]]])
            intpol = rotList[2].upper()
        else:
            stdout.write('\n **ERROR** : mic.rotateModel(): Type %s not known!\n' % rotList[0])
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if echo:
            stdout.write('     -> Rotation matrix\n')
            stdout.write('       %s\n' % Rot[0, :])
            stdout.write('       %s\n' % Rot[1, :])
            stdout.write('       %s\n' % Rot[2, :])
        n1 = numpy.matrix([nx, 0.0, 0.0])
        n2 = numpy.matrix([0.0, ny, 0.0])
        n3 = numpy.matrix([0.0, 0.0, nz])
        n1Ra = numpy.array(numpy.transpose(Rot * n1.T))
        n2Ra = numpy.array(numpy.transpose(Rot * n2.T))
        n3Ra = numpy.array(numpy.transpose(Rot * n3.T))
        n1R = n1Ra[0, :]
        n2R = n2Ra[0, :]
        n3R = n3Ra[0, :]
        n12R = n1R + n2R
        n13R = n1R + n3R
        n23R = n2R + n3R
        n123R = n1R + n2R + n3R
        if debug > 0:
            stdout.write('     -> Rotated basis\n')
            print '        n1R  =', n1R
            print '        n2R  =', n2R
            print '        n3R  =', n3R
            print '        n12R =', n12R
            print '        n13R =', n13R
            print '        n23R =', n23R
            print '        n123R=', n123R
        nxmax = max([n1R[0], n2R[0], n3R[0], n12R[0], n13R[0], n23R[0], n123R[0], 0.0])
        nxmin = min([n1R[0], n2R[0], n3R[0], n12R[0], n13R[0], n23R[0], n123R[0], 0.0])
        nymax = max([n1R[1], n2R[1], n3R[1], n12R[1], n13R[1], n23R[1], n123R[1], 0.0])
        nymin = min([n1R[1], n2R[1], n3R[1], n12R[1], n13R[1], n23R[1], n123R[1], 0.0])
        nzmax = max([n1R[2], n2R[2], n3R[2], n12R[2], n13R[2], n23R[2], n123R[2], 0.0])
        nzmin = min([n1R[2], n2R[2], n3R[2], n12R[2], n13R[2], n23R[2], n123R[2], 0.0])
        if debug > 0:
            stdout.write('     -> Bounding Box (voxels)\n')
            print '        nxmin =', nxmin
            print '        nxmax =', nxmax
            print '        nymin =', nymin
            print '        nymax =', nymax
            print '        nzmin =', nzmin
            print '        nzmax =', nzmax
        sVec = numpy.zeros(3, numpy.float32)
        if nxmin < 0.0:
            sVec[0] = nxmin * -1.0
        if nymin < 0.0:
            sVec[1] = nymin * -1.0
        if nzmin < 0.0:
            sVec[2] = nzmin * -1.0
        if debug > 0:
            stdout.write('     -> Shift vector (voxels)\n')
            print '        x =', sVec[0]
            print '        y =', sVec[1]
            print '        z =', sVec[2]
        offOld = numpy.array(additData['Offset'])
        additData['Offset'] = [offOld[0] - sVec[0] * ldim[0], offOld[1] - sVec[1] * ldim[1], offOld[2] - sVec[2] * ldim[2]]
        additData['TransformMatrix'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        nxN = int(nxmax + sVec[0] + 0.999999999)
        nyN = int(nymax + sVec[1] + 0.999999999)
        nzN = int(nzmax + sVec[2] + 0.999999999)
        if debug > 0:
            stdout.write('     -> Corrected bounding box\n')
            print '        nxN =', nxN
            print '        nyN =', nyN
            print '        nzN =', nzN
        newVoxelModel = numpy.zeros((nzN, nyN, nxN), numpy.float32)
        if intpol == 'YES':
            newVoxelModel = micf77.rotatemodel(voxelModel, newVoxelModel, sVec, Rot, 1, 2, nxN, nyN, nzN, nx, ny, nz)
        else:
            newVoxelModel = micf77.rotatemodel(voxelModel, newVoxelModel, sVec, Rot, 0, 2, nxN, nyN, nzN, nx, ny, nz)
        del voxelModel
        return newVoxelModel

    def rotateModel_old(self, voxelModel, rotList, ldim):
        """
        Rotate the voxel model, currently only isotropic resolutions (by Jaruan)
        
        @param  voxelModel: voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ...u voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param  rotList: list of arithmetic operations
             - TYPE: list[0] = type, list[1]= Option1, list[1]= Option2
        @param ldim: voxel dimensions in x,y,z
               - TYPE: list[0] = lenX, list[1] = lenY, list[2] = lenZ
               - float lenX, lenY, lenZ ... physical voxel dimensions in x,y,z
        
        @return:
            voxelModel: modified voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel, 0..25
        """
        nx = voxelModel.shape[2]
        ny = voxelModel.shape[1]
        nz = voxelModel.shape[0]
        stdout.write(' ... rotate model \n')
        stdout.flush()
        if int(ldim[0] * 1000000) != int(ldim[1] * 1000000) or int(ldim[0] * 1000000) != int(ldim[2] * 1000000):
            stdout.write('\n **ERROR** : rotateModel(): Only isotropic resolutions implemented!\n')
            stdout.flush()
            stdout.write('             -ldim = %s \n' % repr(ldim))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if rotList[0] == 'angle':
            ptN = numpy.zeros(3, numpy.float32)
            ptO = numpy.zeros(3, numpy.float32)
            ROf = numpy.array([int(rotList[3]), int(rotList[2]), int(rotList[1])])
            RAf = numpy.array([-float(rotList[4]), -float(rotList[5]), -float(rotList[6])])
            pi = 3.141592654
            a1 = RAf[0] * pi / 180.0
            a2 = RAf[1] * pi / 180.0
            a3 = RAf[2] * pi / 180.0
            sa1 = numpy.sin(a1)
            sa2 = numpy.sin(a2)
            sa3 = numpy.sin(a3)
            ca1 = numpy.cos(a1)
            ca2 = numpy.cos(a2)
            ca3 = numpy.cos(a3)
            R = numpy.zeros((3, 3, 3), numpy.float32)
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
            Rot = numpy.dot(R[ROf[2] - 1], numpy.dot(R[ROf[1] - 1], R[ROf[0] - 1]))
            n1 = numpy.matrix([nx, 0.0, 0.0])
            n2 = numpy.matrix([0.0, ny, 0.0])
            n3 = numpy.matrix([0.0, 0.0, nz])
            n1Ra = numpy.array(numpy.transpose(Rot * n1.T))
            n2Ra = numpy.array(numpy.transpose(Rot * n2.T))
            n3Ra = numpy.array(numpy.transpose(Rot * n3.T))
            n1R = n1Ra[0, :]
            n2R = n2Ra[0, :]
            n3R = n3Ra[0, :]
            n12R = n1R + n2R
            n13R = n1R + n3R
            n23R = n2R + n3R
            n123R = n1R + n2R + n3R
            nxmax = max([n1R[0], n2R[0], n3R[0], n12R[0], n13R[0], n23R[0], n123R[0]])
            nxmin = min([n1R[0], n2R[0], n3R[0], n12R[0], n13R[0], n23R[0], n123R[0]])
            nymax = max([n1R[1], n2R[1], n3R[1], n12R[1], n13R[1], n23R[1], n123R[1]])
            nymin = min([n1R[1], n2R[1], n3R[1], n12R[1], n13R[1], n23R[1], n123R[1]])
            nzmax = max([n1R[2], n2R[2], n3R[2], n12R[2], n13R[2], n23R[2], n123R[2]])
            nzmin = min([n1R[2], n2R[2], n3R[2], n12R[2], n13R[2], n23R[2], n123R[2]])
            print '---Bounding Box---'
            print nxmin, nxmax, nymin, nymax, nzmin, nzmax
            sVec = numpy.zeros(3, numpy.float32)
            if nxmin < 0.0:
                sVec[0] = nxmin * -1.0
            if nymin < 0.0:
                sVec[1] = nymin * -1.0
            if nzmin < 0.0:
                sVec[2] = nzmin * -1.0
            print '---Shift vector---'
            print sVec
            nxN = int(nxmax + sVec[0] + 0.999999999)
            nyN = int(nymax + sVec[1] + 0.999999999)
            nzN = int(nzmax + sVec[2] + 0.999999999)
            print '---Corr Bounding Box---'
            print nxN, nyN, nzN
            newVoxelModel = numpy.zeros((nzN, nyN, nxN), numpy.float32)
            intpol = rotList[7].upper()
            if intpol == 'YES':
                newVoxelModel = micf77.rotatemodel(voxelModel, newVoxelModel, sVec, Rot, 1, 2, nxN, nyN, nzN, nx, ny, nz)
            else:
                newVoxelModel = micf77.rotatemodel(voxelModel, newVoxelModel, sVec, Rot, 0, 2, nxN, nyN, nzN, nx, ny, nz)
        elif rotList[0] == 'matrix':
            n1 = numpy.matrix([nx, 0.0, 0.0])
            n2 = numpy.matrix([0.0, ny, 0.0])
            n3 = numpy.matrix([0.0, 0.0, nz])
            Rot = numpy.matrix([[float(rotList[1]), float(rotList[2]), float(rotList[3])],
             [
              float(rotList[4]), float(rotList[5]), float(rotList[6])],
             [
              float(rotList[7]), float(rotList[8]), float(rotList[9])]])
            n1Ra = numpy.array(numpy.transpose(Rot * n1.T))
            n2Ra = numpy.array(numpy.transpose(Rot * n2.T))
            n3Ra = numpy.array(numpy.transpose(Rot * n3.T))
            n1R = n1Ra[0, :]
            n2R = n2Ra[0, :]
            n3R = n3Ra[0, :]
            n12R = n1R + n2R
            n13R = n1R + n3R
            n23R = n2R + n3R
            n123R = n1R + n2R + n3R
            nxmax = max([n1R[0], n2R[0], n3R[0], n12R[0], n13R[0], n23R[0], n123R[0]])
            nxmin = min([n1R[0], n2R[0], n3R[0], n12R[0], n13R[0], n23R[0], n123R[0]])
            nymax = max([n1R[1], n2R[1], n3R[1], n12R[1], n13R[1], n23R[1], n123R[1]])
            nymin = min([n1R[1], n2R[1], n3R[1], n12R[1], n13R[1], n23R[1], n123R[1]])
            nzmax = max([n1R[2], n2R[2], n3R[2], n12R[2], n13R[2], n23R[2], n123R[2]])
            nzmin = min([n1R[2], n2R[2], n3R[2], n12R[2], n13R[2], n23R[2], n123R[2]])
            print '---Bounding Box---'
            print nxmin, nxmax, nymin, nymax, nzmin, nzmax
            sVec = numpy.zeros(3, numpy.float32)
            if nxmin < 0.0:
                sVec[0] = nxmin * -1.0
            if nymin < 0.0:
                sVec[1] = nymin * -1.0
            if nzmin < 0.0:
                sVec[2] = nzmin * -1.0
            print '---Shift vector---'
            print sVec
            nxN = int(nxmax + sVec[0] + 0.999999999)
            nyN = int(nymax + sVec[1] + 0.999999999)
            nzN = int(nzmax + sVec[2] + 0.999999999)
            print '---Corr Bounding Box---'
            print nxN, nyN, nzN
            newVoxelModel = numpy.zeros((nzN, nyN, nxN), numpy.float32)
            intpol = rotList[10].upper()
            if intpol == 'YES':
                newVoxelModel = micf77.rotatemodel(voxelModel, newVoxelModel, sVec, Rot, 1, 2, nxN, nyN, nzN, nx, ny, nz)
            else:
                newVoxelModel = micf77.rotatemodel(voxelModel, newVoxelModel, sVec, Rot, 0, 2, nxN, nyN, nzN, nx, ny, nz)
        else:
            stdout.write('\n **ERROR** : mic.rotateModel(): Type %s not known!\n' % rotList[0])
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        del voxelModel
        return newVoxelModel

    def mirrorModel(self, voxelModel, axes):
        """
        Function mirrors the given voxel model and returns a new voxel model.
        
        
        @param  voxelModel: voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel,
        @param  axes: mirror axes 
          - TYPE: list[0] = axis1, list[1] = axis2, ...
          - string direction ... direction of mirror plane (x,y, or z)
        
        @return:
         voxelModel: mirrored voxel model
          - TYPE: numpy.array[kZ, jY, iX] = grayValue
          - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
          - int grayValue  ... gray value of voxel, 0..25
        """
        stdout.write(' ... mirror model   \n')
        stdout.flush()
        nx, ny, nz = self.get_Shape(voxelModel)
        auxVoxelModel = self.createVoxelModel(nx, ny, nz, 'f')
        auxVoxelModel = voxelModel
        for axis in axes:
            if axis.upper() == 'X':
                newVoxelModel = self.createVoxelModel(2 * nx, ny, nz, 'f')
                dpUtils.progressStart('     -> progress mirror X   : ')
                for k in range(nz):
                    progress = float(k + 1) / float(nz) * 10.0
                    for j in range(ny):
                        for i in range(nx):
                            newVoxelModel[k, j, i] = auxVoxelModel[k, j, i]
                            newVoxelModel[k, j, 2 * nx - 1 - i] = auxVoxelModel[k, j, i]

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
            elif axis.upper() == 'Y':
                newVoxelModel = self.createVoxelModel(nx, 2 * ny, nz, 'f')
                dpUtils.progressStart('     -> progress mirror Y   : ')
                for k in range(nz):
                    progress = float(k + 1) / float(nz) * 10.0
                    for j in range(ny):
                        for i in range(nx):
                            newVoxelModel[k, j, i] = auxVoxelModel[k, j, i]
                            newVoxelModel[k, 2 * ny - 1 - j, i] = auxVoxelModel[k, j, i]

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
            elif axis.upper() == 'Z':
                newVoxelModel = self.createVoxelModel(nx, ny, 2 * nz, 'f')
                dpUtils.progressStart('     -> progress mirror Z   : ')
                for k in range(nz):
                    progress = float(k + 1) / float(nz) * 10.0
                    for j in range(ny):
                        for i in range(nx):
                            newVoxelModel[k, j, i] = auxVoxelModel[k, j, i]
                            newVoxelModel[2 * nz - 1 - k, j, i] = auxVoxelModel[k, j, i]

                    dpUtils.progressNext(progress)

                dpUtils.progressEnd()
            else:
                stdout.write('\n **ERROR** : mirrorModel(): axis=%s not supported!\n' % axis)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            nx, ny, nz = self.get_Shape(newVoxelModel)
            del auxVoxelModel
            auxVoxelModel = self.createVoxelModel(nx, ny, nz, 'f')
            auxVoxelModel = newVoxelModel

        del auxVoxelModel
        return newVoxelModel

    def mask(self, voxelModel, maskModel):
        """
        Functions masks a given voxel model by using a mask voxel model.
        
        @param  voxelModel: voxel model with old scale
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param  maskModel: voxel model with mask. Should contain only 2 different
           gray values. Will be scaled automatically to 0/1.
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
               
        @return:
           voxelModel: scaled voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel, 0..25
        """
        stdout.write(' ... mask Model \n')
        stdout.flush()
        minVox, maxVox = self.computeNumpyMinMax(maskModel, 0)
        if int(minVox) == 0 and int(maxVox) == 1:
            pass
        else:
            stdout.write('\n **ERROR** : mask() array min..max is not 0..1, min=%i  max=%i!\n' % (minVox, maxVox))
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        vnx, vny, vnz = self.get_Shape(voxelModel)
        mnx, mny, mnz = self.get_Shape(maskModel)
        if vnx == mnx and vny == mny and vnz == mnz:
            stdout.write('     -> compute\n')
            stdout.flush()
            voxelModel = voxelModel * maskModel
        elif vnx % mnx == 0 and vny % mny == 0 and vnz % mnz == 0 and vnx / mnx == vny / mny and vnx / mnx == vnz / mnz:
            scale = vnx / mnx
            dpUtils.progressStart('     -> Processed Data      : ')
            for k in range(mnz):
                progress = (k + 1) / float(mnz) * 10.0
                for j in range(mny):
                    for i in range(mnx):
                        if maskModel[k, j, i] == 0:
                            iscale = i * scale
                            jscale = j * scale
                            kscale = k * scale
                            for nk in range(scale):
                                for nj in range(scale):
                                    for ni in range(scale):
                                        voxelModel[kscale + nk, jscale + nj, iscale + ni] = 0

                dpUtils.progressNext(progress)

            dpUtils.progressEnd()
        else:
            stdout.write('\n **ERROR** : mask(): Array dimension not compatible!\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return voxelModel

    def combineModels(self, voxelModel, additData, operator, combineFile):
        """
        Functions combines two voxel models with given boolean operator.
        """
        stdout.write(' ... combine models \n')
        stdout.flush()
        vOFF = [
         0.0, 0.0, 0.0]
        if additData.has_key('Offset'):
            vOFF = additData['Offset']
        vnx, vny, vnz = self.get_Shape(voxelModel)
        vlx, vly, vlz = additData['-ldim']
        combineModel, additData2 = myModel.read(combineFile)
        cOFF = [0.0, 0.0, 0.0]
        if additData2.has_key('Offset'):
            cOFF = additData2['Offset']
        cnx, cny, cnz = self.get_Shape(combineModel)
        clx, cly, clz = additData2['-ldim']
        prec = 1e-06
        if abs(vlx - clx) > prec or abs(vly - cly) > prec or abs(vlz - clz) > prec:
            stdout.write('\n **ERROR** : combineModels(): Voxels do not have the same dimensions!\n')
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        dnx = int((cOFF[0] - vOFF[0]) / vlx)
        dny = int((cOFF[1] - vOFF[1]) / vly)
        dnz = int((cOFF[2] - vOFF[2]) / vlz)
        if operator == '+' or operator == '-' or operator == '*' or operator == 'r':
            voxelModel = micf77.combinemodels(voxelModel, combineModel, dnx, dny, dnz, operator, 2)
        else:
            stdout.write("\n **ERROR** : combineModels(): Operator '%s' not known!\n" % operator)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return voxelModel

    def combineMeshModels(self, voxelModel, additData, meshFile, elsetName, grayValue=None, createNew=None):
        """
        Functions combines a mesh and voxel model
        """
        stdout.write(' ... combine mesh+model \n')
        stdout.flush()
        fecModel = fec.fec()
        filename, ext = os.path.splitext(meshFile)
        if ext.upper().find('INP') > -1:
            title, nodes, nsets, elems, elsets, properties = fecModel.readAbaqus(meshFile, props=True)
        elif ext.upper().find('FEM') > -1:
            title, nodes, nsets, elems, elsets, properties = fecModel.readFem(meshFile, props=True)
        elif ext.upper().find('VTK') > -1:
            title, nodes, nsets, elems, elsets, properties = fecModel.readVTK(meshFile, props=True)
        else:
            print ' **ERROR** mic().combineMeshModels: intput file extension of file: "%s" not known!\n\n' % meshFile
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        stdout.write('     -> Processed Nodes     : %10i        \n' % len(nodes))
        stdout.flush()
        stdout.write('     -> Processed Elements  : %10i        \n' % len(elems))
        stdout.flush()
        stdout.write('     -> Processed Node Sets : %10i        \n' % len(nsets))
        stdout.flush()
        stdout.write('     -> Processed Elems Set : %10i        \n' % len(elsets))
        stdout.flush()
        elsetList = []
        elsetList2 = dpUtils.userSplit(elsetName)
        if elsetList2.count('ALL') > 0:
            elsetList = elsets.keys()
        else:
            elsetList = elsetList2
        for curElsetName in elsetList:
            if elsets.has_key(curElsetName) or curElsetName.upper() == 'ALL':
                continue
                stdout.write("\n **ERROR** : combineMeshModels(): ELSET '%s' not found!\n" % curElsetName)
                stdout.write('             elsets=%s\n' % elsets.keys())
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)

        vOFF = numpy.array([0.0, 0.0, 0.0])
        if not createNew:
            if additData.has_key('Offset'):
                vOFF = numpy.array(additData['Offset'])
            nx, ny, nz = self.get_Shape(voxelModel)
            lx, ly, lz = additData['-ldim']
        else:
            lx = float(createNew)
            ly = float(createNew)
            lz = float(createNew)
            xList = []
            yList = []
            zList = []
            for nid in nodes:
                x, y, z = nodes[nid].get_coord()
                xList.append(x)
                yList.append(y)
                zList.append(z)

            minx = numpy.min(xList)
            maxx = numpy.max(xList)
            miny = numpy.min(yList)
            maxy = numpy.max(yList)
            minz = numpy.min(zList)
            maxz = numpy.max(zList)
            radiusList = [
             0.0]
            for curElsetName in elsetList:
                for elid in elsets[curElsetName]:
                    eltype = elems[elid].get_type()
                    if eltype == 'bar2':
                        if len(properties[curElsetName]['geometry']) == 1:
                            radius = properties[curElsetName]['geometry'][0]
                        elif len(properties[curElsetName]['geometry']) == 2:
                            radius1 = properties[curElsetName]['geometry'][0]
                            radius2 = properties[curElsetName]['geometry'][1]
                            radius = numpy.max(properties[curElsetName]['geometry'])
                        else:
                            stdout.write('\n **ERROR** : mic().combineMeshModels(): one or two radius required in case of bar2 elements!\n')
                            stdout.write('\n E N D E D  with ERRORS \n\n')
                            stdout.flush()
                            exit(1)
                        radiusList.append(radius)
                    if eltype == 'tria3' or eltype == 'quad4':
                        thickness = properties[curElsetName]['thickness']
                        radiusList.append(thickness)

            coverVox = int(numpy.max(radiusList) / lx + 0.5)
            nx = int((maxx - minx) / lx + 0.5) + 2 * coverVox
            ny = int((maxy - miny) / ly + 0.5) + 2 * coverVox
            nz = int((maxz - minz) / lz + 0.5) + 2 * coverVox
            additData['Offset'] = [minx - coverVox * lx, miny - coverVox * ly, minz - coverVox * lz]
            additData['TransformMatrix'] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
            additData['-ldim'] = [lx, ly, lz]
            additData['-ndim'] = [nx, ny, nz]
            additData['ElementSpacing'] = [lx, ly, lz]
            additData['DimSize'] = [nx, ny, nz]
            vOFF = numpy.array(additData['Offset'])
            del voxelModel
            voxelModel = self.createVoxelModel(nx, ny, nz, 'f')
        minVox, maxVox = self.computeNumpyMinMax(voxelModel, 0)
        dpUtils.progressStart('     -> Voxelized elements  : ')
        count = 0
        listlen = len(elsetList)
        madId = {}
        matCount = int(maxVox) + 1
        isSmallCount = 0
        for curElsetName in elsetList:
            count += 1
            progress = float(count) / float(listlen) * 10.0
            dpUtils.progressNext(progress)
            solidOut = 0
            beamOut = 0
            for elid in elsets[curElsetName]:
                eltype = elems[elid].get_type()
                pid = properties[curElsetName]['material']
                if pid not in madId:
                    madId[pid] = matCount
                    matCount += 1
                if grayValue == None:
                    curGrayValue = madId[pid]
                else:
                    curGrayValue = int(grayValue)
                if curGrayValue >= minVox and curGrayValue <= maxVox:
                    stdout.write('\n ** WARNING **: The given gray value is within the range of \n ')
                    stdout.write('               gray values of the voxel! \n\n ')
                    stdout.flush()
                if eltype == 'bar2':
                    nlist = elems[elid].get_nodes()
                    n1 = nodes[nlist[0]].get_coord_numpy() - vOFF
                    n2 = nodes[nlist[1]].get_coord_numpy() - vOFF
                    minP = numpy.array([min(n1[0], n2[0]), min(n1[1], n2[1]), min(n1[2], n2[2])])
                    maxP = numpy.array([max(n1[0], n2[0]), max(n1[1], n2[1]), max(n1[2], n2[2])])
                    radius = None
                    radius1 = None
                    radius2 = None
                    if len(properties[curElsetName]['geometry']) == 1:
                        radius = properties[curElsetName]['geometry'][0]
                    elif len(properties[curElsetName]['geometry']) == 2:
                        radius1 = properties[curElsetName]['geometry'][0]
                        radius2 = properties[curElsetName]['geometry'][1]
                        radius = numpy.max(properties[curElsetName]['geometry'])
                    else:
                        stdout.write('\n **ERROR** : mic().combineMeshModels(): one or two radius required in case of bar2 elements!\n')
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                    minP -= numpy.array([radius, radius, radius])
                    maxP += numpy.array([radius, radius, radius])
                    bbox = {}
                    bbox['xmin'] = minP[0]
                    bbox['ymin'] = minP[1]
                    bbox['zmin'] = minP[2]
                    bbox['xmax'] = maxP[0]
                    bbox['ymax'] = maxP[1]
                    bbox['zmax'] = maxP[2]
                    nx0 = int(bbox['xmin'] / lx)
                    ny0 = int(bbox['ymin'] / ly)
                    nz0 = int(bbox['zmin'] / lz)
                    dnx = int((bbox['xmax'] - bbox['xmin']) / lx + 1.5)
                    dny = int((bbox['ymax'] - bbox['ymin']) / ly + 1.5)
                    dnz = int((bbox['zmax'] - bbox['zmin']) / lz + 1.5)
                    outsideFlag = False
                    if nx0 < 0:
                        nx0 = 0
                        outsideFlag = True
                    if ny0 < 0:
                        ny0 = 0
                        outsideFlag = True
                    if nz0 < 0:
                        nz0 = 0
                        outsideFlag = True
                    if nx0 + dnx > nx:
                        dnx = nx - nx0
                        outsideFlag = True
                    if ny0 + dny > ny:
                        dny = ny - ny0
                        outsideFlag = True
                    if nz0 + dnz > nz:
                        dnz = nz - nz0
                        outsideFlag = True
                    if outsideFlag == True:
                        beamOut += 1
                    if radius1 == None and radius2 == None:
                        radius1 = radius
                        radius2 = radius
                    if radius1 > lx or radius2 > lx:
                        voxelModel = miaf77.ispointinsideconecapbbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, radius1 * radius1, radius2 * radius2)
                    else:
                        isSmallCount += 1
                elif eltype == 'tria3' or eltype == 'quad4':
                    nlist = elems[elid].get_nodes()
                    n1h = nodes[nlist[0]].get_coord_numpy() - vOFF
                    n2h = nodes[nlist[1]].get_coord_numpy() - vOFF
                    n3h = nodes[nlist[2]].get_coord_numpy() - vOFF
                    if eltype == 'quad4':
                        n4h = nodes[nlist[3]].get_coord_numpy() - vOFF
                    vec1 = n2h - n1h
                    vec2 = n3h - n1h
                    nvec = dpTensor.UnitVector(dpTensor.CrossProduct(vec1, vec2))
                    thickness = properties[curElsetName]['thickness']
                    nvec = nvec * thickness
                    if eltype == 'tria3':
                        n1 = n1h - nvec
                        n4 = n1h + nvec
                        n2 = n2h - nvec
                        n5 = n2h + nvec
                        n3 = n3h - nvec
                        n6 = n3h + nvec
                        xyzList = [
                         n1, n2, n3, n4, n5, n6]
                    elif eltype == 'quad4':
                        n1 = n1h - nvec
                        n5 = n1h + nvec
                        n2 = n2h - nvec
                        n6 = n2h + nvec
                        n3 = n3h - nvec
                        n7 = n3h + nvec
                        n4 = n4h - nvec
                        n8 = n4h + nvec
                        xyzList = [
                         n1, n2, n3, n4, n5, n6, n7, n8]
                    bbox = {}
                    bbox['xmin'] = 1e+308
                    bbox['ymin'] = 1e+308
                    bbox['zmin'] = 1e+308
                    bbox['xmax'] = -1e+308
                    bbox['ymax'] = -1e+308
                    bbox['zmax'] = -1e+308
                    for xyz in xyzList:
                        x, y, z = xyz
                        if x < bbox['xmin']:
                            bbox['xmin'] = x
                        if x > bbox['xmax']:
                            bbox['xmax'] = x
                        if y < bbox['ymin']:
                            bbox['ymin'] = y
                        if y > bbox['ymax']:
                            bbox['ymax'] = y
                        if z < bbox['zmin']:
                            bbox['zmin'] = z
                        if z > bbox['zmax']:
                            bbox['zmax'] = z

                    nx0 = int(bbox['xmin'] / lx)
                    ny0 = int(bbox['ymin'] / ly)
                    nz0 = int(bbox['zmin'] / lz)
                    dnx = int((bbox['xmax'] - bbox['xmin']) / lx + 0.5)
                    dny = int((bbox['ymax'] - bbox['ymin']) / ly + 0.5)
                    dnz = int((bbox['zmax'] - bbox['zmin']) / lz + 0.5)
                    outsideFlag = False
                    if nx0 < 0:
                        nx0 = 0
                        outsideFlag = True
                    if ny0 < 0:
                        ny0 = 0
                        outsideFlag = True
                    if nz0 < 0:
                        nz0 = 0
                        outsideFlag = True
                    if nx0 + dnx > nx:
                        dnx = nx - nx0
                        outsideFlag = True
                    if ny0 + dny > ny:
                        dny = ny - ny0
                        outsideFlag = True
                    if nz0 + dnz > nz:
                        dnz = nz - nz0
                        outsideFlag = True
                    if outsideFlag == True:
                        solidOut += 1
                    if eltype == 'tria3':
                        voxelModel = miaf77.ispointinsidewedgebbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, n3, n4, n5, n6)
                    if eltype == 'quad4':
                        voxelModel = miaf77.ispointinsidehexbbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, n3, n4, n5, n6, n7, n8)
                elif eltype == 'tetra4' or eltype == 'hexa8' or eltype == 'penta6' or eltype == 'pyra5':
                    nlist = elems[elid].get_nodes()
                    n1 = nodes[nlist[0]].get_coord_numpy() - vOFF
                    n2 = nodes[nlist[1]].get_coord_numpy() - vOFF
                    n3 = nodes[nlist[2]].get_coord_numpy() - vOFF
                    n4 = nodes[nlist[3]].get_coord_numpy() - vOFF
                    if eltype == 'pyra5':
                        n5 = nodes[nlist[4]].get_coord_numpy() - vOFF
                    if eltype == 'hexa8' or eltype == 'penta6':
                        n5 = nodes[nlist[4]].get_coord_numpy() - vOFF
                        n6 = nodes[nlist[5]].get_coord_numpy() - vOFF
                    if eltype == 'hexa8':
                        n7 = nodes[nlist[6]].get_coord_numpy() - vOFF
                        n8 = nodes[nlist[7]].get_coord_numpy() - vOFF
                    bbox = {}
                    bbox['xmin'] = 1e+308
                    bbox['ymin'] = 1e+308
                    bbox['zmin'] = 1e+308
                    bbox['xmax'] = -1e+308
                    bbox['ymax'] = -1e+308
                    bbox['zmax'] = -1e+308
                    for nid in nlist:
                        x, y, z = numpy.array(nodes[nid].get_coord_numpy()) - vOFF
                        if x < bbox['xmin']:
                            bbox['xmin'] = x
                        if x > bbox['xmax']:
                            bbox['xmax'] = x
                        if y < bbox['ymin']:
                            bbox['ymin'] = y
                        if y > bbox['ymax']:
                            bbox['ymax'] = y
                        if z < bbox['zmin']:
                            bbox['zmin'] = z
                        if z > bbox['zmax']:
                            bbox['zmax'] = z

                    nx0 = int(bbox['xmin'] / lx)
                    ny0 = int(bbox['ymin'] / ly)
                    nz0 = int(bbox['zmin'] / lz)
                    dnx = int((bbox['xmax'] - bbox['xmin']) / lx + 0.5)
                    dny = int((bbox['ymax'] - bbox['ymin']) / ly + 0.5)
                    dnz = int((bbox['zmax'] - bbox['zmin']) / lz + 0.5)
                    outsideFlag = False
                    if nx0 < 0:
                        nx0 = 0
                        outsideFlag = True
                    if ny0 < 0:
                        ny0 = 0
                        outsideFlag = True
                    if nz0 < 0:
                        nz0 = 0
                        outsideFlag = True
                    if nx0 + dnx > nx:
                        dnx = nx - nx0
                        outsideFlag = True
                    if ny0 + dny > ny:
                        dny = ny - ny0
                        outsideFlag = True
                    if nz0 + dnz > nz:
                        dnz = nz - nz0
                        outsideFlag = True
                    if outsideFlag == True:
                        solidOut += 1
                    if eltype == 'tetra4':
                        voxelModel = miaf77.ispointinsidetetrahedronbbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, n3, n4)
                    if eltype == 'pyra5':
                        voxelModel = miaf77.ispointinsidepyrabbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, n3, n4, n5)
                    if eltype == 'penta6':
                        voxelModel = miaf77.ispointinsidewedgebbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, n3, n4, n5, n6)
                    if eltype == 'hexa8':
                        voxelModel = miaf77.ispointinsidehexbbox(voxelModel, [nx0, ny0, nz0, dnx, dny, dnz], [lx, ly, lz], curGrayValue, n1, n2, n3, n4, n5, n6, n7, n8)
                else:
                    stdout.write("\n **ERROR** : combineMeshModels(): ELTYPE '%s' not implemented!\n" % eltype)
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

        dpUtils.progressEnd()
        if solidOut > 0:
            stdout.write('\n **WARNING**: mic.combineMeshModels: %i solid elements outside voxel model!\n\n' % solidOut)
        if beamOut > 0:
            stdout.write('\n **WARNING**: mic.combineMeshModels: %i beams outside voxel model!\n\n' % beamOut)
        if isSmallCount > 0:
            stdout.write('\n **WARNING**: mic.combineMeshModels: %i beams have a radius smaller \n              than the voxel resolution and are not printed!\n\n' % isSmallCount)
        return voxelModel

    def boundingBox(self, voxelModel, thresValue):
        """
        Functions returns the bounding box of a given voxel model.
        A threshold value is used to find the bounding box.
        
        @param  voxelModel: voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param  thresValue: All gray values >= this value will be within the
           bounding box.
             - TYPE: int
        
        @return:
            bbox: bounding box
             - TYPE: dict['minX']=minX, dict['minY']=minY, dict['minZ']=minZ
                dict['maxX']=maxX, dict['maxY']=maxY, dict['maxZ']=maxZ
             - minX, minY, minZ: minimum bounding box values starting at 0
             - maxX, maxY, maxZ: maximum bounding box values ending at number of
                voxels minus one
        """
        stdout.write(' ... compute bounding box \n')
        stdout.flush()
        stdout.write('     -> greater ...\n')
        stdout.flush()
        bvoxModel = numpy.greater(voxelModel, thresValue - 1.0)
        nx, ny, nz = self.get_Shape(voxelModel)
        minX = -1
        minY = -1
        minZ = -1
        maxX = nx
        maxY = ny
        maxZ = nz
        stdout.write('     -> find x+             : ')
        stdout.flush()
        sum = 0
        while sum == 0:
            minX += 1
            if minX == maxX:
                stdout.write('NOT ')
                break
            sum = numpy.sum(numpy.sum(bvoxModel[:, :, minX]))

        stdout.write('found\n')
        stdout.flush()
        stdout.write('     -> find y+             : ')
        stdout.flush()
        sum = 0
        while sum == 0:
            minY += 1
            if minY == maxY:
                stdout.write('NOT ')
                break
            sum = numpy.sum(numpy.sum(bvoxModel[:, minY, :]))

        stdout.write('found\n')
        stdout.flush()
        stdout.write('     -> find z+             : ')
        stdout.flush()
        sum = 0
        while sum == 0:
            minZ += 1
            if minZ == maxZ:
                stdout.write('NOT ')
                break
            sum = numpy.sum(numpy.sum(bvoxModel[minZ, :, :]))

        stdout.write('found\n')
        stdout.flush()
        sum = 0
        stdout.write('     -> find x-             : ')
        stdout.flush()
        while sum == 0:
            maxX -= 1
            if maxX == 0:
                stdout.write('NOT ')
                break
            sum = numpy.sum(numpy.sum(bvoxModel[:, :, maxX]))

        stdout.write('found\n')
        stdout.flush()
        sum = 0
        stdout.write('     -> find y-             : ')
        stdout.flush()
        while sum == 0:
            maxY -= 1
            if maxY == 0:
                stdout.write('NOT ')
                break
            sum = numpy.sum(numpy.sum(bvoxModel[:, maxY, :]))

        stdout.write('found\n')
        stdout.flush()
        sum = 0
        stdout.write('     -> find z-             : ')
        stdout.flush()
        while sum == 0:
            maxZ -= 1
            if maxY == 0:
                stdout.write('NOT ')
                break
            sum = numpy.sum(numpy.sum(bvoxModel[maxZ, :, :]))

        stdout.write('found\n')
        stdout.flush()
        bbox = {}
        bbox['minX'] = minX
        bbox['minY'] = minY
        bbox['minZ'] = minZ
        bbox['maxX'] = maxX
        bbox['maxY'] = maxY
        bbox['maxZ'] = maxZ
        stdout.write('     -> xmin ymin zmin      : %i %i %i \n' % (minX, minY, minZ))
        stdout.write('     -> xmax ymax zmax      : %i %i %i \n' % (maxX, maxY, maxZ))
        return bbox

    def cutBlockBoundingBox(self, voxelModel, thresValue, additData=None, extendList=None):
        """
        Crops a model by using a bounding box. Everthing which is
        greater than "value" is inside.
        
        @param  voxelModel: voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param  thresValue: All gray values >= this value will be within the
           bounding box.
             - TYPE: int
        @param  extendList: list of voxel which extends the bbox
             - TYPE: list[0] = px, list[1] = py, , list[1] = pz
        
        @return:
          cropped voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        """
        px = 0
        py = 0
        pz = 0
        if extendList != None:
            px = int(extendList[0])
            py = int(extendList[1])
            pz = int(extendList[2])
        bbox = self.boundingBox(voxelModel, thresValue)
        modVoxelModel = self.cutBlock(voxelModel, [bbox['minX'] - px, bbox['minY'] - py, bbox['minZ'] - pz,
         bbox['maxX'] - bbox['minX'] + 1 + 2 * px,
         bbox['maxY'] - bbox['minY'] + 1 + 2 * py,
         bbox['maxZ'] - bbox['minZ'] + 1 + 2 * pz], additData=additData)
        return modVoxelModel

    def makeEndCaps(self, voxelModel, capList, additData=None):
        """
        Function attach end caps on a voxel model in a given direction.
        The thickness of the end caps is selection by number of voxels.
        Note: Model should be scaled such that no bone voxel has
        the same grayvalue as the embedding material!
        
        @param voxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        
        @param  capList: list of cap parameters
              - TYPE: list[0] = direction, list[1] = thickness, list[2] = grayVal
              - int direction ... direction of end caps 1,2, or 3
              - int thickness ... thickness of endcaps in number of voxels
              - int grayVal   ... gray value which should be applied to these voxels
        
        @return:
           newVoxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        """
        print ' ... make end caps'
        format = 'f'
        direction = capList[0]
        thickness = int(capList[1])
        grayVal = int(capList[2])
        nx, ny, nz = self.get_Shape(voxelModel)
        oxn = oyn = ozn = 0.0
        bottom = False
        top = False
        if direction.find('-') > -1:
            bottom = True
        elif direction.find('+') > -1:
            top = True
        else:
            bottom = True
            top = True
        thickness2 = thickness
        if bottom == False:
            thickness2 = 0
        if direction.find('3') > -1:
            if bottom == True:
                ozn = thickness
            if bottom == True and top == True:
                newVoxelModel = numpy.zeros((nz + 2 * thickness, ny, nx)) + grayVal
            else:
                newVoxelModel = numpy.zeros((nz + thickness, ny, nx)) + grayVal
            newVoxelModel = self.castType(newVoxelModel, format)
            for k in range(nz):
                nK = k + thickness2
                newVoxelModel[nK] = voxelModel[k]

        elif direction.find('2') > -1:
            if bottom == True:
                oyn = thickness
            if bottom == True and top == True:
                newVoxelModel = numpy.zeros((nz, ny + 2 * thickness, nx)) + grayVal
            else:
                newVoxelModel = numpy.zeros((nz, ny + thickness, nx)) + grayVal
            newVoxelModel = self.castType(newVoxelModel, format)
            for k in range(nz):
                for j in range(ny):
                    nJ = j + thickness2
                    for i in range(nx):
                        newVoxelModel[k, nJ, i] = voxelModel[k, j, i]

        elif direction.find('1') > -1:
            if bottom == True:
                oxn = thickness
            if bottom == True and top == True:
                newVoxelModel = numpy.zeros((nz, ny, nx + 2 * thickness)) + grayVal
            else:
                newVoxelModel = numpy.zeros((nz, ny, nx + thickness)) + grayVal
            newVoxelModel = self.castType(newVoxelModel, format)
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        nI = i + thickness2
                        newVoxelModel[k, j, nI] = voxelModel[k, j, i]

        else:
            stdout.write("\n **ERROR** makeEndCaps(): Direction '%s' not implemented!\n\n" % repr(direction))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox - oxn * lx, oy - oyn * ly, oz - ozn * lz]
        return newVoxelModel

    def extendModel(self, voxelModel, extendList, additData=None):
        """
        Function surrounds whole model with a layer of material. 
        The thickness of the end caps is selection by number of voxels.
        
        @param voxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        
        @param  extendList: list of extend parameters
              - TYPE: list[0] = direction, list[1] = thickness
              - int direction ... direction of extension 1/2/3 
              - int thickness ... thickness of endcaps in number of voxels
        
        @return:
           newVoxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        """
        stdout.write(' ... extend model\n')
        stdout.flush()
        direct = extendList[0]
        thick = extendList[1]
        grayval = 0
        if len(extendList) == 3:
            grayval = extendList[2]
        nx, ny, nz = self.get_Shape(voxelModel)
        oxn = oyn = ozn = 0.0
        if direct == 3:
            newVoxelModel = numpy.zeros((nz + thick, ny, nx))
            newVoxelModel = self.castType(newVoxelModel, 'f')
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        newVoxelModel[k, j, i] = voxelModel[k, j, i]
                        if k == nz - 1:
                            for l in range(thick):
                                value = voxelModel[k, j, i]
                                if grayval > 0 and value > 0:
                                    value = grayval
                                newVoxelModel[k + l + 1, j, i] = value

        elif direct == -3:
            ozn = thick
            newVoxelModel = numpy.zeros((nz + thick, ny, nx))
            newVoxelModel = self.castType(newVoxelModel, 'f')
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        newVoxelModel[k + thick, j, i] = voxelModel[k, j, i]
                        if k == 0:
                            for l in range(thick):
                                value = voxelModel[k, j, i]
                                if grayval > 0 and value > 0:
                                    value = grayval
                                newVoxelModel[k + l, j, i] = value

        elif direct == 1:
            newVoxelModel = numpy.zeros((nz, ny, nx + thick))
            newVoxelModel = self.castType(newVoxelModel, 'f')
            for i in range(nx):
                for k in range(nz):
                    for j in range(ny):
                        newVoxelModel[k, j, i] = voxelModel[k, j, i]
                        if i == nx - 1:
                            for l in range(thick):
                                value = voxelModel[k, j, i]
                                if grayval > 0 and value > 0:
                                    value = grayval
                                newVoxelModel[k, j, i + l + 1] = value

        elif direct == -1:
            oxn = thick
            newVoxelModel = numpy.zeros((nz, ny, nx + thick))
            newVoxelModel = self.castType(newVoxelModel, 'f')
            for i in range(nx):
                for k in range(nz):
                    for j in range(ny):
                        newVoxelModel[k, j, i + thick] = voxelModel[k, j, i]
                        if i == 0:
                            for l in range(thick):
                                value = voxelModel[k, j, i]
                                if grayval > 0 and value > 0:
                                    value = grayval
                                newVoxelModel[k, j, i + l] = value

        elif direct == 2:
            newVoxelModel = numpy.zeros((nz, ny + thick, nx))
            newVoxelModel = self.castType(newVoxelModel, 'f')
            for j in range(ny):
                for i in range(nx):
                    for k in range(nz):
                        newVoxelModel[k, j, i] = voxelModel[k, j, i]
                        if j == ny - 1:
                            for l in range(thick):
                                value = voxelModel[k, j, i]
                                if grayval > 0 and value > 0:
                                    value = grayval
                                newVoxelModel[k, j + l + 1, i] = value

        elif direct == -2:
            oyn = thick
            newVoxelModel = numpy.zeros((nz, ny + thick, nx))
            newVoxelModel = self.castType(newVoxelModel, 'f')
            for j in range(ny):
                for i in range(nx):
                    for k in range(nz):
                        newVoxelModel[k, j + thick, i] = voxelModel[k, j, i]
                        if j == 0:
                            for l in range(thick):
                                value = voxelModel[k, j, i]
                                if grayval > 0 and value > 0:
                                    value = grayval
                                newVoxelModel[k, j + l, i] = value

        else:
            stdout.write("\n **ERROR** extendModel(): Direction '%i' not implemented!\n\n" % direct)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox - oxn * lx, oy - oyn * ly, oz - ozn * lz]
        return newVoxelModel

    def closeModel(self, voxelModel, closeList):
        """
        Function closes sliced voxel models by 2D fill
        """
        stdout.write(' ... close model\n')
        stdout.flush()
        direct = closeList[0]
        thres = closeList[1]
        kernel = closeList[2]
        grayVal = closeList[3]
        nx, ny, nz = self.get_Shape(voxelModel)
        inOut = 1
        valid = 3
        minThick = 1
        fbbox = numpy.array([0, nx, 0, ny, 0, nz])
        if direct == 3:
            slice0 = numpy.zeros((1, ny, nx))
            slice1 = numpy.zeros((1, ny, nx))
            for j in range(ny):
                for i in range(nx):
                    slice0[0, j, i] = voxelModel[0, j, i]
                    slice1[0, j, i] = voxelModel[nz - 1, j, i]

            plane2d = 3
            slice0 = micf77.computefill(slice0, thres, valid, inOut, kernel, minThick, plane2d, 0, kernel, fbbox)
            slice0 = micf77.computedilation2(slice0, thres, kernel, 2)
            slice0 = micf77.computeerosion(slice0, thres, kernel, 2)
            slice1 = micf77.computefill(slice1, thres, valid, inOut, kernel, minThick, plane2d, 0, kernel, fbbox)
            slice1 = micf77.computedilation2(slice1, thres, kernel, 2)
            slice1 = micf77.computeerosion(slice1, thres, kernel, 2)
            for j in range(ny):
                for i in range(nx):
                    if slice0[0, j, i] >= thres:
                        voxelModel[0, j, i] = grayVal
                    if slice1[0, j, i] >= thres:
                        voxelModel[nz - 1, j, i] = grayVal

        elif direct == 1:
            slice0 = numpy.zeros((nz, ny, 1))
            slice1 = numpy.zeros((nz, ny, 1))
            for j in range(ny):
                for k in range(nz):
                    slice0[k, j, 0] = voxelModel[k, j, 0]
                    slice1[k, j, 0] = voxelModel[k, j, nx - 1]

            plane2d = 1
            slice0 = micf77.computefill(slice0, thres, valid, inOut, kernel, minThick, plane2d, 0, kernel, fbbox)
            slice0 = micf77.computedilation2(slice0, thres, kernel, 2)
            slice0 = micf77.computeerosion(slice0, thres, kernel, 2)
            slice1 = micf77.computefill(slice1, thres, valid, inOut, kernel, minThick, plane2d, 0, kernel, fbbox)
            slice1 = micf77.computedilation2(slice1, thres, kernel, 2)
            slice1 = micf77.computeerosion(slice1, thres, kernel, 2)
            for j in range(ny):
                for k in range(nz):
                    if slice0[k, j, 0] >= thres:
                        voxelModel[k, j, 0] = grayVal
                    if slice1[k, j, 0] >= thres:
                        voxelModel[k, j, nx - 1] = grayVal

        elif direct == 2:
            slice0 = numpy.zeros((nz, 1, nx))
            slice1 = numpy.zeros((nz, 1, nx))
            for i in range(nx):
                for k in range(nz):
                    slice0[k, 0, i] = voxelModel[k, 0, i]
                    slice1[k, 0, i] = voxelModel[k, ny - 1, i]

            plane2d = 2
            slice0 = micf77.computefill(slice0, thres, valid, inOut, kernel, minThick, plane2d, 0, kernel, fbbox)
            slice0 = micf77.computedilation2(slice0, thres, kernel, 2)
            slice0 = micf77.computeerosion(slice0, thres, kernel, 2)
            slice1 = micf77.computefill(slice1, thres, valid, inOut, kernel, minThick, plane2d, 0, kernel, fbbox)
            slice1 = micf77.computedilation2(slice1, thres, kernel, 2)
            slice1 = micf77.computeerosion(slice1, thres, kernel, 2)
            for i in range(nx):
                for k in range(nz):
                    if slice0[k, 0, i] >= thres:
                        voxelModel[k, 0, i] = grayVal
                    if slice1[k, 0, i] >= thres:
                        voxelModel[k, ny - 1, i] = grayVal

        else:
            stdout.write("\n **ERROR** closeModel(): Direction '%i' not implemented!\n\n" % direct)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return voxelModel

    def coverModel(self, voxelModel, coverList, additData=None):
        """
        Function surrounds whole model with a layer of material. 
        The thickness of the end caps is selection by number of voxels.
        
        @param voxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        
        @param  coverList: list of cover parameters
              - TYPE: list[0] = direction, list[1] = thickness, list[2] = grayVal
              - int thickness ... thickness of endcaps in number of voxels
              - int grayVal   ... gray value which should be applied to these voxels
        
        @return:
           newVoxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        """
        stdout.write(' ... cover model\n')
        stdout.flush()
        thick = coverList[0]
        grayVal = coverList[1]
        nx, ny, nz = self.get_Shape(voxelModel)
        newVoxelModel = numpy.zeros((nz + 2 * thick, ny + 2 * thick, nx + 2 * thick)) + grayVal
        newVoxelModel = self.castType(newVoxelModel, 'f')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    newVoxelModel[k + thick, j + thick, i + thick] = voxelModel[k, j, i]

        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox - thick * lx, oy - thick * ly, oz - thick * lz]
        return newVoxelModel

    def placeInBlockModel(self, voxelModel, blockList, additData=None):
        """
        Function places voxel model inside a bigger block of material
        
        @param voxelModel: voxel model of the RVE
              - TYPE: numpy.array[iZ, jY, kX] = grayValue
              - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
              - int grayValue  ... gray value of voxel,
        
        @param  coverList: list of cover parameters
              - TYPE: list[0] = nx, list[1] = ny,  list[2] = nz, list[3] = grayVal
              - int nx,ny,nz  ... number of voxels in x, y, z of new block
              - int grayVal   ... gray value which should be applied to the new voxels
        
        @return:
           newVoxelModel: voxel model of the RVE
              - TYPE: numpy.array[nz, ny, nx] = grayValue
        """
        stdout.write(' ... place model in block\n')
        stdout.flush()
        NX = int(blockList[0])
        NY = int(blockList[1])
        NZ = int(blockList[2])
        grayVal = int(blockList[3])
        nx, ny, nz = self.get_Shape(voxelModel)
        if nx > NX or ny > NY or nz > NZ:
            stdout.write('\n\n **ERROR** placeInBlockModel(): Old dimensions bigger than new dimensions!\n')
            stdout.write('           OLD: nx,ny,nz=%i,%i,%i  NEW: NX,NY,NZ=%i,%i,%i\n\n' % (nx, ny, nz, NX, NY, NZ))
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        startX = int((NX - nx) / 2.0 + 0.5)
        startY = int((NY - ny) / 2.0 + 0.5)
        startZ = int((NZ - nz) / 2.0 + 0.5)
        newVoxelModel = numpy.zeros((NZ, NY, NX), dtype=numpy.float) + float(grayVal)
        newVoxelModel[startZ:startZ + nz, startY:startY + ny, startX:startX + nx] = voxelModel
        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox - startX * lx, oy - startY * ly, oz - startZ * lz]
        return newVoxelModel

    def embedModel(self, voxelModel, embedList, additData=None):
        """
        Function embeds model in given direction within layers of
        given thickness. The embedding material grayValue is given.
        Note: Model should be scaled such that no bone voxel has
        the same grayvalue as the embedding material!
        
        @param  voxelModel: voxel model
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        @param embedList: List of embedding parameters
               - TYPE: list[0] = direction, list[1] = thickIn, list[2] = thickOut
                   list[3] = grayVal
               - int direction ... direction of end caps 1,2, or 3
               - int thickIn   ... thickness of endcaps in number of voxels which
                                   go inside the bone
               - int thickOut  ... thickness of endcaps in number of voxels which
                                   go outside the bone
               - int grayVal ... gray value which should be applied to these voxels
        
        @return:
          newVoxelModel: new voxel model with embedded end caps
             - TYPE: numpy.array[kZ, jY, iX] = grayValue
             - int iX, jY, kZ ... voxels number ID in x,y,z start a 0, x fastest.
             - int grayValue  ... gray value of voxel,
        """
        stdout.write(' ... embed model \n')
        stdout.flush()
        stdout.write("     -> recast model from '%s' to 'i'\n" % voxelModel.dtype.char)
        stdout.flush()
        curVoxelModel = self.castType(voxelModel, 'i')
        nx, ny, nz = self.get_Shape(voxelModel)
        direct = embedList[0]
        thickIn = int(embedList[1])
        thickOut = int(embedList[2])
        grayVal = float(embedList[3])
        oxn = oyn = ozn = 0.0
        bottom = False
        top = False
        if direct.find('-') > -1:
            bottom = True
        elif direct.find('+') > -1:
            top = True
        else:
            bottom = True
            top = True
        if direct.find('1') > -1:
            if thickIn > nx / 2:
                stdout.write('\n\n **ERROR** embedModel(): ThickIn=%s has to smaller than nx/2=%s!\n\n' % (repr(thickIn), repr(nx / 2)))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            projArea = numpy.zeros((nz, ny))
            projArea = self.castType(projArea, 'i')
            for k in range(nz):
                for j in range(ny):
                    line = voxelModel[k, j, :]
                    if numpy.sum(line) > 0:
                        projArea[k, j] = grayVal
                        if bottom:
                            start = False
                            ip = thickIn + 1
                            while ip > 0:
                                ip -= 1
                                if start == True and line[ip] == 0:
                                    voxelModel[k, j, ip] = grayVal
                                if line[ip] > 0:
                                    start = True

                        if top:
                            start = False
                            ip = nx - 1 - (thickIn + 1)
                            while ip < nx - 1:
                                ip += 1
                                if start == True and line[ip] == 0:
                                    voxelModel[k, j, ip] = grayVal
                                if line[ip] > 0:
                                    start = True

            if thickOut > 0:
                if top == False:
                    oxn = thickOut
                    newVoxelModel = self.createVoxelModel(nx + thickOut, ny, nz, 'i')
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx + thickOut):
                                if i < thickOut:
                                    newVoxelModel[k, j, i] = projArea[k, j]
                                else:
                                    newVoxelModel[k, j, i] = voxelModel[k, j, i - thickOut]

                elif bottom == False:
                    newVoxelModel = self.createVoxelModel(nx + thickOut, ny, nz, 'i')
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx + thickOut):
                                if i < nx:
                                    newVoxelModel[k, j, i] = voxelModel[k, j, i]
                                else:
                                    newVoxelModel[k, j, i] = projArea[k, j]

                else:
                    oxn = thickOut
                    newVoxelModel = self.createVoxelModel(nx + 2 * thickOut, ny, nz, 'i')
                    for k in range(nz):
                        for j in range(ny):
                            for i in range(nx + 2 * thickOut):
                                if i < thickOut or i > nx - 1 + thickOut:
                                    newVoxelModel[k, j, i] = projArea[k, j]
                                else:
                                    newVoxelModel[k, j, i] = voxelModel[k, j, i - thickOut]

        elif direct.find('2') > -1:
            if thickIn > ny / 2:
                stdout.write('\n\n **ERROR** embedModel(): ThickIn=%s has to smaller than ny/2=%s!\n\n' % (repr(thickIn), repr(ny / 2)))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            projArea = numpy.zeros((nz, nx))
            projArea = self.castType(projArea, 'i')
            for k in range(nz):
                for i in range(nx):
                    line = voxelModel[k, :, i]
                    if numpy.sum(line) > 0:
                        projArea[k, i] = grayVal
                        if bottom:
                            start = False
                            jp = thickIn + 1
                            while jp > 0:
                                jp -= 1
                                if start == True and line[jp] == 0:
                                    voxelModel[k, jp, i] = grayVal
                                if line[jp] > 0:
                                    start = True

                        if top:
                            start = False
                            jp = ny - 1 - (thickIn + 1)
                            while jp < ny - 1:
                                jp += 1
                                if start == True and line[jp] == 0:
                                    voxelModel[k, jp, i] = grayVal
                                if line[jp] > 0:
                                    start = True

            if thickOut > 0:
                if top == False:
                    oyn = thickOut
                    newVoxelModel = self.createVoxelModel(nx, ny + thickOut, nz, 'i')
                    for k in range(nz):
                        for j in range(ny + thickOut):
                            for i in range(nx):
                                if j < thickOut:
                                    newVoxelModel[k, j, i] = projArea[k, i]
                                else:
                                    newVoxelModel[k, j, i] = voxelModel[k, j - thickOut, i]

                elif bottom == False:
                    newVoxelModel = self.createVoxelModel(nx, ny + thickOut, nz, 'i')
                    for k in range(nz):
                        for j in range(ny + thickOut):
                            for i in range(nx):
                                if j < ny:
                                    newVoxelModel[k, j, i] = voxelModel[k, j, i]
                                else:
                                    newVoxelModel[k, j, i] = projArea[k, i]

                else:
                    oyn = thickOut
                    newVoxelModel = self.createVoxelModel(nx, ny + 2 * thickOut, nz, 'i')
                    for k in range(nz):
                        for j in range(ny + 2 * thickOut):
                            for i in range(nx):
                                if j < thickOut or j > ny - 1 + thickOut:
                                    newVoxelModel[k, j, i] = projArea[k, i]
                                else:
                                    newVoxelModel[k, j, i] = voxelModel[k, j - thickOut, i]

        elif direct.find('3') > -1:
            if thickIn > nz / 2:
                stdout.write('\n\n **ERROR** embedModel(): ThickIn=%s has to smaller than nz/2=%s!\n\n' % (repr(thickIn), repr(nz / 2)))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            projArea = numpy.zeros((ny, nx))
            projArea = self.castType(projArea, 'i')
            for j in range(ny):
                for i in range(nx):
                    line = voxelModel[:, j, i]
                    if numpy.sum(line) > 0:
                        projArea[j, i] = grayVal
                        if bottom:
                            start = False
                            kp = thickIn + 1
                            while kp > 0:
                                kp -= 1
                                if start == True and line[kp] == 0:
                                    voxelModel[kp, j, i] = grayVal
                                if line[kp] > 0:
                                    start = True

                        if top:
                            start = False
                            kp = nz - 1 - (thickIn + 1)
                            while kp < nz - 1:
                                kp += 1
                                if start == True and line[kp] == 0:
                                    voxelModel[kp, j, i] = grayVal
                                if line[kp] > 0:
                                    start = True

            if thickOut > 0:
                if top == False:
                    ozn = thickOut
                    newVoxelModel = self.createVoxelModel(nx, ny, nz + thickOut, 'i')
                    for k in range(nz + thickOut):
                        if k < thickOut:
                            newVoxelModel[k] = projArea
                        else:
                            newVoxelModel[k] = voxelModel[k - thickOut]

                elif bottom == False:
                    newVoxelModel = self.createVoxelModel(nx, ny, nz + thickOut, 'i')
                    for k in range(nz + thickOut):
                        if k < nz:
                            newVoxelModel[k] = voxelModel[k]
                        else:
                            newVoxelModel[k] = projArea

                else:
                    ozn = thickOut
                    newVoxelModel = self.createVoxelModel(nx, ny, nz + 2 * thickOut, 'i')
                    for k in range(nz + 2 * thickOut):
                        if k < thickOut or k > nz - 1 + thickOut:
                            newVoxelModel[k] = projArea
                        else:
                            newVoxelModel[k] = voxelModel[k - thickOut]

        else:
            stdout.write('\n\n **ERROR** embedModel(): Direction "%s" unknown!\n\n' % repr(direct))
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if additData is not None:
            lx, ly, lz = additData['-ldim']
            ox, oy, oz = additData['Offset']
            additData['Offset'] = [ox - oxn * lx, oy - oyn * ly, oz - ozn * lz]
        if thickOut > 0:
            return self.castType(newVoxelModel, 'f')
        else:
            return voxelModel
            return


#    def smoothVoxelModel--- This code section failed: ---
#
#11242       0  LOAD_GLOBAL           0  'stdout'
#           3  LOAD_ATTR             1  'write'
#           6  LOAD_CONST            1  ' ... smooth model            \n'
#           9  CALL_FUNCTION_1       1 
#          12  POP_TOP          
#          13  LOAD_GLOBAL           0  'stdout'
#          16  LOAD_ATTR             2  'flush'
#          19  CALL_FUNCTION_0       0 
#          22  POP_TOP          
#
#11244      23  LOAD_FAST             2  'dimList'
#          26  LOAD_CONST            2  ''
#          29  BINARY_SUBSCR    
#          30  STORE_FAST            5  'xvox'
#          33  LOAD_FAST             2  'dimList'
#          36  LOAD_CONST            3  1
#          39  BINARY_SUBSCR    
#          40  STORE_FAST            6  'yvox'
#          43  LOAD_FAST             2  'dimList'
#          46  LOAD_CONST            4  2
#          49  BINARY_SUBSCR    
#          50  STORE_FAST            7  'zvox'
#
#11245      53  LOAD_FAST             0  'self'
#          56  LOAD_ATTR             3  'get_Shape'
#          59  LOAD_FAST             1  'curVoxelModel'
#          62  CALL_FUNCTION_1       1 
#          65  UNPACK_SEQUENCE_3     3 
#          68  STORE_FAST            8  'nx'
#          71  STORE_FAST            9  'ny'
#          74  STORE_FAST           10  'nz'
#
#11246      77  LOAD_FAST             8  'nx'
#          80  LOAD_CONST            3  1
#          83  BINARY_ADD       
#          84  STORE_FAST           11  'nx1'
#
#11247      87  LOAD_FAST             9  'ny'
#          90  LOAD_CONST            3  1
#          93  BINARY_ADD       
#          94  LOAD_FAST             8  'nx'
#          97  LOAD_CONST            3  1
#         100  BINARY_ADD       
#         101  BINARY_MULTIPLY  
#         102  STORE_FAST           12  'nxy1'
#
#11249     105  LOAD_GLOBAL           4  'True'
#         108  STORE_FAST           13  'nearIntf'
#
#11250     111  LOAD_GLOBAL           5  'int'
#         114  LOAD_FAST             3  'smoothParam'
#         117  LOAD_CONST            5  3
#         120  BINARY_SUBSCR    
#         121  CALL_FUNCTION_1       1 
#         124  LOAD_CONST            2  ''
#         127  COMPARE_OP            2  '=='
#         130  POP_JUMP_IF_FALSE   142  'to 142'
#         133  LOAD_GLOBAL           6  'False'
#         136  STORE_FAST           13  'nearIntf'
#         139  JUMP_FORWARD          0  'to 142'
#       142_0  COME_FROM                '139'
#
#11251     142  LOAD_GLOBAL           5  'int'
#         145  LOAD_FAST             3  'smoothParam'
#         148  LOAD_CONST            6  4
#         151  BINARY_SUBSCR    
#         152  CALL_FUNCTION_1       1 
#         155  STORE_FAST           14  'bcidParam'
#
#11253     158  BUILD_MAP_0           0 
#         161  STORE_FAST           15  'nodeCoord'
#
#11254     164  BUILD_MAP_0           0 
#         167  STORE_FAST           16  'nodeCoordPre'
#
#11256     170  LOAD_CONST            2  ''
#         173  STORE_FAST           17  'noid'
#
#11257     176  SETUP_LOOP          204  'to 383'
#         179  LOAD_GLOBAL           7  'range'
#         182  LOAD_FAST            10  'nz'
#         185  LOAD_CONST            3  1
#         188  BINARY_ADD       
#         189  CALL_FUNCTION_1       1 
#         192  GET_ITER         
#         193  FOR_ITER            186  'to 382'
#         196  STORE_FAST           18  'k'
#
#11258     199  SETUP_LOOP          177  'to 379'
#         202  LOAD_GLOBAL           7  'range'
#         205  LOAD_FAST             9  'ny'
#         208  LOAD_CONST            3  1
#         211  BINARY_ADD       
#         212  CALL_FUNCTION_1       1 
#         215  GET_ITER         
#         216  FOR_ITER            159  'to 378'
#         219  STORE_FAST           19  'j'
#
#11259     222  SETUP_LOOP          150  'to 375'
#         225  LOAD_GLOBAL           7  'range'
#         228  LOAD_FAST             8  'nx'
#         231  LOAD_CONST            3  1
#         234  BINARY_ADD       
#         235  CALL_FUNCTION_1       1 
#         238  GET_ITER         
#         239  FOR_ITER            132  'to 374'
#         242  STORE_FAST           20  'i'
#
#11260     245  LOAD_FAST            17  'noid'
#         248  LOAD_CONST            3  1
#         251  BINARY_ADD       
#         252  STORE_FAST           17  'noid'
#
#11261     255  LOAD_FAST             4  'activeNodes'
#         258  LOAD_ATTR             8  'has_key'
#         261  LOAD_FAST            17  'noid'
#         264  CALL_FUNCTION_1       1 
#         267  POP_JUMP_IF_FALSE   239  'to 239'
#
#11262     270  LOAD_GLOBAL           9  'float'
#         273  LOAD_FAST             5  'xvox'
#         276  LOAD_FAST            20  'i'
#         279  BINARY_MULTIPLY  
#         280  CALL_FUNCTION_1       1 
#         283  LOAD_GLOBAL           9  'float'
#         286  LOAD_FAST             6  'yvox'
#         289  LOAD_FAST            19  'j'
#         292  BINARY_MULTIPLY  
#         293  CALL_FUNCTION_1       1 
#         296  LOAD_GLOBAL           9  'float'
#         299  LOAD_FAST             7  'zvox'
#         302  LOAD_FAST            18  'k'
#         305  BINARY_MULTIPLY  
#         306  CALL_FUNCTION_1       1 
#         309  BUILD_TUPLE_3         3 
#         312  LOAD_FAST            15  'nodeCoord'
#         315  LOAD_FAST            17  'noid'
#         318  STORE_SUBSCR     
#
#11263     319  LOAD_GLOBAL           9  'float'
#         322  LOAD_FAST             5  'xvox'
#         325  LOAD_FAST            20  'i'
#         328  BINARY_MULTIPLY  
#         329  CALL_FUNCTION_1       1 
#         332  LOAD_GLOBAL           9  'float'
#         335  LOAD_FAST             6  'yvox'
#         338  LOAD_FAST            19  'j'
#         341  BINARY_MULTIPLY  
#         342  CALL_FUNCTION_1       1 
#         345  LOAD_GLOBAL           9  'float'
#         348  LOAD_FAST             7  'zvox'
#         351  LOAD_FAST            18  'k'
#         354  BINARY_MULTIPLY  
#         355  CALL_FUNCTION_1       1 
#         358  BUILD_TUPLE_3         3 
#         361  LOAD_FAST            16  'nodeCoordPre'
#         364  LOAD_FAST            17  'noid'
#         367  STORE_SUBSCR     
#         368  JUMP_BACK           239  'to 239'
#         371  JUMP_BACK           239  'to 239'
#         374  POP_BLOCK        
#       375_0  COME_FROM                '222'
#         375  JUMP_BACK           216  'to 216'
#         378  POP_BLOCK        
#       379_0  COME_FROM                '199'
#         379  JUMP_BACK           193  'to 193'
#         382  POP_BLOCK        
#       383_0  COME_FROM                '176'
#
#11266     383  LOAD_CONST            2  ''
#         386  STORE_FAST           21  'sum'
#
#11267     389  BUILD_MAP_0           0 
#         392  STORE_FAST           22  'intfaceNodes'
#
#11268     395  SETUP_LOOP          749  'to 1147'
#         398  LOAD_GLOBAL           7  'range'
#         401  LOAD_FAST            10  'nz'
#         404  CALL_FUNCTION_1       1 
#         407  GET_ITER         
#         408  FOR_ITER            735  'to 1146'
#         411  STORE_FAST           18  'k'
#
#11269     414  LOAD_FAST            21  'sum'
#         417  LOAD_CONST            3  1
#         420  INPLACE_ADD      
#         421  STORE_FAST           21  'sum'
#
#11270     424  LOAD_GLOBAL           9  'float'
#         427  LOAD_FAST            21  'sum'
#         430  CALL_FUNCTION_1       1 
#         433  LOAD_GLOBAL           9  'float'
#         436  LOAD_FAST            10  'nz'
#         439  CALL_FUNCTION_1       1 
#         442  BINARY_DIVIDE    
#         443  LOAD_CONST            7  100.0
#         446  BINARY_MULTIPLY  
#         447  STORE_FAST           23  'progress'
#
#11271     450  LOAD_GLOBAL           0  'stdout'
#         453  LOAD_ATTR             1  'write'
#         456  LOAD_CONST            8  '     -> Setup Interface Node: %10i %s\r'
#         459  LOAD_FAST            23  'progress'
#         462  LOAD_CONST            9  '%'
#         465  BUILD_TUPLE_2         2 
#         468  BINARY_MODULO    
#         469  CALL_FUNCTION_1       1 
#         472  POP_TOP          
#
#11272     473  LOAD_GLOBAL           0  'stdout'
#         476  LOAD_ATTR             2  'flush'
#         479  CALL_FUNCTION_0       0 
#         482  POP_TOP          
#
#11273     483  SETUP_LOOP          657  'to 1143'
#         486  LOAD_GLOBAL           7  'range'
#         489  LOAD_FAST             9  'ny'
#         492  CALL_FUNCTION_1       1 
#         495  GET_ITER         
#         496  FOR_ITER            643  'to 1142'
#         499  STORE_FAST           19  'j'
#
#11274     502  SETUP_LOOP          634  'to 1139'
#         505  LOAD_GLOBAL           7  'range'
#         508  LOAD_FAST             8  'nx'
#         511  CALL_FUNCTION_1       1 
#         514  GET_ITER         
#         515  FOR_ITER            620  'to 1138'
#         518  STORE_FAST           20  'i'
#
#11275     521  LOAD_FAST             1  'curVoxelModel'
#         524  LOAD_FAST            18  'k'
#         527  LOAD_FAST            19  'j'
#         530  LOAD_FAST            20  'i'
#         533  BUILD_TUPLE_3         3 
#         536  BINARY_SUBSCR    
#         537  STORE_FAST           24  'grayValue'
#
#11277     540  LOAD_FAST            20  'i'
#         543  LOAD_FAST             8  'nx'
#         546  LOAD_CONST            3  1
#         549  BINARY_SUBTRACT  
#         550  COMPARE_OP            2  '=='
#         553  POP_JUMP_IF_FALSE   559  'to 559'
#
#11278     556  JUMP_FORWARD        174  'to 733'
#
#11280     559  LOAD_FAST             1  'curVoxelModel'
#         562  LOAD_FAST            18  'k'
#         565  LOAD_FAST            19  'j'
#         568  LOAD_FAST            20  'i'
#         571  LOAD_CONST            3  1
#         574  BINARY_ADD       
#         575  BUILD_TUPLE_3         3 
#         578  BINARY_SUBSCR    
#         579  STORE_FAST           25  'grayValue2'
#
#11281     582  LOAD_FAST            25  'grayValue2'
#         585  LOAD_FAST            24  'grayValue'
#         588  COMPARE_OP            3  '!='
#         591  POP_JUMP_IF_FALSE   733  'to 733'
#
#11282     594  LOAD_CONST            2  ''
#         597  LOAD_FAST            22  'intfaceNodes'
#         600  LOAD_FAST            12  'nxy1'
#         603  LOAD_FAST            18  'k'
#         606  BINARY_MULTIPLY  
#         607  LOAD_FAST            11  'nx1'
#         610  LOAD_FAST            19  'j'
#         613  BINARY_MULTIPLY  
#         614  BINARY_ADD       
#         615  LOAD_FAST            20  'i'
#         618  LOAD_CONST            4  2
#         621  BINARY_ADD       
#         622  BINARY_ADD       
#         623  STORE_SUBSCR     
#
#11283     624  LOAD_CONST            2  ''
#         627  LOAD_FAST            22  'intfaceNodes'
#         630  LOAD_FAST            12  'nxy1'
#         633  LOAD_FAST            18  'k'
#         636  BINARY_MULTIPLY  
#         637  LOAD_FAST            11  'nx1'
#         640  LOAD_FAST            19  'j'
#         643  LOAD_CONST            3  1
#         646  BINARY_ADD       
#         647  BINARY_MULTIPLY  
#         648  BINARY_ADD       
#         649  LOAD_FAST            20  'i'
#         652  LOAD_CONST            4  2
#         655  BINARY_ADD       
#         656  BINARY_ADD       
#         657  STORE_SUBSCR     
#
#11284     658  LOAD_CONST            2  ''
#         661  LOAD_FAST            22  'intfaceNodes'
#         664  LOAD_FAST            12  'nxy1'
#         667  LOAD_FAST            18  'k'
#         670  LOAD_CONST            3  1
#         673  BINARY_ADD       
#         674  BINARY_MULTIPLY  
#         675  LOAD_FAST            11  'nx1'
#         678  LOAD_FAST            19  'j'
#         681  BINARY_MULTIPLY  
#         682  BINARY_ADD       
#         683  LOAD_FAST            20  'i'
#         686  LOAD_CONST            4  2
#         689  BINARY_ADD       
#         690  BINARY_ADD       
#         691  STORE_SUBSCR     
#
#11285     692  LOAD_CONST            2  ''
#         695  LOAD_FAST            22  'intfaceNodes'
#         698  LOAD_FAST            12  'nxy1'
#         701  LOAD_FAST            18  'k'
#         704  LOAD_CONST            3  1
#         707  BINARY_ADD       
#         708  BINARY_MULTIPLY  
#         709  LOAD_FAST            11  'nx1'
#         712  LOAD_FAST            19  'j'
#         715  LOAD_CONST            3  1
#         718  BINARY_ADD       
#         719  BINARY_MULTIPLY  
#         720  BINARY_ADD       
#         721  LOAD_FAST            20  'i'
#         724  LOAD_CONST            4  2
#         727  BINARY_ADD       
#         728  BINARY_ADD       
#         729  STORE_SUBSCR     
#         730  JUMP_FORWARD          0  'to 733'
#       733_0  COME_FROM                '730'
#       733_1  COME_FROM                '556'
#
#11287     733  LOAD_FAST            19  'j'
#         736  LOAD_FAST             9  'ny'
#         739  LOAD_CONST            3  1
#         742  BINARY_SUBTRACT  
#         743  COMPARE_OP            2  '=='
#         746  POP_JUMP_IF_FALSE   752  'to 752'
#
#11288     749  JUMP_FORWARD        182  'to 934'
#
#11290     752  LOAD_FAST             1  'curVoxelModel'
#         755  LOAD_FAST            18  'k'
#         758  LOAD_FAST            19  'j'
#         761  LOAD_CONST            3  1
#         764  BINARY_ADD       
#         765  LOAD_FAST            20  'i'
#         768  BUILD_TUPLE_3         3 
#         771  BINARY_SUBSCR    
#         772  STORE_FAST           25  'grayValue2'
#
#11291     775  LOAD_FAST            25  'grayValue2'
#         778  LOAD_FAST            24  'grayValue'
#         781  COMPARE_OP            3  '!='
#         784  POP_JUMP_IF_FALSE   934  'to 934'
#
#11292     787  LOAD_CONST            2  ''
#         790  LOAD_FAST            22  'intfaceNodes'
#         793  LOAD_FAST            12  'nxy1'
#         796  LOAD_FAST            18  'k'
#         799  BINARY_MULTIPLY  
#         800  LOAD_FAST            11  'nx1'
#         803  LOAD_FAST            19  'j'
#         806  LOAD_CONST            3  1
#         809  BINARY_ADD       
#         810  BINARY_MULTIPLY  
#         811  BINARY_ADD       
#         812  LOAD_FAST            20  'i'
#         815  LOAD_CONST            4  2
#         818  BINARY_ADD       
#         819  BINARY_ADD       
#         820  STORE_SUBSCR     
#
#11293     821  LOAD_CONST            2  ''
#         824  LOAD_FAST            22  'intfaceNodes'
#         827  LOAD_FAST            12  'nxy1'
#         830  LOAD_FAST            18  'k'
#         833  BINARY_MULTIPLY  
#         834  LOAD_FAST            11  'nx1'
#         837  LOAD_FAST            19  'j'
#         840  LOAD_CONST            3  1
#         843  BINARY_ADD       
#         844  BINARY_MULTIPLY  
#         845  BINARY_ADD       
#         846  LOAD_FAST            20  'i'
#         849  LOAD_CONST            3  1
#         852  BINARY_ADD       
#         853  BINARY_ADD       
#         854  STORE_SUBSCR     
#
#11294     855  LOAD_CONST            2  ''
#         858  LOAD_FAST            22  'intfaceNodes'
#         861  LOAD_FAST            12  'nxy1'
#         864  LOAD_FAST            18  'k'
#         867  LOAD_CONST            3  1
#         870  BINARY_ADD       
#         871  BINARY_MULTIPLY  
#         872  LOAD_FAST            11  'nx1'
#         875  LOAD_FAST            19  'j'
#         878  LOAD_CONST            3  1
#         881  BINARY_ADD       
#         882  BINARY_MULTIPLY  
#         883  BINARY_ADD       
#         884  LOAD_FAST            20  'i'
#         887  LOAD_CONST            4  2
#         890  BINARY_ADD       
#         891  BINARY_ADD       
#         892  STORE_SUBSCR     
#
#11295     893  LOAD_CONST            2  ''
#         896  LOAD_FAST            22  'intfaceNodes'
#         899  LOAD_FAST            12  'nxy1'
#         902  LOAD_FAST            18  'k'
#         905  LOAD_CONST            3  1
#         908  BINARY_ADD       
#         909  BINARY_MULTIPLY  
#         910  LOAD_FAST            11  'nx1'
#         913  LOAD_FAST            19  'j'
#         916  LOAD_CONST            3  1
#         919  BINARY_ADD       
#         920  BINARY_MULTIPLY  
#         921  BINARY_ADD       
#         922  LOAD_FAST            20  'i'
#         925  LOAD_CONST            3  1
#         928  BINARY_ADD       
#         929  BINARY_ADD       
#         930  STORE_SUBSCR     
#         931  JUMP_FORWARD          0  'to 934'
#       934_0  COME_FROM                '931'
#       934_1  COME_FROM                '749'
#
#11297     934  LOAD_FAST            18  'k'
#         937  LOAD_FAST            10  'nz'
#         940  LOAD_CONST            3  1
#         943  BINARY_SUBTRACT  
#         944  COMPARE_OP            2  '=='
#         947  POP_JUMP_IF_FALSE   953  'to 953'
#
#11298     950  CONTINUE            515  'to 515'
#
#11300     953  LOAD_FAST             1  'curVoxelModel'
#         956  LOAD_FAST            18  'k'
#         959  LOAD_CONST            3  1
#         962  BINARY_ADD       
#         963  LOAD_FAST            19  'j'
#         966  LOAD_FAST            20  'i'
#         969  BUILD_TUPLE_3         3 
#         972  BINARY_SUBSCR    
#         973  STORE_FAST           25  'grayValue2'
#
#11301     976  LOAD_FAST            25  'grayValue2'
#         979  LOAD_FAST            24  'grayValue'
#         982  COMPARE_OP            3  '!='
#         985  POP_JUMP_IF_FALSE   515  'to 515'
#
#11302     988  LOAD_CONST            2  ''
#         991  LOAD_FAST            22  'intfaceNodes'
#         994  LOAD_FAST            12  'nxy1'
#         997  LOAD_FAST            18  'k'
#        1000  LOAD_CONST            3  1
#        1003  BINARY_ADD       
#        1004  BINARY_MULTIPLY  
#        1005  LOAD_FAST            11  'nx1'
#        1008  LOAD_FAST            19  'j'
#        1011  BINARY_MULTIPLY  
#        1012  BINARY_ADD       
#        1013  LOAD_FAST            20  'i'
#        1016  LOAD_CONST            3  1
#        1019  BINARY_ADD       
#        1020  BINARY_ADD       
#        1021  STORE_SUBSCR     
#
#11303    1022  LOAD_CONST            2  ''
#        1025  LOAD_FAST            22  'intfaceNodes'
#        1028  LOAD_FAST            12  'nxy1'
#        1031  LOAD_FAST            18  'k'
#        1034  LOAD_CONST            3  1
#        1037  BINARY_ADD       
#        1038  BINARY_MULTIPLY  
#        1039  LOAD_FAST            11  'nx1'
#        1042  LOAD_FAST            19  'j'
#        1045  BINARY_MULTIPLY  
#        1046  BINARY_ADD       
#        1047  LOAD_FAST            20  'i'
#        1050  LOAD_CONST            4  2
#        1053  BINARY_ADD       
#        1054  BINARY_ADD       
#        1055  STORE_SUBSCR     
#
#11304    1056  LOAD_CONST            2  ''
#        1059  LOAD_FAST            22  'intfaceNodes'
#        1062  LOAD_FAST            12  'nxy1'
#        1065  LOAD_FAST            18  'k'
#        1068  LOAD_CONST            3  1
#        1071  BINARY_ADD       
#        1072  BINARY_MULTIPLY  
#        1073  LOAD_FAST            11  'nx1'
#        1076  LOAD_FAST            19  'j'
#        1079  LOAD_CONST            3  1
#        1082  BINARY_ADD       
#        1083  BINARY_MULTIPLY  
#        1084  BINARY_ADD       
#        1085  LOAD_FAST            20  'i'
#        1088  LOAD_CONST            4  2
#        1091  BINARY_ADD       
#        1092  BINARY_ADD       
#        1093  STORE_SUBSCR     
#
#11305    1094  LOAD_CONST            2  ''
#        1097  LOAD_FAST            22  'intfaceNodes'
#        1100  LOAD_FAST            12  'nxy1'
#        1103  LOAD_FAST            18  'k'
#        1106  LOAD_CONST            3  1
#        1109  BINARY_ADD       
#        1110  BINARY_MULTIPLY  
#        1111  LOAD_FAST            11  'nx1'
#        1114  LOAD_FAST            19  'j'
#        1117  LOAD_CONST            3  1
#        1120  BINARY_ADD       
#        1121  BINARY_MULTIPLY  
#        1122  BINARY_ADD       
#        1123  LOAD_FAST            20  'i'
#        1126  LOAD_CONST            3  1
#        1129  BINARY_ADD       
#        1130  BINARY_ADD       
#        1131  STORE_SUBSCR     
#        1132  JUMP_BACK           515  'to 515'
#        1135  JUMP_BACK           515  'to 515'
#        1138  POP_BLOCK        
#      1139_0  COME_FROM                '502'
#        1139  JUMP_BACK           496  'to 496'
#        1142  POP_BLOCK        
#      1143_0  COME_FROM                '483'
#        1143  JUMP_BACK           408  'to 408'
#        1146  POP_BLOCK        
#      1147_0  COME_FROM                '395'
#
#11308    1147  LOAD_FAST            13  'nearIntf'
#        1150  POP_JUMP_IF_FALSE  2940  'to 2940'
#
#11310    1153  LOAD_CONST            2  ''
#        1156  STORE_FAST           21  'sum'
#
#11311    1159  BUILD_MAP_0           0 
#        1162  STORE_FAST           26  'nearIntfaceNodes'
#
#11312    1165  SETUP_LOOP         1772  'to 2940'
#        1168  LOAD_GLOBAL           7  'range'
#        1171  LOAD_FAST            10  'nz'
#        1174  LOAD_CONST            3  1
#        1177  BINARY_SUBTRACT  
#        1178  CALL_FUNCTION_1       1 
#        1181  GET_ITER         
#        1182  FOR_ITER           1751  'to 2936'
#        1185  STORE_FAST           18  'k'
#
#11313    1188  LOAD_FAST            21  'sum'
#        1191  LOAD_CONST            3  1
#        1194  INPLACE_ADD      
#        1195  STORE_FAST           21  'sum'
#
#11314    1198  LOAD_GLOBAL           9  'float'
#        1201  LOAD_FAST            21  'sum'
#        1204  CALL_FUNCTION_1       1 
#        1207  LOAD_GLOBAL           9  'float'
#        1210  LOAD_FAST            10  'nz'
#        1213  LOAD_CONST            3  1
#        1216  BINARY_SUBTRACT  
#        1217  CALL_FUNCTION_1       1 
#        1220  BINARY_DIVIDE    
#        1221  LOAD_CONST            7  100.0
#        1224  BINARY_MULTIPLY  
#        1225  STORE_FAST           23  'progress'
#
#11315    1228  LOAD_GLOBAL           0  'stdout'
#        1231  LOAD_ATTR             1  'write'
#        1234  LOAD_CONST           10  '     -> Setup Near Intf Node: %10i %s\r'
#        1237  LOAD_FAST            23  'progress'
#        1240  LOAD_CONST            9  '%'
#        1243  BUILD_TUPLE_2         2 
#        1246  BINARY_MODULO    
#        1247  CALL_FUNCTION_1       1 
#        1250  POP_TOP          
#
#11316    1251  LOAD_GLOBAL           0  'stdout'
#        1254  LOAD_ATTR             2  'flush'
#        1257  CALL_FUNCTION_0       0 
#        1260  POP_TOP          
#
#11317    1261  SETUP_LOOP         1669  'to 2933'
#        1264  LOAD_GLOBAL           7  'range'
#        1267  LOAD_FAST             9  'ny'
#        1270  LOAD_CONST            3  1
#        1273  BINARY_SUBTRACT  
#        1274  CALL_FUNCTION_1       1 
#        1277  GET_ITER         
#        1278  FOR_ITER           1651  'to 2932'
#        1281  STORE_FAST           19  'j'
#
#11318    1284  SETUP_LOOP         1642  'to 2929'
#        1287  LOAD_GLOBAL           7  'range'
#        1290  LOAD_FAST             8  'nx'
#        1293  LOAD_CONST            3  1
#        1296  BINARY_SUBTRACT  
#        1297  CALL_FUNCTION_1       1 
#        1300  GET_ITER         
#        1301  FOR_ITER           1624  'to 2928'
#        1304  STORE_FAST           20  'i'
#
#11319    1307  LOAD_FAST             1  'curVoxelModel'
#        1310  LOAD_FAST            18  'k'
#        1313  LOAD_FAST            19  'j'
#        1316  LOAD_FAST            20  'i'
#        1319  BUILD_TUPLE_3         3 
#        1322  BINARY_SUBSCR    
#        1323  STORE_FAST           24  'grayValue'
#
#11320    1326  LOAD_FAST            24  'grayValue'
#        1329  LOAD_CONST            2  ''
#        1332  COMPARE_OP            4  '>'
#        1335  POP_JUMP_IF_FALSE  1301  'to 1301'
#
#11321    1338  LOAD_FAST             1  'curVoxelModel'
#        1341  LOAD_FAST            18  'k'
#        1344  LOAD_FAST            19  'j'
#        1347  LOAD_FAST            20  'i'
#        1350  LOAD_CONST            3  1
#        1353  BINARY_ADD       
#        1354  BUILD_TUPLE_3         3 
#        1357  BINARY_SUBSCR    
#        1358  LOAD_FAST            24  'grayValue'
#        1361  COMPARE_OP            3  '!='
#        1364  POP_JUMP_IF_FALSE  1602  'to 1602'
#
#11322    1367  LOAD_FAST            12  'nxy1'
#        1370  LOAD_FAST            18  'k'
#        1373  BINARY_MULTIPLY  
#        1374  LOAD_FAST            11  'nx1'
#        1377  LOAD_FAST            19  'j'
#        1380  BINARY_MULTIPLY  
#        1381  BINARY_ADD       
#        1382  LOAD_FAST            20  'i'
#        1385  LOAD_CONST            3  1
#        1388  BINARY_ADD       
#        1389  BINARY_ADD       
#        1390  STORE_FAST           27  'curNode'
#
#11323    1393  LOAD_FAST            22  'intfaceNodes'
#        1396  LOAD_ATTR             8  'has_key'
#        1399  LOAD_FAST            27  'curNode'
#        1402  CALL_FUNCTION_1       1 
#        1405  POP_JUMP_IF_TRUE   1421  'to 1421'
#
#11324    1408  LOAD_CONST            2  ''
#        1411  LOAD_FAST            26  'nearIntfaceNodes'
#        1414  LOAD_FAST            27  'curNode'
#        1417  STORE_SUBSCR     
#        1418  JUMP_FORWARD          0  'to 1421'
#      1421_0  COME_FROM                '1418'
#
#11325    1421  LOAD_FAST            12  'nxy1'
#        1424  LOAD_FAST            18  'k'
#        1427  BINARY_MULTIPLY  
#        1428  LOAD_FAST            11  'nx1'
#        1431  LOAD_FAST            19  'j'
#        1434  LOAD_CONST            3  1
#        1437  BINARY_ADD       
#        1438  BINARY_MULTIPLY  
#        1439  BINARY_ADD       
#        1440  LOAD_FAST            20  'i'
#        1443  LOAD_CONST            3  1
#        1446  BINARY_ADD       
#        1447  BINARY_ADD       
#        1448  STORE_FAST           27  'curNode'
#
#11326    1451  LOAD_FAST            22  'intfaceNodes'
#        1454  LOAD_ATTR             8  'has_key'
#        1457  LOAD_FAST            27  'curNode'
#        1460  CALL_FUNCTION_1       1 
#        1463  POP_JUMP_IF_TRUE   1479  'to 1479'
#
#11327    1466  LOAD_CONST            2  ''
#        1469  LOAD_FAST            26  'nearIntfaceNodes'
#        1472  LOAD_FAST            27  'curNode'
#        1475  STORE_SUBSCR     
#        1476  JUMP_FORWARD          0  'to 1479'
#      1479_0  COME_FROM                '1476'
#
#11328    1479  LOAD_FAST            12  'nxy1'
#        1482  LOAD_FAST            18  'k'
#        1485  LOAD_CONST            3  1
#        1488  BINARY_ADD       
#        1489  BINARY_MULTIPLY  
#        1490  LOAD_FAST            11  'nx1'
#        1493  LOAD_FAST            19  'j'
#        1496  BINARY_MULTIPLY  
#        1497  BINARY_ADD       
#        1498  LOAD_FAST            20  'i'
#        1501  LOAD_CONST            3  1
#        1504  BINARY_ADD       
#        1505  BINARY_ADD       
#        1506  STORE_FAST           27  'curNode'
#
#11329    1509  LOAD_FAST            22  'intfaceNodes'
#        1512  LOAD_ATTR             8  'has_key'
#        1515  LOAD_FAST            27  'curNode'
#        1518  CALL_FUNCTION_1       1 
#        1521  POP_JUMP_IF_TRUE   1537  'to 1537'
#
#11330    1524  LOAD_CONST            2  ''
#        1527  LOAD_FAST            26  'nearIntfaceNodes'
#        1530  LOAD_FAST            27  'curNode'
#        1533  STORE_SUBSCR     
#        1534  JUMP_FORWARD          0  'to 1537'
#      1537_0  COME_FROM                '1534'
#
#11331    1537  LOAD_FAST            12  'nxy1'
#        1540  LOAD_FAST            18  'k'
#        1543  LOAD_CONST            3  1
#        1546  BINARY_ADD       
#        1547  BINARY_MULTIPLY  
#        1548  LOAD_FAST            11  'nx1'
#        1551  LOAD_FAST            19  'j'
#        1554  LOAD_CONST            3  1
#        1557  BINARY_ADD       
#        1558  BINARY_MULTIPLY  
#        1559  BINARY_ADD       
#        1560  LOAD_FAST            20  'i'
#        1563  LOAD_CONST            3  1
#        1566  BINARY_ADD       
#        1567  BINARY_ADD       
#        1568  STORE_FAST           27  'curNode'
#
#11332    1571  LOAD_FAST            22  'intfaceNodes'
#        1574  LOAD_ATTR             8  'has_key'
#        1577  LOAD_FAST            27  'curNode'
#        1580  CALL_FUNCTION_1       1 
#        1583  POP_JUMP_IF_TRUE   1602  'to 1602'
#
#11333    1586  LOAD_CONST            2  ''
#        1589  LOAD_FAST            26  'nearIntfaceNodes'
#        1592  LOAD_FAST            27  'curNode'
#        1595  STORE_SUBSCR     
#        1596  JUMP_ABSOLUTE      1602  'to 1602'
#        1599  JUMP_FORWARD          0  'to 1602'
#      1602_0  COME_FROM                '1599'
#
#11335    1602  LOAD_FAST             1  'curVoxelModel'
#        1605  LOAD_FAST            18  'k'
#        1608  LOAD_FAST            19  'j'
#        1611  LOAD_FAST            20  'i'
#        1614  LOAD_CONST            3  1
#        1617  BINARY_SUBTRACT  
#        1618  BUILD_TUPLE_3         3 
#        1621  BINARY_SUBSCR    
#        1622  LOAD_FAST            24  'grayValue'
#        1625  COMPARE_OP            3  '!='
#        1628  POP_JUMP_IF_FALSE  1866  'to 1866'
#
#11336    1631  LOAD_FAST            12  'nxy1'
#        1634  LOAD_FAST            18  'k'
#        1637  BINARY_MULTIPLY  
#        1638  LOAD_FAST            11  'nx1'
#        1641  LOAD_FAST            19  'j'
#        1644  BINARY_MULTIPLY  
#        1645  BINARY_ADD       
#        1646  LOAD_FAST            20  'i'
#        1649  LOAD_CONST            4  2
#        1652  BINARY_ADD       
#        1653  BINARY_ADD       
#        1654  STORE_FAST           27  'curNode'
#
#11337    1657  LOAD_FAST            22  'intfaceNodes'
#        1660  LOAD_ATTR             8  'has_key'
#        1663  LOAD_FAST            27  'curNode'
#        1666  CALL_FUNCTION_1       1 
#        1669  POP_JUMP_IF_TRUE   1685  'to 1685'
#
#11338    1672  LOAD_CONST            2  ''
#        1675  LOAD_FAST            26  'nearIntfaceNodes'
#        1678  LOAD_FAST            27  'curNode'
#        1681  STORE_SUBSCR     
#        1682  JUMP_FORWARD          0  'to 1685'
#      1685_0  COME_FROM                '1682'
#
#11339    1685  LOAD_FAST            12  'nxy1'
#        1688  LOAD_FAST            18  'k'
#        1691  BINARY_MULTIPLY  
#        1692  LOAD_FAST            11  'nx1'
#        1695  LOAD_FAST            19  'j'
#        1698  LOAD_CONST            3  1
#        1701  BINARY_ADD       
#        1702  BINARY_MULTIPLY  
#        1703  BINARY_ADD       
#        1704  LOAD_FAST            20  'i'
#        1707  LOAD_CONST            4  2
#        1710  BINARY_ADD       
#        1711  BINARY_ADD       
#        1712  STORE_FAST           27  'curNode'
#
#11340    1715  LOAD_FAST            22  'intfaceNodes'
#        1718  LOAD_ATTR             8  'has_key'
#        1721  LOAD_FAST            27  'curNode'
#        1724  CALL_FUNCTION_1       1 
#        1727  POP_JUMP_IF_TRUE   1743  'to 1743'
#
#11341    1730  LOAD_CONST            2  ''
#        1733  LOAD_FAST            26  'nearIntfaceNodes'
#        1736  LOAD_FAST            27  'curNode'
#        1739  STORE_SUBSCR     
#        1740  JUMP_FORWARD          0  'to 1743'
#      1743_0  COME_FROM                '1740'
#
#11342    1743  LOAD_FAST            12  'nxy1'
#        1746  LOAD_FAST            18  'k'
#        1749  LOAD_CONST            3  1
#        1752  BINARY_ADD       
#        1753  BINARY_MULTIPLY  
#        1754  LOAD_FAST            11  'nx1'
#        1757  LOAD_FAST            19  'j'
#        1760  BINARY_MULTIPLY  
#        1761  BINARY_ADD       
#        1762  LOAD_FAST            20  'i'
#        1765  LOAD_CONST            4  2
#        1768  BINARY_ADD       
#        1769  BINARY_ADD       
#        1770  STORE_FAST           27  'curNode'
#
#11343    1773  LOAD_FAST            22  'intfaceNodes'
#        1776  LOAD_ATTR             8  'has_key'
#        1779  LOAD_FAST            27  'curNode'
#        1782  CALL_FUNCTION_1       1 
#        1785  POP_JUMP_IF_TRUE   1801  'to 1801'
#
#11344    1788  LOAD_CONST            2  ''
#        1791  LOAD_FAST            26  'nearIntfaceNodes'
#        1794  LOAD_FAST            27  'curNode'
#        1797  STORE_SUBSCR     
#        1798  JUMP_FORWARD          0  'to 1801'
#      1801_0  COME_FROM                '1798'
#
#11345    1801  LOAD_FAST            12  'nxy1'
#        1804  LOAD_FAST            18  'k'
#        1807  LOAD_CONST            3  1
#        1810  BINARY_ADD       
#        1811  BINARY_MULTIPLY  
#        1812  LOAD_FAST            11  'nx1'
#        1815  LOAD_FAST            19  'j'
#        1818  LOAD_CONST            3  1
#        1821  BINARY_ADD       
#        1822  BINARY_MULTIPLY  
#        1823  BINARY_ADD       
#        1824  LOAD_FAST            20  'i'
#        1827  LOAD_CONST            4  2
#        1830  BINARY_ADD       
#        1831  BINARY_ADD       
#        1832  STORE_FAST           27  'curNode'
#
#11346    1835  LOAD_FAST            22  'intfaceNodes'
#        1838  LOAD_ATTR             8  'has_key'
#        1841  LOAD_FAST            27  'curNode'
#        1844  CALL_FUNCTION_1       1 
#        1847  POP_JUMP_IF_TRUE   1866  'to 1866'
#
#11347    1850  LOAD_CONST            2  ''
#        1853  LOAD_FAST            26  'nearIntfaceNodes'
#        1856  LOAD_FAST            27  'curNode'
#        1859  STORE_SUBSCR     
#        1860  JUMP_ABSOLUTE      1866  'to 1866'
#        1863  JUMP_FORWARD          0  'to 1866'
#      1866_0  COME_FROM                '1863'
#
#11349    1866  LOAD_FAST             1  'curVoxelModel'
#        1869  LOAD_FAST            18  'k'
#        1872  LOAD_FAST            19  'j'
#        1875  LOAD_CONST            3  1
#        1878  BINARY_ADD       
#        1879  LOAD_FAST            20  'i'
#        1882  BUILD_TUPLE_3         3 
#        1885  BINARY_SUBSCR    
#        1886  LOAD_FAST            24  'grayValue'
#        1889  COMPARE_OP            3  '!='
#        1892  POP_JUMP_IF_FALSE  2122  'to 2122'
#
#11350    1895  LOAD_FAST            12  'nxy1'
#        1898  LOAD_FAST            18  'k'
#        1901  BINARY_MULTIPLY  
#        1902  LOAD_FAST            11  'nx1'
#        1905  LOAD_FAST            19  'j'
#        1908  BINARY_MULTIPLY  
#        1909  BINARY_ADD       
#        1910  LOAD_FAST            20  'i'
#        1913  LOAD_CONST            4  2
#        1916  BINARY_ADD       
#        1917  BINARY_ADD       
#        1918  STORE_FAST           27  'curNode'
#
#11351    1921  LOAD_FAST            22  'intfaceNodes'
#        1924  LOAD_ATTR             8  'has_key'
#        1927  LOAD_FAST            27  'curNode'
#        1930  CALL_FUNCTION_1       1 
#        1933  POP_JUMP_IF_TRUE   1949  'to 1949'
#
#11352    1936  LOAD_CONST            2  ''
#        1939  LOAD_FAST            26  'nearIntfaceNodes'
#        1942  LOAD_FAST            27  'curNode'
#        1945  STORE_SUBSCR     
#        1946  JUMP_FORWARD          0  'to 1949'
#      1949_0  COME_FROM                '1946'
#
#11353    1949  LOAD_FAST            12  'nxy1'
#        1952  LOAD_FAST            18  'k'
#        1955  BINARY_MULTIPLY  
#        1956  LOAD_FAST            11  'nx1'
#        1959  LOAD_FAST            19  'j'
#        1962  BINARY_MULTIPLY  
#        1963  BINARY_ADD       
#        1964  LOAD_FAST            20  'i'
#        1967  LOAD_CONST            3  1
#        1970  BINARY_ADD       
#        1971  BINARY_ADD       
#        1972  STORE_FAST           27  'curNode'
#
#11354    1975  LOAD_FAST            22  'intfaceNodes'
#        1978  LOAD_ATTR             8  'has_key'
#        1981  LOAD_FAST            27  'curNode'
#        1984  CALL_FUNCTION_1       1 
#        1987  POP_JUMP_IF_TRUE   2003  'to 2003'
#
#11355    1990  LOAD_CONST            2  ''
#        1993  LOAD_FAST            26  'nearIntfaceNodes'
#        1996  LOAD_FAST            27  'curNode'
#        1999  STORE_SUBSCR     
#        2000  JUMP_FORWARD          0  'to 2003'
#      2003_0  COME_FROM                '2000'
#
#11356    2003  LOAD_FAST            12  'nxy1'
#        2006  LOAD_FAST            18  'k'
#        2009  LOAD_CONST            3  1
#        2012  BINARY_ADD       
#        2013  BINARY_MULTIPLY  
#        2014  LOAD_FAST            11  'nx1'
#        2017  LOAD_FAST            19  'j'
#        2020  BINARY_MULTIPLY  
#        2021  BINARY_ADD       
#        2022  LOAD_FAST            20  'i'
#        2025  LOAD_CONST            4  2
#        2028  BINARY_ADD       
#        2029  BINARY_ADD       
#        2030  STORE_FAST           27  'curNode'
#
#11357    2033  LOAD_FAST            22  'intfaceNodes'
#        2036  LOAD_ATTR             8  'has_key'
#        2039  LOAD_FAST            27  'curNode'
#        2042  CALL_FUNCTION_1       1 
#        2045  POP_JUMP_IF_TRUE   2061  'to 2061'
#
#11358    2048  LOAD_CONST            2  ''
#        2051  LOAD_FAST            26  'nearIntfaceNodes'
#        2054  LOAD_FAST            27  'curNode'
#        2057  STORE_SUBSCR     
#        2058  JUMP_FORWARD          0  'to 2061'
#      2061_0  COME_FROM                '2058'
#
#11359    2061  LOAD_FAST            12  'nxy1'
#        2064  LOAD_FAST            18  'k'
#        2067  LOAD_CONST            3  1
#        2070  BINARY_ADD       
#        2071  BINARY_MULTIPLY  
#        2072  LOAD_FAST            11  'nx1'
#        2075  LOAD_FAST            19  'j'
#        2078  BINARY_MULTIPLY  
#        2079  BINARY_ADD       
#        2080  LOAD_FAST            20  'i'
#        2083  LOAD_CONST            3  1
#        2086  BINARY_ADD       
#        2087  BINARY_ADD       
#        2088  STORE_FAST           27  'curNode'
#
#11360    2091  LOAD_FAST            22  'intfaceNodes'
#        2094  LOAD_ATTR             8  'has_key'
#        2097  LOAD_FAST            27  'curNode'
#        2100  CALL_FUNCTION_1       1 
#        2103  POP_JUMP_IF_TRUE   2122  'to 2122'
#
#11361    2106  LOAD_CONST            2  ''
#        2109  LOAD_FAST            26  'nearIntfaceNodes'
#        2112  LOAD_FAST            27  'curNode'
#        2115  STORE_SUBSCR     
#        2116  JUMP_ABSOLUTE      2122  'to 2122'
#        2119  JUMP_FORWARD          0  'to 2122'
#      2122_0  COME_FROM                '2119'
#
#11363    2122  LOAD_FAST             1  'curVoxelModel'
#        2125  LOAD_FAST            18  'k'
#        2128  LOAD_FAST            19  'j'
#        2131  LOAD_CONST            3  1
#        2134  BINARY_SUBTRACT  
#        2135  LOAD_FAST            20  'i'
#        2138  BUILD_TUPLE_3         3 
#        2141  BINARY_SUBSCR    
#        2142  LOAD_FAST            24  'grayValue'
#        2145  COMPARE_OP            3  '!='
#        2148  POP_JUMP_IF_FALSE  2394  'to 2394'
#
#11364    2151  LOAD_FAST            12  'nxy1'
#        2154  LOAD_FAST            18  'k'
#        2157  BINARY_MULTIPLY  
#        2158  LOAD_FAST            11  'nx1'
#        2161  LOAD_FAST            19  'j'
#        2164  LOAD_CONST            3  1
#        2167  BINARY_ADD       
#        2168  BINARY_MULTIPLY  
#        2169  BINARY_ADD       
#        2170  LOAD_FAST            20  'i'
#        2173  LOAD_CONST            4  2
#        2176  BINARY_ADD       
#        2177  BINARY_ADD       
#        2178  STORE_FAST           27  'curNode'
#
#11365    2181  LOAD_FAST            22  'intfaceNodes'
#        2184  LOAD_ATTR             8  'has_key'
#        2187  LOAD_FAST            27  'curNode'
#        2190  CALL_FUNCTION_1       1 
#        2193  POP_JUMP_IF_TRUE   2209  'to 2209'
#
#11366    2196  LOAD_CONST            2  ''
#        2199  LOAD_FAST            26  'nearIntfaceNodes'
#        2202  LOAD_FAST            27  'curNode'
#        2205  STORE_SUBSCR     
#        2206  JUMP_FORWARD          0  'to 2209'
#      2209_0  COME_FROM                '2206'
#
#11367    2209  LOAD_FAST            12  'nxy1'
#        2212  LOAD_FAST            18  'k'
#        2215  BINARY_MULTIPLY  
#        2216  LOAD_FAST            11  'nx1'
#        2219  LOAD_FAST            19  'j'
#        2222  LOAD_CONST            3  1
#        2225  BINARY_ADD       
#        2226  BINARY_MULTIPLY  
#        2227  BINARY_ADD       
#        2228  LOAD_FAST            20  'i'
#        2231  LOAD_CONST            3  1
#        2234  BINARY_ADD       
#        2235  BINARY_ADD       
#        2236  STORE_FAST           27  'curNode'
#
#11368    2239  LOAD_FAST            22  'intfaceNodes'
#        2242  LOAD_ATTR             8  'has_key'
#        2245  LOAD_FAST            27  'curNode'
#        2248  CALL_FUNCTION_1       1 
#        2251  POP_JUMP_IF_TRUE   2267  'to 2267'
#
#11369    2254  LOAD_CONST            2  ''
#        2257  LOAD_FAST            26  'nearIntfaceNodes'
#        2260  LOAD_FAST            27  'curNode'
#        2263  STORE_SUBSCR     
#        2264  JUMP_FORWARD          0  'to 2267'
#      2267_0  COME_FROM                '2264'
#
#11370    2267  LOAD_FAST            12  'nxy1'
#        2270  LOAD_FAST            18  'k'
#        2273  LOAD_CONST            3  1
#        2276  BINARY_ADD       
#        2277  BINARY_MULTIPLY  
#        2278  LOAD_FAST            11  'nx1'
#        2281  LOAD_FAST            19  'j'
#        2284  LOAD_CONST            3  1
#        2287  BINARY_ADD       
#        2288  BINARY_MULTIPLY  
#        2289  BINARY_ADD       
#        2290  LOAD_FAST            20  'i'
#        2293  LOAD_CONST            4  2
#        2296  BINARY_ADD       
#        2297  BINARY_ADD       
#        2298  STORE_FAST           27  'curNode'
#
#11371    2301  LOAD_FAST            22  'intfaceNodes'
#        2304  LOAD_ATTR             8  'has_key'
#        2307  LOAD_FAST            27  'curNode'
#        2310  CALL_FUNCTION_1       1 
#        2313  POP_JUMP_IF_TRUE   2329  'to 2329'
#
#11372    2316  LOAD_CONST            2  ''
#        2319  LOAD_FAST            26  'nearIntfaceNodes'
#        2322  LOAD_FAST            27  'curNode'
#        2325  STORE_SUBSCR     
#        2326  JUMP_FORWARD          0  'to 2329'
#      2329_0  COME_FROM                '2326'
#
#11373    2329  LOAD_FAST            12  'nxy1'
#        2332  LOAD_FAST            18  'k'
#        2335  LOAD_CONST            3  1
#        2338  BINARY_ADD       
#        2339  BINARY_MULTIPLY  
#        2340  LOAD_FAST            11  'nx1'
#        2343  LOAD_FAST            19  'j'
#        2346  LOAD_CONST            3  1
#        2349  BINARY_ADD       
#        2350  BINARY_MULTIPLY  
#        2351  BINARY_ADD       
#        2352  LOAD_FAST            20  'i'
#        2355  LOAD_CONST            3  1
#        2358  BINARY_ADD       
#        2359  BINARY_ADD       
#        2360  STORE_FAST           27  'curNode'
#
#11374    2363  LOAD_FAST            22  'intfaceNodes'
#        2366  LOAD_ATTR             8  'has_key'
#        2369  LOAD_FAST            27  'curNode'
#        2372  CALL_FUNCTION_1       1 
#        2375  POP_JUMP_IF_TRUE   2394  'to 2394'
#
#11375    2378  LOAD_CONST            2  ''
#        2381  LOAD_FAST            26  'nearIntfaceNodes'
#        2384  LOAD_FAST            27  'curNode'
#        2387  STORE_SUBSCR     
#        2388  JUMP_ABSOLUTE      2394  'to 2394'
#        2391  JUMP_FORWARD          0  'to 2394'
#      2394_0  COME_FROM                '2391'
#
#11377    2394  LOAD_FAST             1  'curVoxelModel'
#        2397  LOAD_FAST            18  'k'
#        2400  LOAD_CONST            3  1
#        2403  BINARY_ADD       
#        2404  LOAD_FAST            19  'j'
#        2407  LOAD_FAST            20  'i'
#        2410  BUILD_TUPLE_3         3 
#        2413  BINARY_SUBSCR    
#        2414  LOAD_FAST            24  'grayValue'
#        2417  COMPARE_OP            3  '!='
#        2420  POP_JUMP_IF_FALSE  2650  'to 2650'
#
#11378    2423  LOAD_FAST            12  'nxy1'
#        2426  LOAD_FAST            18  'k'
#        2429  BINARY_MULTIPLY  
#        2430  LOAD_FAST            11  'nx1'
#        2433  LOAD_FAST            19  'j'
#        2436  BINARY_MULTIPLY  
#        2437  BINARY_ADD       
#        2438  LOAD_FAST            20  'i'
#        2441  LOAD_CONST            3  1
#        2444  BINARY_ADD       
#        2445  BINARY_ADD       
#        2446  STORE_FAST           27  'curNode'
#
#11379    2449  LOAD_FAST            22  'intfaceNodes'
#        2452  LOAD_ATTR             8  'has_key'
#        2455  LOAD_FAST            27  'curNode'
#        2458  CALL_FUNCTION_1       1 
#        2461  POP_JUMP_IF_TRUE   2477  'to 2477'
#
#11380    2464  LOAD_CONST            2  ''
#        2467  LOAD_FAST            26  'nearIntfaceNodes'
#        2470  LOAD_FAST            27  'curNode'
#        2473  STORE_SUBSCR     
#        2474  JUMP_FORWARD          0  'to 2477'
#      2477_0  COME_FROM                '2474'
#
#11381    2477  LOAD_FAST            12  'nxy1'
#        2480  LOAD_FAST            18  'k'
#        2483  BINARY_MULTIPLY  
#        2484  LOAD_FAST            11  'nx1'
#        2487  LOAD_FAST            19  'j'
#        2490  BINARY_MULTIPLY  
#        2491  BINARY_ADD       
#        2492  LOAD_FAST            20  'i'
#        2495  LOAD_CONST            4  2
#        2498  BINARY_ADD       
#        2499  BINARY_ADD       
#        2500  STORE_FAST           27  'curNode'
#
#11382    2503  LOAD_FAST            22  'intfaceNodes'
#        2506  LOAD_ATTR             8  'has_key'
#        2509  LOAD_FAST            27  'curNode'
#        2512  CALL_FUNCTION_1       1 
#        2515  POP_JUMP_IF_TRUE   2531  'to 2531'
#
#11383    2518  LOAD_CONST            2  ''
#        2521  LOAD_FAST            26  'nearIntfaceNodes'
#        2524  LOAD_FAST            27  'curNode'
#        2527  STORE_SUBSCR     
#        2528  JUMP_FORWARD          0  'to 2531'
#      2531_0  COME_FROM                '2528'
#
#11384    2531  LOAD_FAST            12  'nxy1'
#        2534  LOAD_FAST            18  'k'
#        2537  BINARY_MULTIPLY  
#        2538  LOAD_FAST            11  'nx1'
#        2541  LOAD_FAST            19  'j'
#        2544  LOAD_CONST            3  1
#        2547  BINARY_ADD       
#        2548  BINARY_MULTIPLY  
#        2549  BINARY_ADD       
#        2550  LOAD_FAST            20  'i'
#        2553  LOAD_CONST            4  2
#        2556  BINARY_ADD       
#        2557  BINARY_ADD       
#        2558  STORE_FAST           27  'curNode'
#
#11385    2561  LOAD_FAST            22  'intfaceNodes'
#        2564  LOAD_ATTR             8  'has_key'
#        2567  LOAD_FAST            27  'curNode'
#        2570  CALL_FUNCTION_1       1 
#        2573  POP_JUMP_IF_TRUE   2589  'to 2589'
#
#11386    2576  LOAD_CONST            2  ''
#        2579  LOAD_FAST            26  'nearIntfaceNodes'
#        2582  LOAD_FAST            27  'curNode'
#        2585  STORE_SUBSCR     
#        2586  JUMP_FORWARD          0  'to 2589'
#      2589_0  COME_FROM                '2586'
#
#11387    2589  LOAD_FAST            12  'nxy1'
#        2592  LOAD_FAST            18  'k'
#        2595  BINARY_MULTIPLY  
#        2596  LOAD_FAST            11  'nx1'
#        2599  LOAD_FAST            19  'j'
#        2602  LOAD_CONST            3  1
#        2605  BINARY_ADD       
#        2606  BINARY_MULTIPLY  
#        2607  BINARY_ADD       
#        2608  LOAD_FAST            20  'i'
#        2611  LOAD_CONST            3  1
#        2614  BINARY_ADD       
#        2615  BINARY_ADD       
#        2616  STORE_FAST           27  'curNode'
#
#11388    2619  LOAD_FAST            22  'intfaceNodes'
#        2622  LOAD_ATTR             8  'has_key'
#        2625  LOAD_FAST            27  'curNode'
#        2628  CALL_FUNCTION_1       1 
#        2631  POP_JUMP_IF_TRUE   2650  'to 2650'
#
#11389    2634  LOAD_CONST            2  ''
#        2637  LOAD_FAST            26  'nearIntfaceNodes'
#        2640  LOAD_FAST            27  'curNode'
#        2643  STORE_SUBSCR     
#        2644  JUMP_ABSOLUTE      2650  'to 2650'
#        2647  JUMP_FORWARD          0  'to 2650'
#      2650_0  COME_FROM                '2647'
#
#11391    2650  LOAD_FAST             1  'curVoxelModel'
#        2653  LOAD_FAST            18  'k'
#        2656  LOAD_CONST            3  1
#        2659  BINARY_SUBTRACT  
#        2660  LOAD_FAST            19  'j'
#        2663  LOAD_FAST            20  'i'
#        2666  BUILD_TUPLE_3         3 
#        2669  BINARY_SUBSCR    
#        2670  LOAD_FAST            24  'grayValue'
#        2673  COMPARE_OP            3  '!='
#        2676  POP_JUMP_IF_FALSE  2925  'to 2925'
#
#11392    2679  LOAD_FAST            12  'nxy1'
#        2682  LOAD_FAST            18  'k'
#        2685  LOAD_CONST            3  1
#        2688  BINARY_ADD       
#        2689  BINARY_MULTIPLY  
#        2690  LOAD_FAST            11  'nx1'
#        2693  LOAD_FAST            19  'j'
#        2696  BINARY_MULTIPLY  
#        2697  BINARY_ADD       
#        2698  LOAD_FAST            20  'i'
#        2701  LOAD_CONST            3  1
#        2704  BINARY_ADD       
#        2705  BINARY_ADD       
#        2706  STORE_FAST           27  'curNode'
#
#11393    2709  LOAD_FAST            22  'intfaceNodes'
#        2712  LOAD_ATTR             8  'has_key'
#        2715  LOAD_FAST            27  'curNode'
#        2718  CALL_FUNCTION_1       1 
#        2721  POP_JUMP_IF_TRUE   2737  'to 2737'
#
#11394    2724  LOAD_CONST            2  ''
#        2727  LOAD_FAST            26  'nearIntfaceNodes'
#        2730  LOAD_FAST            27  'curNode'
#        2733  STORE_SUBSCR     
#        2734  JUMP_FORWARD          0  'to 2737'
#      2737_0  COME_FROM                '2734'
#
#11395    2737  LOAD_FAST            12  'nxy1'
#        2740  LOAD_FAST            18  'k'
#        2743  LOAD_CONST            3  1
#        2746  BINARY_ADD       
#        2747  BINARY_MULTIPLY  
#        2748  LOAD_FAST            11  'nx1'
#        2751  LOAD_FAST            19  'j'
#        2754  BINARY_MULTIPLY  
#        2755  BINARY_ADD       
#        2756  LOAD_FAST            20  'i'
#        2759  LOAD_CONST            4  2
#        2762  BINARY_ADD       
#        2763  BINARY_ADD       
#        2764  STORE_FAST           27  'curNode'
#
#11396    2767  LOAD_FAST            22  'intfaceNodes'
#        2770  LOAD_ATTR             8  'has_key'
#        2773  LOAD_FAST            27  'curNode'
#        2776  CALL_FUNCTION_1       1 
#        2779  POP_JUMP_IF_TRUE   2795  'to 2795'
#
#11397    2782  LOAD_CONST            2  ''
#        2785  LOAD_FAST            26  'nearIntfaceNodes'
#        2788  LOAD_FAST            27  'curNode'
#        2791  STORE_SUBSCR     
#        2792  JUMP_FORWARD          0  'to 2795'
#      2795_0  COME_FROM                '2792'
#
#11398    2795  LOAD_FAST            12  'nxy1'
#        2798  LOAD_FAST            18  'k'
#        2801  LOAD_CONST            3  1
#        2804  BINARY_ADD       
#        2805  BINARY_MULTIPLY  
#        2806  LOAD_FAST            11  'nx1'
#        2809  LOAD_FAST            19  'j'
#        2812  LOAD_CONST            3  1
#        2815  BINARY_ADD       
#        2816  BINARY_MULTIPLY  
#        2817  BINARY_ADD       
#        2818  LOAD_FAST            20  'i'
#        2821  LOAD_CONST            4  2
#        2824  BINARY_ADD       
#        2825  BINARY_ADD       
#        2826  STORE_FAST           27  'curNode'
#
#11399    2829  LOAD_FAST            22  'intfaceNodes'
#        2832  LOAD_ATTR             8  'has_key'
#        2835  LOAD_FAST            27  'curNode'
#        2838  CALL_FUNCTION_1       1 
#        2841  POP_JUMP_IF_TRUE   2857  'to 2857'
#
#11400    2844  LOAD_CONST            2  ''
#        2847  LOAD_FAST            26  'nearIntfaceNodes'
#        2850  LOAD_FAST            27  'curNode'
#        2853  STORE_SUBSCR     
#        2854  JUMP_FORWARD          0  'to 2857'
#      2857_0  COME_FROM                '2854'
#
#11401    2857  LOAD_FAST            12  'nxy1'
#        2860  LOAD_FAST            18  'k'
#        2863  LOAD_CONST            3  1
#        2866  BINARY_ADD       
#        2867  BINARY_MULTIPLY  
#        2868  LOAD_FAST            11  'nx1'
#        2871  LOAD_FAST            19  'j'
#        2874  LOAD_CONST            3  1
#        2877  BINARY_ADD       
#        2878  BINARY_MULTIPLY  
#        2879  BINARY_ADD       
#        2880  LOAD_FAST            20  'i'
#        2883  LOAD_CONST            3  1
#        2886  BINARY_ADD       
#        2887  BINARY_ADD       
#        2888  STORE_FAST           27  'curNode'
#
#11402    2891  LOAD_FAST            22  'intfaceNodes'
#        2894  LOAD_ATTR             8  'has_key'
#        2897  LOAD_FAST            27  'curNode'
#        2900  CALL_FUNCTION_1       1 
#        2903  POP_JUMP_IF_TRUE   2922  'to 2922'
#
#11403    2906  LOAD_CONST            2  ''
#        2909  LOAD_FAST            26  'nearIntfaceNodes'
#        2912  LOAD_FAST            27  'curNode'
#        2915  STORE_SUBSCR     
#        2916  JUMP_ABSOLUTE      2922  'to 2922'
#        2919  JUMP_ABSOLUTE      2925  'to 2925'
#        2922  JUMP_BACK          1301  'to 1301'
#        2925  JUMP_BACK          1301  'to 1301'
#        2928  POP_BLOCK        
#      2929_0  COME_FROM                '1284'
#        2929  JUMP_BACK          1278  'to 1278'
#        2932  POP_BLOCK        
#      2933_0  COME_FROM                '1261'
#        2933  JUMP_BACK          1182  'to 1182'
#        2936  POP_BLOCK        
#      2937_0  COME_FROM                '1165'
#        2937  JUMP_FORWARD          0  'to 2940'
#      2940_0  COME_FROM                '1165'
#
#11406    2940  LOAD_CONST            2  ''
#        2943  STORE_FAST           17  'noid'
#
#11407    2946  SETUP_LOOP          646  'to 3595'
#        2949  LOAD_GLOBAL           7  'range'
#        2952  LOAD_FAST            10  'nz'
#        2955  LOAD_CONST            3  1
#        2958  BINARY_ADD       
#        2959  CALL_FUNCTION_1       1 
#        2962  GET_ITER         
#        2963  FOR_ITER            628  'to 3594'
#        2966  STORE_FAST           18  'k'
#
#11408    2969  SETUP_LOOP          619  'to 3591'
#        2972  LOAD_GLOBAL           7  'range'
#        2975  LOAD_FAST             9  'ny'
#        2978  LOAD_CONST            3  1
#        2981  BINARY_ADD       
#        2982  CALL_FUNCTION_1       1 
#        2985  GET_ITER         
#        2986  FOR_ITER            601  'to 3590'
#        2989  STORE_FAST           19  'j'
#
#11409    2992  SETUP_LOOP          592  'to 3587'
#        2995  LOAD_GLOBAL           7  'range'
#        2998  LOAD_FAST             8  'nx'
#        3001  LOAD_CONST            3  1
#        3004  BINARY_ADD       
#        3005  CALL_FUNCTION_1       1 
#        3008  GET_ITER         
#        3009  FOR_ITER            574  'to 3586'
#        3012  STORE_FAST           20  'i'
#
#11410    3015  LOAD_FAST            17  'noid'
#        3018  LOAD_CONST            3  1
#        3021  BINARY_ADD       
#        3022  STORE_FAST           17  'noid'
#
#11411    3025  LOAD_FAST            14  'bcidParam'
#        3028  LOAD_CONST            3  1
#        3031  COMPARE_OP            2  '=='
#        3034  POP_JUMP_IF_FALSE  3226  'to 3226'
#
#11412    3037  LOAD_FAST            22  'intfaceNodes'
#        3040  LOAD_ATTR             8  'has_key'
#        3043  LOAD_FAST            17  'noid'
#        3046  CALL_FUNCTION_1       1 
#        3049  POP_JUMP_IF_FALSE  3226  'to 3226'
#
#11413    3052  LOAD_FAST            18  'k'
#        3055  LOAD_CONST            2  ''
#        3058  COMPARE_OP            2  '=='
#        3061  POP_JUMP_IF_FALSE  3080  'to 3080'
#        3064  LOAD_CONST            5  3
#        3067  LOAD_FAST            22  'intfaceNodes'
#        3070  LOAD_FAST            17  'noid'
#        3073  STORE_SUBSCR     
#        3074  JUMP_BACK          3009  'to 3009'
#        3077  JUMP_FORWARD          0  'to 3080'
#      3080_0  COME_FROM                '3077'
#
#11414    3080  LOAD_FAST            18  'k'
#        3083  LOAD_FAST            10  'nz'
#        3086  COMPARE_OP            2  '=='
#        3089  POP_JUMP_IF_FALSE  3108  'to 3108'
#        3092  LOAD_CONST            5  3
#        3095  LOAD_FAST            22  'intfaceNodes'
#        3098  LOAD_FAST            17  'noid'
#        3101  STORE_SUBSCR     
#        3102  JUMP_BACK          3009  'to 3009'
#        3105  JUMP_FORWARD          0  'to 3108'
#      3108_0  COME_FROM                '3105'
#
#11415    3108  LOAD_FAST            19  'j'
#        3111  LOAD_CONST            2  ''
#        3114  COMPARE_OP            2  '=='
#        3117  POP_JUMP_IF_FALSE  3136  'to 3136'
#        3120  LOAD_CONST            4  2
#        3123  LOAD_FAST            22  'intfaceNodes'
#        3126  LOAD_FAST            17  'noid'
#        3129  STORE_SUBSCR     
#        3130  JUMP_BACK          3009  'to 3009'
#        3133  JUMP_FORWARD          0  'to 3136'
#      3136_0  COME_FROM                '3133'
#
#11416    3136  LOAD_FAST            19  'j'
#        3139  LOAD_FAST             9  'ny'
#        3142  COMPARE_OP            2  '=='
#        3145  POP_JUMP_IF_FALSE  3164  'to 3164'
#        3148  LOAD_CONST            4  2
#        3151  LOAD_FAST            22  'intfaceNodes'
#        3154  LOAD_FAST            17  'noid'
#        3157  STORE_SUBSCR     
#        3158  JUMP_BACK          3009  'to 3009'
#        3161  JUMP_FORWARD          0  'to 3164'
#      3164_0  COME_FROM                '3161'
#
#11417    3164  LOAD_FAST            20  'i'
#        3167  LOAD_CONST            2  ''
#        3170  COMPARE_OP            2  '=='
#        3173  POP_JUMP_IF_FALSE  3192  'to 3192'
#        3176  LOAD_CONST            3  1
#        3179  LOAD_FAST            22  'intfaceNodes'
#        3182  LOAD_FAST            17  'noid'
#        3185  STORE_SUBSCR     
#        3186  JUMP_BACK          3009  'to 3009'
#        3189  JUMP_FORWARD          0  'to 3192'
#      3192_0  COME_FROM                '3189'
#
#11418    3192  LOAD_FAST            20  'i'
#        3195  LOAD_FAST             8  'nx'
#        3198  COMPARE_OP            2  '=='
#        3201  POP_JUMP_IF_FALSE  3223  'to 3223'
#        3204  LOAD_CONST            3  1
#        3207  LOAD_FAST            22  'intfaceNodes'
#        3210  LOAD_FAST            17  'noid'
#        3213  STORE_SUBSCR     
#        3214  JUMP_BACK          3009  'to 3009'
#        3217  JUMP_ABSOLUTE      3223  'to 3223'
#        3220  JUMP_ABSOLUTE      3226  'to 3226'
#        3223  JUMP_FORWARD          0  'to 3226'
#      3226_0  COME_FROM                '3223'
#
#11419    3226  LOAD_FAST            14  'bcidParam'
#        3229  LOAD_CONST            2  ''
#        3232  COMPARE_OP            2  '=='
#        3235  POP_JUMP_IF_FALSE  3409  'to 3409'
#
#11420    3238  LOAD_FAST            22  'intfaceNodes'
#        3241  LOAD_ATTR             8  'has_key'
#        3244  LOAD_FAST            17  'noid'
#        3247  CALL_FUNCTION_1       1 
#        3250  POP_JUMP_IF_FALSE  3409  'to 3409'
#
#11421    3253  LOAD_FAST            18  'k'
#        3256  LOAD_CONST            2  ''
#        3259  COMPARE_OP            2  '=='
#        3262  POP_JUMP_IF_FALSE  3278  'to 3278'
#        3265  LOAD_FAST            22  'intfaceNodes'
#        3268  LOAD_FAST            17  'noid'
#        3271  DELETE_SUBSCR    
#        3272  JUMP_BACK          3009  'to 3009'
#        3275  JUMP_FORWARD          0  'to 3278'
#      3278_0  COME_FROM                '3275'
#
#11422    3278  LOAD_FAST            18  'k'
#        3281  LOAD_FAST            10  'nz'
#        3284  COMPARE_OP            2  '=='
#        3287  POP_JUMP_IF_FALSE  3303  'to 3303'
#        3290  LOAD_FAST            22  'intfaceNodes'
#        3293  LOAD_FAST            17  'noid'
#        3296  DELETE_SUBSCR    
#        3297  JUMP_BACK          3009  'to 3009'
#        3300  JUMP_FORWARD          0  'to 3303'
#      3303_0  COME_FROM                '3300'
#
#11423    3303  LOAD_FAST            19  'j'
#        3306  LOAD_CONST            2  ''
#        3309  COMPARE_OP            2  '=='
#        3312  POP_JUMP_IF_FALSE  3328  'to 3328'
#        3315  LOAD_FAST            22  'intfaceNodes'
#        3318  LOAD_FAST            17  'noid'
#        3321  DELETE_SUBSCR    
#        3322  JUMP_BACK          3009  'to 3009'
#        3325  JUMP_FORWARD          0  'to 3328'
#      3328_0  COME_FROM                '3325'
#
#11424    3328  LOAD_FAST            19  'j'
#        3331  LOAD_FAST             9  'ny'
#        3334  COMPARE_OP            2  '=='
#        3337  POP_JUMP_IF_FALSE  3353  'to 3353'
#        3340  LOAD_FAST            22  'intfaceNodes'
#        3343  LOAD_FAST            17  'noid'
#        3346  DELETE_SUBSCR    
#        3347  JUMP_BACK          3009  'to 3009'
#        3350  JUMP_FORWARD          0  'to 3353'
#      3353_0  COME_FROM                '3350'
#
#11425    3353  LOAD_FAST            20  'i'
#        3356  LOAD_CONST            2  ''
#        3359  COMPARE_OP            2  '=='
#        3362  POP_JUMP_IF_FALSE  3378  'to 3378'
#        3365  LOAD_FAST            22  'intfaceNodes'
#        3368  LOAD_FAST            17  'noid'
#        3371  DELETE_SUBSCR    
#        3372  JUMP_BACK          3009  'to 3009'
#        3375  JUMP_FORWARD          0  'to 3378'
#      3378_0  COME_FROM                '3375'
#
#11426    3378  LOAD_FAST            20  'i'
#        3381  LOAD_FAST             8  'nx'
#        3384  COMPARE_OP            2  '=='
#        3387  POP_JUMP_IF_FALSE  3406  'to 3406'
#        3390  LOAD_FAST            22  'intfaceNodes'
#        3393  LOAD_FAST            17  'noid'
#        3396  DELETE_SUBSCR    
#        3397  JUMP_BACK          3009  'to 3009'
#        3400  JUMP_ABSOLUTE      3406  'to 3406'
#        3403  JUMP_ABSOLUTE      3409  'to 3409'
#        3406  JUMP_FORWARD          0  'to 3409'
#      3409_0  COME_FROM                '3406'
#
#11427    3409  LOAD_FAST            13  'nearIntf'
#        3412  POP_JUMP_IF_FALSE  3009  'to 3009'
#        3415  LOAD_FAST            26  'nearIntfaceNodes'
#        3418  LOAD_ATTR             8  'has_key'
#        3421  LOAD_FAST            17  'noid'
#        3424  CALL_FUNCTION_1       1 
#      3427_0  COME_FROM                '3412'
#        3427  POP_JUMP_IF_FALSE  3009  'to 3009'
#
#11428    3430  LOAD_FAST            18  'k'
#        3433  LOAD_CONST            2  ''
#        3436  COMPARE_OP            2  '=='
#        3439  POP_JUMP_IF_FALSE  3455  'to 3455'
#        3442  LOAD_FAST            26  'nearIntfaceNodes'
#        3445  LOAD_FAST            17  'noid'
#        3448  DELETE_SUBSCR    
#        3449  JUMP_BACK          3009  'to 3009'
#        3452  JUMP_FORWARD          0  'to 3455'
#      3455_0  COME_FROM                '3452'
#
#11429    3455  LOAD_FAST            18  'k'
#        3458  LOAD_FAST            10  'nz'
#        3461  COMPARE_OP            2  '=='
#        3464  POP_JUMP_IF_FALSE  3480  'to 3480'
#        3467  LOAD_FAST            26  'nearIntfaceNodes'
#        3470  LOAD_FAST            17  'noid'
#        3473  DELETE_SUBSCR    
#        3474  JUMP_BACK          3009  'to 3009'
#        3477  JUMP_FORWARD          0  'to 3480'
#      3480_0  COME_FROM                '3477'
#
#11430    3480  LOAD_FAST            19  'j'
#        3483  LOAD_CONST            2  ''
#        3486  COMPARE_OP            2  '=='
#        3489  POP_JUMP_IF_FALSE  3505  'to 3505'
#        3492  LOAD_FAST            26  'nearIntfaceNodes'
#        3495  LOAD_FAST            17  'noid'
#        3498  DELETE_SUBSCR    
#        3499  JUMP_BACK          3009  'to 3009'
#        3502  JUMP_FORWARD          0  'to 3505'
#      3505_0  COME_FROM                '3502'
#
#11431    3505  LOAD_FAST            19  'j'
#        3508  LOAD_FAST             9  'ny'
#        3511  COMPARE_OP            2  '=='
#        3514  POP_JUMP_IF_FALSE  3530  'to 3530'
#        3517  LOAD_FAST            26  'nearIntfaceNodes'
#        3520  LOAD_FAST            17  'noid'
#        3523  DELETE_SUBSCR    
#        3524  JUMP_BACK          3009  'to 3009'
#        3527  JUMP_FORWARD          0  'to 3530'
#      3530_0  COME_FROM                '3527'
#
#11432    3530  LOAD_FAST            20  'i'
#        3533  LOAD_CONST            2  ''
#        3536  COMPARE_OP            2  '=='
#        3539  POP_JUMP_IF_FALSE  3555  'to 3555'
#        3542  LOAD_FAST            26  'nearIntfaceNodes'
#        3545  LOAD_FAST            17  'noid'
#        3548  DELETE_SUBSCR    
#        3549  JUMP_BACK          3009  'to 3009'
#        3552  JUMP_FORWARD          0  'to 3555'
#      3555_0  COME_FROM                '3552'
#
#11433    3555  LOAD_FAST            20  'i'
#        3558  LOAD_FAST             8  'nx'
#        3561  COMPARE_OP            2  '=='
#        3564  POP_JUMP_IF_FALSE  3583  'to 3583'
#        3567  LOAD_FAST            26  'nearIntfaceNodes'
#        3570  LOAD_FAST            17  'noid'
#        3573  DELETE_SUBSCR    
#        3574  JUMP_BACK          3009  'to 3009'
#        3577  JUMP_ABSOLUTE      3583  'to 3583'
#        3580  JUMP_BACK          3009  'to 3009'
#        3583  JUMP_BACK          3009  'to 3009'
#        3586  POP_BLOCK        
#      3587_0  COME_FROM                '2992'
#        3587  JUMP_BACK          2986  'to 2986'
#        3590  POP_BLOCK        
#      3591_0  COME_FROM                '2969'
#        3591  JUMP_BACK          2963  'to 2963'
#        3594  POP_BLOCK        
#      3595_0  COME_FROM                '2946'
#
#11434    3595  LOAD_GLOBAL           0  'stdout'
#        3598  LOAD_ATTR             1  'write'
#        3601  LOAD_CONST           11  '     -> Interface node found: %10i             \n'
#        3604  LOAD_GLOBAL          10  'len'
#        3607  LOAD_FAST            22  'intfaceNodes'
#        3610  CALL_FUNCTION_1       1 
#        3613  BINARY_MODULO    
#        3614  CALL_FUNCTION_1       1 
#        3617  POP_TOP          
#
#11435    3618  LOAD_FAST            13  'nearIntf'
#        3621  POP_JUMP_IF_FALSE  3650  'to 3650'
#
#11436    3624  LOAD_GLOBAL           0  'stdout'
#        3627  LOAD_ATTR             1  'write'
#        3630  LOAD_CONST           12  '     -> Near intf node found: %10i             \n'
#        3633  LOAD_GLOBAL          10  'len'
#        3636  LOAD_FAST            26  'nearIntfaceNodes'
#        3639  CALL_FUNCTION_1       1 
#        3642  BINARY_MODULO    
#        3643  CALL_FUNCTION_1       1 
#        3646  POP_TOP          
#        3647  JUMP_FORWARD          0  'to 3650'
#      3650_0  COME_FROM                '3647'
#
#11490    3650  LOAD_CONST            2  ''
#        3653  STORE_FAST           21  'sum'
#
#11491    3656  LOAD_GLOBAL          10  'len'
#        3659  LOAD_FAST            22  'intfaceNodes'
#        3662  CALL_FUNCTION_1       1 
#        3665  STORE_FAST           28  'nnd'
#
#11492    3668  BUILD_MAP_0           0 
#        3671  STORE_FAST           29  'intfaceNeiNodes'
#
#11493    3674  SETUP_LOOP          394  'to 4071'
#        3677  LOAD_FAST            22  'intfaceNodes'
#        3680  GET_ITER         
#        3681  FOR_ITER            386  'to 4070'
#        3684  STORE_FAST           17  'noid'
#
#11494    3687  LOAD_FAST            21  'sum'
#        3690  LOAD_CONST            3  1
#        3693  INPLACE_ADD      
#        3694  STORE_FAST           21  'sum'
#
#11495    3697  LOAD_FAST            21  'sum'
#        3700  LOAD_CONST           13  10000
#        3703  BINARY_MODULO    
#        3704  LOAD_CONST            2  ''
#        3707  COMPARE_OP            2  '=='
#        3710  POP_JUMP_IF_TRUE   3725  'to 3725'
#        3713  LOAD_FAST            21  'sum'
#        3716  LOAD_CONST            3  1
#        3719  COMPARE_OP            2  '=='
#      3722_0  COME_FROM                '3710'
#        3722  POP_JUMP_IF_FALSE  3787  'to 3787'
#
#11496    3725  LOAD_GLOBAL           9  'float'
#        3728  LOAD_FAST            21  'sum'
#        3731  CALL_FUNCTION_1       1 
#        3734  LOAD_GLOBAL           9  'float'
#        3737  LOAD_FAST            28  'nnd'
#        3740  CALL_FUNCTION_1       1 
#        3743  BINARY_DIVIDE    
#        3744  LOAD_CONST            7  100.0
#        3747  BINARY_MULTIPLY  
#        3748  STORE_FAST           23  'progress'
#
#11497    3751  LOAD_GLOBAL           0  'stdout'
#        3754  LOAD_ATTR             1  'write'
#        3757  LOAD_CONST           14  '     -> Setup Neighbor List : %10i %s\r'
#        3760  LOAD_FAST            23  'progress'
#        3763  LOAD_CONST            9  '%'
#        3766  BUILD_TUPLE_2         2 
#        3769  BINARY_MODULO    
#        3770  CALL_FUNCTION_1       1 
#        3773  POP_TOP          
#
#11498    3774  LOAD_GLOBAL           0  'stdout'
#        3777  LOAD_ATTR             2  'flush'
#        3780  CALL_FUNCTION_0       0 
#        3783  POP_TOP          
#        3784  JUMP_FORWARD          0  'to 3787'
#      3787_0  COME_FROM                '3784'
#
#11499    3787  BUILD_LIST_0          0 
#        3790  LOAD_FAST            29  'intfaceNeiNodes'
#        3793  LOAD_FAST            17  'noid'
#        3796  STORE_SUBSCR     
#
#11500    3797  LOAD_FAST            17  'noid'
#        3800  LOAD_FAST            12  'nxy1'
#        3803  BINARY_SUBTRACT  
#        3804  STORE_FAST           27  'curNode'
#
#11501    3807  LOAD_FAST            22  'intfaceNodes'
#        3810  LOAD_ATTR             8  'has_key'
#        3813  LOAD_FAST            27  'curNode'
#        3816  CALL_FUNCTION_1       1 
#        3819  POP_JUMP_IF_FALSE  3842  'to 3842'
#
#11502    3822  LOAD_FAST            29  'intfaceNeiNodes'
#        3825  LOAD_FAST            17  'noid'
#        3828  BINARY_SUBSCR    
#        3829  LOAD_ATTR            11  'append'
#        3832  LOAD_FAST            27  'curNode'
#        3835  CALL_FUNCTION_1       1 
#        3838  POP_TOP          
#        3839  JUMP_FORWARD          0  'to 3842'
#      3842_0  COME_FROM                '3839'
#
#11503    3842  LOAD_FAST            17  'noid'
#        3845  LOAD_FAST            11  'nx1'
#        3848  BINARY_SUBTRACT  
#        3849  STORE_FAST           27  'curNode'
#
#11504    3852  LOAD_FAST            22  'intfaceNodes'
#        3855  LOAD_ATTR             8  'has_key'
#        3858  LOAD_FAST            27  'curNode'
#        3861  CALL_FUNCTION_1       1 
#        3864  POP_JUMP_IF_FALSE  3887  'to 3887'
#
#11505    3867  LOAD_FAST            29  'intfaceNeiNodes'
#        3870  LOAD_FAST            17  'noid'
#        3873  BINARY_SUBSCR    
#        3874  LOAD_ATTR            11  'append'
#        3877  LOAD_FAST            27  'curNode'
#        3880  CALL_FUNCTION_1       1 
#        3883  POP_TOP          
#        3884  JUMP_FORWARD          0  'to 3887'
#      3887_0  COME_FROM                '3884'
#
#11506    3887  LOAD_FAST            17  'noid'
#        3890  LOAD_CONST            3  1
#        3893  BINARY_SUBTRACT  
#        3894  STORE_FAST           27  'curNode'
#
#11507    3897  LOAD_FAST            22  'intfaceNodes'
#        3900  LOAD_ATTR             8  'has_key'
#        3903  LOAD_FAST            27  'curNode'
#        3906  CALL_FUNCTION_1       1 
#        3909  POP_JUMP_IF_FALSE  3932  'to 3932'
#
#11508    3912  LOAD_FAST            29  'intfaceNeiNodes'
#        3915  LOAD_FAST            17  'noid'
#        3918  BINARY_SUBSCR    
#        3919  LOAD_ATTR            11  'append'
#        3922  LOAD_FAST            27  'curNode'
#        3925  CALL_FUNCTION_1       1 
#        3928  POP_TOP          
#        3929  JUMP_FORWARD          0  'to 3932'
#      3932_0  COME_FROM                '3929'
#
#11509    3932  LOAD_FAST            17  'noid'
#        3935  LOAD_CONST            3  1
#        3938  BINARY_ADD       
#        3939  STORE_FAST           27  'curNode'
#
#11510    3942  LOAD_FAST            22  'intfaceNodes'
#        3945  LOAD_ATTR             8  'has_key'
#        3948  LOAD_FAST            27  'curNode'
#        3951  CALL_FUNCTION_1       1 
#        3954  POP_JUMP_IF_FALSE  3977  'to 3977'
#
#11511    3957  LOAD_FAST            29  'intfaceNeiNodes'
#        3960  LOAD_FAST            17  'noid'
#        3963  BINARY_SUBSCR    
#        3964  LOAD_ATTR            11  'append'
#        3967  LOAD_FAST            27  'curNode'
#        3970  CALL_FUNCTION_1       1 
#        3973  POP_TOP          
#        3974  JUMP_FORWARD          0  'to 3977'
#      3977_0  COME_FROM                '3974'
#
#11512    3977  LOAD_FAST            17  'noid'
#        3980  LOAD_FAST            11  'nx1'
#        3983  BINARY_ADD       
#        3984  STORE_FAST           27  'curNode'
#
#11513    3987  LOAD_FAST            22  'intfaceNodes'
#        3990  LOAD_ATTR             8  'has_key'
#        3993  LOAD_FAST            27  'curNode'
#        3996  CALL_FUNCTION_1       1 
#        3999  POP_JUMP_IF_FALSE  4022  'to 4022'
#
#11514    4002  LOAD_FAST            29  'intfaceNeiNodes'
#        4005  LOAD_FAST            17  'noid'
#        4008  BINARY_SUBSCR    
#        4009  LOAD_ATTR            11  'append'
#        4012  LOAD_FAST            27  'curNode'
#        4015  CALL_FUNCTION_1       1 
#        4018  POP_TOP          
#        4019  JUMP_FORWARD          0  'to 4022'
#      4022_0  COME_FROM                '4019'
#
#11515    4022  LOAD_FAST            17  'noid'
#        4025  LOAD_FAST            12  'nxy1'
#        4028  BINARY_ADD       
#        4029  STORE_FAST           27  'curNode'
#
#11516    4032  LOAD_FAST            22  'intfaceNodes'
#        4035  LOAD_ATTR             8  'has_key'
#        4038  LOAD_FAST            27  'curNode'
#        4041  CALL_FUNCTION_1       1 
#        4044  POP_JUMP_IF_FALSE  3681  'to 3681'
#
#11517    4047  LOAD_FAST            29  'intfaceNeiNodes'
#        4050  LOAD_FAST            17  'noid'
#        4053  BINARY_SUBSCR    
#        4054  LOAD_ATTR            11  'append'
#        4057  LOAD_FAST            27  'curNode'
#        4060  CALL_FUNCTION_1       1 
#        4063  POP_TOP          
#        4064  JUMP_BACK          3681  'to 3681'
#
#11527    4067  JUMP_BACK          3681  'to 3681'
#        4070  POP_BLOCK        
#      4071_0  COME_FROM                '3674'
#
#11531    4071  LOAD_FAST            13  'nearIntf'
#        4074  POP_JUMP_IF_FALSE  4591  'to 4591'
#
#11532    4077  LOAD_CONST            2  ''
#        4080  STORE_FAST           21  'sum'
#
#11533    4083  LOAD_GLOBAL          10  'len'
#        4086  LOAD_FAST            26  'nearIntfaceNodes'
#        4089  CALL_FUNCTION_1       1 
#        4092  STORE_FAST           28  'nnd'
#
#11534    4095  BUILD_MAP_0           0 
#        4098  STORE_FAST           30  'nearIntfaceNeiNodes'
#
#11535    4101  SETUP_LOOP          487  'to 4591'
#        4104  LOAD_FAST            26  'nearIntfaceNodes'
#        4107  GET_ITER         
#        4108  FOR_ITER            476  'to 4587'
#        4111  STORE_FAST           17  'noid'
#
#11536    4114  LOAD_FAST            21  'sum'
#        4117  LOAD_CONST            3  1
#        4120  INPLACE_ADD      
#        4121  STORE_FAST           21  'sum'
#
#11537    4124  LOAD_FAST            21  'sum'
#        4127  LOAD_CONST           13  10000
#        4130  BINARY_MODULO    
#        4131  LOAD_CONST            2  ''
#        4134  COMPARE_OP            2  '=='
#        4137  POP_JUMP_IF_TRUE   4152  'to 4152'
#        4140  LOAD_FAST            21  'sum'
#        4143  LOAD_CONST            3  1
#        4146  COMPARE_OP            2  '=='
#      4149_0  COME_FROM                '4137'
#        4149  POP_JUMP_IF_FALSE  4214  'to 4214'
#
#11538    4152  LOAD_GLOBAL           9  'float'
#        4155  LOAD_FAST            21  'sum'
#        4158  CALL_FUNCTION_1       1 
#        4161  LOAD_GLOBAL           9  'float'
#        4164  LOAD_FAST            28  'nnd'
#        4167  CALL_FUNCTION_1       1 
#        4170  BINARY_DIVIDE    
#        4171  LOAD_CONST            7  100.0
#        4174  BINARY_MULTIPLY  
#        4175  STORE_FAST           23  'progress'
#
#11539    4178  LOAD_GLOBAL           0  'stdout'
#        4181  LOAD_ATTR             1  'write'
#        4184  LOAD_CONST           14  '     -> Setup Neighbor List : %10i %s\r'
#        4187  LOAD_FAST            23  'progress'
#        4190  LOAD_CONST            9  '%'
#        4193  BUILD_TUPLE_2         2 
#        4196  BINARY_MODULO    
#        4197  CALL_FUNCTION_1       1 
#        4200  POP_TOP          
#
#11540    4201  LOAD_GLOBAL           0  'stdout'
#        4204  LOAD_ATTR             2  'flush'
#        4207  CALL_FUNCTION_0       0 
#        4210  POP_TOP          
#        4211  JUMP_FORWARD          0  'to 4214'
#      4214_0  COME_FROM                '4211'
#
#11541    4214  BUILD_LIST_0          0 
#        4217  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4220  LOAD_FAST            17  'noid'
#        4223  STORE_SUBSCR     
#
#11542    4224  LOAD_FAST            17  'noid'
#        4227  LOAD_FAST            12  'nxy1'
#        4230  BINARY_SUBTRACT  
#        4231  STORE_FAST           27  'curNode'
#
#11543    4234  LOAD_FAST            26  'nearIntfaceNodes'
#        4237  LOAD_ATTR             8  'has_key'
#        4240  LOAD_FAST            27  'curNode'
#        4243  CALL_FUNCTION_1       1 
#        4246  POP_JUMP_IF_TRUE   4264  'to 4264'
#        4249  LOAD_FAST            22  'intfaceNodes'
#        4252  LOAD_ATTR             8  'has_key'
#        4255  LOAD_FAST            27  'curNode'
#        4258  CALL_FUNCTION_1       1 
#      4261_0  COME_FROM                '4246'
#        4261  POP_JUMP_IF_FALSE  4284  'to 4284'
#
#11544    4264  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4267  LOAD_FAST            17  'noid'
#        4270  BINARY_SUBSCR    
#        4271  LOAD_ATTR            11  'append'
#        4274  LOAD_FAST            27  'curNode'
#        4277  CALL_FUNCTION_1       1 
#        4280  POP_TOP          
#        4281  JUMP_FORWARD          0  'to 4284'
#      4284_0  COME_FROM                '4281'
#
#11545    4284  LOAD_FAST            17  'noid'
#        4287  LOAD_FAST            11  'nx1'
#        4290  BINARY_SUBTRACT  
#        4291  STORE_FAST           27  'curNode'
#
#11546    4294  LOAD_FAST            26  'nearIntfaceNodes'
#        4297  LOAD_ATTR             8  'has_key'
#        4300  LOAD_FAST            27  'curNode'
#        4303  CALL_FUNCTION_1       1 
#        4306  POP_JUMP_IF_TRUE   4324  'to 4324'
#        4309  LOAD_FAST            22  'intfaceNodes'
#        4312  LOAD_ATTR             8  'has_key'
#        4315  LOAD_FAST            27  'curNode'
#        4318  CALL_FUNCTION_1       1 
#      4321_0  COME_FROM                '4306'
#        4321  POP_JUMP_IF_FALSE  4344  'to 4344'
#
#11547    4324  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4327  LOAD_FAST            17  'noid'
#        4330  BINARY_SUBSCR    
#        4331  LOAD_ATTR            11  'append'
#        4334  LOAD_FAST            27  'curNode'
#        4337  CALL_FUNCTION_1       1 
#        4340  POP_TOP          
#        4341  JUMP_FORWARD          0  'to 4344'
#      4344_0  COME_FROM                '4341'
#
#11548    4344  LOAD_FAST            17  'noid'
#        4347  LOAD_CONST            3  1
#        4350  BINARY_SUBTRACT  
#        4351  STORE_FAST           27  'curNode'
#
#11549    4354  LOAD_FAST            26  'nearIntfaceNodes'
#        4357  LOAD_ATTR             8  'has_key'
#        4360  LOAD_FAST            27  'curNode'
#        4363  CALL_FUNCTION_1       1 
#        4366  POP_JUMP_IF_TRUE   4384  'to 4384'
#        4369  LOAD_FAST            22  'intfaceNodes'
#        4372  LOAD_ATTR             8  'has_key'
#        4375  LOAD_FAST            27  'curNode'
#        4378  CALL_FUNCTION_1       1 
#      4381_0  COME_FROM                '4366'
#        4381  POP_JUMP_IF_FALSE  4404  'to 4404'
#
#11550    4384  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4387  LOAD_FAST            17  'noid'
#        4390  BINARY_SUBSCR    
#        4391  LOAD_ATTR            11  'append'
#        4394  LOAD_FAST            27  'curNode'
#        4397  CALL_FUNCTION_1       1 
#        4400  POP_TOP          
#        4401  JUMP_FORWARD          0  'to 4404'
#      4404_0  COME_FROM                '4401'
#
#11551    4404  LOAD_FAST            17  'noid'
#        4407  LOAD_CONST            3  1
#        4410  BINARY_ADD       
#        4411  STORE_FAST           27  'curNode'
#
#11552    4414  LOAD_FAST            26  'nearIntfaceNodes'
#        4417  LOAD_ATTR             8  'has_key'
#        4420  LOAD_FAST            27  'curNode'
#        4423  CALL_FUNCTION_1       1 
#        4426  POP_JUMP_IF_TRUE   4444  'to 4444'
#        4429  LOAD_FAST            22  'intfaceNodes'
#        4432  LOAD_ATTR             8  'has_key'
#        4435  LOAD_FAST            27  'curNode'
#        4438  CALL_FUNCTION_1       1 
#      4441_0  COME_FROM                '4426'
#        4441  POP_JUMP_IF_FALSE  4464  'to 4464'
#
#11553    4444  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4447  LOAD_FAST            17  'noid'
#        4450  BINARY_SUBSCR    
#        4451  LOAD_ATTR            11  'append'
#        4454  LOAD_FAST            27  'curNode'
#        4457  CALL_FUNCTION_1       1 
#        4460  POP_TOP          
#        4461  JUMP_FORWARD          0  'to 4464'
#      4464_0  COME_FROM                '4461'
#
#11554    4464  LOAD_FAST            17  'noid'
#        4467  LOAD_FAST            11  'nx1'
#        4470  BINARY_ADD       
#        4471  STORE_FAST           27  'curNode'
#
#11555    4474  LOAD_FAST            26  'nearIntfaceNodes'
#        4477  LOAD_ATTR             8  'has_key'
#        4480  LOAD_FAST            27  'curNode'
#        4483  CALL_FUNCTION_1       1 
#        4486  POP_JUMP_IF_TRUE   4504  'to 4504'
#        4489  LOAD_FAST            22  'intfaceNodes'
#        4492  LOAD_ATTR             8  'has_key'
#        4495  LOAD_FAST            27  'curNode'
#        4498  CALL_FUNCTION_1       1 
#      4501_0  COME_FROM                '4486'
#        4501  POP_JUMP_IF_FALSE  4524  'to 4524'
#
#11556    4504  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4507  LOAD_FAST            17  'noid'
#        4510  BINARY_SUBSCR    
#        4511  LOAD_ATTR            11  'append'
#        4514  LOAD_FAST            27  'curNode'
#        4517  CALL_FUNCTION_1       1 
#        4520  POP_TOP          
#        4521  JUMP_FORWARD          0  'to 4524'
#      4524_0  COME_FROM                '4521'
#
#11557    4524  LOAD_FAST            17  'noid'
#        4527  LOAD_FAST            12  'nxy1'
#        4530  BINARY_ADD       
#        4531  STORE_FAST           27  'curNode'
#
#11558    4534  LOAD_FAST            26  'nearIntfaceNodes'
#        4537  LOAD_ATTR             8  'has_key'
#        4540  LOAD_FAST            27  'curNode'
#        4543  CALL_FUNCTION_1       1 
#        4546  POP_JUMP_IF_TRUE   4564  'to 4564'
#        4549  LOAD_FAST            22  'intfaceNodes'
#        4552  LOAD_ATTR             8  'has_key'
#        4555  LOAD_FAST            27  'curNode'
#        4558  CALL_FUNCTION_1       1 
#      4561_0  COME_FROM                '4546'
#        4561  POP_JUMP_IF_FALSE  4108  'to 4108'
#
#11559    4564  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        4567  LOAD_FAST            17  'noid'
#        4570  BINARY_SUBSCR    
#        4571  LOAD_ATTR            11  'append'
#        4574  LOAD_FAST            27  'curNode'
#        4577  CALL_FUNCTION_1       1 
#        4580  POP_TOP          
#        4581  JUMP_BACK          4108  'to 4108'
#        4584  JUMP_BACK          4108  'to 4108'
#        4587  POP_BLOCK        
#      4588_0  COME_FROM                '4101'
#        4588  JUMP_FORWARD          0  'to 4591'
#      4591_0  COME_FROM                '4101'
#
#11563    4591  LOAD_GLOBAL           5  'int'
#        4594  LOAD_FAST             3  'smoothParam'
#        4597  LOAD_CONST            2  ''
#        4600  BINARY_SUBSCR    
#        4601  CALL_FUNCTION_1       1 
#        4604  STORE_FAST           31  'nIt'
#
#11564    4607  LOAD_GLOBAL           9  'float'
#        4610  LOAD_FAST             3  'smoothParam'
#        4613  LOAD_CONST            3  1
#        4616  BINARY_SUBSCR    
#        4617  CALL_FUNCTION_1       1 
#        4620  STORE_FAST           32  'lam'
#
#11565    4623  LOAD_GLOBAL           9  'float'
#        4626  LOAD_FAST             3  'smoothParam'
#        4629  LOAD_CONST            3  1
#        4632  BINARY_SUBSCR    
#        4633  CALL_FUNCTION_1       1 
#        4636  LOAD_CONST           15  2.0
#        4639  BINARY_DIVIDE    
#        4640  STORE_FAST           33  'lam2'
#
#11566    4643  LOAD_CONST           16  1.0
#        4646  LOAD_GLOBAL           9  'float'
#        4649  LOAD_FAST             3  'smoothParam'
#        4652  LOAD_CONST            4  2
#        4655  BINARY_SUBSCR    
#        4656  CALL_FUNCTION_1       1 
#        4659  LOAD_CONST           16  1.0
#        4662  LOAD_FAST            32  'lam'
#        4665  BINARY_DIVIDE    
#        4666  BINARY_SUBTRACT  
#        4667  BINARY_DIVIDE    
#        4668  STORE_FAST           34  'mu'
#
#11567    4671  LOAD_CONST           16  1.0
#        4674  LOAD_GLOBAL           9  'float'
#        4677  LOAD_FAST             3  'smoothParam'
#        4680  LOAD_CONST            4  2
#        4683  BINARY_SUBSCR    
#        4684  CALL_FUNCTION_1       1 
#        4687  LOAD_CONST           16  1.0
#        4690  LOAD_FAST            33  'lam2'
#        4693  BINARY_DIVIDE    
#        4694  BINARY_SUBTRACT  
#        4695  BINARY_DIVIDE    
#        4696  STORE_FAST           35  'mu2'
#
#11568    4699  LOAD_CONST           16  1.0
#        4702  STORE_FAST           36  'weight'
#
#11569    4705  SETUP_LOOP          875  'to 5583'
#        4708  LOAD_GLOBAL           7  'range'
#        4711  LOAD_FAST            31  'nIt'
#        4714  CALL_FUNCTION_1       1 
#        4717  GET_ITER         
#        4718  FOR_ITER            861  'to 5582'
#        4721  STORE_FAST           37  'it'
#
#11570    4724  LOAD_CONST            2  ''
#        4727  STORE_FAST           21  'sum'
#
#11571    4730  LOAD_GLOBAL          10  'len'
#        4733  LOAD_FAST            22  'intfaceNodes'
#        4736  CALL_FUNCTION_1       1 
#        4739  STORE_FAST           28  'nnd'
#
#11572    4742  LOAD_GLOBAL           0  'stdout'
#        4745  LOAD_ATTR             1  'write'
#        4748  LOAD_CONST           17  '     -> Smoothing Iteration :  %4i/%4i  \r'
#        4751  LOAD_FAST            37  'it'
#        4754  LOAD_CONST            3  1
#        4757  BINARY_ADD       
#        4758  LOAD_FAST            31  'nIt'
#        4761  BUILD_TUPLE_2         2 
#        4764  BINARY_MODULO    
#        4765  CALL_FUNCTION_1       1 
#        4768  POP_TOP          
#
#11573    4769  LOAD_GLOBAL           0  'stdout'
#        4772  LOAD_ATTR             2  'flush'
#        4775  CALL_FUNCTION_0       0 
#        4778  POP_TOP          
#
#11575    4779  SETUP_LOOP          389  'to 5171'
#        4782  LOAD_FAST            22  'intfaceNodes'
#        4785  GET_ITER         
#        4786  FOR_ITER            381  'to 5170'
#        4789  STORE_FAST           17  'noid'
#
#11576    4792  LOAD_FAST            22  'intfaceNodes'
#        4795  LOAD_FAST            17  'noid'
#        4798  BINARY_SUBSCR    
#        4799  STORE_FAST           38  'bcid'
#
#11577    4802  LOAD_FAST            16  'nodeCoordPre'
#        4805  LOAD_FAST            17  'noid'
#        4808  BINARY_SUBSCR    
#        4809  LOAD_CONST            2  ''
#        4812  BINARY_SUBSCR    
#        4813  STORE_FAST           39  'xi'
#
#11578    4816  LOAD_FAST            16  'nodeCoordPre'
#        4819  LOAD_FAST            17  'noid'
#        4822  BINARY_SUBSCR    
#        4823  LOAD_CONST            3  1
#        4826  BINARY_SUBSCR    
#        4827  STORE_FAST           40  'yi'
#
#11579    4830  LOAD_FAST            16  'nodeCoordPre'
#        4833  LOAD_FAST            17  'noid'
#        4836  BINARY_SUBSCR    
#        4837  LOAD_CONST            4  2
#        4840  BINARY_SUBSCR    
#        4841  STORE_FAST           41  'zi'
#
#11580    4844  LOAD_CONST           18  ''
#        4847  STORE_FAST           42  'dxi'
#        4850  LOAD_CONST           18  ''
#        4853  STORE_FAST           43  'dyi'
#        4856  LOAD_CONST           18  ''
#        4859  STORE_FAST           44  'dzi'
#
#11581    4862  LOAD_FAST            36  'weight'
#        4865  LOAD_GLOBAL           9  'float'
#        4868  LOAD_GLOBAL          10  'len'
#        4871  LOAD_FAST            29  'intfaceNeiNodes'
#        4874  LOAD_FAST            17  'noid'
#        4877  BINARY_SUBSCR    
#        4878  CALL_FUNCTION_1       1 
#        4881  CALL_FUNCTION_1       1 
#        4884  BINARY_DIVIDE    
#        4885  STORE_FAST           45  'wij'
#
#11582    4888  SETUP_LOOP          159  'to 5050'
#        4891  LOAD_FAST            29  'intfaceNeiNodes'
#        4894  LOAD_FAST            17  'noid'
#        4897  BINARY_SUBSCR    
#        4898  GET_ITER         
#        4899  FOR_ITER            147  'to 5049'
#        4902  STORE_FAST           46  'nnoid'
#
#11583    4905  LOAD_FAST            16  'nodeCoordPre'
#        4908  LOAD_FAST            46  'nnoid'
#        4911  BINARY_SUBSCR    
#        4912  LOAD_CONST            2  ''
#        4915  BINARY_SUBSCR    
#        4916  STORE_FAST           47  'xj'
#
#11584    4919  LOAD_FAST            16  'nodeCoordPre'
#        4922  LOAD_FAST            46  'nnoid'
#        4925  BINARY_SUBSCR    
#        4926  LOAD_CONST            3  1
#        4929  BINARY_SUBSCR    
#        4930  STORE_FAST           48  'yj'
#
#11585    4933  LOAD_FAST            16  'nodeCoordPre'
#        4936  LOAD_FAST            46  'nnoid'
#        4939  BINARY_SUBSCR    
#        4940  LOAD_CONST            4  2
#        4943  BINARY_SUBSCR    
#        4944  STORE_FAST           49  'zj'
#
#11586    4947  LOAD_FAST            38  'bcid'
#        4950  LOAD_CONST            3  1
#        4953  COMPARE_OP            3  '!='
#        4956  POP_JUMP_IF_FALSE  4980  'to 4980'
#        4959  LOAD_FAST            42  'dxi'
#        4962  LOAD_FAST            45  'wij'
#        4965  LOAD_FAST            39  'xi'
#        4968  LOAD_FAST            47  'xj'
#        4971  BINARY_SUBTRACT  
#        4972  BINARY_MULTIPLY  
#        4973  INPLACE_ADD      
#        4974  STORE_FAST           42  'dxi'
#        4977  JUMP_FORWARD          0  'to 4980'
#      4980_0  COME_FROM                '4977'
#
#11587    4980  LOAD_FAST            38  'bcid'
#        4983  LOAD_CONST            4  2
#        4986  COMPARE_OP            3  '!='
#        4989  POP_JUMP_IF_FALSE  5013  'to 5013'
#        4992  LOAD_FAST            43  'dyi'
#        4995  LOAD_FAST            45  'wij'
#        4998  LOAD_FAST            40  'yi'
#        5001  LOAD_FAST            48  'yj'
#        5004  BINARY_SUBTRACT  
#        5005  BINARY_MULTIPLY  
#        5006  INPLACE_ADD      
#        5007  STORE_FAST           43  'dyi'
#        5010  JUMP_FORWARD          0  'to 5013'
#      5013_0  COME_FROM                '5010'
#
#11588    5013  LOAD_FAST            38  'bcid'
#        5016  LOAD_CONST            5  3
#        5019  COMPARE_OP            3  '!='
#        5022  POP_JUMP_IF_FALSE  4899  'to 4899'
#        5025  LOAD_FAST            44  'dzi'
#        5028  LOAD_FAST            45  'wij'
#        5031  LOAD_FAST            41  'zi'
#        5034  LOAD_FAST            49  'zj'
#        5037  BINARY_SUBTRACT  
#        5038  BINARY_MULTIPLY  
#        5039  INPLACE_ADD      
#        5040  STORE_FAST           44  'dzi'
#        5043  JUMP_BACK          4899  'to 4899'
#        5046  JUMP_BACK          4899  'to 4899'
#        5049  POP_BLOCK        
#      5050_0  COME_FROM                '4888'
#
#11589    5050  LOAD_FAST            37  'it'
#        5053  LOAD_CONST            4  2
#        5056  BINARY_MODULO    
#        5057  LOAD_CONST            2  ''
#        5060  COMPARE_OP            2  '=='
#        5063  POP_JUMP_IF_TRUE   5078  'to 5078'
#        5066  LOAD_FAST            37  'it'
#        5069  LOAD_CONST            2  ''
#        5072  COMPARE_OP            2  '=='
#      5075_0  COME_FROM                '5063'
#        5075  POP_JUMP_IF_FALSE  5124  'to 5124'
#
#11590    5078  LOAD_FAST            39  'xi'
#        5081  LOAD_FAST            32  'lam'
#        5084  LOAD_FAST            42  'dxi'
#        5087  BINARY_MULTIPLY  
#        5088  BINARY_ADD       
#        5089  LOAD_FAST            40  'yi'
#        5092  LOAD_FAST            32  'lam'
#        5095  LOAD_FAST            43  'dyi'
#        5098  BINARY_MULTIPLY  
#        5099  BINARY_ADD       
#        5100  LOAD_FAST            41  'zi'
#        5103  LOAD_FAST            32  'lam'
#        5106  LOAD_FAST            44  'dzi'
#        5109  BINARY_MULTIPLY  
#        5110  BINARY_ADD       
#        5111  BUILD_TUPLE_3         3 
#        5114  LOAD_FAST            15  'nodeCoord'
#        5117  LOAD_FAST            17  'noid'
#        5120  STORE_SUBSCR     
#        5121  JUMP_BACK          4786  'to 4786'
#
#11592    5124  LOAD_FAST            39  'xi'
#        5127  LOAD_FAST            34  'mu'
#        5130  LOAD_FAST            42  'dxi'
#        5133  BINARY_MULTIPLY  
#        5134  BINARY_ADD       
#        5135  LOAD_FAST            40  'yi'
#        5138  LOAD_FAST            34  'mu'
#        5141  LOAD_FAST            43  'dyi'
#        5144  BINARY_MULTIPLY  
#        5145  BINARY_ADD       
#        5146  LOAD_FAST            41  'zi'
#        5149  LOAD_FAST            34  'mu'
#        5152  LOAD_FAST            44  'dzi'
#        5155  BINARY_MULTIPLY  
#        5156  BINARY_ADD       
#        5157  BUILD_TUPLE_3         3 
#        5160  LOAD_FAST            15  'nodeCoord'
#        5163  LOAD_FAST            17  'noid'
#        5166  STORE_SUBSCR     
#        5167  JUMP_BACK          4786  'to 4786'
#        5170  POP_BLOCK        
#      5171_0  COME_FROM                '4779'
#
#11594    5171  LOAD_FAST            13  'nearIntf'
#        5174  POP_JUMP_IF_FALSE  5548  'to 5548'
#
#11596    5177  SETUP_LOOP          334  'to 5514'
#        5180  LOAD_FAST            26  'nearIntfaceNodes'
#        5183  GET_ITER         
#        5184  FOR_ITER            326  'to 5513'
#        5187  STORE_FAST           17  'noid'
#
#11597    5190  LOAD_FAST            16  'nodeCoordPre'
#        5193  LOAD_FAST            17  'noid'
#        5196  BINARY_SUBSCR    
#        5197  LOAD_CONST            2  ''
#        5200  BINARY_SUBSCR    
#        5201  STORE_FAST           39  'xi'
#
#11598    5204  LOAD_FAST            16  'nodeCoordPre'
#        5207  LOAD_FAST            17  'noid'
#        5210  BINARY_SUBSCR    
#        5211  LOAD_CONST            3  1
#        5214  BINARY_SUBSCR    
#        5215  STORE_FAST           40  'yi'
#
#11599    5218  LOAD_FAST            16  'nodeCoordPre'
#        5221  LOAD_FAST            17  'noid'
#        5224  BINARY_SUBSCR    
#        5225  LOAD_CONST            4  2
#        5228  BINARY_SUBSCR    
#        5229  STORE_FAST           41  'zi'
#
#11600    5232  LOAD_CONST           18  ''
#        5235  STORE_FAST           42  'dxi'
#        5238  LOAD_CONST           18  ''
#        5241  STORE_FAST           43  'dyi'
#        5244  LOAD_CONST           18  ''
#        5247  STORE_FAST           44  'dzi'
#
#11601    5250  LOAD_FAST            36  'weight'
#        5253  LOAD_GLOBAL           9  'float'
#        5256  LOAD_GLOBAL          10  'len'
#        5259  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        5262  LOAD_FAST            17  'noid'
#        5265  BINARY_SUBSCR    
#        5266  CALL_FUNCTION_1       1 
#        5269  CALL_FUNCTION_1       1 
#        5272  BINARY_DIVIDE    
#        5273  STORE_FAST           45  'wij'
#
#11602    5276  SETUP_LOOP          114  'to 5393'
#        5279  LOAD_FAST            30  'nearIntfaceNeiNodes'
#        5282  LOAD_FAST            17  'noid'
#        5285  BINARY_SUBSCR    
#        5286  GET_ITER         
#        5287  FOR_ITER            102  'to 5392'
#        5290  STORE_FAST           46  'nnoid'
#
#11603    5293  LOAD_FAST            16  'nodeCoordPre'
#        5296  LOAD_FAST            46  'nnoid'
#        5299  BINARY_SUBSCR    
#        5300  LOAD_CONST            2  ''
#        5303  BINARY_SUBSCR    
#        5304  STORE_FAST           47  'xj'
#
#11604    5307  LOAD_FAST            16  'nodeCoordPre'
#        5310  LOAD_FAST            46  'nnoid'
#        5313  BINARY_SUBSCR    
#        5314  LOAD_CONST            3  1
#        5317  BINARY_SUBSCR    
#        5318  STORE_FAST           48  'yj'
#
#11605    5321  LOAD_FAST            16  'nodeCoordPre'
#        5324  LOAD_FAST            46  'nnoid'
#        5327  BINARY_SUBSCR    
#        5328  LOAD_CONST            4  2
#        5331  BINARY_SUBSCR    
#        5332  STORE_FAST           49  'zj'
#
#11606    5335  LOAD_FAST            42  'dxi'
#        5338  LOAD_FAST            45  'wij'
#        5341  LOAD_FAST            39  'xi'
#        5344  LOAD_FAST            47  'xj'
#        5347  BINARY_SUBTRACT  
#        5348  BINARY_MULTIPLY  
#        5349  INPLACE_ADD      
#        5350  STORE_FAST           42  'dxi'
#
#11607    5353  LOAD_FAST            43  'dyi'
#        5356  LOAD_FAST            45  'wij'
#        5359  LOAD_FAST            40  'yi'
#        5362  LOAD_FAST            48  'yj'
#        5365  BINARY_SUBTRACT  
#        5366  BINARY_MULTIPLY  
#        5367  INPLACE_ADD      
#        5368  STORE_FAST           43  'dyi'
#
#11608    5371  LOAD_FAST            44  'dzi'
#        5374  LOAD_FAST            45  'wij'
#        5377  LOAD_FAST            41  'zi'
#        5380  LOAD_FAST            49  'zj'
#        5383  BINARY_SUBTRACT  
#        5384  BINARY_MULTIPLY  
#        5385  INPLACE_ADD      
#        5386  STORE_FAST           44  'dzi'
#        5389  JUMP_BACK          5287  'to 5287'
#        5392  POP_BLOCK        
#      5393_0  COME_FROM                '5276'
#
#11609    5393  LOAD_FAST            37  'it'
#        5396  LOAD_CONST            4  2
#        5399  BINARY_MODULO    
#        5400  LOAD_CONST            2  ''
#        5403  COMPARE_OP            2  '=='
#        5406  POP_JUMP_IF_TRUE   5421  'to 5421'
#        5409  LOAD_FAST            37  'it'
#        5412  LOAD_CONST            2  ''
#        5415  COMPARE_OP            2  '=='
#      5418_0  COME_FROM                '5406'
#        5418  POP_JUMP_IF_FALSE  5467  'to 5467'
#
#11610    5421  LOAD_FAST            39  'xi'
#        5424  LOAD_FAST            33  'lam2'
#        5427  LOAD_FAST            42  'dxi'
#        5430  BINARY_MULTIPLY  
#        5431  BINARY_ADD       
#        5432  LOAD_FAST            40  'yi'
#        5435  LOAD_FAST            33  'lam2'
#        5438  LOAD_FAST            43  'dyi'
#        5441  BINARY_MULTIPLY  
#        5442  BINARY_ADD       
#        5443  LOAD_FAST            41  'zi'
#        5446  LOAD_FAST            33  'lam2'
#        5449  LOAD_FAST            44  'dzi'
#        5452  BINARY_MULTIPLY  
#        5453  BINARY_ADD       
#        5454  BUILD_TUPLE_3         3 
#        5457  LOAD_FAST            15  'nodeCoord'
#        5460  LOAD_FAST            17  'noid'
#        5463  STORE_SUBSCR     
#        5464  JUMP_BACK          5184  'to 5184'
#
#11612    5467  LOAD_FAST            39  'xi'
#        5470  LOAD_FAST            35  'mu2'
#        5473  LOAD_FAST            42  'dxi'
#        5476  BINARY_MULTIPLY  
#        5477  BINARY_ADD       
#        5478  LOAD_FAST            40  'yi'
#        5481  LOAD_FAST            35  'mu2'
#        5484  LOAD_FAST            43  'dyi'
#        5487  BINARY_MULTIPLY  
#        5488  BINARY_ADD       
#        5489  LOAD_FAST            41  'zi'
#        5492  LOAD_FAST            35  'mu2'
#        5495  LOAD_FAST            44  'dzi'
#        5498  BINARY_MULTIPLY  
#        5499  BINARY_ADD       
#        5500  BUILD_TUPLE_3         3 
#        5503  LOAD_FAST            15  'nodeCoord'
#        5506  LOAD_FAST            17  'noid'
#        5509  STORE_SUBSCR     
#        5510  JUMP_BACK          5184  'to 5184'
#        5513  POP_BLOCK        
#      5514_0  COME_FROM                '5177'
#
#11614    5514  SETUP_LOOP           31  'to 5548'
#        5517  LOAD_FAST            26  'nearIntfaceNodes'
#        5520  GET_ITER         
#        5521  FOR_ITER             20  'to 5544'
#        5524  STORE_FAST           17  'noid'
#
#11615    5527  LOAD_FAST            15  'nodeCoord'
#        5530  LOAD_FAST            17  'noid'
#        5533  BINARY_SUBSCR    
#        5534  LOAD_FAST            16  'nodeCoordPre'
#        5537  LOAD_FAST            17  'noid'
#        5540  STORE_SUBSCR     
#        5541  JUMP_BACK          5521  'to 5521'
#        5544  POP_BLOCK        
#      5545_0  COME_FROM                '5514'
#        5545  JUMP_FORWARD          0  'to 5548'
#      5548_0  COME_FROM                '5514'
#
#11618    5548  SETUP_LOOP           28  'to 5579'
#        5551  LOAD_FAST            22  'intfaceNodes'
#        5554  GET_ITER         
#        5555  FOR_ITER             20  'to 5578'
#        5558  STORE_FAST           17  'noid'
#
#11619    5561  LOAD_FAST            15  'nodeCoord'
#        5564  LOAD_FAST            17  'noid'
#        5567  BINARY_SUBSCR    
#        5568  LOAD_FAST            16  'nodeCoordPre'
#        5571  LOAD_FAST            17  'noid'
#        5574  STORE_SUBSCR     
#        5575  JUMP_BACK          5555  'to 5555'
#        5578  POP_BLOCK        
#      5579_0  COME_FROM                '5548'
#        5579  JUMP_BACK          4718  'to 4718'
#        5582  POP_BLOCK        
#      5583_0  COME_FROM                '4705'
#
#11621    5583  LOAD_GLOBAL           0  'stdout'
#        5586  LOAD_ATTR             1  'write'
#        5589  LOAD_CONST           19  '     -> Smoothing Iteration :       %4i  \n'
#        5592  LOAD_FAST            31  'nIt'
#        5595  BINARY_MODULO    
#        5596  CALL_FUNCTION_1       1 
#        5599  POP_TOP          
#
#11623    5600  LOAD_FAST            15  'nodeCoord'
#        5603  RETURN_VALUE     
#          -1  RETURN_LAST      
#
#Parse error at or near `POP_BLOCK' instruction at offset 1138


    def correctSmoothVoxelModel(self, nodeCoord, nodeCoordOrig, elsetNodes, niter):
        """
        Funktion corrects the smooth voxel model using the Jacobian
        
        @param nodeCoord: smoothed nodal coordinates of voxel model 
              - TYPE: dict[key] = (x,y,z)
              - int key     ... node id 
              - float x,y,z ... new nodal coordinates
        
        @param nodeCoordOrig: original (not smoothed) nodal coordinates  
              - TYPE: dict[key] = (x,y,z)
              - int key     ... node id
              - float x,y,z ... nodal coordinates
        
        @param elsetNodes: topological relations of mesh 
              - TYPE: dict[key] = (elid, noid1, ..., noid8)
              - int key  ... elset id
              - int elid ... element id 
              - int noid1-8 ... node ids building the element 
        
        @param niter: number of iterations          
              - int niter
        
        @return:
            nodeCoord: smoothed nodal coordinates of voxel model
              - TYPE: dict[key] = (x,y,z)
              - int key     ... node id 
              - float x,y,z ... new nodal coordinates
        """
        xl = numpy.zeros((3, 8), numpy.float)
        rstField = numpy.zeros((3, 8), numpy.float)
        rstField[(0, 0)] = -1.0
        rstField[(1, 0)] = -1.0
        rstField[(2, 0)] = -1.0
        rstField[(0, 1)] = +1.0
        rstField[(1, 1)] = -1.0
        rstField[(2, 1)] = -1.0
        rstField[(0, 2)] = +1.0
        rstField[(1, 2)] = +1.0
        rstField[(2, 2)] = -1.0
        rstField[(0, 3)] = -1.0
        rstField[(1, 3)] = +1.0
        rstField[(2, 3)] = -1.0
        rstField[(0, 4)] = -1.0
        rstField[(1, 4)] = -1.0
        rstField[(2, 4)] = +1.0
        rstField[(0, 5)] = +1.0
        rstField[(1, 5)] = -1.0
        rstField[(2, 5)] = +1.0
        rstField[(0, 6)] = +1.0
        rstField[(1, 6)] = +1.0
        rstField[(2, 6)] = +1.0
        rstField[(0, 7)] = -1.0
        rstField[(1, 7)] = +1.0
        rstField[(2, 7)] = +1.0
        stdout.write(' ... correct smooth model    \n')
        stdout.flush()
        for jiter in range(niter):
            count = 0
            countJacobi = 0
            ncount = len(elsetNodes)
            for elset in elsetNodes:
                count += 1
                stdout.write('     -> check Jacobian      : %10i %s  \r' % (int(count / float(ncount) * 100.0), '%'))
                if len(elsetNodes[elset]) > 0:
                    for elnds in elsetNodes[elset]:
                        xl[(0, 0)] = nodeCoord[elnds[1][0]][0]
                        xl[(1, 0)] = nodeCoord[elnds[1][0]][1]
                        xl[(2, 0)] = nodeCoord[elnds[1][0]][2]
                        xl[(0, 1)] = nodeCoord[elnds[1][1]][0]
                        xl[(1, 1)] = nodeCoord[elnds[1][1]][1]
                        xl[(2, 1)] = nodeCoord[elnds[1][1]][2]
                        xl[(0, 2)] = nodeCoord[elnds[1][2]][0]
                        xl[(1, 2)] = nodeCoord[elnds[1][2]][1]
                        xl[(2, 2)] = nodeCoord[elnds[1][2]][2]
                        xl[(0, 3)] = nodeCoord[elnds[1][3]][0]
                        xl[(1, 3)] = nodeCoord[elnds[1][3]][1]
                        xl[(2, 3)] = nodeCoord[elnds[1][3]][2]
                        xl[(0, 4)] = nodeCoord[elnds[1][4]][0]
                        xl[(1, 4)] = nodeCoord[elnds[1][4]][1]
                        xl[(2, 4)] = nodeCoord[elnds[1][4]][2]
                        xl[(0, 5)] = nodeCoord[elnds[1][5]][0]
                        xl[(1, 5)] = nodeCoord[elnds[1][5]][1]
                        xl[(2, 5)] = nodeCoord[elnds[1][5]][2]
                        xl[(0, 6)] = nodeCoord[elnds[1][6]][0]
                        xl[(1, 6)] = nodeCoord[elnds[1][6]][1]
                        xl[(2, 6)] = nodeCoord[elnds[1][6]][2]
                        xl[(0, 7)] = nodeCoord[elnds[1][7]][0]
                        xl[(1, 7)] = nodeCoord[elnds[1][7]][1]
                        xl[(2, 7)] = nodeCoord[elnds[1][7]][2]
                        noid = 0
                        for i in range(8):
                            rst = rstField[:, i]
                            noid += 1
                            det = micf77.jacobidet(rst[0], rst[1], rst[2], xl)
                            if det < 0.0:
                                countJacobi += 1
                                xyzOrig = nodeCoordOrig[elnds[1][noid - 1]]
                                xyzCurr = nodeCoord[elnds[1][noid - 1]]
                                nodeCoord[elnds[1][noid - 1]] = (
                                 0.5 * (xyzOrig[0] + xyzCurr[0]),
                                 0.5 * (xyzOrig[1] + xyzCurr[1]),
                                 0.5 * (xyzOrig[2] + xyzCurr[2]))

            stdout.write('     -> neg. Jacobi found %2i: %10i             \n' % (jiter + 1, countJacobi))

        return nodeCoord

    def writeInfo(self, voxelModel, additData, file=None, mode=None):
        stdout.write(' ... write model info                \n')
        nx, ny, nz = self.get_Shape(voxelModel)
        if additData['-ldim'] != None:
            lx, ly, lz = additData['-ldim']
        if additData.has_key('-thres'):
            thresList = additData['-thres']
        else:
            thresList = []
        filename = 'na'
        outfile = 'na'
        filenameShort = 'na'
        path = 'na'
        if additData.has_key('-in'):
            filename = additData['-in']
        if additData.has_key('-out'):
            if additData['-out'] != 'None':
                outfile = additData['-out']
                path, filenameShort = os.path.split(outfile)
                if path == '':
                    path = 'na'
        minVox, maxVox = self.computeNumpyMinMax(voxelModel, 0)
        stdout.write('     -> nx ny nz            : %i %i %i \n' % (nx, ny, nz))
        if additData['-ldim'] != None:
            stdout.write('     -> lx ly lz            : %g %g %g\n' % (lx, ly, lz))
        else:
            stdout.write('     -> lx ly lz            : %s\n' % 'not available')
        stdout.write('     -> minVal maxVal       : %g %g\n' % (minVox, maxVox))
        if 'ElementType' in additData:
            stdout.write('     -> inital data type    : %s \n' % additData['ElementType'])
        stdout.write('     -> final  data type    : %s \n' % voxelModel.dtype.name)
        if additData.has_key('-cogv'):
            xs = repr(additData['-cogv'][0])
            ys = repr(additData['-cogv'][1])
            zs = repr(additData['-cogv'][2])
        else:
            xs = 'na'
            ys = 'na'
            zs = 'na'
        offs = [
         0, 0, 0]
        if additData.has_key('Offset'):
            offs = additData['Offset']
        if file != None:
            voxSum = 0
            graySum = 0.0
            voxSum, graySum = micf77.computesum2(voxelModel, 2)
            BVTVlist = None
            if thresList:
                BVTVlist, sumVox = self.boneStatistics(voxelModel, thresList, echo=False)
            if mode == 'a':
                osfil = open(file, 'a')
            elif mode == 'w':
                osfil = open(file, 'w')
                osfil.write('$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;' % ('VOI_ElementDataFile',
                                                                                                                              'VOI_ElementDataFilePath',
                                                                                                                              'VOI_DimSize_x',
                                                                                                                              'VOI_DimSize_y',
                                                                                                                              'VOI_DimSize_z',
                                                                                                                              'VOI_ElementSpacing_x',
                                                                                                                              'VOI_ElementSpacing_y',
                                                                                                                              'VOI_ElementSpacing_z',
                                                                                                                              'VOI_Size_x',
                                                                                                                              'VOI_Size_y',
                                                                                                                              'VOI_Size_z',
                                                                                                                              'VOI_Offset_x',
                                                                                                                              'VOI_Offset_y',
                                                                                                                              'VOI_Offset_z',
                                                                                                                              'minVal',
                                                                                                                              'maxVal',
                                                                                                                              'voxSum',
                                                                                                                              'graySum',
                                                                                                                              'xsVox',
                                                                                                                              'ysVox',
                                                                                                                              'zsVox',
                                                                                                                              'Date',
                                                                                                                              'Time'))
                if BVTVlist:
                    for region in BVTVlist:
                        if int(region) > 0:
                            osfil.write('$%s;$%s;$%s;' % ('Threshold', 'BVTV100', 'BVTV'))

                osfil.write('\n')
            else:
                stdout.write("\n **ERROR** writeInfo(): Output mode '%s' not known!\n\n" % mode)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if additData['-ldim'] != None:
                osfil.write('%s;%s;%i;%i;%i;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%i;%g;%s;%s;%s;%s;%s;' % (filenameShort, path, nx, ny, nz, lx, ly, lz, nx * lx, ny * ly, nz * lz, offs[0], offs[1], offs[2], minVox, maxVox, voxSum, graySum, xs, ys, zs, time.strftime('%d %b %Y'), time.strftime('%H:%M:%S')))
            else:
                osfil.write('%s;%s;%i;%i;%i;%s;%s;%s;%s;%s;%s;%g;%g;%g;%g;%g;%i;%g;%s;%s;%s;%s;%s;' % (filenameShort, path, nx, ny, nz, 'na', 'na', 'na', 'na', 'na', 'na', offs[0], offs[1], offs[2], minVox, maxVox, voxSum, graySum, xs, ys, zs, time.strftime('%d %b %Y'), time.strftime('%H:%M:%S')))
            if BVTVlist:
                for region in BVTVlist:
                    if int(region) > 0:
                        BVTV = float(BVTVlist[region] / float(sumVox) * 100.0)
                        thres = region
                        osfil.write('%g;%g;%g;' % (thres, BVTV, BVTV / 100.0))

            osfil.write('\n')
            osfil.close()
        return

    def writeInfoOld(self, voxelModel, additData, file=None, mode=None):
        stdout.write(' ... write model info                \n')
        nx, ny, nz = self.get_Shape(voxelModel)
        if additData['-ldim'] != None:
            lx, ly, lz = additData['-ldim']
        if additData.has_key('-thres'):
            thresList = additData['-thres']
        else:
            thresList = []
        filename = ''
        outfile = ''
        if additData.has_key('-in'):
            filename = additData['-in']
        if additData.has_key('-out'):
            outfile = additData['-out']
        minVox, maxVox = self.computeNumpyMinMax(voxelModel, 0)
        stdout.write('     -> nx ny nz            : %i %i %i \n' % (nx, ny, nz))
        if additData['-ldim'] != None:
            stdout.write('     -> lx ly lz            : %g %g %g\n' % (lx, ly, lz))
        else:
            stdout.write('     -> lx ly lz            : %s\n' % 'not available')
        stdout.write('     -> minVal maxVal       : %g %g\n' % (minVox, maxVox))
        stdout.write('     -> inital data type    : %s \n' % additData['ElementType'])
        stdout.write('     -> final  data type    : %s \n' % voxelModel.dtype.name)
        if additData.has_key('-cogv'):
            xs = repr(additData['-cogv'][0])
            ys = repr(additData['-cogv'][1])
            zs = repr(additData['-cogv'][2])
        else:
            xs = 'na'
            ys = 'na'
            zs = 'na'
        if file != None:
            voxSum = 0
            graySum = 0.0
            voxSum, graySum = micf77.computesum2(voxelModel, 2)
            BVTVlist = None
            if thresList:
                BVTVlist, sumVox = self.boneStatistics(voxelModel, thresList, echo=False)
            if mode == 'a':
                osfil = open(file, 'a')
            elif mode == 'w':
                osfil = open(file, 'w')
                osfil.write('$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;$%s;' % ('inFile',
                                                                                              'outFile',
                                                                                              'nx',
                                                                                              'ny',
                                                                                              'nz',
                                                                                              'lx',
                                                                                              'ly',
                                                                                              'lz',
                                                                                              'minVal',
                                                                                              'maxVal',
                                                                                              'voxSum',
                                                                                              'graySum',
                                                                                              'xsVox',
                                                                                              'ysVox',
                                                                                              'zsVox'))
                if BVTVlist:
                    for region in BVTVlist:
                        if int(region) > 0:
                            osfil.write('$%s;$%s;$%s;' % ('Threshold', 'BVTV100', 'BVTV'))

                osfil.write('\n')
            else:
                stdout.write("\n **ERROR** writeInfo(): Output mode '%s' not known!\n\n" % mode)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if additData['-ldim'] != None:
                osfil.write('%s;%s;%i;%i;%i;%g;%g;%g;%g;%g;%i;%g;%s;%s;%s;' % (filename, outfile, nx, ny, nz, lx, ly, lz, minVox, maxVox, voxSum, graySum, xs, ys, zs))
            else:
                osfil.write('%s;%s;%i;%i;%i;%s;%s;%s;%g;%g;%i;%g;%s;%s;%s;' % (filename, outfile, nx, ny, nz, 'na', 'na', 'na', minVox, maxVox, voxSum, graySum, xs, ys, zs))
            if BVTVlist:
                for region in BVTVlist:
                    if int(region) > 0:
                        BVTV = float(BVTVlist[region] / float(sumVox) * 100.0)
                        thres = region
                        osfil.write('%g;%g;%g;' % (thres, BVTV, BVTV / 100.0))

            osfil.write('\n')
            osfil.close()
        return

    def writeBBoxInfo(self, bbox, thres, bbType, filename, file, mode):
        stdout.write(' ... write bounding box info         \n')
        stdout.flush()
        if mode == 'a':
            osfil = open(file, 'a')
        elif mode == 'w':
            osfil = open(file, 'w')
            osfil.write('#%s;%s;%s;%s;%s;%s;%s;%s;%s\n' % ('filename', 'threshold',
                                                           'bbType', 'minX', 'minY',
                                                           'minZ', 'maxX', 'maxY',
                                                           'maxZ'))
        else:
            stdout.write("\n **ERROR** writeInfo(): Output mode '%s' not known!\n\n" % mode)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        osfil.write('%s;%g;%s;%i;%i;%i;%i;%i;%i\n' % (filename, thres, bbType, bbox['minX'], bbox['minY'], bbox['minZ'], bbox['maxX'], bbox['maxY'], bbox['maxZ']))

    def averageRVE(self, voxelModel, ldims, filename, phyWinSize, threshold, maskValue):
        stdout.write(' ... compute RVE averges         \n')
        stdout.flush()
        OS = open(filename, 'w')
        if threshold.upper() == 'NONE':
            OS.write('nx0;ny0;nz0;dnx;dny;dnz;lx;ly;lz;ZV;TV;grayValueAverage\n')
        else:
            OS.write('nx0;ny0;nz0;dnx;dny;dnz;lx;ly;lz;ZV;TV;BVTV\n')
        nx, ny, nz = myModel.get_Shape(voxelModel)
        lx, ly, lz = ldims[0], ldims[1], ldims[2]
        nldims = (nx * lx, ny * ly, nz * lz)
        phyWinSize = float(phyWinSize)
        voxList = [[0], [0], [0]]
        for i in range(3):
            step = int(nldims[i] / phyWinSize)
            for j in range(step):
                nvox = int((j + 1) * phyWinSize / ldims[i] + 0.5)
                voxList[i].append(nvox)

        value = 0.0
        for k in range(len(voxList[2]) - 1):
            for j in range(len(voxList[1]) - 1):
                for i in range(len(voxList[0]) - 1):
                    dnx = voxList[0][i + 1] - voxList[0][i]
                    dny = voxList[1][j + 1] - voxList[1][j]
                    dnz = voxList[2][k + 1] - voxList[2][k]
                    nx0 = voxList[0][i]
                    ny0 = voxList[1][j]
                    nz0 = voxList[2][k]
                    subVoxelModelOrig = self.createVoxelModel(dnx, dny, dnz, 'f')
                    subVoxelModelOrig = micf77.cut(voxelModel, subVoxelModelOrig, nx0, ny0, nz0, 0)
                    thresMask = 0.0
                    ZV = 0
                    if maskValue.upper() != 'NONE':
                        thresMask = float(maskValue)
                        ZV = miaf77.compute_zv(subVoxelModelOrig, thresMask, 0)
                    TV = dnx * dny * dnz
                    if TV == ZV:
                        value = 0.0
                    elif threshold.upper() == 'NONE':
                        voxSum, grayValueSum = micf77.computesum2(subVoxelModelOrig, 0)
                        grayValueSum = grayValueSum - ZV * thresMask
                        grayValueAverage = grayValueSum / float(TV - ZV)
                        value = grayValueAverage
                    else:
                        subVoxelModelThres, BV = miaf77.compute_thres_bv(subVoxelModelOrig, float(threshold), 0)
                        rho = BV / float(TV - ZV)
                        value = rho
                    OS.write('%i;%i;%i;%i;%i;%i;%g;%g;%g;%i;%i;%g\n' % (nx0, ny0, nz0, dnx, dny, dnz, lx, ly, lz, ZV, TV, value))

        OS.close()

    def createFlatVoxelModel(self, dim, format):
        try:
            if format == 'B':
                voxelModel = numpy.zeros(dim, numpy.uint8)
            elif format == 'h':
                voxelModel = numpy.zeros(dim, numpy.int16)
            elif format == 'i':
                voxelModel = numpy.zeros(dim, numpy.int32)
            elif format == 'f':
                voxelModel = numpy.zeros(dim, numpy.float32)
            elif format == 'H':
                voxelModel = numpy.zeros(dim, numpy.uint16)
            else:
                stdout.write(' **ERROR** createFlatVoxelModel(): Binary Type "%s" not supported!\n\n' % format)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        except MemoryError:
            stdout.write(' **ERROR** createFlatVoxelModel: Cannot allocate enough memory!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        return voxelModel

    def createVoxelModel(self, dimX, dimY, dimZ, format):
        try:
            if format == 'B':
                voxelModel = numpy.zeros((dimZ, dimY, dimX), numpy.uint8)
            elif format == 'h':
                voxelModel = numpy.zeros((dimZ, dimY, dimX), numpy.int16)
            elif format == 'i':
                voxelModel = numpy.zeros((dimZ, dimY, dimX), numpy.int32)
            elif format == 'f':
                voxelModel = numpy.zeros((dimZ, dimY, dimX), numpy.float32)
            elif format == 'H':
                voxelModel = numpy.zeros((dimZ, dimY, dimX), numpy.uint16)
            else:
                stdout.write(' **ERROR** createVoxelModel(): Binary Type "%s" not supported!\n\n' % format)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        except MemoryError:
            stdout.write(' **ERROR** createVoxelModel: Cannot allocate enough memory!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

        return voxelModel

    def computeNumpyMinMax(self, curVoxelModel, echo=0):
        return (
         numpy.min(curVoxelModel), numpy.max(curVoxelModel))

    def castType(self, curVoxelModel, format):
        numpyVersion = float(numpy.__version__[0:3])
        minVox = 10000000
        maxVox = 10000000
        if format == 'B' or format == 'H' or format == 'h':
            maxVox = curVoxelModel.max()
            minVox = curVoxelModel.min()
        if format == 'B':
            if int(minVox) < 0 or int(maxVox) > 255:
                stdout.write('\n **ERROR** castType(). min=%s, max=%s, format=%s!\n' % (repr(minVox), repr(maxVox), format))
                stdout.flush()
                stdout.write(' *********** Use "-scale" option to scale your data from 0..255 first. \n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            elif numpyVersion > 1.6:
                curVoxelModel = curVoxelModel.astype(numpy.uint8, order='F')
            else:
                curVoxelModel = curVoxelModel.astype(numpy.uint8)
        elif format == 'H':
            if int(minVox) < 0 or int(maxVox) > 65535:
                stdout.write('\n **ERROR** castType(). min=%s, max=%s, format=%s!\n' % (repr(minVox), repr(maxVox), format))
                stdout.flush()
                stdout.write(' *********** Use "-scale" option to scale your data from 0..65535 first. \n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            elif numpyVersion > 1.6:
                curVoxelModel = curVoxelModel.astype(numpy.uint16, order='F')
            else:
                curVoxelModel = curVoxelModel.astype(numpy.uint16)
        elif format == 'h':
            if int(minVox) < -32768 or int(maxVox) > 32767:
                stdout.write('\n **ERROR** castType(). min=%s, max=%s, format=%s!\n' % (repr(minVox), repr(maxVox), format))
                stdout.flush()
                stdout.write(' *********** Use "-scale" option to scale your data from -32768..+32767 first. \n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            elif numpyVersion > 1.6:
                curVoxelModel = curVoxelModel.astype(numpy.int16, order='F')
            else:
                curVoxelModel = curVoxelModel.astype(numpy.int16)
        elif format == 'i':
            if numpyVersion > 1.6:
                curVoxelModel = curVoxelModel.astype(numpy.int32, order='F')
            else:
                curVoxelModel = curVoxelModel.astype(numpy.int32)
        elif format == 'f':
            if numpyVersion > 1.6:
                curVoxelModel = curVoxelModel.astype(numpy.float32, order='F')
            else:
                curVoxelModel = curVoxelModel.astype(numpy.float32)
        else:
            stdout.write('\n **ERROR** castType(). format=%s! not implemented\n' % format)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return curVoxelModel

    def checkMinMaxDataValue(self, ProcLog, minValue, maxValue):
        """
        Performs check in the case of *.aim L{readAimFile} file reading. Currently not used.
        """
        start = ProcLog.find('Minimum data value') + 19
        end = ProcLog.find('Maximum data value') - 1
        part = ProcLog[start:end]
        part = replace(part, '\n', '')
        part = replace(part, ' ', '')
        minValue2 = int(float(part))
        start = ProcLog.find('Maximum data value') + 19
        end = ProcLog.find('Average data value') - 1
        part = ProcLog[start:end]
        part = replace(part, '\n', '')
        part = replace(part, ' ', '')
        maxValue2 = int(float(part))
        if minValue2 != minValue:
            stdout.write('\n **ERROR** checkMinMaxDataValue(). Minimal data value read from *.aim is: %s\n' % minValue)
            stdout.flush()
            stdout.write('           Minimal data value from process log is: %s\n\n' % minValue2)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        if maxValue2 != maxValue:
            stdout.write('\n **ERROR** checkMinMaxDataValue(). Maximal data value read from *.aim is: %s\n' % maxValue)
            stdout.flush()
            stdout.write('           Maximal data value from process log is: %s\n\n' % maxValue2)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)

    def getDtypeMul(self, rawFormat):
        if rawFormat == 'B':
            dtype = 'uint8'
            mul = 1
        elif rawFormat == 'h':
            dtype = 'int16'
            mul = 2
        elif rawFormat == 'H':
            dtype = 'uint16'
            mul = 2
        elif rawFormat == 'i':
            dtype = 'int32'
            mul = 4
        elif rawFormat == 'f':
            dtype = 'float32'
            mul = 4
        else:
            stdout.write("\n **ERROR** mic.getDtypeMul: binary format '%s' not supported!\n\n" % rawFormat)
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return (dtype, mul)

    def getExtension(self, fileName):
        """ Function returns file extension. """
        parts = fileName.split('.')
        nparts = len(parts)
        ext = parts[nparts - 1]
        return ext

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

    def getShortFilenameAndExtension(self, fileName):
        filename, ext = self.getFilenameAndExtension(fileName)
        filename = filename.split('/')
        return (
         filename[len(filename) - 1], ext)

    def getPath(self, fileName):
        parts = fileName.split('/')
        nparts = len(parts)
        path = ''
        for part in range(nparts - 1):
            if len(parts[part]) > 0:
                path = path + '/' + parts[part]

        return path

    def updateModel(self, partIdElemDict, cleanVoxelModel, thresList, islandSize):
        """
        Auxiliary function for the cleaning function which updates interal 'part' information.
        Used in L{clean()}.
        """
        stdout.write('     -> update model                                   \n')
        stdout.flush()
        nx, ny, nz = self.get_Shape(cleanVoxelModel)
        maxElemBone = 0
        maxElemPartIdBone = 0
        for partId in partIdElemDict.keys():
            ck = partIdElemDict[partId][0][0]
            cj = partIdElemDict[partId][0][1]
            ci = partIdElemDict[partId][0][2]
            if cleanVoxelModel[ck, cj, ci] > 0:
                if len(partIdElemDict[partId]) > maxElemBone:
                    maxElemBone = len(partIdElemDict[partId])
                    maxElemPartIdBone = partId

        noBone = 0
        for partId in partIdElemDict.keys():
            ck = partIdElemDict[partId][0][0]
            cj = partIdElemDict[partId][0][1]
            ci = partIdElemDict[partId][0][2]
            if cleanVoxelModel[ck, cj, ci] > 0:
                noBone = noBone + 1

        noCleanElem = 0
        changedPartIds = []
        for partId in partIdElemDict.keys():
            ck = partIdElemDict[partId][0][0]
            cj = partIdElemDict[partId][0][1]
            ci = partIdElemDict[partId][0][2]
            if cleanVoxelModel[ck, cj, ci] > 0 and partId != maxElemPartIdBone:
                changedPartIds.append(partId)
                for elems in partIdElemDict[partId]:
                    cleanVoxelModel[elems[0], elems[1], elems[2]] = 0
                    noCleanElem = noCleanElem + 1

        for partId in partIdElemDict.keys():
            ck = partIdElemDict[partId][0][0]
            cj = partIdElemDict[partId][0][1]
            ci = partIdElemDict[partId][0][2]
            if cleanVoxelModel[ck, cj, ci] == 0 and len(partIdElemDict[partId]) <= islandSize and partId not in changedPartIds:
                for elems in partIdElemDict[partId]:
                    cleanVoxelModel[elems[0], elems[1], elems[2]] = thresList[0]
                    noCleanElem = noCleanElem - 1

        Vf = (maxElemBone + noCleanElem) / float(nz * ny * nx)
        cleanVf = maxElemBone / float(nz * ny * nx)
        stdout.write('     -> Volume Fraction         :  %5.2f %s\n' % (Vf * 100.0, '%'))
        stdout.flush()
        stdout.write('     -> Volume Fraction cleaned :  %5.2f %s\n' % (cleanVf * 100.0, '%'))
        stdout.flush()
        if noBone > 1:
            return False
        else:
            return True

    def get_Shape(self, voxelModel):
        """
        Returns the shape of a voxel model.
        """
        shapes = voxelModel.shape
        return (
         shapes[2], shapes[1], shapes[0])

    def get_Thresholded_Value(self, curVal):
        """
        Returns the discrete thresholded value for a given continous value between  .
        """
        threshold = 0
        if thresList == None:
            stdout.write('\n **ERROR** get_Threshold(): Threshold not optional when calling this function!\n\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        fullThres = thresList
        fullThres.append(255)
        for no in range(len(fullThres) - 1):
            if curVal > fullThres[no] and curVal <= fullThres[no + 1]:
                threshold = fullThres[no]

        fullThres.pop()
        return threshold

    def make_sphere(self, diameter, threshold):
        """ creates and returns a voxel model of a sphere """
        dn = diameter
        thresholdList = [threshold]
        thresVoxelModel = numpy.zeros((dn, dn, dn))
        for k in range(dn):
            for j in range(dn):
                for i in range(dn):
                    rc = ((i + 0.5 - dn / 2.0) ** 2.0 + (j + 0.5 - dn / 2.0) ** 2.0 + (k + 0.5 - dn / 2.0) ** 2.0) ** 0.5
                    if rc <= (dn - 1) / 2.0 + 1e-06:
                        thresVoxelModel[k, j, i] = threshold

        return thresVoxelModel

    def make_zylinder(self, diameter, height, threshold):
        """ creates and returns a voxel model of a cylinder """
        dn = diameter
        dh = height
        threshold = threshold
        thresholdList = [threshold]
        thresVoxelModel = numpy.zeros((dh, dn, dn))
        for k in range(dh):
            for j in range(dn):
                for i in range(dn):
                    rc = ((i + 0.5 - dn / 2.0) ** 2.0 + (j + 0.5 - dn / 2.0) ** 2.0) ** 0.5
                    if rc <= (dn - 1) / 2.0 + 1e-06:
                        thresVoxelModel[k, j, i] = threshold

        return thresVoxelModel

    def addTriad(self, voxelModel, position, axilen, value):
        dn = int(axilen)
        x = int(position[0])
        y = int(position[1])
        z = int(position[2])
        nx, ny, nz = self.get_Shape(voxelModel)
        for k in range(z, dn):
            if x >= 0 and x < nx and y >= 0 and y < ny and k >= 0 and k < nz:
                voxelModel[k, y, x] = value

        for j in range(y, dn):
            if x >= 0 and x < nx and j >= 0 and j < ny and z >= 0 and z < nz:
                voxelModel[z, j, x] = value

        for i in range(x, dn):
            if i >= 0 and i < nx and y >= 0 and y < ny and z >= 0 and z < nz:
                voxelModel[z, y, i] = value

        for k in range(dn):
            for j in range(dn):
                for i in range(dn):
                    rc = ((i + 0.5 - x) ** 2.0 + (j + 0.5 - y) ** 2.0 + (k + 0.5 - z) ** 2.0) ** 0.5
                    if rc <= dn / 1.6 and rc >= (dn - 1.5) / 1.6:
                        if i >= 0 and i < nx and j >= 0 and j < ny and k >= 0 and k < nz:
                            voxelModel[k, j, i] = value

        return voxelModel

    def writeMidplaneImageFem(self, voxelModel, outFileName, direct, scale=1.0):
        stdout.write(' ... write midplane fem images %s      \n' % direct)
        stdout.flush()
        time1 = time.clock()
        nx, ny, nz = self.get_Shape(voxelModel)
        filename, ext = self.getFilenameAndExtension(outFileName)
        fecModel = fec.fec()
        Nodes = {}
        Elems = {}
        EResults = {}
        nid = 0
        elid = 0
        corners = {}
        uniqueNodes = {}
        Elsets = {'slice': []}
        if direct == 'xm':
            name = filename + '-XM' + '.' + ext
            for k in range(nz + 1):
                for j in range(ny + 1):
                    nid = (ny + 1) * k + j + 1
                    curNode = dpFem.node(nid, nx / 2.0 * scale, j * scale, k * scale)
                    Nodes[nid] = curNode

            for k in range(nz):
                for j in range(ny):
                    elid += 1
                    nlist = [(ny + 1) * k + j + 1, (ny + 1) * k + (j + 1) + 1,
                     (ny + 1) * (k + 1) + (j + 1) + 1, (ny + 1) * (k + 1) + j + 1]
                    curElem = dpFem.element(elid, nlist, 'quad4')
                    Elems[elid] = curElem
                    EResults[elid] = voxelModel[k, j, nx / 2]
                    Elsets['slice'].append(elid)

        if direct == 'zm':
            name = filename + '-ZM' + '.' + ext
            for j in range(ny + 1):
                for i in range(nx + 1):
                    nid = (nx + 1) * j + i + 1
                    curNode = dpFem.node(nid, i * scale, j * scale, nz / 2.0 * scale)
                    Nodes[nid] = curNode

            for j in range(ny):
                for i in range(nx):
                    elid += 1
                    nlist = [(nx + 1) * j + i + 1, (nx + 1) * j + (i + 1) + 1,
                     (nx + 1) * (j + 1) + (i + 1) + 1, (nx + 1) * (j + 1) + i + 1]
                    curElem = dpFem.element(elid, nlist, 'quad4')
                    Elems[elid] = curElem
                    EResults[elid] = voxelModel[nz / 2, j, i]
                    Elsets['slice'].append(elid)

        if direct == 'ym':
            name = filename + '-YM' + '.' + ext
            for k in range(nz + 1):
                for i in range(nx + 1):
                    nid = (nx + 1) * k + i + 1
                    curNode = dpFem.node(nid, i * scale, ny / 2.0 * scale, k * scale)
                    Nodes[nid] = curNode

            for k in range(nz):
                for i in range(nx):
                    elid += 1
                    nlist = [(nx + 1) * k + i + 1, (nx + 1) * k + (i + 1) + 1,
                     (nx + 1) * (k + 1) + (i + 1) + 1, (nx + 1) * (k + 1) + i + 1]
                    curElem = dpFem.element(elid, nlist, 'quad4')
                    Elems[elid] = curElem
                    EResults[elid] = voxelModel[k, ny / 2, i]
                    Elsets['slice'].append(elid)

        fecModel = fec.fec()
        fecModel.write(name, None, Nodes, None, Elems, Elsets, EscaResults=[EResults])
        return

    def kernel3d(self, ktype):
        """
        6    ... 6 neighbors (cross) same as ktype=1.0 
        18   ... 18 neigbors same as ktype=1.732 (r<sqrt(3)) 
        1x1  ... square, 1x1 voxels 
        2x2  ... square, 2x2 voxels
        2.0  ... sphere, radus 2.0 voxels 
        3.5  ... sphere, radus 3.5 voxels     
        """
        if ktype.find('x') > -1:
            x, y = ktype.split('x')
            x = int(x)
            y = int(y)
            if abs(x) == abs(y):
                if x > 1:
                    if x % 2 == 1:
                        return numpy.ones((x, x, x))
                    dpUtils.throwError('mic.kernel3d(): Squared kernel must have an odd size i.e. 3x3, 5x5...!')
                else:
                    dpUtils.throwError('mic.kernel3d(): Squared kernel must be bigger than 1 i.e. 3x3, 5x5...!')
            else:
                string = 'mic.kernel3d(): Kernel not squared! ' + 'x=' + repr(x) + ' y=' + repr(y) + '!'
                dpUtils.throwError(string)
        elif ktype.find('.') > -1:
            r = float(ktype)
            size = int(2.0 * r / 2.0) * 2 + 1
            if size < 3:
                string = 'mic.kernel3d(): Kernel radius to small! ' + 'r=' + repr(r) + ' rmin = 1!'
                dpUtils.throwError(string)
            else:
                kernel = numpy.ones((size, size, size))
                xm = size / 2.0
                for k in range(size):
                    for j in range(size):
                        for i in range(size):
                            rc = numpy.sqrt((i + 0.5 - xm) * (i + 0.5 - xm) + (j + 0.5 - xm) * (j + 0.5 - xm) + (k + 0.5 - xm) * (k + 0.5 - xm))
                            if rc > r:
                                kernel[k, j, i] = 0

                return kernel
        else:
            if ktype.find('6') > -1 and int(ktype) == 6:
                kernel = numpy.zeros((3, 3, 3))
                kernel[(1, 1, 1)] = 1
                kernel[(1, 1, 0)] = 1
                kernel[(1, 1, 2)] = 1
                kernel[(0, 1, 1)] = 1
                kernel[(2, 1, 1)] = 1
                kernel[(1, 0, 1)] = 1
                kernel[(1, 2, 1)] = 1
                return kernel
            if ktype.find('18') > -1 and int(ktype) == 18:
                size = int(ktype)
                kernel = numpy.ones((3, 3, 3))
                kernel[(0, 0, 0)] = 0
                kernel[(2, 0, 0)] = 0
                kernel[(0, 2, 0)] = 0
                kernel[(0, 0, 2)] = 0
                kernel[(2, 2, 0)] = 0
                kernel[(0, 2, 2)] = 0
                kernel[(2, 0, 2)] = 0
                kernel[(2, 2, 2)] = 0
                return kernel
            string = "mic.kernel3d(): Kernel type='" + ktype + "' not implemented!"
            dpUtils.throwError(string)

    def getRawBinForm(self):
        rawBinForm = {}
        rawBinForm['uint8'] = 'B'
        rawBinForm['uint16'] = 'H'
        rawBinForm['int16'] = 'h'
        rawBinForm['int32'] = 'i'
        rawBinForm['float32'] = 'f'
        rawBinForm['B'] = 'B'
        rawBinForm['H'] = 'H'
        rawBinForm['h'] = 'h'
        rawBinForm['i'] = 'i'
        rawBinForm['f'] = 'f'
        rawBinForm[None] = None
        return rawBinForm

    def getRawBinForm2(self):
        rawBinForm = {}
        rawBinForm['B'] = 'uint8'
        rawBinForm['H'] = 'uint16'
        rawBinForm['h'] = 'int16'
        rawBinForm['i'] = 'int32'
        rawBinForm['f'] = 'float32'
        rawBinForm[None] = None
        return rawBinForm


if __name__ == '__main__':
    inName = None
    outName = None
    threshold = None
    fixt = None
    slevt = None
    dlevt = None
    lthres = None
    gthres = None
    ldim = None
    ndim = None
    cleanArg = None
    resolut = None
    resolut2 = None
    resolut3 = None
    resolutf = None
    resolutdf = None
    cut = None
    cap = None
    imr = None
    autot = None
    raw = None
    mask = None
    bool = None
    repl = None
    embed = None
    bbcut = None
    roicut = None
    scale = None
    arith = None
    rot = None
    rota = None
    rotm = None
    rotf = None
    form = None
    muscal = None
    fill = None
    cont = None
    imfil = None
    avg = None
    sobel = None
    grad = None
    lapl = None
    cfill = None
    mean = None
    median = None
    gauss = None
    morph = None
    morpk = None
    laham = None
    histo = None
    shist = None
    mid = None
    sing = None
    cover = None
    block = None
    geom = None
    mcut = None
    mirror = None
    mir2 = None
    flip = None
    temp = None
    voxSmooth = None
    cen = None
    mesh2d3 = None
    mesh2d4 = None
    ifo = None
    bbox = None
    extend = None
    close = None
    cog = None
    mfil = None
    guiInit = "*modulName                  'mic - medical image converter'                                                                                                                                    \n" + "*fileEntryIn      -in       'Image input file name        | filename'                                     slice.mhd                        no      mhd;nhdr;raw;png;jpg;gif;bmp;tif;aim;bin;AIM;bst;isq;ISQ    \n" + "*fileEntryOut     -out      'Image output file name       | filename'                                     test.mhd                         yes     mhd;nhdr;raw;png;jpg;gif;bmp;tif;inp;bin;fem;in;aim;gz;svx  \n" + "*combo            -form     'Raw image write format       | format'                                       uint8                            yes     B;H;h;i;f;uint8;uint16;int16;int32;float32  \n" + "*fileEntryOut     -mid      'Midplanes output file name   | filename'                                     test.png                         yes     png;jpg;gif;case                            \n" + "*subWindowStart             'Read parameters'                                                                                                                                                  \n" + "*splitEntry       -imr      '2D image read parameters     | startId;step;endId'                           1;1;60                           yes     int                                         \n" + "*splitEntry       -raw      'Raw image read parameters    | binFormat;headerSize;endianess;fastestDiret'  B;64;little;x                    yes     str                                         \n" + "*splitEntry       -ldim     'Voxel dimensions in 123      | len1;len2;len3'                      0.05;0.05;0.05                   yes     float                                                \n" + "*splitEntry       -ndim     'Numbers of voxel in 123      | nVox1;nVox2;nVox3'                            60;60;60                         yes     int                                         \n" + "*subWindowEnd               'Read Parameters'                                                                                                                                                  \n" + "*subWindowStart             'Write Parameters'                                                                                                                                                 \n" + "*fileEntryInEdit  -temp     'FE template file name        | filename'                                     main_temp.inp                    yes     inp;txt;mesh;in                             \n" + "*entry            -muscal   'AIM mu scaling               | scaleFactor'                                  8192                             yes     int                                         \n" + "*splitEntry       -smooth   'Mesh smoothing               | niter;lambda;kPB;interface;boundary;jacobian' 5;0.6307;0.1;0;1;0.05            yes     str                                         \n" + "*fileEntryOut     -mesh2d3  'Output 2D contour tria mesh  | filename'                                     test.case                        yes     inp;case;off;msh;geom;stl                   \n" + "*fileEntryOut     -mesh2d4  'Output 2D contour quad mesh  | filename'                                     test.case                        yes     inp;case;off;msh;geom                       \n" + "*splitEntry       -sing     'Output single 2d plane       | filename;direction;sliceId'                   test.png;1+;30                   yes     str                                         \n" + "*entry            -geom     'Output geometric object      | type;grayValue;voxelSize;filename;diameter;height'  cyl;255;0.05;test.mhd;5;10       yes     1                                     \n" + "*subWindowEnd               'Write Parameters'                                                                                                                                                 \n" + "*subWindowStart             'Compute Filters'                                                                                                                                                  \n" + "*entry            -flip     'Flip image                   | axis1;axis2;axis3'                            x;y;z                            yes     str                                         \n" + "*fileEntryIn      -mask     'Mask image                   | filename'                                     mask.mhd                         yes     mhd;nhdr;raw;png;jpg;gif;aim;bin;AIM;bst;isq;ISQ \n" + "*splitEntry       -bool     'Combine images*              | operator;filename'                            +;combine.mhd                    yes     str                                         \n" + "*entry            -repl     'Combine mesh+images*         | filename;elset;grayvalue;resolution'          mesh.inp;SET77;254               yes     str                                         \n" + "*entry            -arith    'Compute equation             | operation1;operation2;operation3;operation4'  ?<0=0;/8;*7;-5;<A;+3;<B;>A;^2;/B yes     str                                         \n" + "*entry            -scale    'Scale gray-value range       | newMin;newMax;[oldMin];[oldMax];[format]'     0;255                            yes     float                                       \n" + "*entry            -rota     'Rotate image by angles*      | axi1;axi2;axi3;angle1;angle2;angle3;interpolate' 3;2;1;1.6;-4.6;-39.8;YES      yes     str                                         \n" + "*entry            -rotm     'Rotate image by matrix*      | R11;R12;R13;R21;R22;R23;R31;R32;R33;interpolate' 0.7071;-0.7071;0;0.7071;0.7071;0;0;0;1;YES   yes     str                          \n" + "*entry            -rotf     'Rotate image by file*        | filename;interpolate'                         test.nhdr;YES                    yes     str                                         \n" + "*entry            -fill     'Fill pors                    | threshold;valid;type;kernel;[minThick];[niter];[kernel2];[bx0];[bx1];[by0];[by1];[bz0];[bz1]' 130;6;out;3   yes     str            \n" + "*entry            -cont     'Extract contour              | threshold'                                    75                               yes     float                                       \n" + "*splitEntry       -avg      'Average regions              | filename;phyWinSize;thres;maskValue'          test.txt;3.0;None;None           yes     str                                         \n" + "*subWindowEnd               'Compute Filters'                                                                                                                                                  \n" + "*subWindowStart             'Neighborhood Filters'                                                                                                                                             \n" + "*splitEntry       -laham    'Laplac Hamming filter        | weight;cutoff;amplitude'                      0.5;0.4;1.0                      yes     float                                       \n" + "*entry            -sobel    'Sobel filter'                                                                None                             yes     str                                         \n" + "*splitEntry       -mean     'Mean filter                  | kernel;thres1;thres2'                         k3x3;0;255                       yes     str                                         \n" + "*splitEntry       -median   'Median filter                | kernel;thres1;thres2'                         k3x3;0;255                       yes     str                                         \n" + "*splitEntry       -gauss    'Gauss filter                 | radius;sigma;thres1;thres2'                   1;0.8;0;75                       yes     str                                         \n" + "*splitEntry       -morph    'Morphological operations     | radius;type;shape;thres'                      3;o;1;140                        yes     str                                         \n" + "*splitEntry       -morpk    'Morphological operations 2   | kernel;type;thres'                            k6;o;140                         yes     str                                         \n" + "*entry            -grad     'Gradient filter'                                                             None                             yes     str                                         \n" + "*entry            -lapl     'Laplacian filter'                                                            None                             yes     str                                         \n" + "*entry            -cfill    'Cube fill filter             |  nlevels'                                     2                                yes     int                                         \n" + "*subWindowEnd               'Neighborhood Filters'                                                                                                                                             \n" + "*subWindowStart             'Resize Filters'                                                                                                                                                   \n" + "*splitEntry       -cut      'Crop image by voxel ROI*     | n0Vox1;n0Vox2;n0Vox3;dnVox1;dnVox2;dnVox3'    0;0;0;30;40;50                   yes     int                                         \n" + "*splitEntry       -bbcut    'Crop image by bounding box*  | threshold;extVox1;extVox2;extVox3'            75;0;0;0                         yes     int                                         \n" + "*splitEntry       -roicut   'Crop image by physical ROI*  | x0;y0;z0;x1;y1;z1'                            0.1;0.2;0.3;1.0;1.5;2.0          yes     float                                       \n" + "*splitEntry       -mcut     'Crop ROI from the middle*    | dnVox1;dnVox2;dnVox3'                         30;40;50                         yes     int                                         \n" + "*splitEntry       -res2     'Change image resolution      | res1;res2;res3'                               20;30;40                         yes     int                                         \n" + "*splitEntry       -res3     'Change resolution by length  | len1;len2;len3'                               0.15;0.17;0.18                   yes     float                                       \n" + "*splitEntry       -resf     'Change resolution by factor  | factor'                                       2                                yes     int                                         \n" + "*splitEntry       -refi     'Refine resolution by factor  | direction;factor'                             1;4                              yes     int                                         \n" + "*entry            -mirr     'Mirror image                 | axis1;axis2;axis3'                            x;y;z                            yes     str                                         \n" + "*splitEntry       -mir2     'Mirror image                 | nVox1;nVox2;nVox3'                            91;92;93                         yes     int                                         \n" + "*subWindowEnd               'Resize Filters'                                                                                                                                                   \n" + "*subWindowStart             'Segmentation Filters'                                                                                                                                             \n" + "*splitEntry       -autot    'Auto threshold with BVTV     | BVTV;error;estimate'                          0.2;0.01;100                     yes     str                                         \n" + "*splitEntry       -fixt     'Fixed threshold by value     | thresRatio;minimum;maximum'                   0.3;minVal;200                   yes     str                                         \n" + "*entry            -slevt    'Auto single level threshold  | threshold'                                    -32768                           yes     float                                       \n" + "*splitEntry       -dlevt    'Auto double level threshold  | thres1;thres2'                                40;80                            yes     float                                       \n" + "*entry            -lthres   'Local adaptive threshold     | type;alpha;LT;UT'                             stddev;1                         yes     str                                         \n" + "*combo            -gthres   'Global adaptove threshold    | type'                                         mean                             yes     mean;median                                 \n" + "*entry            -thres    'Simple threshold by value    | thres1;thres2;thres3'                         75                               yes     float                                       \n" + "*subWindowEnd               'Segmentation Filters'                                                                                                                                             \n" + "*subWindowStart             'Modify Filters'                                                                                                                                                   \n" + "*entry            -clean    'Clean image                  | type'                                         FAST                             yes     str                                         \n" + "*entry            -extend   'Extend image*                | direction;thickVox;[newGrayvalue]'            3;1                              yes     int                                         \n" + "*splitEntry       -close    'Close slice models           | direction;threshold;kernel;newGrayValue'      3;130;3;253                      yes     int                                         \n" + "*splitEntry       -embed    'Embed image*                 | direction;thickVoxIn;thickVoxOut;newGrayValue' 3;20;10;255                     yes     int                                         \n" + "*splitEntry       -cap      'Attach end caps*             | direction;thickVox;newGrayValue'              3;10;254                         yes     int                                         \n" + "*splitEntry       -cover    'Cover image*                 | thickVox;newGrayValue'                        5;253                            yes     int                                         \n" + "*splitEntry       -block    'Place image in block*        | nVox1;nVox2;nVox3;newGrayValue'               100;110;120;0                    yes     int                                         \n" + "*subWindowEnd               'Modify Filters'                                                                                                                                                   \n" + "*subWindowStart             'Multi Run'                                                                                                                                                        \n" + "*table            -mfil     'Multi Filter Run             | Option;Value(s)'                              mcut;;20;30;40;;;mid;;test1.png  yes     *;2                                         \n" + "*subWindowEnd               'Multi Run'                                                                                                                                                        \n" + "*subWindowStart             'Image Infos'                                                                                                                                                      \n" + "*entry            -cen      'Image center                 | threshold;filename;mode'                      75;test.txt;a                    yes     str                                         \n" + "*entry            -bbox     'Bounding box                 | type;threshold;filename;mode;i;j;k'           out;75;test.txt;a                yes     str                                         \n" + "*entry            -cog      'Center of gravity            | newGrayValue'                                 255                              yes     int                                         \n" + "*splitEntry       -ifo      'Output general info as CSV   | filename:mode'                                test.csv;w                       yes     str                                         \n" + "*splitEntry       -histo    'Output histogram             | filename;normalize;nIntervals'                test.dat;n;250                   yes     str                                         \n" + "*entry            -shist    'Show histogram material fit  | type;filename;mode'                           2;test.dat;a                     yes     int                                         \n" + "*subWindowEnd               'Image Infos'                                                                                                                                                      \n"
    argList = argv
    argc = len(argList)
    i = 0
    while i < argc:
        if argList[i][:3] == '-in':
            i += 1
            inName = argList[i]
        elif argList[i][:4] == '-out':
            i += 1
            outName = argList[i]
        elif argList[i][:6] == '-thres':
            i += 1
            threshold = argList[i]
        elif argList[i][:5] == '-fixt':
            i += 1
            fixt = argList[i]
        elif argList[i][:6] == '-slevt':
            i += 1
            slevt = argList[i]
        elif argList[i][:6] == '-dlevt':
            i += 1
            dlevt = argList[i]
        elif argList[i][:7] == '-lthres':
            i += 1
            lthres = argList[i]
        elif argList[i][:7] == '-gthres':
            i += 1
            gthres = argList[i]
        elif argList[i][:5] == '-ldim':
            i += 1
            ldim = argList[i]
        elif argList[i][:5] == '-ndim':
            i += 1
            ndim = argList[i]
        elif argList[i][:6] == '-clean':
            i += 1
            cleanArg = argList[i]
        elif argList[i][:5] == '-res1':
            i += 1
            resolut = argList[i]
        elif argList[i][:5] == '-res2':
            i += 1
            resolut2 = argList[i]
        elif argList[i][:5] == '-res3':
            i += 1
            resolut3 = argList[i]
        elif argList[i][:5] == '-resf':
            i += 1
            resolutf = argList[i]
        elif argList[i][:5] == '-refi':
            i += 1
            resolutdf = argList[i]
        elif argList[i][:4] == '-cut':
            i += 1
            cut = argList[i]
        elif argList[i][:5] == '-mcut':
            i += 1
            mcut = argList[i]
        elif argList[i][:4] == '-cap':
            i += 1
            cap = argList[i]
        elif argList[i][:4] == '-imr':
            i += 1
            imr = argList[i]
        elif argList[i][:6] == '-autot':
            i += 1
            autot = argList[i]
        elif argList[i][:4] == '-raw':
            i += 1
            raw = argList[i]
        elif argList[i][:6] == '-histo':
            i += 1
            histo = argList[i]
        elif argList[i][:6] == '-shist':
            i += 1
            shist = argList[i]
        elif argList[i][:5] == '-mask':
            i += 1
            mask = argList[i]
        elif argList[i][:5] == '-repl':
            i += 1
            repl = argList[i]
        elif argList[i][:5] == '-bool':
            i += 1
            bool = argList[i]
        elif argList[i][:6] == '-embed':
            i += 1
            embed = argList[i]
        elif argList[i][:6] == '-scale':
            i += 1
            scale = argList[i]
        elif argList[i][:6] == '-arith':
            i += 1
            arith = argList[i]
        elif argList[i][:4] == '-rot':
            if argList[i][:5] == '-rota':
                i += 1
                rota = argList[i]
            elif argList[i][:5] == '-rotm':
                i += 1
                rotm = argList[i]
            elif argList[i][:5] == '-rotf':
                i += 1
                rotf = argList[i]
            else:
                i += 1
                rot = argList[i]
        elif argList[i][:6] == '-bbcut':
            i += 1
            bbcut = argList[i]
        elif argList[i][:7] == '-roicut':
            i += 1
            roicut = argList[i]
        elif argList[i][:5] == '-fill':
            i += 1
            fill = argList[i]
        elif argList[i][:5] == '-cont':
            i += 1
            cont = argList[i]
        elif argList[i][:5] == '-form':
            i += 1
            form = argList[i]
        elif argList[i][:7] == '-muscal':
            i += 1
            muscal = argList[i]
        elif argList[i][:6] == '-imfil':
            i += 1
            imfil = argList[i]
        elif argList[i][:4] == '-avg':
            i += 1
            avg = argList[i]
        elif argList[i][:6] == '-sobel':
            i += 1
            sobel = argList[i]
        elif argList[i][:5] == '-mean':
            i += 1
            mean = argList[i]
        elif argList[i][:7] == '-median':
            i += 1
            median = argList[i]
        elif argList[i][:6] == '-gauss':
            i += 1
            gauss = argList[i]
        elif argList[i][:6] == '-morph':
            i += 1
            morph = argList[i]
        elif argList[i][:6] == '-morpk':
            i += 1
            morpk = argList[i]
        elif argList[i][:6] == '-laham':
            i += 1
            laham = argList[i]
        elif argList[i][:5] == '-grad':
            i += 1
            grad = argList[i]
        elif argList[i][:6] == '-lapl':
            i += 1
            lapl = argList[i]
        elif argList[i][:6] == '-cfill':
            i += 1
            cfill = argList[i]
        elif argList[i][:4] == '-cen':
            i += 1
            cen = argList[i]
        elif argList[i][:4] == '-mid':
            i += 1
            mid = argList[i]
        elif argList[i][:5] == '-sing':
            i += 1
            sing = argList[i]
        elif argList[i][:8] == '-mesh2d3':
            i += 1
            mesh2d3 = argList[i]
        elif argList[i][:8] == '-mesh2d4':
            i += 1
            mesh2d4 = argList[i]
        elif argList[i][:7] == '-extend':
            i += 1
            extend = argList[i]
        elif argList[i][:6] == '-close':
            i += 1
            close = argList[i]
        elif argList[i][:6] == '-cover':
            i += 1
            cover = argList[i]
        elif argList[i][:6] == '-block':
            i += 1
            block = argList[i]
        elif argList[i][:5] == '-geom':
            i += 1
            geom = argList[i]
        elif argList[i][:5] == '-mirr':
            i += 1
            mirror = argList[i]
        elif argList[i][:5] == '-mir2':
            i += 1
            mir2 = argList[i]
        elif argList[i][:5] == '-flip':
            i += 1
            flip = argList[i]
        elif argList[i][:5] == '-temp':
            i += 1
            temp = argList[i]
        elif argList[i][:7] == '-smooth':
            i += 1
            voxSmooth = argList[i]
        elif argList[i][:4] == '-ifo':
            i += 1
            ifo = argList[i]
        elif argList[i][:5] == '-bbox':
            i += 1
            bbox = argList[i]
        elif argList[i][:4] == '-cog':
            i += 1
            cog = argList[i]
        elif argList[i][:5] == '-mfil':
            i += 1
            mfil = argList[i]
        elif argList[i][:4] == '-gui':
            stdout.write('%s' % guiInit)
            stdout.flush()
            exit(0)
        elif argList[i][:5] == '-help':
            stdout.write('%s' % __doc__)
            stdout.flush()
            exit(0)
        i += 1

    if not inName:
        stdout.write(__doc__)
        stdout.flush()
        stdout.write('\n **ERROR** input file name not given\n\n')
        stdout.flush()
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    stdout.write('\n S T A R T  %s %s - %s (%ibit)\n\n' % (mic.__name__, mic.__version__, mic.__author__, osId))
    stdout.flush()
    startTime = time.clock()
    ctime1 = time.time()
    myModel = mic()
    ldimensionList = dpUtils.initValues(ldim, valName='-ldim', minVal=3, valType='float')
    ndimensionList = dpUtils.initValues(ndim, valName='-ndim', minVal=3, valType='int')
    thresholdList = dpUtils.initValues(threshold, valName='-thres', minVal=0, valType='int')
    imrList = dpUtils.initValues(imr, valName='-imr', minVal=3, valType='int')
    rawList = dpUtils.initValues(raw, valName='-raw', minVal=4)
    if voxSmooth != None:
        smoothList = dpUtils.initValues(voxSmooth, valName='-smooth', minVal=6, valType='float')
        smoothList.append(0.6)
    else:
        smoothList = None
    template = None
    if temp != None:
        template = temp
    rawBinForm = myModel.getRawBinForm()
    if raw == None:
        rawInfo = None
    else:
        rawInfo = {}
        rawInfo['binFormat'] = rawBinForm[rawList[0]]
        rawInfo['headerSize'] = int(rawList[1])
        if rawList[2] == 'big':
            rawInfo['bigEndian'] = True
        elif rawList[2] == 'little':
            rawInfo['bigEndian'] = False
        else:
            stdout.write('\n **ERROR** Option \'-raw\': Endian format "%s" unknown! \n\n' % rawList[2])
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        rawInfo['fastestDir'] = rawList[3]
    voxelModel, additData = myModel.read(inName, ndim=ndimensionList, imread=imrList, rawInfo=rawInfo)
    if not 'ElementType' in additData.keys():
        stdout.write('\n **INTERNAL ERROR** Data Element Type not given! \n\n')
        stdout.flush()
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    if '-thres' in additData.keys():
        if threshold == None:
            thresholdList = additData['-thres']
            thresVoxelModel = voxelModel
            threshold = 'setFromFile'
        else:
            stdout.write("\n **ERROR** threshold values are read from input file \n           and also given in '-thres' option! \n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    if '-ldim' in additData.keys():
        if ldimensionList == None:
            ldimensionList = additData['-ldim']
        else:
            stdout.write("\n **ERROR** voxel sizes are read from input file \n           and also given in '-ldim' option! \n\n")
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
    additData['-thres'] = thresholdList
    additData['-ldim'] = ldimensionList
    additData['ElementSpacing'] = ldim
    maxRunNo = 1
    if mfil != None:
        if mfil.find(':::') > -1:
            filterList = mfil.split(':::')
        elif mfil.find(';;;') > -1:
            filterList = mfil.split(';;;')
        else:
            filterList = [
             mfil]
        maxRunNo += len(filterList)
    runNo = 0
    while runNo < maxRunNo:
        runNo += 1
        if runNo > 1:
            stdout.write(' === FILTER RUN %i ===\n' % (runNo - 1))
            stdout.flush()
            mesh2d3 = None
            mesh2d4 = None
            mid = None
            sing = None
            outName = None
            flip = None
            mask = None
            bool = None
            repl = None
            arith = None
            scale = None
            rot = None
            rota = None
            rotm = None
            rotf = None
            fill = None
            cont = None
            imfil = None
            avg = None
            sobel = None
            mean = None
            median = None
            gauss = None
            morph = None
            morpk = None
            laham = None
            grad = None
            lapl = None
            cfill = None
            cut = None
            bbcut = None
            roicut = None
            mcut = None
            resolut = None
            resolut2 = None
            resolut3 = None
            resolutf = None
            resolutdf = None
            mirror = None
            mir2 = None
            autot = None
            fixt = None
            threshold = None
            cleanArg = None
            extend = None
            close = None
            embed = None
            cap = None
            cover = None
            block = None
            geom = None
            cen = None
            bbox = None
            cog = None
            ifo = None
            histo = None
            shist = None
            curFilter = filterList[runNo - 2]
            splitList = []
            if curFilter.find('::') > -1:
                splitList = curFilter.split('::')
            elif curFilter.find(';;') > -1:
                splitList = curFilter.split(';;')
            else:
                stdout.write("\n **ERROR** -mfil delimiter not implemented. String = '%s'\n\n" % curFilter)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if len(splitList) != 2:
                stdout.write('\n **ERROR** : -mfil Filter name and parameters has to be given!\n')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            filterName, parameters = splitList
            if filterName == 'mesh2d3':
                mesh2d3 = parameters
            elif filterName == 'mesh2d4':
                mesh2d4 = parameters
            elif filterName == 'out':
                outName = parameters
            elif filterName == 'mid':
                mid = parameters
            elif filterName == 'sing':
                sing = parameters
            elif filterName == 'flip':
                flip = parameters
            elif filterName == 'mask':
                mask = parameters
            elif filterName == 'bool':
                bool = parameters
            elif filterName == 'repl':
                repl = parameters
            elif filterName == 'arith':
                arith = parameters
            elif filterName == 'scale':
                scale = parameters
            elif filterName == 'rot':
                rot = parameters
            elif filterName == 'rota':
                rota = parameters
            elif filterName == 'rotm':
                rotm = parameters
            elif filterName == 'rotf':
                rotf = parameters
            elif filterName == 'fill':
                fill = parameters
            elif filterName == 'cont':
                cont = parameters
            elif filterName == 'imfil':
                imfil = parameters
            elif filterName == 'avg':
                avg = parameters
            elif filterName == 'sobel':
                sobel = parameters
            elif filterName == 'mean':
                mean = parameters
            elif filterName == 'median':
                median = parameters
            elif filterName == 'morph':
                morph = parameters
            elif filterName == 'morpk':
                morpk = parameters
            elif filterName == 'gauss':
                gauss = parameters
            elif filterName == 'laham':
                laham = parameters
            elif filterName == 'grad':
                grad = parameters
            elif filterName == 'lapl':
                lapl = parameters
            elif filterName == 'cfill':
                cfill = parameters
            elif filterName == 'cut':
                cut = parameters
            elif filterName == 'bbcut':
                bbcut = parameters
            elif filterName == 'roicut':
                roicut = parameters
            elif filterName == 'mcut':
                mcut = parameters
            elif filterName == 'res1':
                resolut = parameters
            elif filterName == 'res2':
                resolut2 = parameters
            elif filterName == 'res3':
                resolut3 = parameters
            elif filterName == 'resf':
                resolutf = parameters
            elif filterName == 'refi':
                resolutdf = parameters
            elif filterName == 'mirr':
                mirror = parameters
            elif filterName == 'mir2':
                mir2 = parameters
            elif filterName == 'autot':
                autot = parameters
            elif filterName == 'fixt':
                fixt = parameters
            elif filterName == 'thres':
                threshold = parameters
            elif filterName == 'clean':
                cleanArg = parameters
            elif filterName == 'extend':
                extend = parameters
            elif filterName == 'close':
                close = parameters
            elif filterName == 'embed':
                embed = parameters
            elif filterName == 'cap':
                cap = parameters
            elif filterName == 'cover':
                cover = parameters
            elif filterName == 'block':
                block = parameters
            elif filterName == 'geom':
                geom = parameters
            elif filterName == 'cen':
                cen = parameters
            elif filterName == 'bbox':
                bbox = parameters
            elif filterName == 'cog':
                cog = parameters
            elif filterName == 'ifo':
                ifo = parameters
            elif filterName == 'histo':
                histo = parameters
            elif filterName == 'shist':
                shist = parameters
            else:
                stdout.write("\n **ERROR** : -mfil filter '-%s' is not supported!\n" % filterName)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if threshold != None:
            thresholdList = dpUtils.initValues(threshold, valName='-thres', minVal=0, valType='int')
            additData['-thres'] = thresholdList
        if flip != None:
            flipList = dpUtils.initValues(flip, valName='-flip')
            for axis in flipList:
                if axis.upper() == 'X':
                    voxelModel = micf77.computeflip(voxelModel, 1, 1)
                elif axis.upper() == 'Y':
                    voxelModel = micf77.computeflip(voxelModel, 2, 1)
                elif axis.upper() == 'Z':
                    voxelModel = micf77.computeflip(voxelModel, 3, 1)
                else:
                    stdout.write('\n **ERROR** : -flip axis=%s not supported!\n' % axis)
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)

        if mask != None:
            maskModel, additDataM = myModel.read(mask, ndim=ndimensionList, imread=imrList, rawInfo=rawInfo)
            voxelModel = myModel.mask(voxelModel, maskModel)
        if bool != None:
            boolList = dpUtils.initValues(bool, valName='-bool')
            voxelModel = myModel.combineModels(voxelModel, additData, boolList[0], boolList[1])
        if repl != None:
            replList = dpUtils.initValues(repl, valName='-repl')
            if len(replList) > 3:
                if replList[2].upper() == 'NONE':
                    voxelModel = myModel.combineMeshModels(voxelModel, additData, replList[0], replList[1], grayValue=None, createNew=float(replList[3]))
                else:
                    voxelModel = myModel.combineMeshModels(voxelModel, additData, replList[0], replList[1], grayValue=int(replList[2]), createNew=float(replList[3]))
                ldimensionList = additData['-ldim']
            elif replList[2].upper() == 'NONE':
                voxelModel = myModel.combineMeshModels(voxelModel, additData, replList[0], replList[1], grayValue=None)
            else:
                voxelModel = myModel.combineMeshModels(voxelModel, additData, replList[0], replList[1], grayValue=int(replList[2]))
        if arith != None:
            arithList = dpUtils.initValues(arith, valName='-arith', minVal=0)
            voxelModel = myModel.computeEquation(voxelModel, arithList)
            if thresholdList != None:
                stdout.write('\n ** WARNING **: Function gives no/other segmented image! \n ')
                stdout.write('               Previous threshold is deleted! \n\n ')
                stdout.flush()
                thresholdList == None
                additData['-thres'] = None
        if scale != None:
            scaleList = dpUtils.initValues(scale, valName='-scale')
            newMinValue = scaleList[0]
            newMaxValue = scaleList[1]
            oldMinValue = None
            oldMaxValue = None
            format = None
            if len(scaleList) == 3:
                format = scaleList[2]
            if len(scaleList) == 4:
                oldMinValue = scaleList[2]
                oldMaxValue = scaleList[3]
            if len(scaleList) == 5:
                format = scaleList[4]
            voxelModel = myModel.scaleModel(voxelModel, newMinValue, newMaxValue, oldMinValue, oldMaxValue, format)
            if thresholdList != None:
                stdout.write('\n ** WARNING **: Function gives no segmented image! \n ')
                stdout.write('               Previous threshold is deleted! \n\n ')
                stdout.flush()
                thresholdList == None
                additData['-thres'] = None
        if rot != None:
            rotList = dpUtils.initValues(rot, valName='-rot')
            voxelModel = myModel.rotateModel(voxelModel, rotList, additData)
        if rota != None:
            rotList = dpUtils.initValues(rota, valName='-rota')
            voxelModel = myModel.rotateModel(voxelModel, ['angle'] + rotList, additData)
        if rotm != None:
            rotList = dpUtils.initValues(rotm, valName='-rotm')
            voxelModel = myModel.rotateModel(voxelModel, ['matrix'] + rotList, additData)
        if rotf != None:
            rotList = dpUtils.initValues(rotf, valName='-rotf')
            voxelModel = myModel.rotateModel(voxelModel, ['file'] + rotList, additData)
        if rot != None or rota != None or rotm != None or rotf != None:
            if thresholdList != None:
                stdout.write('\n ** WARNING **: Function gives no segmented image! \n ')
                stdout.write('               Previous threshold is deleted! \n\n ')
                stdout.flush()
                thresholdList == None
                additData['-thres'] = None
        if fill != None:
            fillList = dpUtils.initValues(fill, valName='-fill')
            fThreshold = float(fillList[0])
            valid = float(fillList[1])
            plane2d = 0
            if fillList[2].upper().find('_2D'):
                type_list = fillList[2].upper().split('_2D')
                fillList[2] = type_list[0]
                if len(type_list) > 1:
                    plane2d = int(type_list[1])
            if fillList[2].upper() == 'OUT':
                inOut = 1
            elif fillList[2].upper() == 'IN':
                inOut = 0
            elif fillList[2].upper() == 'T':
                inOut = 2
            elif fillList[2].upper() == 'V':
                inOut = 3
            elif fillList[2].upper() == 'C':
                inOut = 4
            else:
                stdout.write("\n **ERROR** Option '-fill': The third parameter ('%s') has to be IN, OUT, T, C or V! \n\n" % fillList[2])
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            skern = int(fillList[3])
            min_t = 1
            nmax = 0
            nx, ny, nz = myModel.get_Shape(voxelModel)
            fbbox = numpy.array([0, nx, 0, ny, 0, nz])
            if len(fillList) <= 4 and inOut == 4:
                stdout.write("\n **ERROR** Option '-fill': The fifth parameter (min cortex thickness) not set!\n\n")
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if len(fillList) == 5:
                min_t = int(fillList[4])
            if len(fillList) == 6:
                nmax = int(fillList[5])
            ekern = skern
            if len(fillList) == 7:
                nmax = int(fillList[5])
                ekern = int(fillList[6])
            if len(fillList) == 13:
                nmax = int(fillList[5])
                ekern = int(fillList[6])
                fbbox[0] = int(fillList[7])
                fbbox[1] = int(fillList[8])
                fbbox[2] = int(fillList[9])
                fbbox[3] = int(fillList[10])
                fbbox[4] = int(fillList[11])
                fbbox[5] = int(fillList[12])
            voxelModel = micf77.computefill(voxelModel, fThreshold, valid, inOut, skern, min_t, plane2d, nmax, ekern, fbbox)
            if thresholdList != None:
                stdout.write('\n ** WARNING **: Function gives no/other segmented image! \n ')
                stdout.write('               Previous threshold is deleted! \n\n ')
                stdout.flush()
                thresholdList == None
                additData['-thres'] = None
        if cont != None:
            voxelModel = micf77.findcontour(voxelModel, float(cont), 1)
        if imfil != None:
            imfilList = dpUtils.initValues(imfil, valName='-imfil')
            if imfilList[0].upper() == 'SOBEL':
                voxelModel = myModel.scaleModel(voxelModel, 0, 1)
                voxelModel = micf77.computesobel(voxelModel)
                voxelModel = myModel.scaleModel(voxelModel, 0, 255)
            elif imfilList[0].upper() == 'MEAN':
                radius = int(imfilList[1])
                mthres1 = float(imfilList[2])
                mthres2 = float(imfilList[3])
                voxelModel = micf77.computemean(voxelModel, radius, mthres1, mthres2)
            elif imfilList[0].upper() == 'GAUSS':
                radius = int(imfilList[1])
                sigma = float(imfilList[2])
                mthres1 = float(imfilList[3])
                mthres2 = float(imfilList[4])
                voxelModel = micf77.computegauss(voxelModel, radius, sigma, mthres1, mthres2)
            elif imfilList[0].upper() == 'MORPH':
                radius = int(imfilList[1])
                ftype = imfilList[2].upper()
                stype = int(imfilList[3])
                fThreshold = float(imfilList[4])
                if not (stype == 1 or stype == 2):
                    stdout.write("\n **ERROR** Option '-imfil': Kernel shape ('%s') not implemented! \n\n" % imfilList[3])
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                if ftype == 'C':
                    voxelModel = micf77.computedilation(voxelModel, fThreshold, radius, stype)
                    voxelModel = micf77.computeerosion(voxelModel, fThreshold, radius, stype)
                elif ftype == 'O':
                    voxelModel = micf77.computeerosion(voxelModel, fThreshold, radius, stype)
                    voxelModel = micf77.computedilation(voxelModel, fThreshold, radius, stype)
                elif ftype == 'E':
                    voxelModel = micf77.computeerosion(voxelModel, fThreshold, radius, stype)
                elif ftype == 'D':
                    voxelModel = micf77.computedilation(voxelModel, fThreshold, radius, stype)
                else:
                    stdout.write("\n **ERROR** Option '-imfil': Morh. Algorithm ('%s') not implemented! \n\n" % imfilList[2])
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                if thresholdList != None:
                    stdout.write('\n ** WARNING **: Function gives no/other segmented image! \n ')
                    stdout.write('               Previous threshold is deleted! \n\n ')
                    stdout.flush()
                    thresholdList == None
                    additData['-thres'] = None
                additData['-thres'] = [
                 fThreshold]
            else:
                stdout.write("\n **ERROR** Option '-imfil': Filter type ('%s') not implemented! \n\n" % imfilList[0])
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if avg != None:
            filename, phyWinSize, myThres, maskValue = dpUtils.userSplit(avg)
            ldims = (float(additData['-ldim'][0]), float(additData['-ldim'][1]), float(additData['-ldim'][2]))
            myModel.averageRVE(voxelModel, ldims, filename, phyWinSize, myThres, maskValue)
        if laham != None:
            lahamList = dpUtils.initValues(laham, valName='-laham', minVal=3, valType='float')
            voxelModel = myModel.computeLaplaceHamming(voxelModel, lahamList[0], lahamList[1], lahamList[2], ldimensionList, echo=True)
        if sobel != None:
            voxelModel = myModel.scaleModel(voxelModel, 0, 1)
            voxelModel = micf77.computesobel(voxelModel)
            voxelModel = myModel.scaleModel(voxelModel, 0, 255)
        if mean != None:
            imfilList = dpUtils.initValues(mean, valName='-mean')
            mthres1 = float(imfilList[1])
            mthres2 = float(imfilList[2])
            if imfilList[0].find('k') > -1:
                sid = imfilList[0].replace('k', '')
                kernel = myModel.kernel3d(sid)
                voxelModel = micf77.computemeankernel(voxelModel, kernel, mthres1, mthres2)
            else:
                radius = int(imfilList[0])
                voxelModel = micf77.computemean(voxelModel, radius, mthres1, mthres2)
        if median != None:
            imfilList = dpUtils.initValues(median, valName='-mean')
            sid = imfilList[0].replace('k', '')
            mthres1 = float(imfilList[1])
            mthres2 = float(imfilList[2])
            kernel = myModel.kernel3d(sid)
            voxelModel = micf77.computemediankernel(voxelModel, kernel, mthres1, mthres2)
        if gauss != None:
            imfilList = dpUtils.initValues(gauss, valName='-gauss')
            radius = int(imfilList[0])
            sigma = float(imfilList[1])
            mthres1 = float(imfilList[2])
            mthres2 = float(imfilList[3])
            voxelModel = micf77.computegauss(voxelModel, radius, sigma, mthres1, mthres2)
        if morph != None:
            imfilList = dpUtils.initValues(morph, valName='-morph')
            radius = int(imfilList[0])
            ftype = imfilList[1].upper()
            stype = int(imfilList[2])
            fThreshold = float(imfilList[3])
            if not (stype == 1 or stype == 2):
                stdout.write("\n **ERROR** Option '-morph': Kernel shape ('%s') not implemented! \n\n" % imfilList[3])
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if ftype == 'C':
                voxelModel = micf77.computedilation(voxelModel, fThreshold, radius, stype)
                voxelModel = micf77.computeerosion(voxelModel, fThreshold, radius, stype)
            elif ftype == 'O':
                voxelModel = micf77.computeerosion(voxelModel, fThreshold, radius, stype)
                voxelModel = micf77.computedilation(voxelModel, fThreshold, radius, stype)
            elif ftype == 'E':
                voxelModel = micf77.computeerosion(voxelModel, fThreshold, radius, stype)
            elif ftype == 'D':
                voxelModel = micf77.computedilation(voxelModel, fThreshold, radius, stype)
            else:
                stdout.write("\n **ERROR** Option '-morph': Algorithm ('%s') not implemented! \n\n" % imfilList[2])
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if thresholdList != None:
                stdout.write('\n ** WARNING **: Function gives no/other segmented image! \n ')
                stdout.write('               Previous threshold is deleted! \n\n ')
                stdout.flush()
                thresholdList == None
                additData['-thres'] = None
            additData['-thres'] = [
             fThreshold]
        if morpk != None:
            imfilList = dpUtils.initValues(morpk, valName='-morpk')
            sid = imfilList[0].replace('k', '')
            kernel = myModel.kernel3d(sid)
            ftype = imfilList[1].upper()
            fThreshold = float(imfilList[2])
            if ftype == 'C':
                voxelModel = micf77.computedilationkernel(voxelModel, kernel, fThreshold, 1)
                voxelModel = micf77.computeerosionkernel(voxelModel, kernel, fThreshold, 1)
            elif ftype == 'O':
                voxelModel = micf77.computeerosionkernel(voxelModel, kernel, fThreshold, 1)
                voxelModel = micf77.computedilationkernel(voxelModel, kernel, fThreshold, 1)
            elif ftype == 'E':
                voxelModel = micf77.computeerosionkernel(voxelModel, kernel, fThreshold, 1)
            elif ftype == 'D':
                voxelModel = micf77.computedilationkernel(voxelModel, kernel, fThreshold, 1)
            else:
                stdout.write("\n **ERROR** Option '-morpk': Algorithm ('%s') not implemented! \n\n" % imfilList[2])
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if thresholdList != None:
                stdout.write('\n ** WARNING **: Function gives no/other segmented image! \n ')
                stdout.write('               Previous threshold is deleted! \n\n ')
                stdout.flush()
                thresholdList == None
                additData['-thres'] = None
            additData['-thres'] = [
             fThreshold]
        if grad != None:
            voxelModel = micf77.computegradient(voxelModel)
        if lapl != None:
            voxelModel = micf77.computelaplacian(voxelModel)
        if cfill != None:
            voxelModel = myModel.computeCubeFill(voxelModel, int(cfill))
        if cut != None:
            cutList = dpUtils.initValues(cut, valName='-cut', minVal=6, valType='int')
            modVoxelModel = myModel.cutBlock(voxelModel, cutList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
        if bbcut != None:
            bbcutList = dpUtils.initValues(bbcut, valName='-bbcut', valType='int')
            extendList = None
            if len(bbcutList) > 1:
                extendList = [
                 bbcutList[1], bbcutList[2], bbcutList[3]]
            modVoxelModel = myModel.cutBlockBoundingBox(voxelModel, bbcutList[0], extendList=extendList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
        if roicut != None:
            roicutList = dpUtils.initValues(roicut, valName='-roicut', minVal=6, valType='float')
            modVoxelModel = myModel.cutBlockROI(voxelModel, roicutList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
        if mcut != None:
            mcutList = dpUtils.initValues(mcut, valName='-mcut', minVal=3, valType='int')
            modVoxelModel = myModel.mcutBlock(voxelModel, mcutList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
        if resolut != None:
            resList = dpUtils.initValues(resolut, valName='-res1', minVal=5, valType='float')
            stdout.write('\n ** WARNING ** Image enhancement is not supported anymore!\n')
            stdout.write("               Filter calls 'res2' option without enhancement.\n\n")
            stdout.flush()
            resList2 = resList[:3]
            modVoxelModel = myModel.changeResolution2(voxelModel, resList2, ldimensionList)
            del voxelModel
            voxelModel = modVoxelModel
        if resolut2 != None:
            resList2 = dpUtils.initValues(resolut2, valName='-res2', minVal=3, valType='int')
            modVoxelModel = myModel.changeResolution2(voxelModel, resList2, ldimensionList)
            del voxelModel
            voxelModel = modVoxelModel
        if resolut3 != None:
            resList3 = dpUtils.initValues(resolut3, valName='-res3', minVal=3, valType='float')
            modVoxelModel = myModel.changeResolution3(voxelModel, resList3, ldimensionList)
            del voxelModel
            voxelModel = modVoxelModel
        if resolutf != None:
            resf = int(resolutf)
            modVoxelModel = myModel.changeResolutionF(voxelModel, resf, ldimensionList)
            del voxelModel
            voxelModel = modVoxelModel
        if resolutdf != None:
            dir, factor = dpUtils.initValues(resolutdf, valName='-refi', minVal=2, valType='int')
            modVoxelModel = myModel.changeResolutionDF(voxelModel, dir, factor, ldimensionList)
            del voxelModel
            voxelModel = modVoxelModel
        if mirror != None:
            mirrorList = dpUtils.initValues(mirror, valName='-mirr')
            mirrVoxelModel = myModel.mirrorModel(voxelModel, mirrorList)
            del voxelModel
            voxelModel = mirrVoxelModel
        if mir2 != None:
            mir2List = dpUtils.initValues(mir2, valName='-mir2', minVal=3, valType='int')
            mirrVoxelModel = myModel.createVoxelModel(mir2List[0], mir2List[1], mir2List[2], 'f')
            voxelModel, mirrVoxelModel = micf77.computemirror2(voxelModel, mirrVoxelModel, [0, 0, 0], 1)
            del voxelModel
            voxelModel = mirrVoxelModel
        optionNo = 0
        if autot != None:
            optionNo += 1
        if fixt != None:
            optionNo += 1
        if slevt != None:
            optionNo += 1
        if lthres != None:
            optionNo += 1
        if gthres != None:
            optionNo += 1
        if dlevt != None:
            optionNo += 1
        if threshold != None and threshold != 'setFromFile':
            optionNo += 1
        if optionNo > 1:
            stdout.write("\n ** ERROR **: Only '-thres' or '-autot' or '-fixt' or '-slevt'  or '-dlevt' \n")
            stdout.write("                or '-lthres' or '-gthres' option can be used simultaneously! \n ")
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        else:
            boneStatistics = None
            if autot != None:
                autotList = dpUtils.initValues(autot, valName='-autot', minVal=3, valType='float')
                newThres = myModel.findThreshold(autotList, voxelModel)
                newThresList = [newThres]
                voxelModel = micf77.computethreshold(voxelModel, newThres, 2)
                thresholdList = newThresList
                additData['-thres'] = newThresList
            if fixt != None:
                fixtList = dpUtils.initValues(fixt, valName='-fixt', minVal=3)
                voxelModel, newThreshold = myModel.computeFixThreshold(voxelModel, fixtList[0], fixtList[1], fixtList[2], echo=True)
                additData['-thres'] = [newThreshold]
            if slevt != None:
                bgThresList = dpUtils.initValues(slevt, valName='-slevt', minVal=1, valType='float')
                newThres = micf77.computesinglelevelthreshold(voxelModel, bgThresList[0], 1)
                newThresList = [
                 newThres]
                voxelModel = micf77.computethreshold(voxelModel, newThres, 2)
                thresholdList = newThresList
                additData['-thres'] = newThresList
            if dlevt != None:
                dLevThresList = dpUtils.initValues(dlevt, valName='-dlevt', minVal=2, valType='float')
                newThres1, newThres2 = micf77.computedoublelevelthreshold(voxelModel, dLevThresList[0], dLevThresList[1], 1)
                newThresList = [newThres1, newThres2]
                voxelModel = myModel.threshold(voxelModel, newThresList)
                thresholdList = newThresList
                additData['-thres'] = newThresList
            if lthres != None:
                lthresThresList = dpUtils.initValues(lthres, valName='-lthres')
                type = lthresThresList[0]
                if type.lower() != 'stddev':
                    stdout.write("\n ** ERROR **: option -lthres : type '%s' not implemented  \n" % type)
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                alpha = float(lthresThresList[1])
                LT = 0.0
                UT = 0.0
                if len(lthresThresList) == 4:
                    LT = float(lthresThresList[2])
                    UT = float(lthresThresList[3])
                elif len(lthresThresList) == 3:
                    LT = float(lthresThresList[2])
                    UT = micf77.computesinglelevelthreshold(voxelModel, LT, 1)
                elif len(lthresThresList) == 2:
                    voxSum, graySum = micf77.computesum2(voxelModel, 1)
                    LT = graySum / float(voxSum)
                    UT = micf77.computesinglelevelthreshold(voxelModel, LT, 1)
                else:
                    stdout.write('\n ** ERROR **: option -lthres : 2 or 4 parameters allowed !\n')
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
                voxelModel = micf77.computelatstddev(voxelModel, LT, UT, alpha)
                additData['-thres'] = [1]
            if gthres != None:
                gthresThresList = dpUtils.initValues(gthres, valName='-gthres', minVal=1)
                type = gthresThresList[0]
                stdout.write(' ... compute glob adapt thres %s  \n' % type.upper())
                stdout.flush()
                if type.lower() == 'mean':
                    voxSum, graySum = micf77.computesum2(voxelModel, 0)
                    mean = graySum / float(voxSum)
                    stdout.write('     -> compute mean : %g  \n' % mean)
                    stdout.flush()
                    voxelModel = micf77.computethreshold0(voxelModel, mean, 2)
                    newThresList = [mean]
                    additData['-thres'] = newThresList
                elif type.lower() == 'median':
                    median = numpy.median(voxelModel)
                    stdout.write('     -> compute median : %g  \n' % median)
                    stdout.flush()
                    voxelModel = micf77.computethreshold0(voxelModel, median, 2)
                    newThresList = [median]
                    additData['-thres'] = newThresList
                else:
                    stdout.write("\n ** ERROR **: option -gthres : type '%s' not implemented  \n" % type)
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
            if threshold != None and threshold != 'setFromFile':
                voxelModel = myModel.threshold(voxelModel, thresholdList)
        if cleanArg != None:
            if cleanArg.upper() == 'FAST':
                stdout.write(' ... compute clean BONE \n')
                stdout.flush()
                maxMid = micf77.checkcomputeclean(voxelModel, 1)
                if maxMid < 2147483647:
                    voxelModel = micf77.computeclean(voxelModel, maxMid, 1)
                else:
                    stdout.write('\n ** ERROR **: option -clean : Input array to big  \n')
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
            elif cleanArg.upper() == 'NUMPY':
                stdout.write(' ... compute clean NUMPY \n')
                stdout.flush()
                voxelModel = myModel.clean_numpy(voxelModel)
            elif cleanArg.upper() == 'NUMPY1':
                stdout.write(' ... compute clean NUMPY1 \n')
                stdout.flush()
                voxelModel = myModel.clean_numpy(voxelModel, invert=True)
            elif cleanArg.upper() == 'FAST1' or cleanArg.upper() == 'FAST2':
                if additData['-thres'] != None:
                    if len(additData['-thres']) == 1:
                        if cleanArg.upper() == 'FAST2':
                            stdout.write(' ... compute clean BONE \n')
                            stdout.flush()
                            maxMid = micf77.checkcomputeclean(voxelModel, 1)
                            voxelModel = micf77.computeclean(voxelModel, maxMid, 1)
                        curthres = additData['-thres'][0]
                        arithList = ['*-1', '+' + repr(curthres)]
                        stdout.write(' ... compute clean AIR\n')
                        stdout.flush()
                        voxelModel = myModel.computeEquation(voxelModel, arithList, 0)
                        maxMid = micf77.checkcomputeclean(voxelModel, 1)
                        voxelModel = micf77.computeclean(voxelModel, maxMid, 1)
                        voxelModel = myModel.computeEquation(voxelModel, arithList, 0)
                    else:
                        stdout.write('\n ** ERROR **: Fast inverse cleaning not performed. Exactly one threshold needed! \n\n ')
                        stdout.flush()
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                else:
                    stdout.write('\n ** ERROR **: Fast inverse cleaning not performed. No segmented image provided! \n\n ')
                    stdout.flush()
                    stdout.write('\n E N D E D  with ERRORS \n\n')
                    stdout.flush()
                    exit(1)
            elif len(cleanArg) == 2:
                if additData['-thres'] != None:
                    cleanList = dpUtils.initValues(cleanArg, valName='-clean', minVal=2, valType='int')
                    modVoxelModel = myModel.clean(voxelModel, thresholdList, cleanList)
                    del voxelModel
                    voxelModel = modVoxelModel
                else:
                    stdout.write('\n ** WARNING **: Cleaning not performed. Threshold needed! \n\n ')
                    stdout.flush()
            else:
                stdout.write('\n ** ERROR **: Clean options: 1 value = FAST / FAST1 / FAST2 or 2 values required! \n\n ')
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if embed != None:
            embedList = dpUtils.initValues(embed, valName='-embed', minVal=4)
            modVoxelModel = myModel.embedModel(voxelModel, embedList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
            if threshold != None or autot != None:
                if threshold != 'setFromFile':
                    thresholdList.append(int(embedList[3]))
        if cap != None:
            capList = dpUtils.initValues(cap, valName='-cap', minVal=3)
            modVoxelModel = myModel.makeEndCaps(voxelModel, capList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
            if threshold != None or autot != None:
                if threshold != 'setFromFile':
                    thresholdList.append(int(capList[2]))
        if extend != None:
            extendList = dpUtils.initValues(extend, valName='-extend', valType='int')
            modVoxelModel = myModel.extendModel(voxelModel, extendList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
        if close != None:
            closeList = dpUtils.initValues(close, valName='-close', minVal=4, valType='int')
            voxelModel = myModel.closeModel(voxelModel, closeList)
        if cover != None:
            coverList = dpUtils.initValues(cover, valName='-cover', minVal=2, valType='int')
            modVoxelModel = myModel.coverModel(voxelModel, coverList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
            if threshold != None or autot != None:
                if threshold != 'setFromFile':
                    thresholdList.append(int(coverList[1]))
        if block != None:
            blockList = dpUtils.initValues(block, valName='-block', minVal=4, valType='int')
            modVoxelModel = myModel.placeInBlockModel(voxelModel, blockList, additData=additData)
            del voxelModel
            voxelModel = modVoxelModel
            if threshold != None or autot != None:
                if threshold != 'setFromFile':
                    thresholdList.append(int(blockList[3]))
        if histo != None:
            histoList = dpUtils.initValues(histo, valName='-histo')
            if len(histoList) < 2:
                stdout.write('\n ** ERROR **: -histo option needs at least 2 arguments! \n ')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            elif len(histoList) == 2:
                histoList.append(None)
            elif len(histoList) > 3:
                stdout.write('\n ** ERROR **: -histo option needs maximum 3 arguments! \n ')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            myModel.writeHistogram(histoList[0], voxelModel, normalize=histoList[1], noInterval=histoList[2])
        if shist != None:
            shistList = dpUtils.initValues(shist, valName='-shist')
            if len(shistList) == 1:
                myModel.showHistogramFit(voxelModel, int(shist))
            elif len(shistList) == 3:
                path, inName2 = os.path.split(inName)
                myModel.showHistogramFit(voxelModel, int(shistList[0]), infile=inName2, outfile=shistList[1], mode=shistList[2])
            else:
                stdout.write('\n ** ERROR **: -shist option needs 1 or 3 arguments! \n ')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if cen != None:
            cenList = dpUtils.initValues(cen, valName='-cen')
            if len(cenList) == 1:
                cenRes = myModel.computeCgalCenter(voxelModel, float(cenList[0]), ldimensionList)
            elif len(cenList) == 3:
                cenRes = myModel.computeCgalCenter(voxelModel, float(cenList[0]), ldimensionList, file=cenList[1], mode=cenList[2])
            else:
                stdout.write('\n ** ERROR **: -cen option needs 1 or 3 arguments! \n ')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if bbox != None:
            bboxList = dpUtils.initValues(bbox, valName='-bbox')
            bbType = bboxList[0]
            thres = float(bboxList[1])
            bboxComp = {}
            if bbType == 'out':
                bboxComp = myModel.boundingBox(voxelModel, thres)
            elif bbType == 'in':
                nx, ny, nz = myModel.get_Shape(voxelModel)
                xs = int(nx / 2)
                ys = int(ny / 2)
                zs = int(nz / 2)
                if len(bboxList) == 7:
                    xs = int(bboxList[4])
                    ys = int(bboxList[5])
                    zs = int(bboxList[6])
                bboxAux = micf77.computeinnerbb(voxelModel, thres, xs, ys, zs, 1)
                bboxComp['minX'] = bboxAux[0, 0]
                bboxComp['minY'] = bboxAux[1, 0]
                bboxComp['minZ'] = bboxAux[2, 0]
                bboxComp['maxX'] = bboxAux[0, 1]
                bboxComp['maxY'] = bboxAux[1, 1]
                bboxComp['maxZ'] = bboxAux[2, 1]
            else:
                stdout.write("\n ** ERROR **: -bbox type '%s' no implemented! \n" % bbType)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            if len(bboxList) >= 4:
                file = bboxList[2]
                mode = bboxList[3]
                parts = inName.split('/')
                nParts = len(parts)
                myModel.writeBBoxInfo(bboxComp, thres, bbType, parts[nParts - 1], file, mode)
        if cog != None:
            xcog, ycog, zcog = micf77.computecog(voxelModel, 1)
            stdout.write('     -> xs, ys, zs (Voxel)  : %g %g %g \n' % (xcog, ycog, zcog))
            if additData.has_key('-ldim'):
                lx = float(additData['-ldim'][0])
                ly = float(additData['-ldim'][1])
                lz = float(additData['-ldim'][2])
                stdout.write('     -> xs, ys, zs (Real)   : %g %g %g \n' % (xcog * lx, ycog * ly, zcog * lz))
            additData['-cogv'] = [
             xcog, ycog, zcog]
            if repr(cog).upper() != 'OFF':
                grayValue = float(cog)
                voxelModel = micf77.addplane(voxelModel, 1, int(xcog + 0.5), grayValue, 1)
                voxelModel = micf77.addplane(voxelModel, 2, int(ycog + 0.5), grayValue, 1)
                voxelModel = micf77.addplane(voxelModel, 3, int(zcog + 0.5), grayValue, 1)
                nx, ny, nz = myModel.get_Shape(voxelModel)
                axilen = int(0.2 * (nx + ny + nz) / 3.0)
                position = [0, 0, 0]
                voxelModel = myModel.addTriad(voxelModel, position, axilen, grayValue)
        if ifo != None:
            parts = inName.split('/')
            nParts = len(parts)
            additData['-in'] = parts[nParts - 1]
            additData['-out'] = 'None'
            if outName:
                parts = outName.split('/')
                nParts = len(parts)
                additData['-out'] = parts[nParts - 1]
            ifoList = dpUtils.initValues(ifo, valName='-ifo', minVal=2)
            myModel.writeInfo(voxelModel, additData, file=ifoList[0], mode=ifoList[1])
        else:
            myModel.writeInfo(voxelModel, additData)
        i = 1
        dt = time.localtime()
        history = {}
        history['Session'] = str(dt[0]) + '/' + str(dt[1]) + '/' + str(dt[2]) + '-' + str(dt[3]) + ':' + str(dt[4])
        while i < argc:
            curStr = argv[i + 1]
            curStr = curStr.replace('<', '.LT.')
            curStr = curStr.replace('>', '.GT.')
            history[argv[i]] = curStr
            i = i + 2

        if outName != None:
            geomList = None
            if additData.has_key('TransformMatrix') and additData.has_key('Offset') and additData.has_key('CenterOfRotation') and additData.has_key('AnatomicalOrientation'):
                geomList = {}
                geomList['TransformMatrix'] = additData['TransformMatrix']
                geomList['Offset'] = additData['Offset']
                geomList['CenterOfRotation'] = additData['CenterOfRotation']
                geomList['AnatomicalOrientation'] = additData['AnatomicalOrientation']
            if outName.upper() != 'NO_OUTPUT':
                myModel.write(outName, voxelModel, thresList=thresholdList, dimList=ldimensionList, format=rawBinForm[form], history=history, template=template, smooth=smoothList, geomList=geomList, muscal=muscal)
        if mesh2d3 != None:
            filename, ext = os.path.splitext(mesh2d3)
            if ext.lower() == '.stl':
                myModel.writeMeshStlBin(mesh2d3, voxelModel, ldimensionList, smoothList, additData=additData)
            else:
                myModel.writeMesh2d(mesh2d3, voxelModel, ldimensionList, smoothList, 'tria3')
        if mesh2d4 != None:
            myModel.writeMesh2d(mesh2d4, voxelModel, ldimensionList, smoothList, 'quad4')
        if mid != None:
            midList = dpUtils.initValues(mid, valName='-mid', minVal=1)
            filename, ext = myModel.getFilenameAndExtension(midList[0])
            uext = ext.upper()
            if uext == 'CASE':
                scaleX = 1.0
                scaleY = 1.0
                scaleZ = 1.0
                if ldimensionList != None:
                    scaleX = ldimensionList[0]
                    scaleY = ldimensionList[1]
                    scaleZ = ldimensionList[2]
                myModel.writeMidplaneImageFem(voxelModel, filename + '.case', 'xm', scale=scaleX)
                myModel.writeMidplaneImageFem(voxelModel, filename + '.case', 'ym', scale=scaleY)
                myModel.writeMidplaneImageFem(voxelModel, filename + '.case', 'zm', scale=scaleZ)
            elif uext == 'JPG' or uext == 'PNG' or uext == 'GIF' or uext == 'TIF' or uext == 'BMP':
                myModel.writeMidplaneImages(filename + '.' + ext, voxelModel)
            else:
                stdout.write('\n ** ERROR **: -mid option: Filetype %1 not supported! \n ' % ext)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        if sing != None:
            singList = dpUtils.initValues(sing, valName='-sing', minVal=3)
            filename = singList[0]
            direct = singList[1]
            layer = int(singList[2])
            myModel.writeSingleImages(voxelModel, direct, layer, filename)
        if geom != None:
            geomModelList = dpUtils.initValues(geom, valName='-geom')
            myModel.writeGeomObject(geomModelList)

    stdout.flush()
    endTime = time.clock()
    ctime2 = time.time()
    stdout.write('\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime - startTime, ctime2 - ctime1))
    stdout.flush()

