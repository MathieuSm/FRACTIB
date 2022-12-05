$! DOWNSCALE ISQ
$! Mathieu Simon
$! December 2022
$!
$! Define Parameters
$! ---------------------------------------------
$! P1 = sampleID
$! P2 = meas_nb 
$! P3 = file name
$! P4 = scale_factor
$!
$! ---------------------------------------------
$  path           := DK0:[MICROCT.DATA.'P1'.'P2']
$  isq             := 'path''P3'.ISQ
$  aim            := 'path''P3'.AIM
$  gobj           := 'path''P3'.GOBJ
$  downscale     := 'path''P3'_DOWNSCALED.AIM
$  gobjdown     := 'path''P3'_GOBJ.AIM
$  scale_factor   := 'P4'
$!
$! ---------------------------------------------
$  show symbol path
$  show symbol name_ISQ
$  show symbol name_AIM
$  show symbol name_GOBJ
$  show symbol scale_factor
$  show symbol name_down_scale
$!
$ debug	:= true
$! Debug = true -> safety files are written out
$!
$!Start IPL:
$ ipl_batch

/db_scanco_activate true

!! ------------------------------------------------------------------------------
!! Read ISQ and convert to AIM
!! ------------------------------------------------------------------------------
/isq_to_aim
  -aim_name                  aim
  -isq_filename              "isq
  -pos                       0 0 0
  -dim                       -1 -1 -1

!! ------------------------------------------------------------------------------
!! Mask ISQ with GOBJ (set values outside GOBJ to 0
!! ------------------------------------------------------------------------------

/gobj_maskaimpeel_ow
 -input_output              aim
 -gobj_filename             "gobj
 -peel_iter                 0

!! ------------------------------------------------------------------------------
!! Crop AIM to limits imposed by GOBJ
!! ------------------------------------------------------------------------------

/bounding_box_cut
  -input                     aim
  -output                    crop
  -z_only                    false
  -border                    0 0 0

!! ------------------------------------------------------------------------------
!! Downscale AIM
!! ------------------------------------------------------------------------------

/noip_scale_elsize
  -input                     aim
  -output                    rescaled
  -down_scale                "scale_factor
  -up_scale                  1.000
  -average                   true

!! ------------------------------------------------------------------------------
!! Write to disk
!! ------------------------------------------------------------------------------

/write_v020
  -name                      rescaled
  -filename                  "name_DOWNSCALE
  -compress_type             none
  -version_020               true

!! ------------------------------------------------------------------------------
!! Read GOBJ and convert to aim
!! ------------------------------------------------------------------------------
/gobj_to_aim
  -gobj_filename            "name_GOBJ
  -output                   aim
  -peel_iter                0

!! ------------------------------------------------------------------------------
!! Rescale GOBJ
!! ------------------------------------------------------------------------------

/noip_scale_elsize
  -input                     aim
  -output                    rescaled
  -down_scale                "scale_factor
  -up_scale                  1.000
  -average                   false

!! ------------------------------------------------------------------------------
!! Write to disk
!! ------------------------------------------------------------------------------

/write_v020
  -name                      rescaled
  -filename                  "name_GOBJDOWN
  -compress_type             none
  -version_020               true

 
..
$
$ exit
