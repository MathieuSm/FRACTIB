$! IPL script to generate Input images for hFE
$! Mereo Project
$! Authot: Michael Indermaur
$! Date: 15.07.2020
$! Bern, Switzerland
$!
$!------------------------------------------------------------------------------------------------
$!
$! Submission of the script is
$! @IPL_common_image_size  432_L_77_F  C0003103 
$!
$! Definition of files:
$  Folder              := 'P1' ! folder name
$  meas_ID             := 'P2' ! file name (e.g.CXXXXXXX)
$!
$! Define Parameters
$! -----------------
$ path1             :=DK0:[MICROCT.DATA.MATHIEU.FRACTIB.'P1']
$ aim               :='path1''P2'.AIM
$ aim_final         :='path1''P2'_UNCOMP.AIM
$ seg               :='path1''P2'_SEG.AIM
$ seg_final         :='path1''P2'_SEG_UNCOMP.AIM
$ trab_mask         :='path1''P2'_TRAB_MASK.AIM
$ cort_mask         :='path1''P2'_CORT_MASK.AIM
$ trab_mask_final   :='path1''P2'_TRAB_MASK_UNCOMP.AIM
$ cort_mask_final   :='path1''P2'_CORT_MASK_UNCOMP.AIM
$!
$!
$!
$ ipl_batch

/db_scanco_activate true

!! ------------------------------------------------------------------------------
!! Create whitebox from SEG.AIM
!! ------------------------------------------------------------------------------

/read
  -name                      seg
  -filename                  "seg
  
/copy
  -in                        seg
  -out                       whitebox
  
/set_value
  -input                     whitebox
  -value_object              1
  -value_background          1

/list true false

!! ------------------------------------------------------------------------------
!! Create SEG FINAL.AIM
!! ------------------------------------------------------------------------------

/read
  -name                      seg
  -filename                  "seg

/multiply_volumes
  -input1                    seg
  -input2                    whitebox
  -output                    seg_final
  -common_region_only        true
  -multiply_zero_pixels      true
  
/write_v020
  -name                      seg_final
  -filename                  "seg_final
  -compress_type             none
  -version_020               true

/list true false

!! ------------------------------------------------------------------------------
!! Create BMD FINAL.AIM
!! ------------------------------------------------------------------------------

/read
  -name                      aim
  -filename                  "aim

/multiply_volumes
  -input1                    aim
  -input2                    whitebox
  -output                    aim_final
  -common_region_only        true
  -multiply_zero_pixels      true
  
/write_v020
  -name                      aim_final
  -filename                  "aim_final
  -compress_type             none
  -version_020               true

/list true false

!! ------------------------------------------------------------------------------
!! Create CORT_MASK_FINAL.AIM
!! ------------------------------------------------------------------------------  

/read
  -name                      cort_mask
  -filename                  "cort_mask

/multiply_volumes
  -input1                    cort_mask
  -input2                    whitebox
  -output                    cort_mask_final
  -common_region_only        true
  -multiply_zero_pixels      true

/write_v020
  -name                      cort_mask_final
  -filename                  "cort_mask_final
  -compress_type             none
  -version_020               true

/list true false

!! ------------------------------------------------------------------------------
!! Create TRAB_MASK_FINAL.AIM
!! ------------------------------------------------------------------------------  

/read
  -name                      trab_mask
  -filename                  "trab_mask

/multiply_volumes
  -input1                    trab_mask
  -input2                    whitebox
  -output                    trab_mask_final
  -common_region_only        true
  -multiply_zero_pixels      true

/write_v020
  -name                      trab_mask_final
  -filename                  "trab_mask_final
  -compress_type             none
  -version_020               true
  
/list true false
..
$
$ exit
