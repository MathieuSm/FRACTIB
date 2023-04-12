$! *************************************************************
$! Autor:     Mathieu Simon, from original of Denis Elia Schenk
$! Project:   Standard Evaluation according to Scanco
$! Date:      December 2022
$! Location:  ISTB University of Bern
$! *************************************************************
$!
$! This IPL script runs the same commands as standard evaluation
$! of Scanco in order to obtain cortical and trabecular masks for
$! subsequent hFE
$!
$! *************************************************************
$!
$! userdefined inputs
$! ------------------
$! P1 = sample ID                                     [00000269]
$! P2 = measurement folder		                        [00002159]
$! P3 = measurement number		                        [C0001901]
$! P4 = unit selection for segmentation               [2 (normal) and 5 if density]
$! P5 = Peel iteration 1 (minimal cortical thickness) [default = 6]
$!
$! Define Parameters
$! -----------------
$ dir
$! Paths:
$ path              :=DK0:[MICROCT.DATA.'P1'.'P2']
$ aim               :='path''P3'.AIM
$ gobj              :='path''P3'.GOBJ
$!
$ cort_mask         :='path''P3'_cort_mask.AIM
$ trab_mask         :='path''P3'_trab_mask.AIM
$ cort_gobj         :='path''P3'_cort_mask.GOBJ
$ trab_gobj         :='path''P3'_trab_mask.GOBJ
$ seg               :='path''P3'_seg.AIM
$!
$! Parameters:
$ unit              := 'P4'
$ peel_iter_1       := 'P5'
$!
$!unit = 2 (Density) for normal images | unit = 5 for Hosseini (if images are already in BMD)
$!
$!
$!
$ loop              :=trab_close
$!(if := trab_loop, this is the long loop)
$!(if := trab_close, this is the short loop)
$!
$!
$ ipl_batch

/db_scanco_activate true

!! ------------------------------------------------------------------------------
!! Read AIM and GOBJ files
!! ------------------------------------------------------------------------------

/read
  -name                      aim
  -filename                  "aim

/gobj_to_aim
  -gobj_filename             "gobj
  -output                    gobj
  -peel_iter                 0

/list true false

!!-----------------------------------------------------------------------------
!! Gross gaussian segmentation to eliminate thin (mostly trabeculae) structures
!! Threshold is high to start with cortical structure as a basis
!!------------------------------------------------------------------------------

/gobj_maskaimpeel_ow
  -input_output              aim
  -gobj_filename             "gobj
  -peel_iter         	     0

!! To not loose any slices through the seg_gauss,
!! one mirrored layer on top and on the bottom is added

! This requires to add a bounding box with empty boarder
! in z-direction...

/bounding_box_cut
  -input                     aim
  -output                    aim_bbc
  -z_only                    false
  -border                    0 0 1

/list true false

! ..and apply an offset and fill it with mirrored pixels

/offset_add aim_bbc
  -add_offset                2

/fill_offset_mirror aim_bbc

!! gauss filtering and segmentation (unit 2 = Dens (mg HA/ccm)

/seg_gauss
  -input                     aim_bbc
  -output                    cort
  -sigma                     2.000000
  -support                   3
  -lower_in_perm_aut_al      500.000000
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            127
  -unit                      "unit

/list true false

!!-----------------------------------------------------------------------------
!! Start of cortex mask_off alternative sequence
!! Start working on trabecular mask
!!------------------------------------------------------------------------------

! by subtracting the filtered aim file (mostly cortical region) from the
! gobj mask, a somehow "negative" trabecular mask is created


/subtract_aims
  -input1                    gobj
  -input2                    cort
  -output                    trab

/list true false

! the newly generated internal file "trab" has a couple of layers at the boundary.
! these layers are now deleated (6 voxels)
!! ensure a minimal cortical thickness of 6 voxels (61um -->  366um)
!! for FRACTIB project, voxel size of 72.5um minimal cortical thickness of 5 voxels (362.5um)

/gobj_maskaimpeel_ow
  -input_output              trab
  -gobj_filename             "gobj
  -peel_iter                 "peel_iter_1

! the values are now inverted (the file now contains value 0 for background and value
! 127 for segmented bone
!!Inversion

/set_value
  -input                     trab
  -value_object              0
  -value_background          127

/rename
  -old_name                  trab
  -new_name                  cort

! the internal file "trab" is renamed to "cort" and a rank operation is applied, which
! delete all the objects not connected to the biggest object. The maskaimpeel removes
! the background, so that at the end the cortical mask is produced.
!! ranking operation to delete unconected volumes

/cl_ow_rank_extract
  -input_output              cort
  -first_rank                1
  -last_rank                 1
  -connect_boundary          true
  -value_in_range            127
  -topology                  6

!!Remove the background
/gobj_maskaimpeel_ow
  -input_output              cort
  -gobj_filename             "gobj
  -peel_iter                 0

/list true false

!!------------------------------------------------------------------------------
!! End of cortex_mask_off alternative sequence
!! Start of trabecular masking
!!------------------------------------------------------------------------------

! We end up with only trabecular structure
!! Creating only trabecular structure

/subtract_aims
  -input1                    gobj
  -input2                    cort
  -output                    trab

/list true false

! Because we need a filled trabecular structure, we remove all the objects, not
! connected to the biggest object

/cl_ow_rank_extract
  -input_output              trab
  -first_rank                1
  -last_rank                 1
  -connect_boundary          true
  -value_in_range            127
  -topology                  6

!! 1st Erosion-CL-Dilation loop on Trab mask to remove Tb speckles not connected to it.
!! Erosion distance is small
!!------------------------------------------------------------------------------------
! metric 11 = 3-4-5 champfer metric
! use_previous_margin TRUE: If PREVIOUS Dilation was performed: use the existing margin
! Output aim was NOT shrunk and is still as big as original aim box

!/bounding_box_cut
!  -input                     trab
!  -output                    trab_bbc
!  -z_only                    false
!  -border                    30 30 0

/list true false

/erosion
  -input                     trab
  -output                    trab_ero
  -erode_distance            3
  -use_previous_margin       false
  -metric                    11

/offset_add
  -input                     trab_ero
  -add_offset                -4 -4 -4

/list true false

! Comment of Scanco: For Erosion, continuous_at_boundary setting is ignored for this step.
! Info for Erosion: if object touches the very border of the AIM box:
!   this artificial 'box cut face' is not eroded.
!   however: an empty pixel boundary prevents this 'erosion suppression'
!   So: object touches box --> treated 'continuous' at that place.
!   So: object had > 1 empty pixel border --> treated 'isolated' there.
!   Hint: you may want to use the /bounding_box operation beforehand
!   to control the above conditions.

/cl_ow_rank_extract
  -input_output              trab_ero
  -first_rank                1
  -last_rank                 1
  -connect_boundary          true
  -value_in_range            127
  -topology                  6


! IMPORTANT: If the only the dilation is done, without a /offsett_add afterwards, the aim box
! would get widend by [4 4 4] pixels

/dilation
  -input                     trab_ero
  -output                    trab_dil
  -dilate_distance           3
  -continuous_at_boundary    0 0 1
  -use_previous_margin       true
  -metric                    11

/list true false

! Large close/open sequence to smooth the contour
! probably close_distance of 10 pixels would be enough! --> verify!
! metric 11 = 3-4-5 champfer metric
! continuous_at_boundary: 0 = isolated object, empty border added; 1 = object is mirrored
! into border. Set this value to 1 only if the object is not surrounded by void space in its
! aim box in this direction.

/close
  -input                     trab_dil
  -output                    trab_closing_op
  -close_distance            15
  -continuous_at_boundary    0 0 1
  -metric                    11

/list true false

/delete
  -name                      trab_dil

/delete
  -name                      trab_ero

/delete
  -name                      cort

/gobj_maskaimpeel_ow
  -input_output              trab_closing_op
  -gobj_filename             "gobj
  -peel_iter                 6

/bounding_box_cut
  -input                     trab_closing_op
  -output                    "loop
  -z_only                    false
  -border                    0 0 0

/list true false

/delete trab_close_op


!! *************************************************************************
!! This loop will only be executed if variable loop := trab_close
!! If this variable is set to trab_short, the processing will continue
!! directly at the end of this loop!
!! *************************************************************************


! At this stage, the Tb mask (trab_bbc) is close to be correct,
! but sharp corners may have been deleted by the open process
! From here we put those corners back.

! The input of the open function is trab_loop. This will only be executed,
! if loop variable is set to this (not trab_close)

/open
  -input                     trab_loop
  -output                    trab_open
  -open_distance             15
  -metric                    11

/list true false

/set_value
  -input                     trab_open
  -value_object              127
  -value_background          0

/set_value
  -input                     trab_loop
  -value_object              127
  -value_background          0

/subtract_aims
  -input1                    trab_loop
  -input2                    trab_open
  -output                    corners

/list true false

!! 2nd Erosion-CL-Dilation loop, for removing small speckles and keep only sharp
!!corner lost by the large open process
!!------------------------------------------------------------------------------------

/erosion
  -input                     corners
  -output                    corn_ero
  -erode_distance            3
  -use_previous_margin       false
  -metric                    11

/list true false

/cl_nr_extract
  -input                     corn_ero
  -output                    corn_cl
  -min_number                800
  -max_number                0
  -value_in_range            127
  -topology                  6

/dilation
  -input                     corn_cl
  -output                    corn_major
  -dilate_distance           3
  -continuous_at_boundary    0 0 0
  -use_previous_margin       false
  -metric                    11

/list true false

/cl_nr_extract
  -input                     corn_major
  -output                    corn_cl
  -min_number                800
  -max_number                0
  -value_in_range            127
  -topology                  6

/cl_nr_extract
  -input                     corn_cl
  -output                    corners
  -min_number                1
  -max_number                500000
  -value_in_range            127
  -topology                  6

/delete
  -name                      corn_cl

/delete
  -name                      corn_major

! Adding the corners to the almost finished trabecular mask

/add_aims
  -input1                    corners
  -input2                    trab_open
  -output                    trab2_loop

/list true false

/close
  -input                     trab2_loop
  -output                    trab_loop
  -close_distance            30
  -continuous_at_boundary    0 0 0
  -metric                    11

/list true false

/delete
  -name                      corners

/delete
  -name                      trab_open

!/delete
!  -name                     trab_bbc

!! *****************************************************************
!! END of the loop!
!! *****************************************************************
!!------------------------------------------------------------------------------
!! Heuristics: Ensuring a minimum cortical thickness (6 voxels), slicewise CL clean up
!! and fixing up the box size of the contour object
!!------------------------------------------------------------------------------

/gobj_maskaimpeel_ow
  -input_output              trab_close
  -gobj_filename             "gobj
  -peel_iter                 6

/cl_slicewise_extractow
  -input_output              trab_close
  -lo_vol_fract_in_perc      50.000000
  -up_vol_fract_in_perc      100.000000
  -value_in_range            127
  -topology                  6

/list true false

!former bbc --> deleted

/rename
  -old_name                  trab_close_bbc
  -new_name                  trab_close


/subtract_aims
  -input1                    gobj
  -input2                    trab_close
  -output                    cort

/list true false

/delete
  -name                      trab_close

/cl_slicewise_extractow
  -input_output              cort
  -lo_vol_fract_in_perc      50.000000
  -up_vol_fract_in_perc      100.000000
  -value_in_range            127
  -topology                  6


/rename cort cort_final

/list true false


!! End of Heuristics
!!------------------------------------------------------------------------------
!! Heuristics: TRAB_MASK.GOBJ = GOBJ - CORT_MASK.GOBJ
!!------------------------------------------------------------------------------

/subtract_aims
  -input1                    gobj
  -input2                    cort_final
  -output                    trab_final

/list true false

!! ------------------------------------------------------------------------------
!! Write masks
!! ------------------------------------------------------------------------------

/write_v020
  -name                      cort_final
  -filename                  "cort_mask
  -compress_type             none
  -version_020               true

/write_v020
  -name                      trab_final
  -filename                  "trab_mask
  -compress_type             none
  -version_020               true

/list true false

!! ------------------------------------------------------------------------------
!! Mask to GOBJ for later segmentation
!! ------------------------------------------------------------------------------

!! Curvature smooth set to 5 = no smoothing applied
/togobj_from_aim
  -input                     cort_final
  -gobj_filename             "cort_gobj
  -min_elements              0
  -max_elements              0
  -curvature_smooth          5

!! Curvature smooth set to 5 = no smoothing applied
/togobj_from_aim
  -input                     trab_final
  -gobj_filename             "trab_gobj
  -min_elements              0
  -max_elements              0
  -curvature_smooth          5

!! ------------------------------------------------------------------------------
!! The two masks (cortical and trabecular) are now generated.
!! Next step is to generate the segmented cortical shell and trabecular structure
!! ------------------------------------------------------------------------------

!! ------------------------------------------------------------------------------
!! Writing out the aim file with equal size as other files
!! ------------------------------------------------------------------------------

/copy
  -in                        gobj
  -out                       gobj_copy

/set_value
  -input                     gobj_copy
  -value_object              1
  -value_background          1

/multiply_volumes
  -input1                    aim
  -input2                    gobj_copy
  -output                    aim_final
  -common_region_only        true
  -multiply_zero_pixels      true

/write_v020
  -name                      aim_final
  -filename                  "aim
  -compress_type             none
  -version_020               true

/list true false

!! =================================================================================================================
!! STEP 2 - SEGMENTATION
!! =================================================================================================================

!! ------------------------------------------------------------------------------
!! Start of cortical segmentation
!! ------------------------------------------------------------------------------

/copy
  -in                        aim_final
  -out                       cort

!! To not loose layers on top and bottom through seg_gauss

/bounding_box_cut
  -input                     cort
  -output                    cort_bbc
  -z_only                    false
  -border                    1 1 1

/list true false

/offset_add
  -input                     cort_bbc
  -new_offset                1 1 1


/fill_offset_mirror cort_bbc

/seg_gauss
  -input                     cort_bbc
  -output                    cort_gauss
  -sigma                     0.8
  -support                   1
  -lower_in_perm_aut_al      450.00
  -upper_in_perm_aut_al      3000.000000
  -value_in_range            127
  -unit                      "unit

/list true false

/gobj_maskaimpeel_ow
  -input_output              cort_gauss
  -gobj_filename             "cort_gobj
  -peel_iter                 0

!/bounding_box_cut
!  -input                     cort_gauss
!  -output                    cort_seg
!  -z_only                    false
!  -border                    0 0 0

/multiply_volumes
  -input1                    cort_gauss
  -input2                    gobj_copy
  -output                    cort_seg
  -common_region_only        false
  -multiply_zero_pixels      true

!!-------------------------------------------------------------------------------
!! End of cortical segmentation
!!-------------------------------------------------------------------------------

!! ------------------------------------------------------------------------------
!! Start of trabecular segmentation
!! ------------------------------------------------------------------------------

/copy
  -in                        aim_final
  -out                       trab

!! To not loose layers on top and bottom through seg_gauss

/bounding_box_cut
  -input                     trab
  -output                    trab_bbc
  -z_only                    true
  -border                    1 1 1

/list true false

/offset_add
  -input                     trab_bbc
  -new_offset                1 1 1

!/delete trab

/fill_offset_mirror trab_bbc

/seg_gauss
  -input                     trab_bbc
  -output                    trab_gauss
  -sigma                     0.8
  -support                   1
  -lower_in_perm_aut_al      320
  -upper_in_perm_aut_al      3000
  -value_in_range            127
  -unit                      "unit

/list true false

/gobj_maskaimpeel_ow
  -input_output              trab_gauss
  -gobj_filename             "trab_gobj
  -peel_iter                 0

!/bounding_box_cut
!  -input                     trab_gauss
!  -output                    trab_seg
!  -z_only                    false
!  -border                    0 0 0

/multiply_volumes
  -input1                    trab_gauss
  -input2                    gobj_copy
  -output                    trab_seg
  -common_region_only        false
  -multiply_zero_pixels      true

!!-----------------------------------------------
!! setting the values for both segmentation files
!!-----------------------------------------------

/set_value
  -input                     cort_seg
  -value_object              2
  -value_background          0

/set_value
  -input                     trab_seg
  -value_object              1
  -value_background          0

/offset_set
  -input                     cort_seg
  -new_offset                0 0 0

/offset_set
  -input                     trab_seg
  -new_offset                0 0 0

/add_aims
  -input1                    cort_seg
  -input2                    trab_seg
  -output                    full

/list true false

!! Write out full segmentation (Cort Seg + Trab Seg)
/write_v020
  -name                      full
  -filename                  "seg
  -compress_type             none
  -version_020               true

/list true false

!! =================================================================================================================
!! END OF SCRIPT
!! =================================================================================================================


!! ------------------------------------------------------------------------------
!! End of Standard Evaluation
!! ------------------------------------------------------------------------------
!! Files produced:
!! - AIM file (grey values)
!! - cortical mask
!! - trabecular mask
!! - full segmentation (0 = background, 1= trabecular bone 2= cortical bone
!! ------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------

..
$
$ exit
