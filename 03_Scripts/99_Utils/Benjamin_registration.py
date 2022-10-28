import numpy as np
import matplotlib.pyplot as plt
import SimpleITK as sitk
from importlib.machinery import SourceFileLoader
import scipy.ndimage as ndimage
import cv2
import transformations as trans
import os
import transformations as trans
import time
from skimage.measure import find_contours
import pandas as pd
from numba import njit
from multiprocessing import Pool
# import PID
from simple_pid import PID
from tqdm import tqdm
from scipy import stats

pro = SourceFileLoader("image_processing",
                       "/home/ben/kDrive/Artorg/PhD/FemHals_image_processing/02_Python/image_processing.py").load_module()
read = SourceFileLoader("image_reader",
                        "/home/ben/kDrive/Artorg/PhD/FemHals_image_processing/02_Python/image_reader.py").load_module()


def center_of_image(img_sitk):
    return (int(0.5 * img_sitk.GetSpacing()[0] * img_sitk.GetSize()[0] - 1),
            int(0.5 * img_sitk.GetSpacing()[1] * img_sitk.GetSize()[1] - 1),
            int(0.5 * img_sitk.GetSpacing()[2] * img_sitk.GetSize()[2] - 1))


def center_of_gravity(img, resolution=np.array([1, 1, 1])):
    import scipy.ndimage as ndimage
    if type(img) == np.ndarray:
        pass
    else:
        img = sitk.GetArrayFromImage(img).transpose(2, 1, 0)
    CG_temp = ndimage.measurements.center_of_mass(img)
    return [CG_temp[0] * resolution[0], CG_temp[1] * resolution[1], CG_temp[2] * resolution[2]]


def add_layers(img, add_layer):
    if type(img) != np.ndarray:
        img_np = sitk.GetArrayFromImage(img).transpose(2, 1, 0)
        img_sitk = img
    else:
        img_np = img
        img_sitk = sitk.GetImageFromArray(img.transpose(2, 1, 0))
        print
        '         !! be careful: resolution lost because input = np image !!'
    img_bigger_np = np.zeros(())
    dim_x, dim_y, dim_z = int(img_sitk.GetSize()[0] * add_layer[0]), int(img_sitk.GetSize()[1] * add_layer[1]), int(
        img_sitk.GetSize()[2] * add_layer[2])
    dimension_moving = img_sitk.GetSize()
    del img_sitk
    img_bigger_np = np.zeros((dim_x, dim_y, dim_z))
    img_bigger_np[
    int(dim_x / 2 - dimension_moving[0] / 2):int(dim_x / 2 - dimension_moving[0] / 2) + dimension_moving[0],
    int(dim_y / 2 - dimension_moving[1] / 2):int(dim_y / 2 - dimension_moving[1] / 2) + dimension_moving[1],
    int(dim_z / 2 - dimension_moving[2] / 2):int(dim_z / 2 - dimension_moving[2] / 2) + dimension_moving[2]] = img_np
    del img_np
    img_bigger_sitk = sitk.GetImageFromArray(img_bigger_np.transpose(2, 1, 0))
    img_bigger_sitk = sitk.Cast(img_bigger_sitk, sitk.sitkFloat32)
    try:
        img_bigger_sitk.SetSpacing(img.GetSpacing())
    except:
        pass
    return img_bigger_sitk


def transform_R_t(img_sitk, R, t, center='gravity', add_layer=[1.2, 1.2, 1], background=0):
    if all(np.array(add_layer).astype(float) != np.array([1.0, 1., 1.0])):
        img_sitk = add_layers(img_sitk, add_layer)
    affine = sitk.AffineTransform(3)
    affine.SetMatrix(R.ravel())
    if center == 'gravity':
        COG_voxel = np.array(center_of_gravity(img_sitk))
        COG_mm = (img_sitk.GetSpacing()[0] * COG_voxel[0], img_sitk.GetSpacing()[1] * COG_voxel[1],
                  img_sitk.GetSpacing()[2] * COG_voxel[2])
        affine.SetCenter(COG_mm)
    if center == 'center':
        center_image = (np.rint(np.array(img_sitk.GetSize())[0] / 2).astype('int'),
                        np.rint(np.array(img_sitk.GetSize())[1] / 2).astype('int'),
                        np.rint(np.array(img_sitk.GetSize())[2] / 2).astype('int'))
        center_image_mm = np.array(center_image) * img_sitk.GetSpacing()[0]
        affine.SetCenter(center_image_mm)
    if not isinstance(center, str):
        affine.SetCenter(center)
    if background != 0:
        if background.upper() == 'MEAN':
            statistics_image_filter = sitk.StatisticsImageFilter()
            statistics_image_filter.Execute(img_sitk)
            background = statistics_image_filter.GetMean()
    print(f'   ... transformation background = {background:.2f}')

    resampled_temp = sitk.Resample(img_sitk, img_sitk, affine, sitk.sitkLinear, background)
    affine = sitk.AffineTransform(3)
    affine.SetTranslation(-1 * np.array(t))
    resampled = sitk.Resample(resampled_temp, resampled_temp, affine, sitk.sitkLinear, background)
    del resampled_temp
    resampled.SetSpacing(img_sitk.GetSpacing())
    pro.img_description(resampled)
    return resampled


def get_parametermap(map_type='rigid', background=0.0):
    """
    input:
    """
    if map_type.upper() == 'RIGID':
        # rigid:
        ParameterMap_rigid = sitk.ParameterMap()
        ParameterMap_rigid['FixedInternalImagePixelType'] = ("short",)
        ParameterMap_rigid['MovingInternalImagePixelType'] = ("short",)

        ParameterMap_rigid['AutomaticParameterEstimation'] = ('true',)
        ParameterMap_rigid['AutomaticScalesEstimation'] = ('true',)
        ParameterMap_rigid['AutomaticTransformInitialization'] = ('true',)
        ParameterMap_rigid['HowToCombineTransforms'] = ('Compose',)
        ParameterMap_rigid['NumberOfHistogramBins'] = ('4', '8', '16', '32', '64', '128')  ###

        ParameterMap_rigid['CheckNumberOfSamples'] = ('true',)
        ParameterMap_rigid['DefaultPixelValue'] = (f'{background}',)
        ParameterMap_rigid['FinalBSplineInterpolationOrder'] = ('3',)
        ParameterMap_rigid['BSplineInterpolationOrder'] = ('3',)
        # ParameterMap_rigid['FixedImagePyramid'] = ('FixedSmoothingImagePyramid',)
        ParameterMap_rigid['FixedImagePyramid'] = ('FixedRecursiveImagePyramid',)
        ParameterMap_rigid['ImageSampler'] = ('RandomCoordinate',)
        ParameterMap_rigid['Interpolator'] = ('LinearInterpolator',)
        ParameterMap_rigid['MaximumNumberOfIterations'] = ('3500',)  ###
        ParameterMap_rigid['MaximumNumberOfSamplingAttempts'] = ('6',)
        ParameterMap_rigid['Metric'] = ('AdvancedMattesMutualInformation',)
        # ParameterMap_rigid['MovingImagePyramid'] = ('MovingSmoothingImagePyramid',)
        ParameterMap_rigid['MovingImagePyramid'] = ('MovingRecursiveImagePyramid',)
        ParameterMap_rigid['NewSamplesEveryIteration'] = ('true',)
        ParameterMap_rigid['NumberOfResolutions'] = ('7',)
        # ParameterMap_rigid['GridSpacingSchedule'] = (
        # '8.0', '8.0', '8.0', '4.0', '4.0', '4.0', '2.0', '2.0', '2.0', '1.0', '1.0', '1.0',)
        ParameterMap_rigid['NumberOfSamplesForExactGradient'] = ('4096',)
        ParameterMap_rigid['NumberOfSpatialSamples'] = ('2048',)
        ParameterMap_rigid['Optimizer'] = ('AdaptiveStochasticGradientDescent',)
        ParameterMap_rigid['Registration'] = ('MultiResolutionRegistration',)
        ParameterMap_rigid['ResampleInterpolator'] = ('FinalBSplineInterpolator',)
        ParameterMap_rigid['Resampler'] = ('DefaultResampler',)
        ParameterMap_rigid['ResultImageFormat'] = ('mhd',)
        ParameterMap_rigid['Transform'] = ('EulerTransform',)
        ParameterMap_rigid['WriteIterationInfo'] = ('false',)
        ParameterMap_rigid['WriteResultImage'] = ('true',)
        return ParameterMap_rigid

    if map_type.upper() == 'BSPLINE':
        ParameterMap_bspline = sitk.ParameterMap()
        ParameterMap_bspline['FixedInternalImagePixelType'] = ("short",)
        ParameterMap_bspline['MovingInternalImagePixelType'] = ("short",)
        ParameterMap_bspline['FixedImagePyramid'] = ("FixedRecursiveImagePyramid",)
        ParameterMap_bspline['MovingImagePyramid'] = ("MovingRecursiveImagePyramid",)
        ParameterMap_bspline['AutomaticParameterEstimation'] = ('true',)
        ParameterMap_bspline['CheckNumberOfSamples'] = ('true',)
        ParameterMap_bspline['DefaultPixelValue'] = (f'{background}',)
        ParameterMap_bspline['BSplineInterpolationOrder'] = ('3',)
        ParameterMap_bspline['FinalBSplineInterpolationOrder'] = ('3',)
        ParameterMap_bspline['FinalGridSpacingInVoxels'] = (
            '24',)  # ParameterMap_bspline['FinalGridSpacingInPhysicalUnits'] = ('8.000000',)
        # ParameterMap_bspline['FixedImagePyramid'] = ('FixedSmoothingImagePyramid',)
        # ParameterMap_bspline['GridSpacingSchedule'] = ('4.0', '4.0', '4.0', '2.0', '2.0', '2.0', '1.0', '1.0','1.0',)  # ('2.803221', '1.988100', '1.410000', '1.000000')

        ParameterMap_bspline['NumberOfHistogramBins'] = ('16', '32', '64')  # , '128')   ###

        ParameterMap_bspline['ImageSampler'] = ('RandomCoordinate',)
        ParameterMap_bspline['Interpolator'] = ('LinearInterpolator',)
        ParameterMap_bspline['MaximumNumberOfIterations'] = ('2000',)  # ('256',)
        ParameterMap_bspline['MaximumNumberOfSamplingAttempts'] = ('5',)  ###
        ParameterMap_bspline['Metric'] = ('AdvancedMattesMutualInformation')
        ParameterMap_bspline['Metric0Weight'] = ('1.0',)
        ParameterMap_bspline['Metric1Weight'] = ('1.0',)
        ParameterMap_bspline['MovingImagePyramid'] = ('MovingSmoothingImagePyramid',)
        ParameterMap_bspline['ImagePyramidSchedule'] = ('8', '8', '8', '4', '4', '4', '2', '2', '2', '1', '1', '1',)
        ParameterMap_bspline['NewSamplesEveryIteration'] = ('true',)
        ParameterMap_bspline['NumberOfResolutions'] = ('6',)
        ParameterMap_bspline['NumberOfSamplesForExactGradient'] = ('4096',)
        ParameterMap_bspline['NumberOfSpatialSamples'] = ('2048',)  ###
        ParameterMap_bspline['Optimizer'] = ('AdaptiveStochasticGradientDescent',)
        ParameterMap_bspline['Registration'] = ('MultiResolutionRegistration',)
        # ParameterMap_bspline['Registration'] = ("MultiResolutionRegistration",)
        ParameterMap_bspline['Interpolator'] = ("BSplineInterpolator",)
        ParameterMap_bspline['Metric'] = ('AdvancedMattesMutualInformation',)
        ParameterMap_bspline['ResampleInterpolator'] = ('FinalBSplineInterpolator',)
        ParameterMap_bspline['Resampler'] = ('DefaultResampler',)
        ParameterMap_bspline['HowToCombineTransforms'] = ('Compose',)

        ParameterMap_bspline['ResultImageFormat'] = ('mhd',)
        # ParameterMap_bspline['SP_A'] = ('100',)
        # ParameterMap_bspline['SP_alpha'] = ('0.6',)
        ParameterMap_bspline['Transform'] = ('BSplineTransform',)
        ParameterMap_bspline['WriteIterationInfo'] = ('false',)
        ParameterMap_bspline['WriteResultImage'] = ('true',)
        return ParameterMap_bspline


def get_transformation(ElastixImageFilter):
    # https://lightrun.com/answers/fepegar-torchio-add-high-level-brain-mri-transforms-for-preprocessing
    sitk.WriteParameterFile(ElastixImageFilter.GetTransformParameterMap()[0], '/tmp/param0.txt')

    with open('/tmp/param0.txt', "r") as f:
        for line in f:
            if "(TransformParameters" in line:
                transformation_params = np.array(
                    line.split("ers ")[1].split(")")[0].split(" "), np.float64
                )
            if "(CenterOfRotationPoint" in line:
                center_rotation = np.array(
                    line.split("int ")[1].split(")")[0].split(" "), np.float64
                )
    os.remove('/tmp/param0.txt')
    # print(f'center_rotation = {center_rotation}')
    trn = sitk.Euler3DTransform()
    trn.SetParameters((transformation_params[0], transformation_params[1], transformation_params[2], 0, 0, 0))
    C = np.matrix(np.array([[center_rotation[0]], [center_rotation[1]], [center_rotation[2]]]))
    T = np.eye(4, dtype=np.float64)
    T[:3, :3] = np.matrix(np.reshape(np.array(trn.GetMatrix(), dtype=np.float64), [3, 3], order='F'))
    T[:3, 3:4] = C - T[:3, :3] * C - T[:3, :3] * np.matrix(
        np.array([[transformation_params[3]], [transformation_params[4]], [transformation_params[5]]]))
    T_ = trans.euler_matrix(transformation_params[0], transformation_params[1], transformation_params[2], 'rxyz')
    print(np.array([transformation_params[3], transformation_params[4], transformation_params[5]]))
    T_[3, :3] = np.array([transformation_params[3], transformation_params[4], transformation_params[5]])
    # print(transformation_params)
    return T, T_


def TransformixTransformations(MovingImage, TransformParameterMap, ResultsDirectory=False):
    ## Compute jacobian of deformation field using transformix
    TransformixImageFilter = sitk.TransformixImageFilter()
    # TransformixImageFilter.ComputeDeformationFieldOn()
    TransformixImageFilter.ComputeSpatialJacobianOn()
    # TransformixImageFilter.ComputeDeterminantOfSpatialJacobianOn()
    TransformixImageFilter.SetMovingImage(sitk.GetImageFromArray(MovingImage.transpose(2, 1, 0)))
    TransformixImageFilter.SetTransformParameterMap(TransformParameterMap)
    TransformixImageFilter.SetOutputDirectory(ResultsDirectory)

    TransformixImageFilter.Execute()


def dice(fixed_img, moving_img):
    if type(fixed_img) != np.ndarray:
        fixed_img = sitk.GetArrayFromImage(fixed_img).transpose(2, 1, 0)
    if type(moving_img) != np.ndarray:
        moving_img = sitk.GetArrayFromImage(moving_img).transpose(2, 1, 0)
    DSC = 2 * np.sum(fixed_img * moving_img) / np.sum(fixed_img + moving_img)
    return DSC


# https://towardsdatascience.com/metrics-to-evaluate-your-semantic-segmentation-model-6bcb99639aa2
def iou(y_true, y_pred, smooth=1):
    intersection = np.sum(np.abs(y_true * y_pred))
    union = np.sum(y_true) + np.sum(y_pred) - intersection
    iou = np.mean((intersection + smooth) / (union + smooth))
    return iou


def dice_(y_true, y_pred, smooth=1):
    intersection = np.sum(np.abs(y_true * y_pred))
    total_area = np.sum(y_true) + np.sum(y_pred)
    dice = np.mean((2. * intersection + smooth) / (total_area + smooth))
    return dice


def split(arys, sections, axis=[0, 1, 2], return_index=False):
    if sections == 0:
        sections = 1
    if not isinstance(arys, list):
        arys = [arys]
    for ax in axis:
        arys = [np.array_split(ary, sections, axis=ax) for ary in arys]
        arys = [ary for aa in arys for ary in aa]  # Flatten
    if return_index:
        order = np.arange(sections ** 3).reshape(2, 2, 2, order='F').flatten()

        index_list = [ary.shape for ary in arys]
        return arys, index_list
    else:
        return arys


def display_time(delta_tictoc, title=''):
    m, s = divmod(delta_tictoc, 60)
    h, m = divmod(m, 60)
    print(f"duration {title}: {h:.0f} h, {m:.0f} min {s:.2f} sec")


start_time_1 = time.time()

# dim = (100, 100, 100)
# img_index = np.empty(dim).astype(
#     f'<U{int(np.ceil(np.log10(dim[0])) + np.ceil(np.log10(dim[1])) + np.ceil(np.log10(dim[1])) + 2)}')
# for x in np.arange(img_index.shape[0]):
#     for y in np.arange(img_index.shape[1]):
#         for z in np.arange(img_index.shape[2]):
#             img_index[x, y, z] = f'{x}.{y}.{z}'


# pre_post_pairs_list = [
#     ['C0002153', 'C0002402'],
#     ['C0002156', 'C0002404'],
#     ['C0002159', 'C0002406'],
#     ['C0002161', 'C0002411'],
#     ['C0002167', 'C0002408'],
# ]

pre_post_pairs_list = [
    ['C0002076', 'C0002374'],  # 1 -> 391_L
    ['C0002085', 'C0002376'],  # 2 -> 409_L [::-1] in Z and flipped in X  post
    ['C0002086', 'C0002375'],  # 3 -> 409_R
    ['C0002090', 'C0002377'],  # 4 -> 419_L
    ['C0002091', 'C0002378'],  # 5 -> 419_R
    ['C0002097', 'C0002379'],  # 6 -> 433_R
    ['C0002114', 'C0002380'],  # 7 -> 400_L
    ['C0002115', 'C0002381'],  # 8 -> 404_R
    ['C0002121', 'C0002382'],  # 9 -> 402_R
    ['C0002122', 'C0002384'],  # 10 -> 401_R
    ['C0002123', 'C0002385'],  # 11 -> 415_R
    ['C0002125', 'C0002386'],  # 12 -> 403_L
    ['C0002129', 'C0002387'],  # 13 -> 421_L
    ['C0002130', 'C0002388'],  # 14 -> 392_R
    ['C0002137', 'C0002389'],  # 15 -> 392_L
    ['C0002140', 'C0002390'],  # 16 -> 394_R
    ['C0002141', 'C0002392'],  # 17 -> 394_L
    ['C0002142', 'C0002391'],  # 18 -> 397_R
    ['C0002143', 'C0002393'],  # 19 -> 389_L
    ['C0002144', 'C0002396'],  # 20 -> 397_L
    ['C0002145', 'C0002397'],  # 21 -> 389_R
    ['C0002146', 'C0002398'],  # 22 -> 385_L
    ['C0002148', 'C0002399'],  # 23 -> 393_R
    ['C0002150', 'C0002400'],  # 24 -> 386_L
    ['C0002151', 'C0002401'],  # 25 -> 390_L
    ['C0002153', 'C0002402'],  # 26 -> 399_R
    ['C0002155', 'C0002403'],  # 27 -> 402_L
    ['C0002156', 'C0002404'],  # 28 -> 390_R
    ['C0002157', 'C0002405'],  # 29 -> 399_L
    ['C0002158', 'C0002410'],  # 30 -> 388_R transpose(2, 0, 1) post
    ['C0002159', 'C0002406'],  # 31 -> 387_L
    ['C0002161', 'C0002411'],  # 32 -> 395_L
    ['C0002167', 'C0002408'],  # 33 -> 388_L
    ['C0002168', 'C0002409'],  # 34 -> 385_R
    ['C0002169', 'C0002412'],  # 35 -> 395_R
    ['C0002170', 'C0002413'],  # 36 -> 398_R
    ['C0002172', 'C0002414'],  # 37 -> 403_R
    ['C0002173', 'C0002415'],  # 38 -> 404_L
    ['C0002174', 'C0002416'],  # 39 -> 400_R
    ['C0002176', 'C0002417'],  # 40 -> 416_R
    ['C0002177', 'C0002418'],  # 41 -> 401_L
    ['C0002365', 'C0002419'],  # 42 -> 437_L

]

parent_dir = '/home/ben/kDrive/Artorg/PhD/FemHals_fracture_line/01_Data/'
child_pre_gray_dir = '02_Pre_Fracture_uCT_GREY_DOWNSCALED_0_065_mm_BBCUT_/'
child_post_gray_dir = '02_Post_Fracture_uCT_GREY_DOWNSCALED_0_065_mm_BBCUT_/'
child_pre_seg_dir = '03_Pre_Fracture_uCT_SEG_DOWNSCALED_0_065_mm_BBCUT_/'
child_post_seg_dir = '03_Post_Fracture_uCT_SEG_DOWNSCALED_0_065_mm_BBCUT_/'
child_pre_mask_dir = '04_Pre_Fracture_uCT_MASK_DOWNSCALED_0_065_mm_BBCUT_/'
child_post_mask_dir = '04_Post_Fracture_uCT_MASK_DOWNSCALED_0_065_mm_BBCUT_/'

ResultsDirectory = '/home/ben/kDrive/Artorg/PhD/FemHals_fracture_line/01_Data/'
# ResultsDirectory = '/home/ben/Documents/18_Femhals_Fracture_line/03_Results/'
# ResultsDirectory_Jac = ResultsDirectory + 'Jacobian'
# ResultsDirectory_Det = ResultsDirectory + 'Determinant'
# ResultsDirectory_Def = ResultsDirectory + 'Deformation'

# todo: start the child position from the translated/rotated parent
# todo: check the geometry of the subdivisions (cubes / indexes)

code_dic = {}
counter = 0
# debug = True
debug = False
# plot = False
plot = True
write = False
for pair in pre_post_pairs_list[25:26]:
    filename_pre, filename_post = pair
    outputImageFileName = parent_dir + 'deformation.mhd'

    # read:_____________________________________________________________________________________________________
    img_gray_pre_np, header_pre = read.read(parent_dir + child_pre_gray_dir + filename_pre + '.mhd')
    img_gray_post_np, header_post = read.read(parent_dir + child_post_gray_dir + filename_post + '.mhd')
    img_seg_pre_np, header_pre = read.read(parent_dir + child_pre_seg_dir + filename_pre + '.mhd')
    img_seg_post_np, header_post = read.read(parent_dir + child_post_seg_dir + filename_post + '.mhd')
    # img_mask_pre_np, header_pre = read.read(parent_dir + child_pre_mask_dir + filename_pre + '.mhd')
    # img_mask_post_np, header_post = read.read(parent_dir + child_post_mask_dir + filename_post + '.mhd')
    if pair[1] == 'C0002376':
        img_gray_post_np = np.flip(img_gray_post_np[:, :, ::-1], 0)
        img_seg_post_np = np.flip(img_seg_post_np[:, :, ::-1], 0)
    if pair[1] == 'C0002410':
        img_gray_post_np = img_gray_post_np.transpose(2, 0, 1)
        img_seg_post_np = img_seg_post_np.transpose(2, 0, 1)

    # read.plot_stack(img_gray_pre_np, f'{pair[0]} img_gray_pre_np')
    # read.plot_stack(img_gray_post_np, f'{pair[1]} img_gray_post_np')
    # read.plot_stack(img_seg_pre_np, f'{pair[0]} img_seg_pre_np')
    # read.plot_stack(img_seg_post_np, f'{pair[1]} img_seg_post_np')
    # read.plot_stack(img_mask_pre_np, 'img_mask_pre_np')
    # read.plot_stack(img_mask_post_np, 'img_mask_post_np')

    # img_seg_pre_np, threshold_pre = pro.segmentation(img_gray_pre_np, header_post)
    # img_seg_post_np, threshold_post = pro.segmentation(img_gray_post_np, header_post)

    add_layer = 120
    img_gray_pre_np = pro.add_layers(img_gray_pre_np, add_layer, background=np.median(img_gray_pre_np))
    img_gray_post_np = pro.add_layers(img_gray_post_np, add_layer, background=np.median(img_gray_post_np))
    img_seg_pre_np = pro.add_layers(img_seg_pre_np, add_layer)
    img_seg_post_np = pro.add_layers(img_seg_post_np, add_layer)
    # img_mask_pre_np = pro.add_layers(img_mask_pre_np, add_layer)
    # img_mask_post_np = pro.add_layers(img_mask_post_np, add_layer)

    # read.plot_stack(img_gray_pre_np, 'img_gray_pre_np')
    # read.plot_stack(img_gray_post_np, 'img_gray_post_np')
    # read.plot_stack(img_seg_pre_np, 'img_seg_pre_np')
    # read.plot_stack(img_seg_post_np, 'img_seg_post_np')
    # read.plot_stack(img_mask_pre_np, 'img_mask_pre_np')
    # read.plot_stack(img_mask_post_np, 'img_mask_post_np')

    # center image at center of gravity with same shape:
    pre_seg_cog = np.rint(ndimage.measurements.center_of_mass(img_seg_pre_np)).astype(int)
    post_seg_cog = np.rint(ndimage.measurements.center_of_mass(img_seg_post_np)).astype(int)

    x_m = np.max([pre_seg_cog[0], post_seg_cog[0]]) - add_layer
    x_p = np.max(
        [img_seg_pre_np.shape[0] - pre_seg_cog[0] - add_layer, img_seg_post_np.shape[0] - post_seg_cog[0] - add_layer])
    y_m = np.max([pre_seg_cog[1], post_seg_cog[1]]) - add_layer
    y_p = np.max(
        [img_seg_pre_np.shape[1] - pre_seg_cog[1] - add_layer, img_seg_post_np.shape[1] - post_seg_cog[1] - add_layer])
    z_m = np.max([pre_seg_cog[2], post_seg_cog[2]]) - add_layer
    z_p = np.max(
        [img_seg_pre_np.shape[2] - pre_seg_cog[2] - add_layer, img_seg_post_np.shape[2] - post_seg_cog[2] - add_layer])

    img_gray_pre_np = img_gray_pre_np[pre_seg_cog[0] - x_m:pre_seg_cog[0] + x_p,
                      pre_seg_cog[1] - y_m:pre_seg_cog[1] + y_p, pre_seg_cog[2] - z_m:pre_seg_cog[2] + z_p]
    img_gray_post_np = img_gray_post_np[post_seg_cog[0] - x_m:post_seg_cog[0] + x_p,
                       post_seg_cog[1] - y_m:post_seg_cog[1] + y_p, post_seg_cog[2] - z_m:post_seg_cog[2] + z_p]
    img_seg_pre_np = img_seg_pre_np[pre_seg_cog[0] - x_m:pre_seg_cog[0] + x_p,
                     pre_seg_cog[1] - y_m:pre_seg_cog[1] + y_p, pre_seg_cog[2] - z_m:pre_seg_cog[2] + z_p]
    img_seg_post_np = img_seg_post_np[post_seg_cog[0] - x_m:post_seg_cog[0] + x_p,
                      post_seg_cog[1] - y_m:post_seg_cog[1] + y_p, post_seg_cog[2] - z_m:post_seg_cog[2] + z_p]
    # img_mask_pre_np = img_mask_pre_np[pre_seg_cog[0] - x_m:pre_seg_cog[0] + x_p,
    #                   pre_seg_cog[1] - y_m:pre_seg_cog[1] + y_p, pre_seg_cog[2] - z_m:pre_seg_cog[2] + z_p]
    # img_mask_post_np = img_mask_post_np[post_seg_cog[0] - x_m:post_seg_cog[0] + x_p,
    #                    post_seg_cog[1] - y_m:post_seg_cog[1] + y_p, post_seg_cog[2] - z_m:post_seg_cog[2] + z_p]

    # add_layer = 10
    # img_gray_pre_np = pro.add_layers(img_gray_pre_np, add_layer, background=np.median(img_gray_pre_np))
    # img_gray_post_np = pro.add_layers(img_gray_post_np, add_layer, background=np.median(img_gray_post_np))
    # img_seg_pre_np = pro.add_layers(img_seg_pre_np, add_layer)
    # img_seg_post_np = pro.add_layers(img_seg_post_np, add_layer)
    # img_mask_pre_np = pro.add_layers(img_mask_pre_np, add_layer)
    # img_mask_post_np = pro.add_layers(img_mask_post_np, add_layer)

    img_seg_fixed_sitk = sitk.GetImageFromArray(img_seg_pre_np.transpose(2, 1, 0))
    img_gray_fixed_sitk = sitk.GetImageFromArray(img_gray_pre_np.transpose(2, 1, 0))
    img_seg_moving_sitk = sitk.GetImageFromArray(img_seg_post_np.transpose(2, 1, 0))
    img_gray_moving_sitk = sitk.GetImageFromArray(img_gray_post_np.transpose(2, 1, 0))


    # ______________________________________________________________
    # rigid registration
    moving_sitk = img_gray_moving_sitk
    fixed_sitk = img_gray_fixed_sitk

    header_results = header_pre
    seg_fixed_np, thr = pro.segmentation(sitk.GetArrayFromImage(fixed_sitk).transpose(2, 1, 0), header_results)

    fixed_sitk.SetSpacing(header_pre['resolution'])
    moving_sitk.SetSpacing(header_pre['resolution'])
    fixed_sitk.SetOrigin((0, 0, 0))
    moving_sitk.SetOrigin((0, 0, 0))
    moving_np = sitk.GetArrayFromImage(moving_sitk).transpose(2, 1, 0)
    fixed_np = sitk.GetArrayFromImage(fixed_sitk).transpose(2, 1, 0)

    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetFixedImage(fixed_sitk)
    ElastixImageFilter.SetMovingImage(moving_sitk)
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)

    ParameterMap_rigid_global = get_parametermap('rigid', background=np.median(img_gray_post_np))
    ElastixImageFilter.SetParameterMap(ParameterMap_rigid_global)
    ElastixImageFilter.Execute()

    img_gray_rigid_registered_sitk = ElastixImageFilter.GetResultImage()

    img_gray_rigid_registered_sitk = sitk.Cast(img_gray_rigid_registered_sitk, sitk.sitkFloat64)
    img_gray_rigid_registered_np = sitk.GetArrayFromImage(img_gray_rigid_registered_sitk).transpose(2, 1, 0)
    img_seg_rigid_registered_np, threshold_result = pro.segmentation(img_gray_rigid_registered_np, header_post)

    if plot:
        read.plot([img_seg_pre_np, img_seg_rigid_registered_np], f'{pair[0]} rigid')

    # ______________________________________________________________
    # nonrigid registration
    moving_sitk = img_gray_rigid_registered_sitk
    fixed_sitk = img_gray_fixed_sitk

    header_results = header_pre
    seg_fixed_np, thr = pro.segmentation(sitk.GetArrayFromImage(fixed_sitk).transpose(2, 1, 0), header_results)

    fixed_sitk.SetSpacing(header_pre['resolution'])
    moving_sitk.SetSpacing(header_pre['resolution'])
    fixed_sitk.SetOrigin((0, 0, 0))
    moving_sitk.SetOrigin((0, 0, 0))
    moving_np = sitk.GetArrayFromImage(moving_sitk).transpose(2, 1, 0)
    fixed_np = sitk.GetArrayFromImage(fixed_sitk).transpose(2, 1, 0)

    ElastixImageFilter = sitk.ElastixImageFilter()
    ElastixImageFilter.SetFixedImage(fixed_sitk)
    ElastixImageFilter.SetMovingImage(moving_sitk)
    ElastixImageFilter.SetOutputDirectory(ResultsDirectory)

    statistics_image_filter = sitk.StatisticsImageFilter()
    statistics_image_filter.Execute(img_gray_fixed_sitk)

    bckd = statistics_image_filter.GetMean()

    ParameterMap_nonrigid = get_parametermap('bspline', background=bckd)
    ElastixImageFilter.SetParameterMap(ParameterMap_nonrigid)
    ElastixImageFilter.Execute()

    TransformParameterMap = ElastixImageFilter.GetTransformParameterMap()

    img_gray_nonrigid_registered_sitk = ElastixImageFilter.GetResultImage()
    img_gray_nonrigid_registered_sitk = sitk.Cast(img_gray_nonrigid_registered_sitk, sitk.sitkFloat64)
    img_gray_nonrigid_registered_np = sitk.GetArrayFromImage(img_gray_nonrigid_registered_sitk).transpose(2, 1, 0)
    img_seg_nonrigid_registered_np, threshold_result = pro.segmentation(img_gray_nonrigid_registered_np, header_post)

    # dice:
    dice_coef = dice(img_seg_pre_np, img_seg_nonrigid_registered_np)
    TransformixTransformations(moving_np, TransformParameterMap, ResultsDirectory)

    JacobianImage = sitk.ReadImage(ResultsDirectory + f'/fullSpatialJacobian.mhd')
    JacobianImage_np = sitk.GetArrayFromImage(JacobianImage).transpose(2, 1, 0, 3)

    # Unimodular decomposition:
    # unimodular decompostion of             F: F = J_ F_hat  with J_ = J¹/³I  and  F_hat = J⁻¹/³F
    # unimodular decompostion of C (C = F.T F): C = J_ C_hat  with J_ = J²/³I  and  C_hat = J⁻²/³ F.T F
    # volumic -> J²/³ - 1 or ||J_|| - 1 -> ||J²/³ I|| - 1 or J - 1 ?
    # deviatoric -> ||E_hat|| = ||0.5(C_hat - I)|| = ||0.5(J⁻²/³ F.T F  - I)||

    tic = time.time()
    F_np = JacobianImage_np.reshape(-1, 9).reshape(-1, 3, 3)
    display_time(time.time() - tic, title='F_np')

    tic = time.time()
    J_np = np.linalg.det(F_np)
    display_time(time.time() - tic, title='J_np')
    #
    # tic = time.time()
    # J__np = np.einsum("i,ijk->ijk", J_np ** (2 / 3), np.tile(np.eye(3).flatten(), len(J_np)).reshape(-1, 3, 3))
    # display_time(time.time() - tic, title='J__np')

    t = time.time()
    C_hat_np = np.einsum("i,ijk->ijk", J_np ** (-2 / 3), np.einsum("ilj,ilk->ijk", F_np, F_np))
    display_time(time.time() - t, title='C_hat')

    t = time.time()
    E_hat_norm_np = np.linalg.norm(0.5 * (C_hat_np - np.eye(3)), axis=(1, 2))
    display_time(time.time() - t, title='E_hat_norm')

    img_J_np = J_np.reshape(JacobianImage_np.shape[:3])
    img_E_hat_norm_np = E_hat_norm_np.reshape(JacobianImage_np.shape[:3])

    img_summary = np.hstack([np.vstack([0.2 * (1 + img_seg_pre_np - img_seg_rigid_registered_np),
                                        0.2 * (1 + img_seg_pre_np - img_seg_nonrigid_registered_np)]),
                             np.vstack([.5 * img_J_np + .1 * img_seg_rigid_registered_np,
                                        img_E_hat_norm_np + img_seg_rigid_registered_np * .2])])

    if write:
        read.write_MHD(img_E_hat_norm_np, ResultsDirectory + f'{pair[0]}_{pair[1]}_E_hat_norm_post2pre.mhd', write_float=True)
        read.write_MHD(img_E_hat_norm_np + img_seg_rigid_registered_np * .2,
                       ResultsDirectory + f'{pair[0]}_{pair[1]}_E_hat_norm_superposed_post2pre.mhd',
                       write_float=True)
        read.write_MHD(img_summary, ResultsDirectory + f'{pair[0]}_{pair[1]}_summary_post2pre.mhd', write_float=True)
    if plot:
        read.plot(img_summary)
        read.plot([img_seg_pre_np, img_seg_nonrigid_registered_np, img_E_hat_norm_np])
        read.plot([img_seg_pre_np, img_E_hat_norm_np, img_seg_nonrigid_registered_np])

duration_sec = time.time() - start_time_1
m, s = divmod(duration_sec, 60)
h, m = divmod(m, 60)
print("duration: %d h %02d min %.2f sec" % (h, m, s))
