"""
This file is used to read ABaqus ODB files for the pipeline PSL_Accurate. It differs a bit from the readODB.py of
the pipeline PSL_Fast, as in the accurate pipeline we use distinct phases for cortex and trabecular bone.

"""

import sys
from odbAccess import *
from abaqusConstants import *
import numpy as np
import csv
import os

np.set_printoptions(precision=None, suppress=None)

# read sys args
# -------------------------------------------
odbSample = sys.argv[1]
odbPath = sys.argv[2]
fieldVar = sys.argv[3]
target_frame_number = int(sys.argv[4])

odbOptimPath = odbPath + "/" + "Optimization"

try:
    os.mkdir(odbOptimPath)
except:
    pass

odbName = odbSample

print("... 1) Reading ODB file " + odbPath + odbName + ".odb")

odb_file = odbPath + odbName + '.odb'
odb_file_upgraded = odbPath + odbName + '_new.odb'

# It might happen, that odbs simulated on the server must be upgraded before reading locally. The following code block
# updates older odb formats to the new 2021 release
if isUpgradeRequiredForOdb(upgradeRequiredOdbPath=odb_file):
    upgradeOdb(existingOdbPath=odb_file, upgradedOdbPath=odb_file_upgraded)
    odb = openOdb(odb_file_upgraded, readOnly=True)
    print("\n...converting odb file to 2021 version")
else:
    odb = openOdb(odb_file, readOnly=True)

# Define Step and Frame, Field Output and Element Set
# --------------------------------------------
# 1) Step (only one step: STEP-1)
step = odb.steps['Step-1']
print('... 2) Selected step: ', step.name)

# 2) Frame (last frame index -1)
first_frame = step.frames[0]
second_frame = step.frames[1]
target_frame = step.frames[target_frame_number]

# Field names of field outputs in specific frame of a step
# for fieldName in target_frame.fieldOutputs.keys():
#    print('available field names = ', fieldName)

# FieldOutput
# E tensor = [E11, E22, E33, E12, E13, E23]
elres_fieldVar = target_frame.fieldOutputs[fieldVar]
elres_COORD = first_frame.fieldOutputs['COORD']

elres_BVTVc = second_frame.fieldOutputs['SDV_BVTVc']
elres_BVTVt = second_frame.fieldOutputs['SDV_BVTVt']
elres_PBVc = second_frame.fieldOutputs['SDV_PBVc']
elres_PBVt = second_frame.fieldOutputs['SDV_PBVt']

# Select element set
# print 'Element sets BONE available =', 'BONE' in odb.rootAssembly.instances['PART-1-1'].elementSets.keys()

# Create subset with elset BONE and position CENTROID
# ---------------------------------------------------
# position = CENTROID requires *ELEMENT OUTPUT, POSITION=CENTROIDAL in input file. --> only values integrated to the
# centroid of each element
print('... 3) create subset with elset BONE')
elres_fieldVar_bone = elres_fieldVar.getSubset(region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'],
                                               elementType='C3D8', position=CENTROID)
elres_COORD_bone = elres_COORD.getSubset(region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'],
                                         elementType='C3D8', position=CENTROID)
elres_BVTV_cort = elres_BVTVc.getSubset(region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'],
                                           elementType='C3D8', position=CENTROID)
elres_BVTV_trab = elres_BVTVt.getSubset(region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'],
                                           elementType='C3D8', position=CENTROID)

elres_PBV_cort = elres_PBVc.getSubset(region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'],
                                         elementType='C3D8', position=CENTROID)
elres_PBV_trab = elres_PBVt.getSubset(region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'],
                                         elementType='C3D8', position=CENTROID)

# Store centroid coordinates and element values to csv file
# ---------------------------------------
print('... 4) store centroid coordinates and element values to _res.csv')
results = ['centroid_x', 'centroid_y', 'centroid_z',
           str(fieldVar) + '11', str(fieldVar) + '22',
           str(fieldVar) + '33', str(fieldVar) + '12',
           str(fieldVar) + '13', str(fieldVar) + '23', 'BVTVpbv', 'BVTVpbv_cort', 'BVTVpbv_trab']

results_mandel = ['centroid_x', 'centroid_y', 'centroid_z',
                  str(fieldVar) + '11', str(fieldVar) + '22',
                  str(fieldVar) + '33', str(fieldVar) + '12_mandel',
                  str(fieldVar) + '13_mandel', str(fieldVar) + '23_mandel', 'BVTVpbv', 'BVTVpbv_cort', 'BVTVpbv_trab']

for i in range(len(elres_COORD_bone.values)):
    c = elres_COORD_bone.values[i]
    e = elres_fieldVar_bone.values[i]
    bvtvpbv = elres_BVTV_cort.values[i] * elres_PBV_cort.values[i] + elres_BVTV_trab.values[i] * elres_PBV_trab.values[i]
    bvtv_cort = elres_BVTV_cort.values[i] * elres_PBV_cort.values[i]
    bvtv_trab = elres_BVTV_trab.values[i] * elres_PBV_trab.values[i]
    output = [c.data[0],
              c.data[1],
              c.data[2],
              e.data[0],
              e.data[1],
              e.data[2],
              e.data[3],
              e.data[4],
              e.data[5],
              bvtvpbv.data,
              bvtv_cort.data,
              bvtv_trab.data]
    output_mandel = [c.data[0],
               c.data[1],
               c.data[2],
               e.data[0],
               e.data[1],
               e.data[2],
               np.sqrt(2) / 2 * e.data[3],
               np.sqrt(2) / 2 * e.data[4],
               np.sqrt(2) / 2 * e.data[5],
               bvtvpbv.data,
               bvtv_cort.data,
               bvtv_trab.data]
    results.extend(output)
    results_mandel.extend(output_mandel)
print(output[11:])

# We have to cut the header again in order to convert the lists to numerical arrays. Otherwise they will be strings...
np_results = np.asarray(results[12:])
np_results_shape = np.reshape(np_results, (-1, 12))
np_results_mandel = np.asanyarray(results_mandel[12:])
np_results_mandel_shape = np.reshape(np_results_mandel, (-1, 12))

# for i in range(50):
#     print(np_results_mandel_shape[i])
outfile = odbOptimPath + "/" + odbName + '_' + fieldVar + '_res.csv'

with open(outfile, 'w') as f:
    writer = csv.writer(f)
    for row in np_results_shape:
        writer.writerow(row)

print('... 5) Create np array in shape of mesh with respective element values')
# Create 3d array with sorted elements
# ------------------------------------
# convert np_results to float values
centroid_x = np_results_shape[0:, 0].astype(np.float)
centroid_y = np_results_shape[0:, 1].astype(np.float)
centroid_z = np_results_shape[0:, 2].astype(np.float)

# Search for maximum coordinates in x,y,z for setting up array length
centroid_x_max = np.amax(centroid_x)
centroid_y_max = np.amax(centroid_y)
centroid_z_max = np.amax(centroid_z)
centroid_x_min = np.amin(centroid_x)
centroid_y_min = np.amin(centroid_y)
centroid_z_min = np.amin(centroid_z)

# Check if min coordinates are the same. If TRUE, compute real FE element size (this is skipped, as it caused some
# problems with specific samples, even the coordinates were exactly the same)
if round(centroid_x_min, 5) == round(centroid_y_min, 5) and round(round(centroid_z_min, 5) % round(centroid_x_min, 5),
                                                                  5) == 0:
    elsize_real = centroid_x_min * 2
else:
    print(round(centroid_z_min % centroid_x_min, 5))
    print("Centroid_z_min {} is not equal to multiple of centroid_x_min {}", centroid_z_min, centroid_x_min)
    elsize_real = centroid_x_min * 2

# Compute number of elements in each direction
n_elem_x = int(round((centroid_x_max + centroid_x_min) / elsize_real, 5))
n_elem_y = int(round((centroid_y_max + centroid_y_min) / elsize_real, 5))
n_elem_z = int(round((centroid_z_max + centroid_z_min) / elsize_real, 5))

print(n_elem_x, n_elem_y, n_elem_z)

# create 3d array with shape n_elem with zeros
np_3d_e11 = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_3d_e22 = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_3d_e33 = np.zeros((n_elem_x, n_elem_y, n_elem_z))

np_3d_e12 = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_3d_e13 = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_3d_e23 = np.zeros((n_elem_x, n_elem_y, n_elem_z))

np_3d_mandel = np.zeros((n_elem_x, n_elem_y, n_elem_z, 1, 6))

np_3d_bvtv = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_3d_bvtv_cort = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_3d_bvtv_trab = np.zeros((n_elem_x, n_elem_y, n_elem_z))
np_bone_elem = np.zeros((n_elem_x, n_elem_y, n_elem_z))

# Compute position of values in np_3d_results from centroid coordinates
np_pos = np_results_shape[1:, ].astype(np.float)

# Compute equivalent strain
eqstrain = (3 / 2 * (np_3d_e11 ** 2 + np_3d_e22 ** 2 + np_3d_e33 ** 2)) + (
        3 / 2 * (np_3d_e12 ** 2 + np_3d_e13 ** 2 + np_3d_e23 ** 2))
np_eq_vM_strain = 2 / 3 * np.sqrt(eqstrain)

# Compute positions from centroid coordinates
for i in range(np_pos.shape[0]):  # starts at 0
    # change centroid_x to position index (0:n_elem_x-1)
    np_pos[i, 0] = int(round((np_pos[i, 0] - (elsize_real / 2)) / elsize_real))
    np_pos[i, 1] = int(round((np_pos[i, 1] - (elsize_real / 2)) / elsize_real))
    np_pos[i, 2] = int(round((np_pos[i, 2] - (elsize_real / 2)) / elsize_real))

# Debug write out
# --------------------------------------------
outfile_np_pos = odbOptimPath + "/" + odbName + '_np_pos.txt'

with open(outfile_np_pos, 'w') as f:
    writer = csv.writer(f)
    for row in np_pos:
        writer.writerow(row)
# --------------------------------------------

for i in range(np_pos.shape[0]):
    # compute position index from centroid
    pos_x = int(np.round(np_pos[i, 0], 0))
    pos_y = int(np.round(np_pos[i, 1], 0))
    pos_z = int(np.round(np_pos[i, 2], 0))
    # pos_x = int(round((np_pos[i, 0] - (elsize_real / 2)) / elsize_real))
    # pos_y = int(round((np_pos[i, 1] - (elsize_real / 2)) / elsize_real))
    # pos_z = int(round((np_pos[i, 2] - (elsize_real / 2)) / elsize_real))
    # np_resuts_shape[i,3] = E11
    np_3d_e11[pos_x, pos_y, pos_z] = np_results_shape[i, 3]
    np_3d_e22[pos_x, pos_y, pos_z] = np_results_shape[i, 4]
    np_3d_e33[pos_x, pos_y, pos_z] = np_results_shape[i, 5]

    np_3d_e12[pos_x, pos_y, pos_z] = np_results_shape[i, 6]
    np_3d_e13[pos_x, pos_y, pos_z] = np_results_shape[i, 7]
    np_3d_e23[pos_x, pos_y, pos_z] = np_results_shape[i, 8]

    # Mandel notation
    np_3d_mandel[pos_x, pos_y, pos_z] = [np_results_mandel_shape[i, 3],
                                         np_results_mandel_shape[i, 4],
                                         np_results_mandel_shape[i, 5],
                                         np_results_mandel_shape[i, 8],
                                         np_results_mandel_shape[i, 7],
                                         np_results_mandel_shape[i, 6]]

    np_3d_bvtv[pos_x, pos_y, pos_z] = np_results_shape[i, 9]
    np_3d_bvtv_cort[pos_x, pos_y, pos_z] = np_results_shape[i, 10]
    np_3d_bvtv_trab[pos_x, pos_y, pos_z] = np_results_shape[i, 11]

np_bone_elem[np_3d_e11 != 0.0] = 1

# Save numpy array containing computed variables
np.save(odbOptimPath + "/" + odbName + '_bone_elem.npy', np_bone_elem)
np.save(odbOptimPath + "/" + odbName + '_bvtv_pbv.npy', np_3d_bvtv)
np.save(odbOptimPath + "/" + odbName + '_bvtv_cort.npy', np_3d_bvtv_cort)
np.save(odbOptimPath + "/" + odbName + '_bvtv_trab.npy', np_3d_bvtv_trab)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '11.npy', np_3d_e11)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '22.npy', np_3d_e22)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '33.npy', np_3d_e33)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '12.npy', np_3d_e12)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '13.npy', np_3d_e13)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '23.npy', np_3d_e23)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '_vect6.npy', np_3d_mandel)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '_mandel_shape.npy', np_results_mandel_shape)
np.save(odbOptimPath + "/" + odbName + '_' + fieldVar + '_equivstrain.npy', np_eq_vM_strain)

odb.close()