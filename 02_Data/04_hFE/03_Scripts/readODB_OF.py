"""
Read ODB OF values
------------------
This script reads the specified ODB files and extracts OF values as an array of
shape (number of increments, number of elements). The respective array is then stored in
odbOptimPath as OF_value_array.npy.
"""

# Imports
import sys
from odbAccess import *
from abaqusConstants import *
# import session
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

# Create optimization folder if not yet existing
odbOptimPath = odbPath + "/" + "Optimization"
try:
    os.mkdir(odbOptimPath)
except:
    pass

# Read ODB and update if necessary
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


# Read OF values from ODBsample
# --------------------------------------------
# 1) Step (only one step: STEP-1)
step = odb.steps['Step-1']
print('... 2) Selected step: ', step.name)

full_OF = []
mean_OF = []
sd_OF = []

for i, frame in enumerate(step.frames):
    OF_field_bone = frame.fieldOutputs['SDV_OFvalue'].getSubset(
        region=odb.rootAssembly.instances['PART-1-1'].elementSets['BONE'], elementType='C3D8', position=CENTROID)

    # Read OF values of specified frame
    OF_values = []
    for i in range(len(OF_field_bone.values)):
        OF_values.append(OF_field_bone.values[i].data)
    OF_values = np.asarray(OF_values)

    # Store mean OF value
    full_OF.append(OF_values)
    mean_OF.append(OF_values.mean())
    sd_OF.append(OF_values.std())

full_OF = np.asarray(full_OF)
sd_OF = np.asarray(sd_OF)
mean_OF = np.asarray(mean_OF)

print(full_OF.shape)

# Save numpy array containing computed variables
np.save(odbOptimPath + "/" + odbName + '_OF_values.npy', full_OF)

odb.close()
