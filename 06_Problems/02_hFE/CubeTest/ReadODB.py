# odbRead.py
# A script to read the Abaqus/CAE Visualization module tutorial
# output database and read displacement data from the node at 
# the center of the hemispherical punch.

from odbAccess import *
import numpy as np

odb = openOdb(path='Cube_UMAT.odb')

# Create a variable that refers to the
# last frame of the first step.

LastFrame = odb.steps['Compression'].frames[-1]

# Create a variable that refers to the displacement 'U'
# in the last frame of the first step.

F11 = LastFrame.fieldOutputs['SDV_F11']
F22 = LastFrame.fieldOutputs['SDV_F22']
F33 = LastFrame.fieldOutputs['SDV_F33']

for Instance in odb.rootAssembly.instances:
    for Element in Instance.elements:
        for IntegrationPoint in range(8):
            print 'Element label: ', F11.getSubset(region=Element).values[IntegrationPoint].elementLabel
            print 'Integration point: ', F11.getSubset(region=Element).values[IntegrationPoint].integrationPoint
            print 'F_11: ', F11.getSubset(region=Element).values[IntegrationPoint].data
            print 'F_22: ', F22.getSubset(region=Element).values[IntegrationPoint].data
            print 'F_33: ', F33.getSubset(region=Element).values[IntegrationPoint].data

Nodes = odb.rootAssembly.instances['CUBEINSTANCE'].nodeSets['NODESET']

F11.values[0].instance.elements[0].connectivity


for Instance in odb.rootAssembly.instances:
    for Element in Instance.elements:
        for IntegrationPoint in range(8):
            F11.getSubset(region=Element).values[IntegrationPoint].data



odb.close()
