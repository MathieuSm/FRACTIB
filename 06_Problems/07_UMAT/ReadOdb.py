# ReadODB.py
# A script to read deformation gradient from odb file from ABAQUS.

import csv
import sys
import numpy as np
from odbAccess import *

#print 'Open odb file'

Odb = openOdb('PP.odb')

# Create variable that refers to the last frame of the step.
Steps = Odb.steps.keys()
Step = Steps[len(Steps)-1]

Frames = Odb.steps[Step].frames
Shift = 1
while len(Frames) == 0:
    Step = Steps[len(Steps)-(1+Shift)]
    Frames = Odb.steps[Step].frames
    Shift += 1
Frame = Odb.steps[Step].frames[len(Frames)-1]

# Create variable refering to model instance
Instances = Odb.rootAssembly.instances.keys()
Instance = Odb.rootAssembly.instances[Instances[0]]

F11 = Frame.fieldOutputs['SDV_F11']
F12 = Frame.fieldOutputs['SDV_F12']
F13 = Frame.fieldOutputs['SDV_F13']
F21 = Frame.fieldOutputs['SDV_F21']
F22 = Frame.fieldOutputs['SDV_F22']
F23 = Frame.fieldOutputs['SDV_F23']
F31 = Frame.fieldOutputs['SDV_F31']
F32 = Frame.fieldOutputs['SDV_F32']
F33 = Frame.fieldOutputs['SDV_F33']

F = [F11,F12,F13,F21,F22,F23,F31,F32,F33]
F_Names = ['F11','F12','F13','F21','F22','F23','F31','F32','F33']
Nodes = Instance.nodes
NodeLabels = np.array([])
for iNode, Node in enumerate(Nodes):
    NodeLabels = np.append(NodeLabels,Node.label)

# Initialize loop
N = len(Instance.elements)
ElementsDG = np.zeros((12,N))

for ElementNumber, Element in enumerate(Instance.elements[:N]):
    sys.stdout.write("\r" + "El. " + str(ElementNumber+1) + "/" + str(N))
    sys.stdout.flush()

    # Compute element central position
    XYZ = np.zeros((len(Element.connectivity),3))
    NodeNumber = 0
    for NodeLabel in Element.connectivity:
        Node = np.where(NodeLabels==NodeLabel)[0][0]
        for Axis in range(3):
            XYZ[NodeNumber,Axis] = Nodes[Node].coordinates[Axis]
        NodeNumber += 1

    # Get element mean deformation gradient
    F_IntegrationPoints = np.zeros((1,9))

    for F_Component in range(9):
        F_Value = F[F_Component].getSubset(region=Element).values[0].data
        F_IntegrationPoints[0,F_Component] = F_Value

    # Add data to arrays and increment
    ElementsDG[:3,ElementNumber] = np.mean(XYZ,axis=0)
    ElementsDG[3:,ElementNumber] = F_IntegrationPoints

with open('Elements_DG.csv', 'w') as File:
    Writer = csv.writer(File)
    for Row in ElementsDG.T:
        Writer.writerow(Row)            
Odb.close()
