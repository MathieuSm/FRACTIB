# ReadODB.py
# A script to read deformation gradient from odb file from ABAQUS.

import csv
import sys
import numpy as np
from odbAccess import *

#print 'Open odb file'

Odb = openOdb(sys.argv[1])

# Create variable refering to model instance
Instances = Odb.rootAssembly.instances.keys()
Instance = Odb.rootAssembly.instances[Instances[0]]
N = len(Instance.elements)

# Get nodes labels
Node = Instance.nodes[-1]

# Iterate over each step and frame
UFList = []
Steps = Odb.steps.keys()
for S, Step in enumerate(Steps):
    Frames = Odb.steps[Step].frames
    for F in range(len(Frames)):

        sys.stdout.write("\r" + "Step " + str(S+1) + "/" + str(len(Steps)) + " Frame " + str(F+1) + "/" + str(len(Frames)) + "  ")
        sys.stdout.flush()

        # Get fields output
        U = Frames[F].fieldOutputs['U']
        RF = Frames[F].fieldOutputs['RF']
        
        # Get nodes displacements/rotations and reaction forces/moments
        NodesURFM = np.zeros(8)
        NodesURFM[0] = S
        NodesURFM[1] = F
        NodesURFM[2:5] = U.getSubset(region=Node).values[0].data
        NodesURFM[5:] = RF.getSubset(region=Node).values[0].data
        
        UFList.append(NodesURFM)
        
sys.stdout.write("\n")
sys.stdout.flush()

with open('RefNode.csv', 'w') as File:
    Writer = csv.writer(File)
    Writer.writerow(['Step','Frame',
                     'Ux', 'Uy', 'Uz',
                     'Fx', 'Fy', 'Fz'])            
    for Row in UFList:
        Writer.writerow(Row)

Odb.close()
