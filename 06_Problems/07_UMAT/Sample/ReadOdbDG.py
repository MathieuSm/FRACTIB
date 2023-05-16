# ReadODB.py
# A script to read deformation gradient from odb file from ABAQUS.

import csv
import sys
import numpy as np
from odbAccess import *

#print 'Open odb file'
NFrames = sys.argv[2]
Odb = openOdb(sys.argv[1])

# Create variable refering to model instance
Instances = Odb.rootAssembly.instances.keys()
Instance = Odb.rootAssembly.instances[Instances[0]]
N = len(Instance.elements)

# Get nodes labels
Nodes = Instance.nodes
NodeLabels = np.array([])
for iNode, Node in enumerate(Nodes):
    NodeLabels = np.append(NodeLabels,Node.label)

# Iterate over each step and frame
DGList = []
UFList = []
Steps = Odb.steps.keys()
for S, Step in enumerate(Steps):
    Frames = Odb.steps[Step].frames
    if len(Frames) > 0:
        for F in np.round(np.linspace(0,len(Frames)-1,NFrames)).astype(int):

            # Get fields output
            F11 = Frames[F].fieldOutputs['SDV_F11']
            F12 = Frames[F].fieldOutputs['SDV_F12']
            F13 = Frames[F].fieldOutputs['SDV_F13']
            F21 = Frames[F].fieldOutputs['SDV_F21']
            F22 = Frames[F].fieldOutputs['SDV_F22']
            F23 = Frames[F].fieldOutputs['SDV_F23']
            F31 = Frames[F].fieldOutputs['SDV_F31']
            F32 = Frames[F].fieldOutputs['SDV_F32']
            F33 = Frames[F].fieldOutputs['SDV_F33']
            DG = [F11,F12,F13,F21,F22,F23,F31,F32,F33]

            # Get elements deformation gradient
            ElementsDG = np.zeros((N,22))
            for iElement, Element in enumerate(Instance.elements[:N]):

                Line = "Step " + str(S+1) + "/" + str(len(Steps))
                Line +=  " Frame " + str(F+1) + "/" + str(len(Frames))
                Line +=  " Element " + str(iElement+1) + "/" + str(N)
                sys.stdout.write("\r" + Line + "          ")
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
                    F_Value = DG[F_Component].getSubset(region=Element).values[0].data
                    F_IntegrationPoints[0,F_Component] = F_Value

                # Add data to array
                ElementsDG[iElement,0] = S
                ElementsDG[iElement,1] = F
                ElementsDG[iElement,2:10] = Element.connectivity
                ElementsDG[iElement,10:13] = np.mean(XYZ,axis=0)
                for i in range(9):
                    ElementsDG[iElement,13+i] = DG[i].getSubset(region=Element).values[0].data

            DGList.append(ElementsDG)
            
sys.stdout.write("\n")
sys.stdout.flush()

with open(sys.argv[1][:-4] + '_DG.csv', 'w') as File:
    Writer = csv.writer(File)
    Writer.writerow(['Step','Frame', 'N1', 'N2','N3','N4', 'N5', 'N6','N7','N8', 'X', 'Y', 'Z',
                     'F11','F12','F13','F21','F22','F23','F31','F32','F33'])  
    for Frame in DGList:
        for Row in Frame:
            Writer.writerow(Row)

Odb.close()
