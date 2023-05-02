# ReadODB.py
# A script to read deformation gradient from odb file from ABAQUS.

import csv
import sys
import numpy as np
from odbAccess import *

#print 'Open odb file'

Odb = openOdb('PP.odb')

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
    for F in range(len(Frames)):

        sys.stdout.write("\r" + "Step " + str(S+1) + "/" + str(len(Steps)) + " Frame " + str(F+1) + "/" + str(len(Frames)) + " ")
        sys.stdout.flush()

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

        U = Frames[F].fieldOutputs['U']
        RF = Frames[F].fieldOutputs['RF']

        # Get elements deformation gradient
        ElementsDG = np.zeros((N,22))
        for iElement, Element in enumerate(Instance.elements[:N]):

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
            ElementsDG[iElement,13:] = F_IntegrationPoints
            DGList.append(ElementsDG)
        
        # Get nodes displacements/rotations and reaction forces/moments
        NodesURFM = np.zeros((len(Nodes),12))
        for iNode, Node in enumerate(Nodes):
            NodesURFM[iNode,0] = S
            NodesURFM[iNode,1] = F
            NodesURFM[iNode,2] = Node.label
            NodesURFM[iNode,3:6] = Node.coordinates
            NodesURFM[iNode,6:9] = U.getSubset(region=Node).values[0].data
            NodesURFM[iNode,9:] = RF.getSubset(region=Node).values[0].data
        
        UFList.append(NodesURFM)
        
sys.stdout.write("\n")
sys.stdout.flush()

with open('Elements.csv', 'w') as File:
    Writer = csv.writer(File)
    Writer.writerow(['Step','Frame','Node 1','Node 2','Node 3','Node 4',
                     'Node 5','Node 6','Node 7','Node 8', 'X', 'Y', 'Z',
                     'F11','F12','F13','F21','F22','F23','F31','F32','F33'])            
    for Frame in DGList:
        for Row in Frame:
            Writer.writerow(Row)

with open('Nodes.csv', 'w') as File:
    Writer = csv.writer(File)
    Writer.writerow(['Step','Frame','Label',
                     'X', 'Y', 'Z',
                     'Ux', 'Uy', 'Uz',
                     'Fx', 'Fy', 'Fz'])            
    for Frame in UFList:
        for Row in Frame:
            Writer.writerow(Row)

Odb.close()
