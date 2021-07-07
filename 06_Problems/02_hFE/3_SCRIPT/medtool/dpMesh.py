# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/misc/dpMesh.py
# Compiled at: 2015-04-09 20:17:41
"""
#########################################################################
File   : dpMesh.py
Author : D.H.Pahr
Date   : 6.09.2007
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
generates  a triangulated mesh and supports a library for many mesh 
related functions. 

Usage: python mesh.py ::
                -power    densityParameter
                -out      OutfileName   
                [-help]   (optional)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REQUIRED Arguments::
 -power : Parameter for the computation of the number of directions   
          on a unit sphere (=8*4^power). The direction algorithm projects 
          triangles on a unit sphere surface and refines this mesh. 
          The refinement is done by dividing each triangle into four 
          triangles. Default = 2 which gives in sum 128 triangles.
 -out  :  Output file: file where stiffness information is written 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OPTIONAL Arguments::
 -help :  Print usage
 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
from sys import argv,exit, stdout
from time import *
from string import split, replace, atof
import fec
import dpFem
import dpTensor
import dpMesh
import numpy
import numpy.linalg as linalg
import dpUtils

class dpMesh():

    def __init__(self, nodes=None, elems=None, modelName='default'):
        self.modelName = modelName
        self.nodes = nodes
        self.elems = elems
        self.thick = {}
        self.normals = {}
        self.skipNodes = {}
        self.phyStartPos = {}
        self.phyEndPos = {}

    def getElems(self):
        return self.elems

    def getNodes(self):
        return self.nodes

    def getThick(self):
        return self.thick

    def getNormals(self):
        return self.normals

    def setElems(self, elems):
        self.elems = elems

    def setNodes(self, nodes):
        self.nodes = nodes

    def setThick(self, thick):
        self.thick = thick

    def getThick(self, correctFkt=None):
        if correctFkt != None:
            for noid in self.nodes:
                var = {}
                var['thickOld'] = self.thick[noid]
                execfile(correctFkt, var, var)
                self.thick[noid] = var['thickNew']

        return self.thick

    def getElementsAroundNode(self):
        nodes = self.nodes
        elems = self.elems
        nodeList = {}
        for node in nodes:
            nodeList[node] = []

        for elemId in elems:
            elemNodes = elems[elemId].get_nodes()
            for elNode in elemNodes:
                nodeList[elNode].append(elemId)

        return nodeList

    def getElementsAroundElement(self, elid):
        noList = self.elems[elid].get_nodes()
        n1 = noList[0]
        n2 = noList[1]
        n3 = noList[2]
        elemsAroundNode = self.getElementsAroundNode()
        curElemList1 = elemsAroundNode[n1]
        curElemList2 = elemsAroundNode[n2]
        curElemList3 = elemsAroundNode[n3]
        curElemListH = curElemList1 + curElemList2 + curElemList3
        curElemList = []
        check = {}
        for el in curElemListH:
            check[el] = False

        for el in curElemListH:
            if curElemListH.count(el) == 1:
                curElemList.append(el)
            elif check[el] == False:
                curElemList.append(el)
                check[el] = True

        return curElemList

    def getNodesAroundNode(self):
        nodes = self.nodes
        elems = self.elems
        elemsAroundNode = self.getElementsAroundNode()
        nodeList = {}
        for node in nodes:
            nodeList[node] = []
            elemList = elemsAroundNode[node]
            for elemId in elemList:
                elemNodes = elems[elemId].get_nodes()
                for elNode in elemNodes:
                    if nodeList[node].count(elNode) == 0 and elNode != node:
                        nodeList[node].append(elNode)

        return nodeList

    def sortElementsAroundNode(self):
        elems = self.elems
        elemsAroundNode = self.getElementsAroundNode()
        newList = {}
        for node in elemsAroundNode:
            newElementOrder = []
            consideredNodes = [node]
            newList[node] = []
            elemList = elemsAroundNode[node]
            curElemList = elemsAroundNode[node]
            foundNeighbourElem = 0
            for k in range(len(curElemList) - 1):
                if k == 0:
                    curElemId = curElemList.pop(0)
                    newElementOrder.append(curElemId)
                else:
                    curElemId = foundNeighbourElem
                    curElemList.remove(curElemId)
                curElemIdNodes = elems[curElemId].get_nodes()
                for secondNode in curElemIdNodes:
                    if consideredNodes.count(secondNode) == 0:
                        break

                curPairNodes = [
                 node, secondNode]
                consideredNodes.append(secondNode)
                foundNeighbourElem = self.getNeighborElement(curPairNodes, curElemId, curElemList)
                newElementOrder.append(foundNeighbourElem)

            newList[node] = newElementOrder

        return newList

    def getNeighborElement(self, nodePair, firstElem, elemList=None):
        if elemList == None:
            elemsAroundNode = self.getElementsAroundNode()
            elemList = elemsAroundNode[nodePair[0]]
        elems = self.elems
        neighbour = None
        found = False
        for elemId in elemList:
            nodesOfElement = elems[elemId].get_nodes()
            if nodesOfElement.count(nodePair[0]) == 1 and nodesOfElement.count(nodePair[1]) == 1:
                if elemId != firstElem:
                    neighbour = elemId
                    found = True
                    break

        if not found:
            return False
        else:
            return neighbour
            return

    def getCOG(self):
        nodes = self.nodes
        elems = self.elems
        cog = {}
        for elemId in elems:
            nL = elems[elemId].get_nodes()
            if len(nL) == 3:
                p1 = numpy.array(nodes[nL[0]].get_coord())
                p2 = numpy.array(nodes[nL[1]].get_coord())
                p3 = numpy.array(nodes[nL[2]].get_coord())
                xs = (p1[0] + p2[0] + p3[0]) / 3.0
                ys = (p1[1] + p2[1] + p3[1]) / 3.0
                zs = (p1[2] + p2[2] + p3[2]) / 3.0
            elif len(nL) == 4:
                p1 = numpy.array(nodes[nL[0]].get_coord())
                p2 = numpy.array(nodes[nL[1]].get_coord())
                p3 = numpy.array(nodes[nL[2]].get_coord())
                p4 = numpy.array(nodes[nL[3]].get_coord())
                xs1 = (p1[0] + p2[0] + p3[0]) / 3.0
                ys1 = (p1[1] + p2[1] + p3[1]) / 3.0
                zs1 = (p1[2] + p2[2] + p3[2]) / 3.0
                xs2 = (p1[0] + p4[0] + p3[0]) / 3.0
                ys2 = (p1[1] + p4[1] + p3[1]) / 3.0
                zs2 = (p1[2] + p4[2] + p3[2]) / 3.0
                xs = (xs1 + xs2) / 2.0
                ys = (ys1 + ys2) / 2.0
                zs = (zs1 + zs2) / 2.0
            else:
                stdout.write('\n **ERROR** getCOG(): Exactly three nodes required!\n\n')
                stdout.flush()
            mul = 1.0
            cog[elemId] = (xs * mul, ys * mul, zs * mul)

        return cog

    def createDualSimplexMesh(self):
        nodes = {}
        elems = {}
        elemInfo = myMesh.sortElementsAroundNode()
        nodeInfo = myMesh.getCOG()
        for elemId in elemInfo:
            elems[elemId] = []
            curElem = dpFem.element(elemId, elemInfo[elemId], 'simplex6')
            elems[elemId] = curElem

        for nodeId in nodeInfo:
            x = float(nodeInfo[nodeId][0])
            y = float(nodeInfo[nodeId][1])
            z = float(nodeInfo[nodeId][2])
            curNode = dpFem.node(nodeId, x, y, z)
            nodes[nodeId] = curNode

        return (
         nodes, elems)

    def createIcosahedronMesh(self):
        nodes = {}
        elems = {}
        r = 2.0 / 5.0 * numpy.sqrt(5.0)
        hB = -1.0 / 5.0 * numpy.sqrt(5.0)
        hT = 1.0 / 5.0 * numpy.sqrt(5.0)
        pi = numpy.pi
        nB = dpFem.node(1, 0.0, 0.0, -1.0)
        nodes[1] = nB
        nT = dpFem.node(12, 0.0, 0.0, 1.0)
        nodes[12] = nT
        firstRow = []
        nodeId = 1
        for i in range(5):
            alpha = i * 2.0 * pi / 5.0
            x = r * numpy.cos(alpha)
            y = r * numpy.sin(alpha)
            z = hB
            curNode = dpFem.node(nodeId, x, y, z)
            nodeId += 1
            nodes[nodeId] = curNode
            firstRow.append(nodeId)

        secondRow = []
        for i in range(5):
            alpha = i * 2.0 * pi / 5.0 + pi / 5.0
            x = r * numpy.cos(alpha)
            y = r * numpy.sin(alpha)
            z = hT
            curNode = dpFem.node(nodeId, x, y, z)
            nodeId += 1
            nodes[nodeId] = curNode
            secondRow.append(nodeId)

        elemId = 0
        for nodeId in firstRow:
            elemId += 1
            elems[elemId] = []
            if nodeId < 6:
                nextNodeId = nodeId + 1
            else:
                nextNodeId = 2
            curElem = dpFem.element(elemId, [1, nextNodeId, nodeId], 'tria3')
            elems[elemId] = curElem

        for nodeId in secondRow:
            elemId += 1
            elems[elemId] = []
            if nodeId < 11:
                nextNodeId = nodeId + 1
            else:
                nextNodeId = 7
            curElem = dpFem.element(elemId, [12, nodeId, nextNodeId], 'tria3')
            elems[elemId] = curElem

        for nodeId in firstRow:
            upperId = nodeId + 5
            elemId += 1
            elems[elemId] = []
            if nodeId < 6:
                nextNodeId = nodeId + 1
            else:
                nextNodeId = 2
            curElem = dpFem.element(elemId, [nodeId, nextNodeId, upperId], 'tria3')
            elems[elemId] = curElem

        for nodeId in secondRow:
            if nodeId < 11:
                lowerId = nodeId - 4
            else:
                lowerId = 2
            elemId += 1
            elems[elemId] = []
            if nodeId < 11:
                nextNodeId = nodeId + 1
            else:
                nextNodeId = 7
            curElem = dpFem.element(elemId, [nodeId, lowerId, nextNodeId], 'tria3')
            elems[elemId] = curElem

        self.elems = elems
        self.nodes = nodes

    def refineTriaMesh(self):
        nodes = self.nodes
        elems = self.elems
        newElems = {}
        newNodes = self.nodes
        newNode = {}
        noidMax = 0
        for noid in nodes:
            if noid > noidMax:
                noidMax = noid

        newNoid = noidMax
        elemId = 0
        for elid in elems:
            noList = elems[elid].get_nodes()
            if len(noList) > 3:
                stdout.write("\n **ERROR** dpMesh.refineTriaMesh(): Elements have '%s' nodes! \n\n" % repr(len(noList)))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            n1 = noList[0]
            n2 = noList[1]
            n3 = noList[2]
            n1x, n1y, n1z = nodes[n1].get_coord()
            n2x, n2y, n2z = nodes[n2].get_coord()
            n3x, n3y, n3z = nodes[n3].get_coord()
            tria = ((n1x, n1y, n1z), (n2x, n2y, n2z), (n3x, n3y, n3z))
            p4, p5, p6 = self.splitTriangle(tria)
            n4 = 0
            n5 = 0
            n6 = 0
            if newNode.has_key((n1, n2)) == False:
                newNoid += 1
                newNode[n1, n2] = newNoid
                newNode[n2, n1] = newNoid
                n4 = newNoid
                curNode = dpFem.node(newNoid, p4[0], p4[1], p4[2])
                newNodes[newNoid] = curNode
            else:
                n4 = newNode[n1, n2]
            if newNode.has_key((n2, n3)) == False:
                newNoid += 1
                newNode[n2, n3] = newNoid
                newNode[n3, n2] = newNoid
                n5 = newNoid
                curNode = dpFem.node(newNoid, p5[0], p5[1], p5[2])
                newNodes[newNoid] = curNode
            else:
                n5 = newNode[n2, n3]
            if newNode.has_key((n3, n1)) == False:
                newNoid += 1
                newNode[n3, n1] = newNoid
                newNode[n1, n3] = newNoid
                n6 = newNoid
                curNode = dpFem.node(newNoid, p6[0], p6[1], p6[2])
                newNodes[newNoid] = curNode
            else:
                n6 = newNode[n3, n1]
            elemNodeList = [
             [
              n1, n4, n6], [n4, n2, n5], [n4, n5, n6], [n5, n3, n6]]
            for triaNodes in elemNodeList:
                elemId += 1
                curElem = dpFem.element(elemId, triaNodes, 'tria3')
                newElems[elemId] = curElem

        self.nodes = newNodes
        self.elems = newElems

    def refine2DMeshNormalThickness(self):
        nodes = self.nodes
        elems = self.elems
        newElems = {}
        newNodes = self.nodes
        newNode = {}
        newThick = self.thick
        newNorm = self.normals
        nodePairs = {}
        noidMax = 0
        for noid in nodes:
            if noid > noidMax:
                noidMax = noid
            self.skipNodes[noid] = True

        newNoid = noidMax
        elemId = 0
        for elid in elems:
            noList = elems[elid].get_nodes()
            if len(noList) == 3:
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                n1x, n1y, n1z = nodes[n1].get_coord()
                n2x, n2y, n2z = nodes[n2].get_coord()
                n3x, n3y, n3z = nodes[n3].get_coord()
                tria = ((n1x, n1y, n1z), (n2x, n2y, n2z), (n3x, n3y, n3z))
                p4, p5, p6 = self.splitTriangle(tria)
                n4 = 0
                n5 = 0
                n6 = 0
                if newNode.has_key((n1, n2)) == False:
                    newNoid += 1
                    newNode[n1, n2] = newNoid
                    newNode[n2, n1] = newNoid
                    nodePairs[newNoid] = (
                     n1, n2)
                    n4 = newNoid
                    curNode = dpFem.node(newNoid, p4[0], p4[1], p4[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n1] + self.normals[n2]) / 2.0
                    if self.thick[n1] > 0.0 and self.thick[n2] > 0.0:
                        newThick[newNoid] = (self.thick[n1] + self.thick[n2]) / 2
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n4 = newNode[n1, n2]
                if newNode.has_key((n2, n3)) == False:
                    newNoid += 1
                    newNode[n2, n3] = newNoid
                    newNode[n3, n2] = newNoid
                    nodePairs[newNoid] = (
                     n2, n3)
                    n5 = newNoid
                    curNode = dpFem.node(newNoid, p5[0], p5[1], p5[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n3] + self.normals[n2]) / 2.0
                    if self.thick[n3] > 0.0 and self.thick[n2] > 0.0:
                        newThick[newNoid] = (self.thick[n3] + self.thick[n2]) / 2.0
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n5 = newNode[n2, n3]
                if newNode.has_key((n3, n1)) == False:
                    newNoid += 1
                    newNode[n3, n1] = newNoid
                    newNode[n1, n3] = newNoid
                    nodePairs[newNoid] = (
                     n3, n1)
                    n6 = newNoid
                    curNode = dpFem.node(newNoid, p6[0], p6[1], p6[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n3] + self.normals[n1]) / 2.0
                    if self.thick[n3] > 0.0 and self.thick[n1] > 0.0:
                        newThick[newNoid] = (self.thick[n3] + self.thick[n1]) / 2.0
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n6 = newNode[n3, n1]
                elemNodeList = [[n1, n4, n6], [n4, n2, n5], [n4, n5, n6], [n5, n3, n6]]
                for triaNodes in elemNodeList:
                    elemId += 1
                    curElem = dpFem.element(elemId, triaNodes, 'tria3')
                    newElems[elemId] = curElem

            elif len(noList) == 4:
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                n4 = noList[3]
                n1x, n1y, n1z = nodes[n1].get_coord()
                n2x, n2y, n2z = nodes[n2].get_coord()
                n3x, n3y, n3z = nodes[n3].get_coord()
                n4x, n4y, n4z = nodes[n4].get_coord()
                tria1 = ((n1x, n1y, n1z), (n2x, n2y, n2z), (n3x, n3y, n3z))
                tria2 = ((n1x, n1y, n1z), (n3x, n3y, n3z), (n4x, n4y, n4z))
                p5, p6, du = self.splitTriangle(tria1)
                du, p7, p8 = self.splitTriangle(tria2)
                n5 = 0
                n6 = 0
                n7 = 0
                n8 = 0
                n9 = 0
                if newNode.has_key((n1, n2)) == False:
                    newNoid += 1
                    newNode[n1, n2] = newNoid
                    newNode[n2, n1] = newNoid
                    nodePairs[newNoid] = (
                     n1, n2)
                    n5 = newNoid
                    curNode = dpFem.node(newNoid, p5[0], p5[1], p5[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n1] + self.normals[n2]) / 2.0
                    if self.thick[n1] > 0.0 and self.thick[n2] > 0.0:
                        newThick[newNoid] = (self.thick[n1] + self.thick[n2]) / 2.0
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n5 = newNode[n1, n2]
                if newNode.has_key((n2, n3)) == False:
                    newNoid += 1
                    newNode[n2, n3] = newNoid
                    newNode[n3, n2] = newNoid
                    nodePairs[newNoid] = (
                     n2, n3)
                    n6 = newNoid
                    curNode = dpFem.node(newNoid, p6[0], p6[1], p6[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n3] + self.normals[n2]) / 2.0
                    if self.thick[n3] > 0.0 and self.thick[n2] > 0.0:
                        newThick[newNoid] = (self.thick[n3] + self.thick[n2]) / 2.0
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n6 = newNode[n2, n3]
                if newNode.has_key((n3, n4)) == False:
                    newNoid += 1
                    newNode[n3, n4] = newNoid
                    newNode[n4, n3] = newNoid
                    nodePairs[newNoid] = (
                     n3, n4)
                    n7 = newNoid
                    curNode = dpFem.node(newNoid, p7[0], p7[1], p7[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n3] + self.normals[n4]) / 2.0
                    if self.thick[n3] > 0.0 and self.thick[n4] > 0.0:
                        newThick[newNoid] = (self.thick[n3] + self.thick[n4]) / 2.0
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n7 = newNode[n3, n4]
                if newNode.has_key((n4, n1)) == False:
                    newNoid += 1
                    newNode[n4, n1] = newNoid
                    newNode[n1, n4] = newNoid
                    nodePairs[newNoid] = (
                     n4, n1)
                    n8 = newNoid
                    curNode = dpFem.node(newNoid, p8[0], p8[1], p8[2])
                    newNodes[newNoid] = curNode
                    newNorm[newNoid] = (self.normals[n4] + self.normals[n1]) / 2.0
                    if self.thick[n4] > 0.0 and self.thick[n1] > 0.0:
                        newThick[newNoid] = (self.thick[n4] + self.thick[n1]) / 2.0
                    else:
                        newThick[newNoid] = 0.0
                        self.skipNodes[newNoid] = True
                else:
                    n8 = newNode[n4, n1]
                newNoid += 1
                n9 = newNoid
                x = (p5[0] + p7[0]) / 2.0
                y = (p5[1] + p7[1]) / 2.0
                z = (p5[2] + p7[2]) / 2.0
                curNode = dpFem.node(newNoid, x, y, z)
                newNodes[newNoid] = curNode
                newNorm[newNoid] = (self.normals[n1] + self.normals[n2] + self.normals[n3] + self.normals[n4]) / 4.0
                if self.thick[n1] > 0.0 and self.thick[n2] > 0.0 and self.thick[n3] > 0.0 and self.thick[n4] > 0.0:
                    newThick[newNoid] = self.thick[n1] + self.thick[n2] + self.thick[n3] + self.thick[n4] / 4.0
                else:
                    newThick[newNoid] = 0.0
                    self.skipNodes[newNoid] = True
                elemNodeList = [
                 [
                  n1, n5, n9, n8], [n5, n2, n6, n9], [n9, n6, n3, n7], [n8, n9, n7, n4]]
                for triaNodes in elemNodeList:
                    elemId += 1
                    curElem = dpFem.element(elemId, triaNodes, 'quad4')
                    newElems[elemId] = curElem

            else:
                stdout.write("\n **ERROR** dpMesh.refine2DMeshNormalThickness(): Elements have '%s' nodes! \n\n" % repr(len(noList)))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)

        self.nodes = newNodes
        self.elems = newElems
        self.thick = newThick
        self.normals = newNorm
        return nodePairs

    def checkBoundaryEdgeMidnodes(self, boundTol, nlx, nly, nlz):
        elemsAroundNode = self.getElementsAroundNode()
        if boundTol > 0.0:
            for noid in self.nodes:
                curElemList = elemsAroundNode[noid]
                elCount = 0
                for elid in curElemList:
                    pt = self.computeCOG(elid)
                    if abs(pt[0]) < boundTol or abs(pt[0] - nlx) < boundTol or abs(pt[1]) < boundTol or abs(pt[1] - nly) < boundTol or abs(pt[2]) < boundTol or abs(pt[2] - nlz) < boundTol:
                        elCount += 1

                if elCount == len(curElemList) and self.thick[noid] > 0.0:
                    self.skipNodes[noid] = True
                    self.thick[noid] = 0.0

    def checkBoundaryEdgeMidnodes2(self, boundTol, nlx, nly, nlz, oldNodePos=None):
        stdout.write(' ... check boundary elements    \n')
        stdout.flush()
        elemsAroundNode = self.getElementsAroundNode()
        faces = [[], [], [], [], [], []]
        faceN = [
         numpy.array([-1, 0, 0]), numpy.array([0, -1, 0]), numpy.array([0, 0, -1]), numpy.array([1, 0, 0]), numpy.array([0, 1, 0]), numpy.array([0, 0, 1])]
        err = {}
        eCount = 0
        if boundTol > 0.0:
            for elid in self.elems:
                err[elid] = 0
                pt = self.computeCOG(elid)
                if abs(pt[0]) < boundTol:
                    faces[0].append(elid)
                if abs(pt[1]) < boundTol:
                    faces[1].append(elid)
                if abs(pt[2]) < boundTol:
                    faces[2].append(elid)
                if abs(pt[0] - nlx) < boundTol:
                    faces[3].append(elid)
                if abs(pt[1] - nly) < boundTol:
                    faces[4].append(elid)
                if abs(pt[2] - nlz) < boundTol:
                    faces[5].append(elid)

            for faceId in range(6):
                fNorm = faceN[faceId]
                for elid in faces[faceId]:
                    eNorm = self.getNormalVector(elid)
                    if numpy.dot(fNorm, eNorm) < 0.0:
                        eCount += 1
                        err[elid] = 1

        if eCount > 0:
            stdout.write('     !! inflated elements    : %i !!\n' % eCount)
            stdout.flush()
            fecModel = fec.fec()
            fecModel.write('test-inflated-elems.case', 'waring', self.nodes, None, self.elems, None, EscaResults=[err])
            if oldNodePos:
                for elid in err:
                    if err[elid] == 1:
                        nList = self.elems[elid].get_nodes()
                        for noid in nList:
                            x = oldNodePos[noid][0]
                            y = oldNodePos[noid][1]
                            z = oldNodePos[noid][2]
                            self.nodes[noid].set_coord(x, y, z)

        return

    def refineTriaMeshLocal(self, meanDist, rho, echo=False):
        if echo == True:
            stdout.flush()
            stdout.write('\n ... local refine: meltT/meltE/flip/divide ')
            stdout.flush()
        elemIds = []
        for elid in self.elems:
            elemIds.append(elid)

        checkedElems = {}
        countMelt = 0
        countMelt2 = 0
        for elid in elemIds:
            if not checkedElems.has_key(elid):
                elemsAroundNode = self.getElementsAroundNode()
                noList = self.elems[elid].get_nodes()
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                le1 = self.computeDistance(n1, n2)
                le2 = self.computeDistance(n3, n2)
                le3 = self.computeDistance(n1, n3)
                nodePairDiv = []
                edgeLength = []
                oppNode = []
                if le1 < meanDist / rho:
                    nodePairDiv.append((n1, n2))
                    edgeLength.append(le1)
                    oppNode.append(n3)
                if le2 < meanDist / rho:
                    nodePairDiv.append((n3, n2))
                    edgeLength.append(le2)
                    oppNode.append(n1)
                if le3 < meanDist / rho:
                    nodePairDiv.append((n1, n3))
                    edgeLength.append(le3)
                    oppNode.append(n2)
                if len(nodePairDiv) == 3:
                    countMelt2 += 1
                    xs, ys, zs = self.computeCOG(elid)
                    nodeList = self.elems[elid].get_nodes()
                    n1 = nodeList[0]
                    curNode = self.nodes[n1]
                    curNode.set_x(xs)
                    curNode.set_y(ys)
                    curNode.set_z(zs)
                    curElemList = self.getElementsAroundElement(elid)
                    nEl1 = self.getNeighborElement([n1, n2], elid, curElemList)
                    nEl2 = self.getNeighborElement([n1, n3], elid, curElemList)
                    nEl3 = self.getNeighborElement([n2, n3], elid, curElemList)
                    del self.elems[elid]
                    del self.elems[nEl1]
                    del self.elems[nEl2]
                    del self.elems[nEl3]
                    del self.nodes[n2]
                    del self.nodes[n3]
                    curElemList.remove(elid)
                    curElemList.remove(nEl1)
                    curElemList.remove(nEl2)
                    curElemList.remove(nEl3)
                    for curElid in curElemList:
                        nodeList = self.elems[curElid].get_nodes()
                        newNodeList = []
                        for curNoid in nodeList:
                            if curNoid != n2 and curNoid != n3:
                                newNodeList.append(curNoid)
                            else:
                                newNodeList.append(n1)

                        self.elems[curElid].set_nodes(newNodeList)
                        self.elems[curElid].set_type('tria3')

                    checkedElems[nEl1] = True
                    checkedElems[nEl2] = True
                    checkedElems[nEl3] = True
                if len(nodePairDiv) > 1 and len(nodePairDiv) < 3:
                    i = -1
                    if len(nodePairDiv) == 1:
                        i = 0
                    elif len(nodePairDiv) == 2:
                        if edgeLength[0] < edgeLength[1]:
                            i = 0
                        else:
                            i = 1
                    countMelt += 1
                    n3a = nodePairDiv[i][0]
                    n3b = nodePairDiv[i][1]
                    nD = oppNode[i]
                    curElemList = elemsAroundNode[n3a]
                    curElemList2 = elemsAroundNode[n3b]
                    nElid = self.getNeighborElement([n3a, n3b], elid, curElemList)
                    noList = self.elems[nElid].get_nodes()
                    for node in noList:
                        if node != n3a and node != n3b:
                            nC = node
                            break

                    if len(elemsAroundNode[nD]) > 3 and len(elemsAroundNode[nC]) > 3 and len(curElemList) > 3 and len(curElemList2) > 3:
                        n1x, n1y, n1z = self.nodes[n3a].get_coord()
                        n2x, n2y, n2z = self.nodes[n3b].get_coord()
                        curNode = self.nodes[n3a]
                        curNode.set_x((n1x + n2x) / 2.0)
                        curNode.set_y((n1y + n2y) / 2.0)
                        curNode.set_z((n1z + n2z) / 2.0)
                        for curElid in curElemList2:
                            nodeList = self.elems[curElid].get_nodes()
                            newNodeList = []
                            for curNoid in nodeList:
                                if curNoid != n3b:
                                    newNodeList.append(curNoid)
                                else:
                                    newNodeList.append(n3a)

                            self.elems[curElid].set_nodes(newNodeList)
                            self.elems[curElid].set_type('tria3')

                        del self.nodes[n3b]
                        del self.elems[elid]
                        del self.elems[nElid]
                    checkedElems[nElid] = True

        if echo == True:
            stdout.flush()
            stdout.write(' %s/%s' % (repr(countMelt2), repr(countMelt)))
            stdout.flush()
        elemIds = []
        for elid in self.elems:
            elemIds.append(elid)

        checkedElems = {}
        countFlip = 0
        for elid in elemIds:
            if not checkedElems.has_key(elid):
                elemsAroundNode = self.getElementsAroundNode()
                noList = self.elems[elid].get_nodes()
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                le1 = self.computeDistance(n1, n2)
                le2 = self.computeDistance(n3, n2)
                le3 = self.computeDistance(n1, n3)
                if le1 > le2 and le1 > le3:
                    nB = n1
                    nC = n2
                    nA = n3
                    lAB = le2
                    lAC = le3
                if le2 > le1 and le2 > le3:
                    nB = n3
                    nC = n2
                    nA = n1
                    lAB = le1
                    lAC = le3
                if le3 > le1 and le3 > le2:
                    nB = n1
                    nC = n3
                    nA = n2
                    lAB = le2
                    lAC = le1
                curElemList = elemsAroundNode[nB]
                nElid = self.getNeighborElement([nB, nC], elid, curElemList)
                noList = self.elems[nElid].get_nodes()
                for nd in noList:
                    if nd != nB and nd != nC:
                        nD = nd
                        break

                lBC = self.computeDistance(nB, nC)
                lDC = self.computeDistance(nD, nC)
                lBD = self.computeDistance(nB, nD)
                lAD = self.computeDistance(nA, nD)
                if lBC > meanDist * rho and lDC < meanDist * rho and lBD < meanDist * rho and lAB < meanDist * rho and lAC < meanDist * rho and lBC > lAD:
                    self.elems[elid].set_nodes([nB, nA, nD])
                    self.elems[nElid].set_nodes([nA, nD, nC])
                    countFlip += 1
                checkedElems[nElid] = True

        if echo == True:
            stdout.flush()
            stdout.write('/%s' % repr(countFlip))
            stdout.flush()
        nodeIds = []
        elemIds = []
        for elid in self.elems:
            elemIds.append(elid)

        for nid in self.nodes:
            nodeIds.append(nid)

        maxNodeId = max(nodeIds)
        maxElemId = max(elemIds)
        checkedElems = {}
        countSubDiv = 0
        for elid in elemIds:
            if not checkedElems.has_key(elid):
                elemsAroundNode = self.getElementsAroundNode()
                noList = self.elems[elid].get_nodes()
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                le1 = self.computeDistance(n1, n2)
                le2 = self.computeDistance(n3, n2)
                le3 = self.computeDistance(n1, n3)
                nodePairDiv = []
                oppNode = []
                edgeLength = []
                if le1 > meanDist * rho:
                    nodePairDiv.append((n1, n2))
                    oppNode.append(n3)
                    edgeLength.append(le1)
                if le2 > meanDist * rho:
                    nodePairDiv.append((n3, n2))
                    oppNode.append(n1)
                    edgeLength.append(le2)
                if le3 > meanDist * rho:
                    nodePairDiv.append((n1, n3))
                    oppNode.append(n2)
                    edgeLength.append(le3)
                if len(nodePairDiv) > 1:
                    i = -1
                    if len(nodePairDiv) == 1:
                        i = 0
                    elif len(nodePairDiv) == 2:
                        if edgeLength[0] > edgeLength[1]:
                            i = 0
                        else:
                            i = 1
                    elif len(nodePairDiv) == 3:
                        if edgeLength[0] > edgeLength[1] and edgeLength[0] > edgeLength[2]:
                            i = 0
                        if edgeLength[1] > edgeLength[0] and edgeLength[1] > edgeLength[2]:
                            i = 1
                        if edgeLength[2] > edgeLength[0] and edgeLength[2] > edgeLength[1]:
                            i = 2
                    countSubDiv += 1
                    n3b = 0
                    n3a = oppNode[i]
                    n1a = nodePairDiv[i][0]
                    n2a = nodePairDiv[i][1]
                    curElemList = elemsAroundNode[n1a]
                    nElid = self.getNeighborElement([n1a, n2a], elid, curElemList)
                    noList = self.elems[nElid].get_nodes()
                    for node in noList:
                        if node != n1a:
                            if node != n2a:
                                n3b = node
                                break

                    le1b = self.computeDistance(n1a, n2a)
                    le2b = self.computeDistance(n3b, n2a)
                    le3b = self.computeDistance(n1a, n3b)
                    if le1b > edgeLength[i] or le2b > edgeLength[i] or le3b > edgeLength[i]:
                        pass
                    else:
                        n1x, n1y, n1z = self.nodes[n1a].get_coord()
                        n2x, n2y, n2z = self.nodes[n2a].get_coord()
                        maxNodeId += 1
                        n4 = maxNodeId
                        x = (n1x + n2x) / 2.0
                        y = (n1y + n2y) / 2.0
                        z = (n1z + n2z) / 2.0
                        curNode = dpFem.node(n4, x, y, z)
                        self.nodes[n4] = curNode
                        self.elems[elid].set_nodes([n1a, n3a, n4])
                        self.elems[elid].set_type('tria3')
                        self.elems[nElid].set_nodes([n1a, n4, n3b])
                        self.elems[nElid].set_type('tria3')
                        checkedElems[nElid] = True
                        maxElemId += 1
                        elid1 = maxElemId
                        maxElemId += 1
                        elid2 = maxElemId
                        curElem1 = dpFem.element(elid1, [n4, n3a, n2a], 'tria3')
                        self.elems[elid1] = curElem1
                        curElem2 = dpFem.element(elid2, [n4, n2a, n3b], 'tria3')
                        self.elems[elid2] = curElem2

        if echo == True:
            stdout.write('/%s\n' % repr(countSubDiv))
            stdout.flush()

    def flipEdge(self, featureAngle, quality=1.0):
        stdout.write(' ... flip tria edges     \n')
        elemIds = []
        for elid in self.elems:
            noList = self.elems[elid].get_nodes()
            Q = self.computeTriaQuality(noList[0], noList[1], noList[2])
            if Q > quality:
                elemIds.append(elid)

        checkedElems = {}
        countFlip = 0
        maxElids = len(elemIds)
        eCount = 0
        oldprogress = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for elid in elemIds:
            eCount += 1
            progress = float(eCount) / float(maxElids) * 10.0
            if progress > oldprogress:
                dpUtils.progressNext(progress)
            oldprogress = progress
            if not checkedElems.has_key(elid):
                elemsAroundNode = self.getElementsAroundNode()
                noList = self.elems[elid].get_nodes()
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                le1 = self.computeDistance(n1, n2)
                le2 = self.computeDistance(n3, n2)
                le3 = self.computeDistance(n1, n3)
                if le1 > le2 and le1 > le3:
                    nB = n1
                    nC = n2
                    nA = n3
                if le2 > le1 and le2 > le3:
                    nB = n3
                    nC = n2
                    nA = n1
                if le3 > le1 and le3 > le2:
                    nB = n1
                    nC = n3
                    nA = n2
                curElemList = elemsAroundNode[nB]
                nElid = self.getNeighborElement([nB, nC], elid, curElemList)
                if nElid:
                    noList = self.elems[nElid].get_nodes()
                    for nd in noList:
                        if nd != nB and nd != nC:
                            nD = nd
                            break

                    elidQ1 = self.computeTriaQuality(nA, nB, nC)
                    nElidQ1 = self.computeTriaQuality(nB, nC, nD)
                    elidQ2 = self.computeTriaQuality(nA, nB, nD)
                    nElidQ2 = self.computeTriaQuality(nA, nC, nD)
                    if max(elidQ1, nElidQ1) > max(elidQ2, nElidQ2):
                        oldN1 = self.getNormalVector(elid)
                        oldN2 = self.getNormalVector(nElid)
                        dotprod = numpy.dot(oldN1, oldN2)
                        if dotprod > -1.0 and dotprod < 1.0:
                            angle = numpy.arccos(dotprod) * 180.0 / numpy.pi
                        else:
                            angle = 0.0
                        if angle < featureAngle:
                            self.elems[elid].set_nodes([nA, nB, nD])
                            self.elems[nElid].set_nodes([nA, nC, nD])
                            newN1 = self.getNormalVector(elid)
                            newN2 = self.getNormalVector(nElid)
                            if numpy.dot(oldN1, newN1) < 0.0:
                                self.elems[elid].set_nodes([nA, nD, nB])
                            if numpy.dot(oldN2, newN2) < 0.0:
                                self.elems[nElid].set_nodes([nA, nD, nC])
                            countFlip += 1
                    checkedElems[nElid] = True

        dpUtils.progressEnd()
        stdout.write('     -> edges flipped        : %i\n' % countFlip)
        stdout.flush()

    def computeTriaQuality(self, n1, n2, n3):
        x1, y1, z1 = self.nodes[n1].get_coord()
        x2, y2, z2 = self.nodes[n2].get_coord()
        x3, y3, z3 = self.nodes[n3].get_coord()
        a = numpy.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
        b = numpy.sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3))
        c = numpy.sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2))
        h = max(a, b, c)
        rho = 0.5 * numpy.sqrt((b + c - a) * (c + a - b) * (a + b - c) / (a + b + c))
        alpha = 1.0 / 2.0 / numpy.sqrt(3.0)
        return h / rho * alpha

    def computeMinTriaAngle(self, n0, n1, n2):
        l1 = self.computeDistance(n0, n1)
        l2 = self.computeDistance(n1, n2)
        l3 = self.computeDistance(n2, n0)
        a = 0.0
        b = 0.0
        c = 0.0
        if l1 > l2 and l1 > l3:
            c = l1
            b = l2
            a = l3
        if l2 > l1 and l2 > l3:
            c = l2
            b = l1
            a = l3
        if l3 > l2 and l3 > l1:
            c = l3
            b = l2
            a = l1
        val = (b * b + c * c - a * a) / 2.0 / b / c
        print val, numpy.arccos(val)
        alp = numpy.arccos(val)
        bet = b / a * numpy.sin(alp)
        gam = c / a * numpy.sin(alp)
        return min(alp * 180.0 / numpy.pi, bet * 180.0 / numpy.pi, gam * 180.0 / numpy.pi)

    def computeDistance(self, n1, n2):
        n1x, n1y, n1z = self.nodes[n1].get_coord()
        n2x, n2y, n2z = self.nodes[n2].get_coord()
        return numpy.sqrt((n1x - n2x) * (n1x - n2x) + (n1y - n2y) * (n1y - n2y) + (n1z - n2z) * (n1z - n2z))

    def computeCOG(self, elid):
        nodeList = self.elems[elid].get_nodes()
        xs = 0
        ys = 0
        zs = 0
        for nd in nodeList:
            xs += self.nodes[nd].get_x() / 3.0
            ys += self.nodes[nd].get_y() / 3.0
            zs += self.nodes[nd].get_z() / 3.0

        return (xs, ys, zs)

    def scaleMesh(self, xm, ym, zm, sx, sy, sz):
        nodes = self.nodes
        for node in nodes:
            nodes[node].set_x(xm + nodes[node].get_x() * sx)
            nodes[node].set_y(ym + nodes[node].get_y() * sy)
            nodes[node].set_z(zm + nodes[node].get_z() * sz)

    def bisectTriangle(self, Tri):
        """ 
        Splits one triange into four triangles. 
        """
        P1 = Tri[0]
        P2 = Tri[1]
        P3 = Tri[2]
        P1x = P1[0]
        P1y = P1[1]
        P1z = P1[2]
        P2x = P2[0]
        P2y = P2[1]
        P2z = P2[2]
        P3x = P3[0]
        P3y = P3[1]
        P3z = P3[2]
        P4 = (
         (P1x + P2x) / 2, (P1y + P2y) / 2, (P1z + P2z) / 2)
        P5 = ((P3x + P2x) / 2, (P3y + P2y) / 2, (P3z + P2z) / 2)
        P6 = ((P1x + P3x) / 2, (P1y + P3y) / 2, (P1z + P3z) / 2)
        return (
         P4, P5, P6)

    def splitTriangle(self, Tri):
        """ 
        Splits one triange into four triangles. 
        """
        P1 = Tri[0]
        P2 = Tri[1]
        P3 = Tri[2]
        P1x = P1[0]
        P1y = P1[1]
        P1z = P1[2]
        P2x = P2[0]
        P2y = P2[1]
        P2z = P2[2]
        P3x = P3[0]
        P3y = P3[1]
        P3z = P3[2]
        P4 = (
         (P1x + P2x) / 2, (P1y + P2y) / 2, (P1z + P2z) / 2)
        P5 = ((P3x + P2x) / 2, (P3y + P2y) / 2, (P3z + P2z) / 2)
        P6 = ((P1x + P3x) / 2, (P1y + P3y) / 2, (P1z + P3z) / 2)
        return (
         P4, P5, P6)

    def projectToUnitSphere(self):
        nodes = self.nodes
        elems = self.elems
        for noid in nodes:
            node = nodes[noid]
            x, y, z = node.get_coord()
            fac = 1.0 / numpy.sqrt(x * x + y * y + z * z)
            node.set_x(fac * x)
            node.set_y(fac * y)
            node.set_z(fac * z)

    def flipElementsNormal(self, refPoint):
        stdout.flush()
        stdout.write(' ... flip element normal: \n')
        stdout.flush()
        nodeObjs = self.nodes
        elemObjs = self.elems
        rp = numpy.array(refPoint)
        for elem in elemObjs:
            nodeList = elemObjs[elem].get_nodes()
            n1 = numpy.array(nodeObjs[nodeList[0]].get_coord())
            n2 = numpy.array(nodeObjs[nodeList[1]].get_coord())
            n3 = numpy.array(nodeObjs[nodeList[2]].get_coord())
            vecR = n1 - rp
            vec1 = n2 - n1
            vec2 = n3 - n1
            vecN = dpTensor.CrossProduct(vec1, vec2)
            ori = dpTensor.SimpleContraction(vecR, vecN)
            if ori < 0:
                nodeList.reverse()
            elemObjs[elem].set_nodes(nodeList)
            elemObjs[elem].set_type('tria3')

        self.nodes = nodeObjs
        self.elems = elemObjs

    def flipElementsNormal2(self):
        stdout.write(' ... flip element normal 2:  \n')
        stdout.flush()
        prec = 1e-06
        for elem in self.elems:
            elid = elem
            break

        nn, vcog = self.getNormalVectorCOG(elid)
        r0 = vcog + nn * prec
        icount = 0
        dedectIntersect = {}
        stdout.write('  -> find intersection    \n')
        stdout.flush()
        for elem in self.elems:
            if not dedectIntersect.has_key(elem):
                nodeList = self.elems[elem].get_nodes()
                n1 = numpy.array(self.nodes[nodeList[0]].get_coord())
                n2 = numpy.array(self.nodes[nodeList[1]].get_coord())
                n3 = numpy.array(self.nodes[nodeList[2]].get_coord())
                nr = n2 - n1
                ns = n3 - n1
                a11 = nr[0]
                a12 = ns[0]
                a13 = -nn[0]
                a21 = nr[1]
                a22 = ns[1]
                a23 = -nn[1]
                a31 = nr[2]
                a32 = ns[2]
                a33 = -nn[2]
                b = r0 - n1
                DET = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13)
                r1 = -1.0
                r2 = -1.0
                s1 = -1.0
                s2 = -1.0
                t1 = -1.0
                t2 = -1.0
                if abs(DET) > prec:
                    r = 1.0 / DET * ((a33 * a22 - a32 * a23) * b[0] - (a33 * a12 - a32 * a13) * b[1] + (a23 * a12 - a22 * a13) * b[2])
                    s = 1.0 / DET * (-(a33 * a21 - a31 * a23) * b[0] + (a33 * a11 - a31 * a13) * b[1] - (a23 * a11 - a21 * a13) * b[2])
                    t = 1.0 / DET * ((a32 * a21 - a31 * a22) * b[0] - (a32 * a11 - a31 * a12) * b[1] + (a22 * a11 - a21 * a12) * b[2])
                else:
                    stdout.write('\n **ERROR** (): dpMesh.flipElementsNormal2: DETERMING during intersection computation is zero!\n\n')
                    stdout.flush()
                    exit(1)
                if r >= 0 and s >= 0 and r <= 1 and s <= 1 and t > 0 and r + s <= 1:
                    icount += 1
                    if abs(r) < prec or abs(r - 1.0) < prec or abs(s) < prec or abs(s - 1.0) < prec or abs(r + s - 1) < prec:
                        curElemList = self.getElementsAroundElement(elem)
                        for el in curElemList:
                            dedectIntersect[el] = True

                        stdout.write(' **WARNING** (): dpMesh.flipElementsNormal2: Intersection is on node or edge!\n')
                        stdout.flush()
                        print 'el=', elid, ' r=', r, ' s=', s, ' t=', t

        if icount % 2 == 0:
            pass
        else:
            nList = self.elems[elid].get_nodes()
            nList.reverse()
            self.elems[elid].set_nodes(nList)
            self.elems[elid].set_type('tria3')
        stdout.write('  -> flip                 \n')
        stdout.flush()
        searchList = []
        searchList.append(elid)
        doneList = {}
        while len(searchList) > 0:
            elid = searchList.pop(0)
            refNodes = self.elems[elid].get_nodes()
            n1 = refNodes[0]
            n2 = refNodes[1]
            n3 = refNodes[2]
            curElemList = self.getElementsAroundElement(elid)
            nodePairs = [(n1, n2), (n2, n3), (n3, n1)]
            for pair in nodePairs:
                elidN = self.getNeighborElement(pair, elid, curElemList)
                if not doneList.has_key(elidN):
                    searchList.append(elidN)
                    if elidN == None:
                        stdout.write('\n **ERROR** dpMesh.flipElementsNormal2(): No Neighbor Element found!\n')
                        stdout.flush()
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                    nList = self.elems[elidN].get_nodes()
                    ndList = nList + nList
                    orient = False
                    for i in range(5):
                        if ndList[i] == pair[1] and ndList[i + 1] == pair[0]:
                            orient = True

                    if orient == False:
                        nList.reverse()
                        self.elems[elidN].set_nodes(nList)
                        self.elems[elidN].set_type('tria3')
                    doneList[elidN] = True

            doneList[elid] = True

        return

    def getNormalVectorCOG(self, elid):
        nodeList = self.elems[elid].get_nodes()
        n1 = numpy.array(self.nodes[nodeList[0]].get_coord())
        n2 = numpy.array(self.nodes[nodeList[1]].get_coord())
        n3 = numpy.array(self.nodes[nodeList[2]].get_coord())
        vec1 = n2 - n1
        vec2 = n3 - n1
        nn = dpTensor.UnitVector(dpTensor.CrossProduct(vec1, vec2))
        elCOG = self.getCOG()
        vcog = numpy.array(elCOG[elid])
        return (
         nn, vcog)

    def getNormalVector(self, elid):
        nodeList = self.elems[elid].get_nodes()
        n1 = numpy.array(self.nodes[nodeList[0]].get_coord())
        n2 = numpy.array(self.nodes[nodeList[1]].get_coord())
        n3 = numpy.array(self.nodes[nodeList[2]].get_coord())
        vec1 = n2 - n1
        vec2 = n3 - n1
        return dpTensor.UnitVector(dpTensor.CrossProduct(vec1, vec2))

    def correctThicknessAndNormalAtBoundary(self, lenX, lenY, lenZ, tol):
        stdout.write(' ... correct boundary normals \n')
        stdout.write('     -> tolerance (voxels)   : %g\n' % tol)
        stdout.flush()
        for noid in self.nodes:
            pt = numpy.array(self.nodes[noid].get_coord())
            direct = []
            if abs(pt[0]) < tol:
                direct.append(1)
            if abs(pt[0] - lenX) < tol:
                direct.append(-1)
            if abs(pt[1]) < tol:
                direct.append(2)
            if abs(pt[1] - lenY) < tol:
                direct.append(-2)
            if abs(pt[2]) < tol:
                direct.append(3)
            if abs(pt[2] - lenZ) < tol:
                direct.append(-3)
            if len(direct) > 0:
                for cdir in direct:
                    if self.thick.has_key(noid):
                        self.normals[noid] = dpTensor.UnitVector(self.normals[noid])
                        tNormal = self.normals[noid] * self.thick[noid]
                        tNormal[abs(cdir) - 1] = 0.0
                        self.thick[noid] = dpTensor.Length(tNormal)
                    if cdir > 0:
                        self.normals[noid][abs(cdir) - 1] = 1e-06
                    else:
                        self.normals[noid][abs(cdir) - 1] = -1e-06
                    self.normals[noid] = dpTensor.UnitVector(self.normals[noid])

    def computeThicknessAndNormal(self, thickVoxelModel, dims, delta, direct=1, distortParam=None, smoothParam=None, correctParam=None, refinParam=None):
        edgeLengthRatio = 3.0
        thicknessRatio = 2.0
        cosA = 0.5
        if distortParam != None:
            edgeLengthRatio = float(distortParam['edgeLengthRatio'])
            thicknessRatio = float(distortParam['thicknessRatio'])
            cosA = float(distortParam['cosA'])
        itMax = 15
        thickWeight = 0.5
        normalWeight = 2.0
        surfWeight = 0.5
        shapes = thickVoxelModel.shape
        nx = shapes[2]
        ny = shapes[1]
        nz = shapes[0]
        if smoothParam != None:
            itMax = int(smoothParam['itMax'])
            thickWeight = float(smoothParam['thickWeight'])
            normalWeight = float(smoothParam['normalWeight'])
            surfWeight = float(smoothParam['surfWeight'])
            jumpThick = float(smoothParam['jumpThick'])
            maxStepDist = float(smoothParam['maxStepDist'])
        thickCorrect = 0.0
        minThickness = 0.5
        minThickness = 5.0
        if correctParam != None:
            thickCorrect = float(correctParam['thickCorrect'])
            minThickness = float(correctParam['minThickness'])
            maxThickness = float(correctParam['maxThickness'])
            boundCorrect = float(correctParam['boundCorrect'])
        refinIter = 0
        if correctParam != None:
            refinIter = int(refinParam['refinIter'])
        stdout.write('\n ===   M E S H I N G    L E V E L   %i  ===\n' % 0)
        stdout.flush()
        self.computeNormal(direct, boundTol=boundCorrect, ndim=[nx, ny, nz])
        if boundCorrect > 0.0:
            self.correctThicknessAndNormalAtBoundary(nx * dims[0], ny * dims[1], nz * dims[2], float(boundCorrect))
        self.computeThickness(thickVoxelModel, dims, delta, minThickness=minThickness, maxThickness=maxThickness, thickCorrect=thickCorrect, jumpThick=jumpThick, boundTol=boundCorrect)
        if itMax == 0:
            stdout.write('     -> no smoothing               \n')
            stdout.flush()
        else:
            stdout.write('     -> smooth surface             \n')
            stdout.flush()
            nodesAroundNode = self.getNodesAroundNode()
            n1b = numpy.array([0.0, 0.0, 0.0])
            n2b = numpy.array([0.0, 0.0, 0.0])
            checkElements = {}
            meanEdgeLength = self.getMeanNodalDistance()
            for noid in self.nodes:
                checkElements[noid] = False

        for noid in self.nodes:
            if self.thick[noid] > 0.0:
                n = self.normals[noid]
                n1a = numpy.array(self.nodes[noid].get_coord())
                n1b = n1a - n * self.thick[noid]
                xi = 0.0
                yi = 0.0
                zi = 0.0
                for nnoid in nodesAroundNode[noid]:
                    nn = self.normals[nnoid]
                    n2a = numpy.array(self.nodes[nnoid].get_coord())
                    n2b = n2a - nn * self.thick[nnoid]
                    r2 = dpTensor.Length(n2b - n1b)
                    if meanEdgeLength / r2 > edgeLengthRatio or r2 / meanEdgeLength > edgeLengthRatio:
                        checkElements[noid] = True
                        checkElements[nnoid] = True
                    r1 = dpTensor.Length(n2a - n1a)
                    r2 = dpTensor.Length(n2b - n1b)
                    if r1 / r2 > edgeLengthRatio:
                        checkElements[noid] = True
                        checkElements[nnoid] = True
                    ra = dpTensor.UnitVector(n2a - n1a)
                    rb = dpTensor.UnitVector(n2b - n1b)
                    cosAc = ra[0] * rb[0] + ra[1] * rb[1] + ra[2] * rb[2]
                    if cosAc < cosA:
                        checkElements[noid] = True
                        checkElements[nnoid] = True
                    if self.thick[noid] > thicknessRatio * r1:
                        checkElements[noid] = True
                        checkElements[nnoid] = True
                    if self.thick[nnoid] / self.thick[noid] > thicknessRatio:
                        checkElements[noid] = True

        nSmooth = 0
        checkElementsOld = {}
        for i in range(1):
            for noid in self.nodes:
                checkElementsOld[noid] = checkElements[noid]
                if checkElementsOld[noid] == True:
                    nSmooth += 1

            for noid in self.nodes:
                if checkElementsOld[noid] == True:
                    for nnoid in nodesAroundNode[noid]:
                        checkElements[nnoid] = True

        stdout.write('     -> involved nodes       : %i / %g %s\n' % (nSmooth, nSmooth / float(len(checkElements)) * 100.0, '%'))
        stdout.flush()
        nodes, elems, nOutNodes = self.getInnerSurface()
        n1 = numpy.array([0.0, 0.0, 0.0])
        n2 = numpy.array([0.0, 0.0, 0.0])
        newNormal = {}
        newPos = {}
        maxDist = 0
        stepLen = maxStepDist / float(itMax)
        oldprogress = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for it in range(itMax):
            progress = float(it) / float(itMax) * 10.0
            if progress > oldprogress:
                dpUtils.progressNext(progress)
            oldprogress = progress
            if stepLen > 0.0:
                maxDist += stepLen
                nodes, elems, nOutNodes = self.getInnerSurface(maxDist=maxDist)
            for noid in nodes:
                if checkElements[noid] == True:
                    xi = 0.0
                    yi = 0.0
                    zi = 0.0
                    for nnoid in nodesAroundNode[noid]:
                        newx, newy, newz = nodes[nnoid].get_coord()
                        xi += newx
                        yi += newy
                        zi += newz

                    nk = float(len(nodesAroundNode[noid]))
                    oldx, oldy, oldz = nodes[noid].get_coord()
                    newPos[noid] = (
                     (surfWeight * xi / nk + oldx) / (1 + surfWeight),
                     (surfWeight * yi / nk + oldy) / (1 + surfWeight),
                     (surfWeight * zi / nk + oldz) / (1 + surfWeight))

            dpUtils.progressEnd()
            for noid in nodes:
                if newPos.has_key(noid):
                    nodes[noid].set_x(newPos[noid][0])
                    nodes[noid].set_y(newPos[noid][1])
                    nodes[noid].set_z(newPos[noid][2])

        stdout.write('     -> iteration            : %i / %i       \n' % (itMax, itMax))
        stdout.write('     -> nodes outside step   : %i            \n' % nOutNodes)
        stdout.flush()
        for noid in nodes:
            if checkElements[noid] == True and self.thick[noid] > 0.0:
                n1[0] = self.nodes[noid].get_x()
                n1[1] = self.nodes[noid].get_y()
                n1[2] = self.nodes[noid].get_z()
                n2[0] = nodes[noid].get_x()
                n2[1] = nodes[noid].get_y()
                n2[2] = nodes[noid].get_z()
                oldNormal = self.normals[noid]
                oldThick = self.thick[noid]
                newNormal[noid] = dpTensor.UnitVector(n1 - n2)
                if dpTensor.SimpleContraction(newNormal[noid], self.normals[noid]) < 0.0:
                    self.normals[noid] = -newNormal[noid]
                else:
                    self.normals[noid] = newNormal[noid]

        oldNormal = {}
        for noid in nodes:
            oldNormal[noid] = self.normals[noid]

        for noid in nodes:
            mNormal = numpy.array([0.0, 0.0, 0.0])
            nnids = 0
            for nnoid in nodesAroundNode[noid]:
                if self.thick[nnoid] > 0.0:
                    mNormal[0] = mNormal[0] + oldNormal[nnoid][0]
                    mNormal[1] = mNormal[1] + oldNormal[nnoid][1]
                    mNormal[2] = mNormal[2] + oldNormal[nnoid][2]
                    nnids += 1

            if nnids > 0:
                mNormal = mNormal / float(nnids)
                self.normals[noid] = dpTensor.UnitVector((oldNormal[noid] + mNormal * normalWeight) / (1.0 + normalWeight))

        self.computeThickness(thickVoxelModel, dims, delta, minThickness=minThickness, maxThickness=maxThickness, thickCorrect=thickCorrect, jumpThick=jumpThick, boundTol=boundCorrect)
        oldThick = {}
        for noid in nodes:
            oldThick[noid] = self.thick[noid]

        for noid in nodes:
            mthick = 0.0
            nnids = 0
            for nnoid in nodesAroundNode[noid]:
                mthick += oldThick[nnoid]
                nnids += 1

            if nnids > 0:
                mthick = mthick / float(nnids)
                self.thick[noid] = (oldThick[noid] + mthick * thickWeight) / (1.0 + thickWeight)
                if stepLen > 0.0:
                    if self.thick[noid] > maxDist:
                        self.thick[noid] = maxDist

        if boundCorrect > 0.0:
            self.correctThicknessAndNormalAtBoundary(nx * dims[0], ny * dims[1], nz * dims[2], float(boundCorrect))
        for i in range(refinIter):
            stdout.write('\n ===   M E S H I N G    L E V E L   %i  ===\n' % (i + 1))
            stdout.flush()
            nodePairs = self.refine2DMeshNormalThickness()
            oldNodePos = {}
            if boundCorrect > 0.0:
                self.checkBoundaryEdgeMidnodes(boundCorrect, nx * dims[0], ny * dims[1], nz * dims[2])
                for noid in self.nodes:
                    oldNodePos[noid] = self.nodes[noid].get_coord()

            self.computeThickness(thickVoxelModel, dims, delta, minThickness=minThickness, maxThickness=maxThickness, thickCorrect=thickCorrect, jumpThick=jumpThick, boundTol=boundCorrect)
            if boundCorrect > 0.0:
                self.checkBoundaryEdgeMidnodes2(boundCorrect, nx * dims[0], ny * dims[1], nz * dims[2], oldNodePos=oldNodePos)
            stdout.write(' ========================================== \n\n')
            stdout.flush()

        return

    def computeNormal(self, direct, boundTol=0.0, ndim=None):
        stdout.write(' ... compute normal    \n')
        stdout.flush()
        vNormals = {}
        elemsAroundNode = self.getElementsAroundNode()
        for elid in self.elems:
            vNormals[elid] = self.getNormalVector(elid)

        mul = 0.0
        if direct > 0:
            mul = 1.0
        elif direct < 0:
            mul = -1.0
        else:
            stdout.write('\n **ERROR** dpMesh.computeNormal(): Incorrect direction (-1...+1!\n')
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        elemsCOG = {}
        if boundTol > 0.0:
            elemsCOG = self.getCOG()
        nx = 0
        ny = 0
        nz = 0
        if boundTol > 0.0:
            nx = ndim[0]
            ny = ndim[1]
            nz = ndim[2]
        for noid in self.nodes:
            elemsNodeList = elemsAroundNode[noid]
            noElArNo = len(elemsNodeList)
            nSum = numpy.array([0.0, 0.0, 0.0])
            nSumX = numpy.array([0.0, 0.0, 0.0])
            nSumY = numpy.array([0.0, 0.0, 0.0])
            nSumZ = numpy.array([0.0, 0.0, 0.0])
            noX = 0
            noY = 0
            noZ = 0
            for elid in elemsNodeList:
                if boundTol > 0.0:
                    pt = elemsCOG[elid]
                    if abs(pt[0]) < boundTol or abs(pt[0] - nx) < boundTol:
                        nSumX = nSumX + vNormals[elid]
                        noX += 1
                    elif abs(pt[1]) < boundTol or abs(pt[1] - ny) < boundTol:
                        nSumY = nSumY + vNormals[elid]
                        noY += 1
                    elif abs(pt[2]) < boundTol or abs(pt[2] - nz) < boundTol:
                        nSumZ = nSumZ + vNormals[elid]
                        noZ += 1
                    else:
                        nSum = nSum + vNormals[elid]
                else:
                    nSum = nSum + vNormals[elid]

            if noX == noElArNo:
                nSum = nSumX
            if noY == noElArNo:
                nSum = nSumY
            if noZ == noElArNo:
                nSum = nSumZ
            l = numpy.sqrt(nSum[0] * nSum[0] + nSum[1] * nSum[1] + nSum[2] * nSum[2])
            if abs(l) < 1e-06:
                stdout.write('\n **ERROR** dpMesh.computeNormal(): No valid normal found for Node %i' % noid)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            n = dpTensor.UnitVector(nSum)
            self.normals[noid] = n * mul

    def computeThickness(self, thickVoxelModel, dims, delta, minThickness=0.5, maxThickness=5.0, thickCorrect=0.0, jumpThick=1.5, boundTol=0.0):
        shapes = thickVoxelModel.shape
        nx = shapes[2]
        ny = shapes[1]
        nz = shapes[0]
        maxNodes = len(self.nodes)
        stdout.write(' ... compute thickness \n')
        stdout.flush()
        prec = 1e-06
        lx = dims[0]
        ly = dims[1]
        lz = dims[2]
        lm = (lx + ly + lz) / 3.0
        deltaNid = {}
        nInside = 0
        maxDelta = 0.0
        nOutsideCorr = 0
        nCortexCorr = 0
        for noid in self.nodes:
            deltaNid[noid] = delta
            n = self.normals[noid]
            x, y, z = self.nodes[noid].get_coord()
            dx = n[0] * deltaNid[noid]
            dy = n[1] * deltaNid[noid]
            dz = n[2] * deltaNid[noid]
            voxX = int(x / lx + dx / lx)
            voxY = int(y / ly + dy / ly)
            voxZ = int(z / lz + dz / lz)
            pvoxX = x + dx
            pvoxY = y + dy
            pvoxZ = z + dz
            if not (pvoxX <= nx * lx and pvoxY <= ny * ly and pvoxZ <= nz * lz and pvoxX >= 0.0 and pvoxY >= 0.0 and pvoxZ >= 0.0):
                nOutsideCorr += 1
                while 1:
                    dx = n[0] * deltaNid[noid]
                    dy = n[1] * deltaNid[noid]
                    dz = n[2] * deltaNid[noid]
                    voxX = int(x / lx + dx / lx)
                    voxY = int(y / ly + dy / ly)
                    voxZ = int(z / lz + dz / lz)
                    pvoxX = x + dx
                    pvoxY = y + dy
                    pvoxZ = z + dz
                    if pvoxX <= nx * lx and pvoxY <= ny * ly and pvoxZ <= nz * lz and pvoxX >= 0.0 and pvoxY >= 0.0 and pvoxZ >= 0.0:
                        break
                    if deltaNid[noid] <= lm:
                        stdout.write('\n **ERROR** : dpMesh.computeThickness(): Start point of node %i is outside!' % noid)
                        stdout.write('\n             of voxel model. Automatic correction failed!')
                        stdout.write('\n             Reduce -start value, try to cover model or search from inside!\n\n')
                        stdout.write('\n E N D E D  with ERRORS \n\n')
                        stdout.flush()
                        exit(1)
                    deltaNid[noid] = deltaNid[noid] - lm

            else:
                grayValue = thickVoxelModel[voxZ, voxY, voxX]
                if grayValue > 0:
                    nCortexCorr += 1
                    while 1:
                        dx = n[0] * deltaNid[noid]
                        dy = n[1] * deltaNid[noid]
                        dz = n[2] * deltaNid[noid]
                        voxX = int(x / lx + dx / lx)
                        voxY = int(y / ly + dy / ly)
                        voxZ = int(z / lz + dz / lz)
                        pvoxX = x + dx
                        pvoxY = y + dy
                        pvoxZ = z + dz
                        if not (pvoxX <= nx * lx and pvoxY <= ny * ly and pvoxZ <= nz * lz and pvoxX >= 0.0 and pvoxY >= 0.0 and pvoxZ >= 0.0):
                            stdout.write('\n **ERROR** : dpMesh.computeThickness(): Start point of node %i is in cortex' % noid)
                            stdout.write('\n             of voxel model. Automatic correction failed!')
                            stdout.write('\n             Increase the -start value and check voxel model/surface!\n\n')
                            stdout.write('\n E N D E D  with ERRORS \n\n')
                            stdout.flush()
                            exit(1)
                            break
                        else:
                            grayValue = thickVoxelModel[voxZ, voxY, voxX]
                            if not grayValue > 0:
                                break
                        deltaNid[noid] = deltaNid[noid] + lm

        if nOutsideCorr > 0:
            stdout.write('     !! nodes outside bbox   : %i !!\n' % nOutsideCorr)
        if nCortexCorr > 0:
            stdout.write('     !! nodes in cortex      : %i !!\n' % nCortexCorr)
        stdout.flush()
        count = 0
        outCount = 0
        minThickCount = 0
        maxThickCount = 0
        checkid = 0
        rayEntry = ''
        fecModel = fec.fec()
        allCheckNodes = {}
        allCheckElems = {}
        curNoid = 1
        curElid = 1
        nscale = 1.5
        escale = 3.0
        for noid in self.nodes:
            n = self.normals[noid]
            x, y, z = self.nodes[noid].get_coord()
            dx = n[0] * deltaNid[noid]
            dy = n[1] * deltaNid[noid]
            dz = n[2] * deltaNid[noid]
            rgx = x + dx
            rgy = y + dy
            rgz = z + dz
            rex = rgx - nscale * n[0]
            rey = rgy - nscale * n[1]
            rez = rgz - nscale * n[2]
            allCheckNodes[curNoid] = dpFem.node(curNoid, rgx, rgy, rgz)
            curNoid += 1
            allCheckNodes[curNoid] = dpFem.node(curNoid, rex, rey, rez)
            allCheckElems[curElid] = dpFem.element(curElid, [curNoid - 1, curNoid], 'bar2')
            curNoid += 1
            curElid += 1
            allCheckNodes[curNoid] = dpFem.node(curNoid, rgx - 0.5 * lx * escale, rgy - 0.5 * ly * escale, rgz - 0.5 * lz * escale)
            curNoid += 1
            allCheckNodes[curNoid] = dpFem.node(curNoid, rgx + 1.0 * lx * escale, rgy - 0.5 * ly * escale, rgz - 0.5 * lz * escale)
            curNoid += 1
            allCheckNodes[curNoid] = dpFem.node(curNoid, rgx - 0.5 * lx * escale, rgy + 1.0 * ly * escale, rgz - 0.5 * lz * escale)
            curNoid += 1
            allCheckNodes[curNoid] = dpFem.node(curNoid, rgx, rgy, rgz + 1.0 * lz)
            allCheckElems[curElid] = dpFem.element(curElid, [curNoid - 3, curNoid - 2, curNoid - 1, curNoid], 'tetra4')
            curNoid += 1
            curElid += 1

        fecModel.write('test-start_pt-normal' + '.case', 'ray', allCheckNodes, None, allCheckElems, None)
        noSkipNodes = 0
        for noid in self.nodes:
            if self.skipNodes.has_key(noid):
                noSkipNodes += 1

        oldprogress = 0
        dpUtils.progressStart('     -> Processed Data      : ')
        for noid in self.nodes:
            if not self.skipNodes.has_key(noid):
                count += 1
                progress = float(count) / float(maxNodes - noSkipNodes) * 10.0
                if progress > oldprogress:
                    dpUtils.progressNext(progress)
                oldprogress = progress
                n = self.normals[noid]
                x, y, z = self.nodes[noid].get_coord()
                if noid == checkid:
                    print '=========== START NODE ======================'
                    self.nodes[noid].show()
                dx = n[0] * deltaNid[noid]
                dy = n[1] * deltaNid[noid]
                dz = n[2] * deltaNid[noid]
                voxX = int(x / lx + dx / lx)
                voxY = int(y / ly + dy / ly)
                voxZ = int(z / lz + dz / lz)
                isAtBoundary = False
                bvoxX = int((x + boundTol) / lx)
                bvoxY = int((y + boundTol) / ly)
                bvoxZ = int((z + boundTol) / lz)
                if bvoxX == 0 or bvoxY == 0 or bvoxZ == 0 or bvoxX == nx or bvoxY == ny or bvoxZ == nz:
                    isAtBoundary = True
                rgx = x + dx
                rgy = y + dy
                rgz = z + dz
                ngx = -n[0]
                ngy = -n[1]
                ngz = -n[2]
                if abs(ngx) < prec:
                    ngx = prec
                if abs(ngy) < prec:
                    ngy = prec
                if abs(ngz) < prec:
                    ngz = prec
                rex = (voxX + 1) * lx
                rey = (voxY + 1) * ly
                rez = (voxZ + 1) * lz
                stepX = 1
                stepY = 1
                stepZ = 1
                if ngx < 0:
                    stepX = -1
                    rex -= lx
                if ngy < 0:
                    stepY = -1
                    rey -= ly
                if ngz < 0:
                    stepZ = -1
                    rez -= lz
                startVox = {}
                endVox = {}
                phyEndVox = {}
                runNo = 0
                runHole = True
                firstStartFound = False
                while runHole:
                    runNo += 1
                    findBone = False
                    run = True
                    while run:
                        grayValue = -1
                        if voxX < nx and voxY < ny and voxZ < nz and voxX >= 0 and voxY >= 0 and voxZ >= 0:
                            grayValue = thickVoxelModel[voxZ, voxY, voxX]
                        elif boundTol > 0.0 and abs(voxX - nx / 2) < 2 * nx and abs(voxY - ny / 2) < 2 * ny and abs(voxZ - nz / 2) < 2 * nz:
                            if voxX >= nx and voxY < ny and voxZ < nz and voxY >= 0 and voxZ >= 0 and isAtBoundary == True:
                                grayValue = thickVoxelModel[voxZ, voxY, nx - 1]
                            elif voxY >= ny and voxX < nx and voxZ < nz and voxX >= 0 and voxZ >= 0 and isAtBoundary == True:
                                grayValue = thickVoxelModel[voxZ, ny - 1, voxX]
                            elif voxZ >= nz and voxX < nx and voxY < ny and voxX >= 0 and voxY >= 0 and isAtBoundary == True:
                                grayValue = thickVoxelModel[nz - 1, voxY, voxX]
                            elif voxX < 0 and voxY < ny and voxZ < nz and voxY >= 0 and voxZ >= 0 and isAtBoundary == True:
                                grayValue = thickVoxelModel[voxZ, voxY, 0]
                            elif voxY < 0 and voxX < nx and voxZ < nz and voxX >= 0 and voxZ >= 0 and isAtBoundary == True:
                                grayValue = thickVoxelModel[voxZ, 0, voxX]
                            elif voxZ < 0 and voxX < nx and voxY < ny and voxX >= 0 and voxY >= 0 and isAtBoundary == True:
                                grayValue = thickVoxelModel[0, voxY, voxX]
                            else:
                                runHole = False
                                break
                        else:
                            runHole = False
                            break
                        if noid == checkid:
                            print noid, 'grayValue', grayValue, 'voxX, voxY, voxZ', voxX, voxY, voxZ, 'nx ny nz=', nx, ny, nz
                        if findBone == False and grayValue > 0:
                            findBone = True
                            startVox[runNo] = (
                             voxX, voxY, voxZ)
                            if firstStartFound == False:
                                if rayEntry == 'x':
                                    self.phyStartPos[noid] = numpy.array([rgx + tx * ngx, rgy + tx * ngy, rgz + tx * ngz])
                                elif rayEntry == 'y':
                                    self.phyStartPos[noid] = numpy.array([rgx + ty * ngx, rgy + ty * ngy, rgz + ty * ngz])
                                elif rayEntry == 'z':
                                    self.phyStartPos[noid] = numpy.array([rgx + tz * ngx, rgy + tz * ngy, rgz + tz * ngz])
                                else:
                                    stdout.write('\n **INTERNAL ERROR** computeThickness rayEntry not defined! \n\n')
                                    stdout.write('\n E N D E D  with ERRORS \n\n')
                                    stdout.flush()
                                    exit(1)
                                if isAtBoundary:
                                    lenX = nx * lx
                                    lenY = ny * ly
                                    lenZ = nz * lz
                                    pt = numpy.array(self.nodes[noid].get_coord())
                                    dirX = False
                                    dirY = False
                                    dirZ = False
                                    if abs(pt[0]) < boundTol or abs(pt[0] - lenX) < boundTol:
                                        dirX = True
                                    if abs(pt[1]) < boundTol or abs(pt[1] - lenY) < boundTol:
                                        dirY = True
                                    if abs(pt[2]) < boundTol or abs(pt[2] - lenZ) < boundTol:
                                        dirZ = True
                                    if not dirX:
                                        self.nodes[noid].set_x(self.phyStartPos[noid][0])
                                    if not dirY:
                                        self.nodes[noid].set_y(self.phyStartPos[noid][1])
                                    if not dirZ:
                                        self.nodes[noid].set_z(self.phyStartPos[noid][2])
                                else:
                                    self.nodes[noid].set_x(self.phyStartPos[noid][0])
                                    self.nodes[noid].set_y(self.phyStartPos[noid][1])
                                    self.nodes[noid].set_z(self.phyStartPos[noid][2])
                            if firstStartFound == False:
                                firstStartFound = True
                            if noid == checkid:
                                print '\n .....  start found ........'
                                print 'phyStart = ', self.phyStartPos[noid], 'noid ', noid, ' plane', rayEntry
                                print 'tx ty tz Plane', tx, ty, tz, rayEntry
                                print 'vox  ', voxX, voxY, voxZ
                                print 'rg   ', rgx, rgy, rgz
                                print 'ng   ', ngx, ngy, ngz
                        if findBone == True and grayValue == 0:
                            endVox[runNo] = (
                             preVoxX, preVoxY, preVoxZ)
                            if rayEntry == 'x':
                                phyEndVox[runNo] = numpy.array([rgx + tx * ngx, rgy + tx * ngy, rgz + tx * ngz])
                            elif rayEntry == 'y':
                                phyEndVox[runNo] = numpy.array([rgx + ty * ngx, rgy + ty * ngy, rgz + ty * ngz])
                            elif rayEntry == 'z':
                                phyEndVox[runNo] = numpy.array([rgx + tz * ngx, rgy + tz * ngy, rgz + tz * ngz])
                            else:
                                stdout.write('\n **INTERNAL ERROR** computeThickness rayEntry not defined! \n\n')
                                stdout.write('\n E N D E D  with ERRORS \n\n')
                                stdout.flush()
                                exit(1)
                            if noid == checkid:
                                print '.....  end found ........'
                                print 'phyEnd = ', phyEndVox[runNo], 'noid ', noid, ' plane', rayEntry
                                print 'curDist=', dpTensor.Length(self.phyStartPos[noid] - phyEndVox[runNo])
                            break
                        preVoxX = voxX
                        preVoxY = voxY
                        preVoxZ = voxZ
                        tx = (rex - rgx) / ngx
                        ty = (rey - rgy) / ngy
                        tz = (rez - rgz) / ngz
                        rayEntry = ''
                        if abs(tx) < abs(ty):
                            if abs(tx) < abs(tz):
                                voxX += stepX
                                rex += stepX * lx
                                rayEntry = 'x'
                            else:
                                voxZ += stepZ
                                rez += stepZ * lz
                                rayEntry = 'z'
                        elif abs(ty) < abs(tz):
                            voxY += stepY
                            rey += stepY * ly
                            rayEntry = 'y'
                        else:
                            voxZ += stepZ
                            rez += stepZ * lz
                            rayEntry = 'z'

                if noid == checkid:
                    print ' .....  voxel summary ........'
                    print 'start=', startVox, '  end=', endVox
                thick = maxThickness
                thickList = []
                nThick = 0
                runNo = 0
                if len(startVox) > 0:
                    for runNo in endVox:
                        Len = dpTensor.Length(phyEndVox[runNo] - self.phyStartPos[noid]) + thickCorrect
                        thickList.append(Len)

                    if noid == checkid:
                        print 'thickList', thickList
                        stdout.flush()
                    nThick = len(thickList)
                    if nThick == 0:
                        if not isAtBoundary:
                            stdout.write(' **WARNING* : Start but no end found. No valid thickness dedected for node %i!\n' % noid)
                        thick = minThickness
                    for thickNo in range(nThick):
                        if thickList[nThick - 1 - thickNo] < maxThickness:
                            thick = thickList[nThick - 1 - thickNo]
                            break

                if thick < minThickness:
                    thick = minThickness
                if thick == minThickness:
                    minThickCount += 1
                if thick == maxThickness and len(startVox) == 0:
                    thick = minThickness
                    if isAtBoundary == False:
                        minThickCount += 1
                if thick == maxThickness:
                    maxThickCount += 1
                self.thick[noid] = thick
                if noid == checkid:
                    print ' .....  thickness summary ........'
                    print 'minThickness=', minThickness
                    print 'maxThickness=', maxThickness
                    print 'thickList=', thickList
                    print 'nThick=', nThick
                    print 'thick    =', thick
                    print '========================================'
                    stdout.flush()

        dpUtils.progressEnd()
        if minThickCount > 0:
            stdout.write('     !! nodes with min Thick : %i !!\n' % minThickCount)
        if maxThickCount > 0:
            stdout.write('     !! nodes with max Thick : %i !!\n' % maxThickCount)
        stdout.flush()
        if jumpThick > 0.0:
            parMeanThick = jumpThick
            nodesAroundNode = self.getNodesAroundNode()
            nIter = 3
            stdout.write('     -> local correction 1                     \n')
            stdout.flush()
            for i in range(nIter):
                for noid in self.nodes:
                    meanThick = 0
                    nodeThick = self.thick[noid]
                    nNo = 0
                    for nnoid in nodesAroundNode[noid]:
                        meanThick += self.thick[nnoid]
                        nNo += 1

                    if nNo > 0:
                        meanThick = meanThick / float(nNo)
                        if nodeThick > parMeanThick * meanThick:
                            self.thick[noid] = meanThicknoid

            stdout.write('     -> local correction 2                     \n')
            for i in range(nIter):
                for noid in self.nodes:
                    meanThick = 0
                    nodeThick = self.thick[noid]
                    nNo = 0
                    for nnoid in nodesAroundNode[noid]:
                        meanThick += self.thick[nnoid]
                        nNo += 1

                    if nNo > 0:
                        meanThick = meanThick / float(nNo)
                        if meanThick > parMeanThick * nodeThick or nodeThick > parMeanThick * meanThick:
                            self.thick[noid] = meanThick

        stdout.write('     -> processed nodes      : %i          \n' % len(self.nodes))
        stdout.flush()
        return

    def getInnerSurface(self, maxDist=None):
        nodes = {}
        nThick = 0
        for noid in self.nodes:
            n = self.normals[noid]
            t = self.thick[noid]
            if maxDist != None:
                if t > maxDist:
                    t = maxDist
                    nThick += 1
            x, y, z = self.nodes[noid].get_coord()
            r = numpy.array([x, y, z])
            rn = r - n * t
            curNode = dpFem.node(noid, rn[0], rn[1], rn[2])
            nodes[noid] = curNode

        return (
         nodes, self.elems, nThick)

    def smoothLaplacian(self, itMax, noidsToSmooth=None):
        stdout.write(' ... Laplace smoothing \n')
        stdout.flush()
        nodesAroundNode = self.getNodesAroundNode()
        if not noidsToSmooth:
            noidsToSmooth = self.nodes
        for it in range(itMax):
            for noid in noidsToSmooth:
                xi = 0.0
                yi = 0.0
                zi = 0.0
                for nnoid in nodesAroundNode[noid]:
                    x, y, z = self.nodes[nnoid].get_coord()
                    xi += x
                    yi += y
                    zi += z

                nk = float(len(nodesAroundNode[noid]))
                self.nodes[noid].set_x(xi / nk)
                self.nodes[noid].set_y(yi / nk)
                self.nodes[noid].set_z(zi / nk)

    def smoothTaubin(self, itMax, lam, kPB, noidsToSmooth=None):
        stdout.write(' ... Taubin smoothing \n')
        stdout.flush()
        weight = 1.0
        mu = 1.0 / (kPB - 1.0 / lam)
        nodesAroundNode = self.getNodesAroundNode()
        if not noidsToSmooth:
            noidsToSmooth = self.nodes
        for it in range(itMax):
            for noid in noidsToSmooth:
                xi, yi, zi = self.nodes[noid].get_coord()
                dxi = 0.0
                dyi = 0.0
                dzi = 0.0
                wij = weight / float(len(nodesAroundNode[noid]))
                for nnoid in nodesAroundNode[noid]:
                    xj, yj, zj = self.nodes[nnoid].get_coord()
                    dxi += wij * (xi - xj)
                    dyi += wij * (yi - yj)
                    dzi += wij * (zi - zj)

                if it % 2 == 0 or it == 0:
                    self.nodes[noid].set_x(xi + lam * dxi)
                    self.nodes[noid].set_y(yi + lam * dyi)
                    self.nodes[noid].set_z(zi + lam * dzi)
                else:
                    self.nodes[noid].set_x(xi + mu * dxi)
                    self.nodes[noid].set_y(yi + mu * dyi)
                    self.nodes[noid].set_z(zi + mu * dzi)

    def smoothEdge(self, itMax, noidsToSmooth=None):
        stdout.write(' ... Edge smoothing \n')
        stdout.flush()
        nodesAroundNode = self.getNodesAroundNode()
        if not noidsToSmooth:
            noidsToSmooth = self.nodes
        for it in range(itMax):
            for noid in noidsToSmooth:
                x1, y1, z1 = self.nodes[noid].get_coord()
                dx = 0.0
                dy = 0.0
                dz = 0.0
                meanDist = 0.0
                nDist = 0
                for nnoid in nodesAroundNode[noid]:
                    x2, y2, z2 = self.nodes[nnoid].get_coord()
                    dist = numpy.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
                    meanDist += dist
                    nDist += 1

                meanDist = meanDist / float(nDist)
                for nnoid in nodesAroundNode[noid]:
                    x2, y2, z2 = self.nodes[nnoid].get_coord()
                    dist = numpy.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
                    dx += (dist - meanDist) / meanDist * (x2 - x1) / dist
                    dy += (dist - meanDist) / meanDist * (y2 - y1) / dist
                    dz += (dist - meanDist) / meanDist * (z2 - z1) / dist

                self.nodes[noid].set_x(x1 + dx)
                self.nodes[noid].set_y(y1 + dy)
                self.nodes[noid].set_z(z1 + dz)

    def smoothShape(self, itMax, noidsToSmooth=None):
        stdout.write(' ... shape smoothing \n')
        stdout.flush()
        Q = {}
        for elid in self.elems:
            noList = self.elems[elid].get_nodes()
            if len(noList) != 3:
                stdout.write('\n **ERROR** dpMesh.smoothShape: only triangles implemented\n\n')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            x1, y1, z1 = self.nodes[noList[0]].get_coord()
            x2, y2, z2 = self.nodes[noList[1]].get_coord()
            x3, y3, z3 = self.nodes[noList[2]].get_coord()
            a = numpy.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
            b = numpy.sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) + (z1 - z3) * (z1 - z3))
            c = numpy.sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2))
            s = (a + b + c) / 2
            h = max(a, b, c)
            rho = numpy.sqrt((s - a) * (s - b) * (s - c) / s)
            Q[elid] = h / rho

        fecModel = fec.fec()
        fecModel.write('test-Quality.case', 'trias', self.nodes, None, self.elems, None, EscaResults=[Q])
        elemsAroundNode = self.getElementsAroundNode()
        if not noidsToSmooth:
            noidsToSmooth = self.nodes
        for it in range(itMax):
            for noid in noidsToSmooth:
                x1, y1, z1 = self.nodes[noid].get_coord()
                r1 = numpy.array([x1, y1, z1])
                sum_Q_Ij = numpy.array([0.0, 0.0, 0.0])
                sum_Q = 0.0
                for elid in elemsAroundNode[noid]:
                    noList = self.elems[elid].get_nodes()
                    i = 0
                    j = 0
                    if noList[0] == noid:
                        i = 1
                        j = 2
                    if noList[1] == noid:
                        i = 0
                        j = 2
                    if noList[2] == noid:
                        i = 0
                        j = 1
                    x2, y2, z2 = self.nodes[noList[i]].get_coord()
                    x3, y3, z3 = self.nodes[noList[j]].get_coord()
                    r2 = numpy.array([x2, y2, z2])
                    r3 = numpy.array([x3, y3, z3])
                    r32 = r2 - r3
                    r31 = r1 - r3
                    rex1 = dpTensor.CrossProduct(r31, r32)
                    rex2 = dpTensor.CrossProduct(r32, rex1)
                    l32 = dpTensor.Length(r32)
                    rex2n = l32 * dpTensor.UnitVector(rex2)
                    Ij = r3 + 0.5 * r32 + rex2n
                    sum_Q_Ij += Ij * Q[elid]
                    sum_Q += Q[elid]

                r1n = sum_Q_Ij / sum_Q
                self.nodes[noid].set_x(r1n[0])
                self.nodes[noid].set_y(r1n[1])
                self.nodes[noid].set_z(r1n[2])

        return

    def writeMeshQuality(self, filename):
        stdout.write(' ... write mesh quality \n')
        stdout.flush()
        Q = {}
        for elid in self.elems:
            noList = self.elems[elid].get_nodes()
            if len(noList) != 3:
                stdout.write('\n **ERROR** dpMesh.smoothShape: only triangles implemented\n\n')
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            Q[elid] = self.computeTriaQuality(noList[0], noList[1], noList[2])

        fecModel = fec.fec()
        fecModel.write(filename, 'quality', self.nodes, None, self.elems, None, EscaResults=[Q])
        return

    def getMeanNodalDistance(self):
        nodesAroundNode = self.getNodesAroundNode()
        checkElems = {}
        l = 0.0
        count = 0
        for noid in self.nodes:
            checkElems[noid] = True
            for nnoid in nodesAroundNode[noid]:
                if not checkElems.has_key(nnoid):
                    count += 1
                    l += self.computeDistance(noid, nnoid)

        return l / float(count)

    def getMeanThickness(self):
        tmean = 0.0
        for noid in self.nodes:
            tmean += self.thick[noid]

        return tmean / float(len(self.thick))

    def getMinNodalDistance(self):
        nodesAroundNode = self.getNodesAroundNode()
        checkElems = {}
        lmin = 1e+19
        for noid in self.nodes:
            checkElems[noid] = True
            for nnoid in nodesAroundNode[noid]:
                if not checkElems.has_key(nnoid):
                    dist = self.computeDistance(noid, nnoid)
                    if dist < lmin:
                        lmin = dist

        return lmin

    def getTetraVolume(self, n1, n2, n3, n4):
        V = 0.0
        a = numpy.array(n1.get_coord())
        b = numpy.array(n2.get_coord())
        c = numpy.array(n3.get_coord())
        d = numpy.array(n4.get_coord())
        h = numpy.cross(d - b, d - c)
        V = abs(numpy.dot(d - a, h)) / 6.0
        return V

    def getPentaVolume(self, n1, n2, n3, n4, n5, n6):
        V = 0.0
        V1 = self.getTetraVolume(n1, n2, n3, n5)
        V2 = self.getTetraVolume(n3, n6, n5, n1)
        V3 = self.getTetraVolume(n1, n5, n4, n6)
        return V1 + V2 + V3

    def getHexaVolume(self, n1, n2, n3, n4, n5, n6, n7, n8):
        V = 0.0
        V1 = self.getPentaVolume(n1, n2, n4, n5, n6, n8)
        V2 = self.getPentaVolume(n2, n3, n4, n6, n7, n8)
        return V1 + V2

    def getTriaArea(self, n1, n2, n3):
        A = 0.0
        p1 = numpy.array(n1.get_coord())
        p2 = numpy.array(n2.get_coord())
        p3 = numpy.array(n3.get_coord())
        p21 = p2 - p1
        p31 = p3 - p1
        A = 0.5 * dpTensor.Length(numpy.cross(p21, p31))
        return A

    def determinat4d(self, m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33):
        value = m03 * m12 * m21 * m30 - m02 * m13 * m21 * m30 - m03 * m11 * m22 * m30 + m01 * m13 * m22 * m30 + m02 * m11 * m23 * m30 - m01 * m12 * m23 * m30 - m03 * m12 * m20 * m31 + m02 * m13 * m20 * m31 + m03 * m10 * m22 * m31 - m00 * m13 * m22 * m31 - m02 * m10 * m23 * m31 + m00 * m12 * m23 * m31 + m03 * m11 * m20 * m32 - m01 * m13 * m20 * m32 - m03 * m10 * m21 * m32 + m00 * m13 * m21 * m32 + m01 * m10 * m23 * m32 - m00 * m11 * m23 * m32 - m02 * m11 * m20 * m33 + m01 * m12 * m20 * m33 + m02 * m10 * m21 * m33 - m00 * m12 * m21 * m33 - m01 * m10 * m22 * m33 + m00 * m11 * m22 * m33
        return value

    def isPointInsideBbox(self, pt, bbox):
        if pt[0] > bbox['xmin'] and pt[0] < bbox['xmax']:
            if pt[1] > bbox['ymin'] and pt[1] < bbox['ymax']:
                if pt[2] > bbox['zmin'] and pt[2] < bbox['zmax']:
                    return True
        return False

    def isPointInsideWedge(self, pt, n1, n2, n3, n4, n5, n6):
        if self.isPointInsideTetrahedron(pt, n1, n2, n3, n4):
            return True
        else:
            if self.isPointInsideTetrahedron(pt, n4, n5, n6, n2):
                return True
            if self.isPointInsideTetrahedron(pt, n2, n3, n6, n4):
                return True
            return False

    def isPointInsideTetrahedron(self, pt, n1, n2, n3, n4):
        pt1 = numpy.array(n1.get_coord())
        pt2 = numpy.array(n2.get_coord())
        pt3 = numpy.array(n3.get_coord())
        pt4 = numpy.array(n4.get_coord())
        D0 = self.determinat4d(pt1[0], pt1[1], pt1[2], 1.0, pt2[0], pt2[1], pt2[2], 1.0, pt3[0], pt3[1], pt3[2], 1.0, pt4[0], pt4[1], pt4[2], 1.0)
        D1 = self.determinat4d(pt[0], pt[1], pt[2], 1.0, pt2[0], pt2[1], pt2[2], 1.0, pt3[0], pt3[1], pt3[2], 1.0, pt4[0], pt4[1], pt4[2], 1.0)
        D2 = self.determinat4d(pt1[0], pt1[1], pt1[2], 1.0, pt[0], pt[1], pt[2], 1.0, pt3[0], pt3[1], pt3[2], 1.0, pt4[0], pt4[1], pt4[2], 1.0)
        D3 = self.determinat4d(pt1[0], pt1[1], pt1[2], 1.0, pt2[0], pt2[1], pt2[2], 1.0, pt[0], pt[1], pt[2], 1.0, pt4[0], pt4[1], pt4[2], 1.0)
        D4 = self.determinat4d(pt1[0], pt1[1], pt1[2], 1.0, pt2[0], pt2[1], pt2[2], 1.0, pt3[0], pt3[1], pt3[2], 1.0, pt[0], pt[1], pt[2], 1.0)
        prec = 1e-06
        if D0 > -prec and D0 < prec:
            stdout.write('\n **WARNING **: isInsideTetrahedron() Tetrahedron is degenerated!')
        if not (D0 > D1 + D2 + D3 + D4 - prec and D0 < D1 + D2 + D3 + D4 + prec):
            stdout.write('\n **WARNING **: isInsideTetrahedron() DO != D1+D2+D3+D4!')
        if D0 > 0.0:
            if D1 > 0.0 and D2 > 0.0 and D3 > 0.0 and D4 > 0.0:
                return True
            else:
                return False

        else:
            if D1 < 0.0 and D2 < 0.0 and D3 < 0.0 and D4 < 0.0:
                return True
            return False

    def getVolume(self, echo=True):
        if echo == True:
            stdout.flush()
            stdout.write(' ... compute volume           \n')
            stdout.flush()
        keys = self.elems.keys()
        elType = self.elems[keys[0]].get_type()
        volume = 0.0
        if elType == 'penta6':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                V = self.getPentaVolume(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]], self.nodes[noList[3]], self.nodes[noList[4]], self.nodes[noList[5]])
                volume += V

        elif elType == 'tetra4':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                V = self.getTetraVolume(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]], self.nodes[noList[3]])
                volume += V

        elif elType == 'hexa8':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                V = self.getHexaVolume(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]], self.nodes[noList[3]], self.nodes[noList[4]], self.nodes[noList[5]], self.nodes[noList[6]], self.nodes[noList[7]])
                volume += V

        else:
            stdout.write("\n **ERROR** dpMesh.getVolume: element type '%s' not implemented\n\n" % elType)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return volume

    def deleteTempNodes(self, echo=False):
        if echo == True:
            stdout.write(' ... delete temp nodes    \n')
            stdout.flush()
        elemsAroundNode = self.getElementsAroundNode()
        ndel = 0
        for noid in elemsAroundNode:
            if len(elemsAroundNode[noid]) == 0:
                ndel += 1
                del self.nodes[noid]

        if echo == True:
            stdout.write('     -> deleted nodes        : %i \n' % ndel)
            stdout.flush()

    def checkDetJacobi(self, echo=True):
        if echo == True:
            stdout.flush()
            stdout.write(' ... check Jacobi determinant   \n')
            stdout.flush()
        keys = self.elems.keys()
        elType = self.elems[keys[0]].get_type()
        detJ = {}
        if elType == 'penta6':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                detJ[elid] = self.getPentaDetJacobi(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]], self.nodes[noList[3]], self.nodes[noList[4]], self.nodes[noList[5]])

        elif elType == 'hexa8':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                detJ[elid] = self.getHexaDetJacobi(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]], self.nodes[noList[3]], self.nodes[noList[4]], self.nodes[noList[5]], self.nodes[noList[6]], self.nodes[noList[7]])

        else:
            stdout.write("\n **ERROR** dpMesh.getDetJacobi: element type '%s' not implemented\n\n" % elType)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        keys = detJ.keys()
        maxDetJ = detJ[keys[0]]
        minDetJ = detJ[keys[0]]
        critElem = []
        for elid in detJ:
            if detJ[elid] > maxDetJ:
                maxDetJ = detJ[elid]
            if detJ[elid] < minDetJ:
                minDetJ = detJ[elid]
            if detJ[elid] < 0.0:
                if echo == True:
                    stdout.write('     -> **Warning** neg Jacobian elid=%i\n' % elid)
                    stdout.flush()
                critElem.append(elid)

        if echo == True:
            stdout.write('     -> min |J| orig/scaled  : %g / %g\n' % (minDetJ, minDetJ / maxDetJ))
            stdout.flush()
        return (minDetJ, critElem)

    def repairMesh(self, critElem, niter, echo=True):
        if echo == True:
            stdout.write(' ... try repair mesh  \n')
            stdout.flush()
        nodeDist = {}
        for elid in self.elems:
            if self.elems[elid].get_type() == 'penta6':
                noList = self.elems[elid].get_nodes()
                p1 = numpy.array(self.nodes[noList[0]].get_coord())
                p2 = numpy.array(self.nodes[noList[1]].get_coord())
                p3 = numpy.array(self.nodes[noList[2]].get_coord())
                p4 = numpy.array(self.nodes[noList[3]].get_coord())
                p5 = numpy.array(self.nodes[noList[4]].get_coord())
                p6 = numpy.array(self.nodes[noList[5]].get_coord())
                nodeDist[noList[3]] = -p1 + p4
                nodeDist[noList[4]] = -p2 + p5
                nodeDist[noList[5]] = -p3 + p6
            else:
                stdout.write("\n **ERROR** dpMesh.repairMesh: element type '%s' not implemented\n\n" % elType)
                stdout.flush()
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)

        minDetJ = -1.0
        for it in range(niter):
            if echo == True:
                stdout.write('     -> iteration            : %i / %i \n' % (it + 1, niter))
                stdout.flush()
            critNodes = []
            for elid in critElem:
                noList = self.elems[elid].get_nodes()
                critNodes.append(noList[3])
                critNodes.append(noList[4])
                critNodes.append(noList[5])

            for noid in critNodes:
                pA = numpy.array(self.nodes[noid].get_coord())
                pN = pA - nodeDist[noid] / float(niter)
                self.nodes[noid].set_coord(pN[0], pN[1], pN[2])

            minDetJ, critElem = self.checkDetJacobi(echo=False)
            if minDetJ > 0.0:
                break

        if minDetJ < 0.0:
            stdout.write('     -> **Warning** repair not successful \n')
            stdout.flush()
        else:
            stdout.write('     -> positive |J| reached\n')
            stdout.flush()

    def getPentaDetJacobi(self, n1, n2, n3, n4, n5, n6):
        xl = numpy.zeros((3, 6), numpy.float)
        X1, Y1, Z1 = n1.get_coord()
        X2, Y2, Z2 = n2.get_coord()
        X3, Y3, Z3 = n3.get_coord()
        X4, Y4, Z4 = n4.get_coord()
        X5, Y5, Z5 = n5.get_coord()
        X6, Y6, Z6 = n6.get_coord()
        xl[(0, 0)] = X1
        xl[(1, 0)] = Y1
        xl[(2, 0)] = Z1
        xl[(0, 1)] = X2
        xl[(1, 1)] = Y2
        xl[(2, 1)] = Z2
        xl[(0, 2)] = X3
        xl[(1, 2)] = Y3
        xl[(2, 2)] = Z3
        xl[(0, 3)] = X4
        xl[(1, 3)] = Y4
        xl[(2, 3)] = Z4
        xl[(0, 4)] = X5
        xl[(1, 4)] = Y5
        xl[(2, 4)] = Z5
        xl[(0, 5)] = X6
        xl[(1, 5)] = Y6
        xl[(2, 5)] = Z6
        Nrst = []
        sq35 = numpy.sqrt(3.0 / 5.0)
        sq13 = numpy.sqrt(1.0 / 3.0)
        Nrst.append([1.0 / 6.0, 1.0 / 6.0, -sq35])
        Nrst.append([1.0 / 6.0, 4.0 / 6.0, -sq35])
        Nrst.append([4.0 / 6.0, 1.0 / 6.0, -sq35])
        Nrst.append([1.0 / 6.0, 1.0 / 6.0, sq35])
        Nrst.append([1.0 / 6.0, 4.0 / 6.0, sq35])
        Nrst.append([4.0 / 6.0, 1.0 / 6.0, sq35])
        Nrst.append([1.0 / 3.0, 1.0 / 3.0, -sq13])
        Nrst.append([1.0 / 3.0, 1.0 / 3.0, sq13])
        noid = 0
        detJList = []
        for rst in Nrst:
            noid += 1
            xi = rst[0]
            et = rst[1]
            ze = rst[2]
            a = 1.0 - xi - et
            shp = numpy.zeros((4, 6), numpy.float)
            shp[(3, 0)] = 0.5 * a * (1.0 - ze)
            shp[(3, 1)] = 0.5 * xi * (1.0 - ze)
            shp[(3, 2)] = 0.5 * et * (1.0 - ze)
            shp[(3, 3)] = 0.5 * a * (1.0 + ze)
            shp[(3, 4)] = 0.5 * xi * (1.0 + ze)
            shp[(3, 5)] = 0.5 * et * (1.0 + ze)
            shp[(0, 0)] = -0.5 * (1.0 - ze)
            shp[(0, 1)] = 0.5 * (1.0 - ze)
            shp[(0, 2)] = 0.0
            shp[(0, 3)] = -0.5 * (1.0 + ze)
            shp[(0, 4)] = 0.5 * (1.0 + ze)
            shp[(0, 5)] = 0.0
            shp[(1, 0)] = -0.5 * (1.0 - ze)
            shp[(1, 1)] = 0.0
            shp[(1, 2)] = 0.5 * (1.0 - ze)
            shp[(1, 3)] = -0.5 * (1.0 + ze)
            shp[(1, 4)] = 0.0
            shp[(1, 5)] = 0.5 * (1.0 + ze)
            shp[(2, 0)] = -0.5 * a
            shp[(2, 1)] = -0.5 * xi
            shp[(2, 2)] = -0.5 * et
            shp[(2, 3)] = 0.5 * a
            shp[(2, 4)] = 0.5 * xi
            shp[(2, 5)] = 0.5 * et
            xs = numpy.zeros((3, 3), numpy.float)
            for i in range(3):
                for j in range(3):
                    for k in range(6):
                        xs[i, j] = xs[i, j] + xl[i, k] * shp[j, k]

            detJ = xs[(0, 0)] * (xs[(1, 1)] * xs[(2, 2)] - xs[(1, 2)] * xs[(2, 1)]) - xs[(0,
                                                                                          1)] * (xs[(1,
                                                                                                     0)] * xs[(2,
                                                                                                               2)] - xs[(1,
                                                                                                                         2)] * xs[(2,
                                                                                                                                   0)]) + xs[(0,
                                                                                                                                              2)] * (xs[(1,
                                                                                                                                                         0)] * xs[(2,
                                                                                                                                                                   1)] - xs[(1,
                                                                                                                                                                             1)] * xs[(2,
                                                                                                                                                                                       0)])
            detJList.append(detJ)

        return min(detJList)

    def getHexaDetJacobi(self, n1, n2, n3, n4, n5, n6, n7, n8):
        xl = numpy.zeros((3, 8), numpy.float)
        X1, Y1, Z1 = n1.get_coord()
        X2, Y2, Z2 = n2.get_coord()
        X3, Y3, Z3 = n3.get_coord()
        X4, Y4, Z4 = n4.get_coord()
        X5, Y5, Z5 = n5.get_coord()
        X6, Y6, Z6 = n6.get_coord()
        X7, Y7, Z7 = n7.get_coord()
        X8, Y8, Z8 = n8.get_coord()
        xl[(0, 0)] = X1
        xl[(1, 0)] = Y1
        xl[(2, 0)] = Z1
        xl[(0, 1)] = X2
        xl[(1, 1)] = Y2
        xl[(2, 1)] = Z2
        xl[(0, 2)] = X3
        xl[(1, 2)] = Y3
        xl[(2, 2)] = Z3
        xl[(0, 3)] = X4
        xl[(1, 3)] = Y4
        xl[(2, 3)] = Z4
        xl[(0, 4)] = X5
        xl[(1, 4)] = Y5
        xl[(2, 4)] = Z5
        xl[(0, 5)] = X6
        xl[(1, 5)] = Y6
        xl[(2, 5)] = Z6
        xl[(0, 6)] = X7
        xl[(1, 6)] = Y7
        xl[(2, 6)] = Z7
        xl[(0, 7)] = X8
        xl[(1, 7)] = Y8
        xl[(2, 7)] = Z8
        Nrst = []
        Nrst.append([0.0, 0.0, 0.0])
        Nrst.append([1.0, 0.0, 0.0])
        Nrst.append([1.0, 1.0, 0.0])
        Nrst.append([0.0, 1.0, 0.0])
        Nrst.append([0.0, 0.0, 1.0])
        Nrst.append([1.0, 0.0, 1.0])
        Nrst.append([1.0, 1.0, 1.0])
        Nrst.append([0.0, 1.0, 1.0])
        Nrst.append([0.5, 0.5, 0.5])
        noid = 0
        detJList = []
        for rst in Nrst:
            noid += 1
            xi = rst[0]
            et = rst[1]
            ze = rst[2]
            a = 1.0 - xi - et
            shp = numpy.zeros((4, 8), numpy.float)
            rm = 1.0 - xi
            sm = 1.0 - et
            tm = 1.0 - ze
            shp[(3, 0)] = rm * sm * tm
            shp[(3, 1)] = xi * sm * tm
            shp[(3, 2)] = xi * et * tm
            shp[(3, 3)] = rm * et * tm
            shp[(3, 4)] = rm * sm * ze
            shp[(3, 5)] = xi * sm * ze
            shp[(3, 6)] = xi * et * ze
            shp[(3, 7)] = rm * et * ze
            shp[(0, 0)] = -sm * tm
            shp[(0, 1)] = sm * tm
            shp[(0, 2)] = et * tm
            shp[(0, 3)] = -et * tm
            shp[(0, 4)] = -sm * ze
            shp[(0, 5)] = sm * ze
            shp[(0, 6)] = et * ze
            shp[(0, 7)] = -et * ze
            shp[(1, 0)] = -rm * tm
            shp[(1, 1)] = -xi * tm
            shp[(1, 2)] = xi * tm
            shp[(1, 3)] = rm * tm
            shp[(1, 4)] = -rm * ze
            shp[(1, 5)] = -xi * ze
            shp[(1, 6)] = xi * ze
            shp[(1, 7)] = rm * ze
            shp[(2, 0)] = -rm * sm
            shp[(2, 1)] = -xi * sm
            shp[(2, 2)] = -xi * et
            shp[(2, 3)] = -rm * et
            shp[(2, 4)] = rm * sm
            shp[(2, 5)] = xi * sm
            shp[(2, 6)] = xi * et
            shp[(2, 7)] = rm * et
            xs = numpy.zeros((3, 3), numpy.float)
            for i in range(3):
                for j in range(3):
                    for k in range(8):
                        xs[i, j] = xs[i, j] + xl[i, k] * shp[j, k]

            detJ = xs[(0, 0)] * (xs[(1, 1)] * xs[(2, 2)] - xs[(1, 2)] * xs[(2, 1)]) - xs[(0,
                                                                                          1)] * (xs[(1,
                                                                                                     0)] * xs[(2,
                                                                                                               2)] - xs[(1,
                                                                                                                         2)] * xs[(2,
                                                                                                                                   0)]) + xs[(0,
                                                                                                                                              2)] * (xs[(1,
                                                                                                                                                         0)] * xs[(2,
                                                                                                                                                                   1)] - xs[(1,
                                                                                                                                                                             1)] * xs[(2,
                                                                                                                                                                                       0)])
            detJList.append(detJ)

        return min(detJList)

    def getHexaDetJacobi2(self, n1, n2, n3, n4, n5, n6, n7, n8):
        xl = numpy.zeros((3, 8), numpy.float)
        X1, Y1, Z1 = n1.get_coord()
        X2, Y2, Z2 = n2.get_coord()
        X3, Y3, Z3 = n3.get_coord()
        X4, Y4, Z4 = n4.get_coord()
        X5, Y5, Z5 = n5.get_coord()
        X6, Y6, Z6 = n6.get_coord()
        X7, Y7, Z7 = n7.get_coord()
        X8, Y8, Z8 = n8.get_coord()
        xl[(0, 0)] = X1
        xl[(1, 0)] = Y1
        xl[(2, 0)] = Z1
        xl[(0, 1)] = X2
        xl[(1, 1)] = Y2
        xl[(2, 1)] = Z2
        xl[(0, 2)] = X3
        xl[(1, 2)] = Y3
        xl[(2, 2)] = Z3
        xl[(0, 3)] = X4
        xl[(1, 3)] = Y4
        xl[(2, 3)] = Z4
        xl[(0, 4)] = X5
        xl[(1, 4)] = Y5
        xl[(2, 4)] = Z5
        xl[(0, 5)] = X6
        xl[(1, 5)] = Y6
        xl[(2, 5)] = Z6
        xl[(0, 6)] = X7
        xl[(1, 6)] = Y7
        xl[(2, 6)] = Z7
        xl[(0, 7)] = X8
        xl[(1, 7)] = Y8
        xl[(2, 7)] = Z8
        Nrst = []
        sq13i = 1.0
        Nrst.append([+sq13i, +sq13i, -sq13i])
        Nrst.append([+sq13i, -sq13i, -sq13i])
        Nrst.append([-sq13i, +sq13i, -sq13i])
        Nrst.append([-sq13i, -sq13i, -sq13i])
        Nrst.append([+sq13i, +sq13i, +sq13i])
        Nrst.append([+sq13i, -sq13i, +sq13i])
        Nrst.append([-sq13i, +sq13i, +sq13i])
        Nrst.append([-sq13i, -sq13i, +sq13i])
        Nrst.append([0.0, 0.0, 0.0])
        noid = 0
        detJList = []
        for rst in Nrst:
            noid += 1
            xi = rst[0]
            et = rst[1]
            ze = rst[2]
            a = 1.0 - xi - et
            shp = numpy.zeros((4, 8), numpy.float)
            shp[(3, 0)] = (1.0 - xi) * (1.0 - et) * (1.0 - ze) / 8.0
            shp[(3, 1)] = (1.0 + xi) * (1.0 - et) * (1.0 - ze) / 8.0
            shp[(3, 2)] = (1.0 + xi) * (1.0 + et) * (1.0 - ze) / 8.0
            shp[(3, 3)] = (1.0 - xi) * (1.0 + et) * (1.0 - ze) / 8.0
            shp[(3, 4)] = (1.0 - xi) * (1.0 - et) * (1.0 + ze) / 8.0
            shp[(3, 5)] = (1.0 + xi) * (1.0 - et) * (1.0 + ze) / 8.0
            shp[(3, 6)] = (1.0 + xi) * (1.0 + et) * (1.0 + ze) / 8.0
            shp[(3, 7)] = (1.0 - xi) * (1.0 + et) * (1.0 + ze) / 8.0
            shp[(0, 0)] = -(1.0 - et) * (1.0 - ze) / 8.0
            shp[(0, 1)] = (1.0 - et) * (1.0 - ze) / 8.0
            shp[(0, 2)] = (1.0 + et) * (1.0 - ze) / 8.0
            shp[(0, 3)] = -(1.0 + et) * (1.0 - ze) / 8.0
            shp[(0, 4)] = -(1.0 - et) * (1.0 + ze) / 8.0
            shp[(0, 5)] = (1.0 - et) * (1.0 + ze) / 8.0
            shp[(0, 6)] = (1.0 + et) * (1.0 + ze) / 8.0
            shp[(0, 7)] = -(1.0 + et) * (1.0 + ze) / 8.0
            shp[(1, 0)] = -(1.0 - xi) * (1.0 - ze) / 8.0
            shp[(1, 1)] = -(1.0 + xi) * (1.0 - ze) / 8.0
            shp[(1, 2)] = (1.0 + xi) * (1.0 - ze) / 8.0
            shp[(1, 3)] = (1.0 - xi) * (1.0 - ze) / 8.0
            shp[(1, 4)] = -(1.0 - xi) * (1.0 + ze) / 8.0
            shp[(1, 5)] = -(1.0 + xi) * (1.0 + ze) / 8.0
            shp[(1, 6)] = (1.0 + xi) * (1.0 + ze) / 8.0
            shp[(1, 7)] = (1.0 - xi) * (1.0 + ze) / 8.0
            shp[(2, 0)] = -(1.0 - xi) * (1.0 - et) / 8.0
            shp[(2, 1)] = -(1.0 + xi) * (1.0 - et) / 8.0
            shp[(2, 2)] = -(1.0 + xi) * (1.0 + et) / 8.0
            shp[(2, 3)] = -(1.0 - xi) * (1.0 + et) / 8.0
            shp[(2, 4)] = (1.0 - xi) * (1.0 - et) / 8.0
            shp[(2, 5)] = (1.0 + xi) * (1.0 - et) / 8.0
            shp[(2, 6)] = (1.0 + xi) * (1.0 + et) / 8.0
            shp[(2, 7)] = (1.0 - xi) * (1.0 + et) / 8.0
            xs = numpy.zeros((3, 3), numpy.float)
            for i in range(3):
                for j in range(3):
                    for k in range(8):
                        xs[i, j] = xs[i, j] + xl[i, k] * shp[j, k]

            detJ = xs[(0, 0)] * (xs[(1, 1)] * xs[(2, 2)] - xs[(1, 2)] * xs[(2, 1)]) - xs[(0,
                                                                                          1)] * (xs[(1,
                                                                                                     0)] * xs[(2,
                                                                                                               2)] - xs[(1,
                                                                                                                         2)] * xs[(2,
                                                                                                                                   0)]) + xs[(0,
                                                                                                                                              2)] * (xs[(1,
                                                                                                                                                         0)] * xs[(2,
                                                                                                                                                                   1)] - xs[(1,
                                                                                                                                                                             1)] * xs[(2,
                                                                                                                                                                                       0)])
            print '  noid=', noid, '  detJ=', detJ
            detJList.append(detJ)

        return min(detJList)

    def getArea(self, echo=True):
        if echo == True:
            stdout.flush()
            stdout.write(' ... compute area           \n')
            stdout.flush()
        keys = self.elems.keys()
        elType = self.elems[keys[0]].get_type()
        area = 0.0
        if elType == 'tria3':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                A = self.getTriaArea(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]])
                area += A

        elif elType == 'quad4':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                A1 = self.getTriaArea(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]])
                A2 = self.getTriaArea(self.nodes[noList[0]], self.nodes[noList[2]], self.nodes[noList[3]])
                A = A1 + A2
                area += A

        else:
            stdout.write("\n **ERROR** dpMesh.getArea: element type '%s' not implemented\n\n" % elType)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return area

    def getElementArea(self, echo=True):
        if echo == True:
            stdout.flush()
            stdout.write(' ... compute element area       \n')
            stdout.flush()
        keys = self.elems.keys()
        elType = self.elems[keys[0]].get_type()
        area = {}
        if elType == 'tria3':
            for elid in self.elems:
                noList = self.elems[elid].get_nodes()
                A = self.getTriaArea(self.nodes[noList[0]], self.nodes[noList[1]], self.nodes[noList[2]])
                area[elid] = A

        else:
            stdout.write("\n **ERROR** dpMesh.getArea: element type '%s' not implemented\n\n" % elType)
            stdout.flush()
            stdout.write('\n E N D E D  with ERRORS \n\n')
            stdout.flush()
            exit(1)
        return area

    def getBoundingBox(self, echo=True):
        if echo == True:
            stdout.flush()
            stdout.write(' ... compute bounding box           \n')
            stdout.flush()
        xmin = 100000000000000.0
        xmax = -100000000000000.0
        ymin = 100000000000000.0
        ymax = -100000000000000.0
        zmin = 100000000000000.0
        zmax = -100000000000000.0
        for elid in self.elems:
            noList = self.elems[elid].get_nodes()
            for noid in noList:
                x = self.nodes[noid].get_x()
                y = self.nodes[noid].get_y()
                z = self.nodes[noid].get_z()
                if x > xmax:
                    xmax = x
                if x < xmin:
                    xmin = x
                if y > ymax:
                    ymax = y
                if y < ymin:
                    ymin = y
                if z > zmax:
                    zmax = z
                if z < zmin:
                    zmin = z

        return (
         xmin, xmax, ymin, ymax, zmin, zmax)

    def checkMeshForHoles(self, echo=True):
        if echo == True:
            stdout.flush()
            stdout.write(' ... check mesh for holes        \n')
            stdout.flush()
        elemsAroundNode = self.getElementsAroundNode()
        nodePairsDict = {}
        for elid in self.elems:
            noList = self.elems[elid].get_nodes()
            checkPairs = []
            if len(noList) == 3:
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                checkPairs = [
                 (
                  n1, n2), (n2, n3), (n1, n3)]
            elif len(noList) == 4:
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                n4 = noList[3]
                checkPairs = [
                 (
                  n1, n2), (n2, n3), (n3, n4), (n1, n4)]
            else:
                stdout.write('\n **ERROR** dpMesh.checkMeshForHoles(): Elements have %i nodes! \n\n' % len(noList))
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
            for pairs in checkPairs:
                if nodePairsDict.has_key((pairs[0], pairs[1])) or nodePairsDict.has_key((pairs[1], pairs[0])):
                    continue
                    nodePairsDict[pairs[0], pairs[1]] = True
                    curElemList = elemsAroundNode[pairs[0]]
                    if not self.getNeighborElement((pairs[0], pairs[1]), elid, curElemList):
                        return False

        return True

    def equivalenceNodes(self, tolerance=0.001, echo=True):
        if echo == True:
            stdout.write(' ... equivalence nodes           \n')
            stdout.write('     -> tolerance           : %10g     \n' % tolerance)
            stdout.flush()
        posVects = {}
        doubleNodes = 0
        elemsAroundNode = self.getElementsAroundNode()
        for noid in self.nodes:
            x, y, z = self.nodes[noid].get_coord()
            posVec = numpy.sqrt(x * x + y * y + z * z)
            posVecInt = posVec
            if posVects.has_key(posVecInt):
                noid1 = posVects[posVecInt][0]
                noid2 = noid
                x1, y1, z1 = self.nodes[noid1].get_coord()
                x2, y2, z2 = self.nodes[noid2].get_coord()
                dist = numpy.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
                if dist < tolerance:
                    posVects[posVecInt].append(noid)
                    doubleNodes += 1
            else:
                posVects[posVecInt] = []
                posVects[posVecInt].append(noid)

        if echo == True:
            stdout.write('     -> nodes found         : %10i   \n' % doubleNodes)
        for posVec in posVects:
            if len(posVects[posVec]) > 1:
                origNoid = posVects[posVec][0]
                for i in range(len(posVects[posVec])):
                    if i > 0:
                        oldNoid = posVects[posVec][i]
                        elemList = elemsAroundNode[oldNoid]
                        for elid in elemList:
                            nodeList = self.elems[elid].get_nodes()
                            nid = nodeList.index(oldNoid)
                            nodeList[nid] = origNoid
                            self.elems[elid].set_nodes(nodeList)

                        del self.nodes[oldNoid]

    def equivalenceNodes2(self, tolerance=0.001, echo=True):
        if echo == True:
            stdout.write(' ... equivalence nodes           \n')
            stdout.write('     -> tolerance           : %10g     \n' % tolerance)
            stdout.flush()
        doubleNodes = 0
        equivNodes = {}
        nodeMap = {}
        for noid in self.nodes:
            if not equivNodes.has_key(noid):
                nodeMap[noid] = noid
                x1, y1, z1 = self.nodes[noid].get_coord()
                for noid2 in self.nodes:
                    if noid2 != noid:
                        x2, y2, z2 = self.nodes[noid2].get_coord()
                        dist = numpy.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2))
                        if dist < tolerance:
                            equivNodes[noid2] = True
                            doubleNodes += 1
                            nodeMap[noid2] = noid

        if echo == True:
            stdout.write('     -> nodes found         : %10i   \n' % doubleNodes)
            stdout.flush()
        for elid in self.elems:
            nodeList = self.elems[elid].get_nodes()
            newNodeList = []
            for i in range(len(nodeList)):
                newNodeList.append(nodeMap[nodeList[i]])

            self.elems[elid].set_nodes(newNodeList)

        for noid in equivNodes:
            del self.nodes[noid]

    def flipElements(self, rho, echo=True):
        if echo == True:
            stdout.write(' ... flip element egdes \n')
        meanDist = self.getMeanNodalDistance()
        elemIds = []
        for elid in self.elems:
            elemIds.append(elid)

        checkedElems = {}
        countFlip = 0
        for elid in elemIds:
            if not checkedElems.has_key(elid):
                elemsAroundNode = self.getElementsAroundNode()
                noList = self.elems[elid].get_nodes()
                n1 = noList[0]
                n2 = noList[1]
                n3 = noList[2]
                le1 = self.computeDistance(n1, n2)
                le2 = self.computeDistance(n3, n2)
                le3 = self.computeDistance(n1, n3)
                if le1 > le2 and le1 > le3:
                    nB = n1
                    nC = n2
                    nA = n3
                    lAB = le2
                    lAC = le3
                if le2 > le1 and le2 > le3:
                    nB = n3
                    nC = n2
                    nA = n1
                    lAB = le1
                    lAC = le3
                if le3 > le1 and le3 > le2:
                    nB = n1
                    nC = n3
                    nA = n2
                    lAB = le2
                    lAC = le1
                curElemList = elemsAroundNode[nB]
                nElid = self.getNeighborElement([nB, nC], elid, curElemList)
                noList = self.elems[nElid].get_nodes()
                for nd in noList:
                    if nd != nB and nd != nC:
                        nD = nd
                        break

                lBC = self.computeDistance(nB, nC)
                lDC = self.computeDistance(nD, nC)
                lBD = self.computeDistance(nB, nD)
                lAD = self.computeDistance(nA, nD)
                if lBC > meanDist * rho and lDC < meanDist * rho and lBD < meanDist * rho and lAB < meanDist * rho and lAC < meanDist * rho and lBC > lAD:
                    self.elems[elid].set_nodes([nB, nA, nD])
                    self.elems[nElid].set_nodes([nA, nD, nC])
                    countFlip += 1
                checkedElems[nElid] = True

        if echo == True:
            stdout.write('     -> edges flipped       : %i   \n' % countFlip)
            stdout.flush()


if __name__ == '__main__':
    import psyco
    psyco.full()
    power = None
    outName = None
    lout = None
    guiInit = "*entry        -power 'Power for Directions = 8*4^n'       2           yes  1      \n" + "*fileEntryOut -out   'Paraview Output File Name'          test.case   no   case       \n"
    argList = argv
    argc = len(argList)
    i = 0
    while i < argc:
        if argList[i][:3] == '-po':
            i += 1
            power = argList[i]
        elif argList[i][:3] == '-ou':
            i += 1
            outName = argList[i]
        elif argList[i][:2] == '-h':
            print __doc__
            exit(0)
        i += 1

    if not outName:
        stdout.write(__doc__)
        stdout.flush()
        stdout.write('\n **ERROR** output file name not given\n\n')
        stdout.flush()
        stdout.write('\n E N D E D  with ERRORS \n\n')
        stdout.flush()
        exit(1)
    stdout.write('\n S T A R T  dpMesh.py - by D.H.PAHR 2007\n\n')
    stdout.flush()
    startTime = clock()
    ctime1 = time()
    import mia
    myMiaModel = mia.mia()
    power2 = int(power)
    triaMesh = myMiaModel.setupSphereTriangles(power2, False)
    fecModel = fec.fec()
    allNodes, allElems, results_E = myMiaModel.setupSphereFem(triaMesh)
    myMesh = dpMesh(allNodes, allElems)
    elemsAroundNode = myMesh.getElementsAroundNode()
    nodesAroundNode = myMesh.getNodesAroundNode()
    sortElemsAroundNode = myMesh.sortElementsAroundNode()
    elemsCOG = myMesh.getCOG()
    allNewNodes, allNewElems = myMesh.createDualSimplexMesh()
    myMesh2 = dpMesh(allNewNodes, allNewElems)
    refPoint = (0.0, 0.0, 0.0)
    myMesh2.flipElementsNormal(refPoint)
    allNewNodes = myMesh2.getNodes()
    allNewElems = myMesh2.getElems()
    name, ext = outName.split('.')
    fecModel.write('simx_' + name + '.off', 'Sphere', allNewNodes, None, allNewElems, None)
    myMesh3 = dpMesh()
    myMesh3.createIcosahedronMesh()
    for n in range(power2):
        myMesh3.refineTriaMesh()

    myMesh3.projectToUnitSphere()
    allNewNodes = myMesh3.getNodes()
    allNewElems = myMesh3.getElems()
    name, ext = outName.split('.')
    fecModel.write('ico_' + name + '.off', 'Sphere', allNewNodes, None, allNewElems, None)
    sResults = [
     results_E]
    fecModel.write(outName, 'Sphere', allNodes, None, allElems, None, NscaResults=sResults)
    endTime = clock()
    ctime2 = time()
    stdout.write('\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime - startTime, ctime2 - ctime1))
    stdout.flush()
# okay decompiling dpMesh.pyc
