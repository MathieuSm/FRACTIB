# uncompyle6 version 2.13.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.12 (default, Nov 20 2017, 18:23:56) 
# [GCC 5.4.0 20160609]
# Embedded file name: /usr2/pahr/entwicklung/medtool/scripts/misc/dpFem.py
# Compiled at: 2015-02-15 17:45:31
"""
#########################################################################
File   : DPFEM.py
Author : D.H.Pahr
Date   : 1.1.2005
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage: python dpFem.py   
                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requirements:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Description: 
  This program is an fem library.  
  
#########################################################################
"""
from sys import exit, stdout
import numpy

class element:

    def __init__(self, id, nodes, type):
        self.id = id
        self.nodes = nodes
        self.type = type
        self.part = None
        self.mat = None
        self.elems = {}
        self.bbox = {}
        return

    def get_id(self):
        return self.id

    def set_id(self, _id):
        self.id = _id

    def get_type(self):
        return self.type

    def set_type(self, _type):
        self.type = _type

    def get_nodes(self):
        return self.nodes

    def set_nodes(self, nlist):
        self.nodes = nlist

    def append_node(self, node):
        self.nodes.append(node)

    def get_part(self):
        return self.part

    def set_part(self, _id):
        self.part = _id

    def get_mat(self):
        return self.mat

    def set_mat(self, _mat):
        self.mat = _mat

    def get_elems(self):
        return self.elems

    def set_elems(self, elDict):
        self.elems = elDict

    def show(self):
        print '\nELEMENT Info:'
        print 'id     =', self.id
        print 'nodes  =', self.nodes
        print 'type   =', self.type

    def get_center(self):
        x = 0.0
        y = 0.0
        z = 0.0
        for noid in self.nodes:
            x += self.nodes[noid].get_x()
            y += self.nodes[noid].get_y()
            z += self.nodes[noid].get_z()

        return (x, y, z)


class node:

    def __init__(self, id, x, y=None, z=None):
        self.id = id
        self.x = x
        self.y = y
        self.z = z
        self.elemList = []

    def get_coord(self):
        return (
         self.x, self.y, self.z)

    def get_coord_numpy(self):
        return numpy.array([self.x, self.y, self.z])

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_id(self):
        return self.id

    def set_id(self, _id):
        self.id = _id

    def set_coord(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def set_coord_numpy(self, arr):
        self.x = arr[0]
        self.y = arr[1]
        self.z = arr[2]

    def set_x(self, x):
        self.x = x

    def set_y(self, y):
        self.y = y

    def set_z(self, z):
        self.z = z

    def append_to_elemList(self, elementID):
        self.elemList.append(elementID)

    def get_elemList(self):
        return self.elemList

    def get_dimension(self):
        if self.y == None:
            return 1
        else:
            if self.z == None:
                return 2
            return 3
            return

    def show(self):
        print '\nNODE Info:'
        print 'id       =', self.id
        print 'x,y,z    =', self.x, self.y, self.z
        print 'elemList =', self.elemList


class element_type:

    def get_order(self):
        order = [
         'bar2',
         'tria3', 'tria6',
         'quad4', 'quad8',
         'hexa8', 'hexa20',
         'tetra4', 'tetra10',
         'penta6', 'penta15']
        return order

    def get_nnodes(self):
        nnodes = {'bar2': 2,'tria3': 3,
           'tria6': 6,'quad4': 4,
           'quad8': 8,'hexa8': 8,
           'hexa20': 20,'tetra4': 4,
           'tetra10': 10,'penta6': 6,
           'penta15': 15}
        return nnodes

    def get_self2xmf(self):
        elemType = {'bar2': 'Polyline','tria3': 'Triangle',
           'tria6': 'Tri_6','quad4': 'Quadrilateral',
           'quad8': 'Quad_8','hexa8': 'Hexahedron',
           'hexa20': 'Hex_20','tetra4': 'Tetrahedron',
           'tetra10': 'Tet_10','penta6': 'Wedge',
           'penta15': 'Wedge_15'}
        return elemType

    def get_frd2self(self):
        elemType = {11: 'bar2',7: 'tria3',
           8: 'tria6',9: 'quad4',
           10: 'quad8',1: 'hexa8',
           4: 'hexa20',3: 'tetra4',
           6: 'tetra10',2: 'penta6',
           5: 'penta15'}
        return elemType

    def check_type(self, _type):
        if _type in self.get_self2xmf():
            return True
        else:
            return False


class arr_solution:
    """ solution array object """

    def __init__(self, _name, _type, _values, _time, _step):
        self.dataName = _name
        self.dataType = _type
        self.dataValues = _values
        self.time = _time
        self.step = _step


class domain:
    """ mesh + results of a corresponding domain for heavy data """

    def __init__(self, _name):
        self.name = _name
        self.struct = {}
        self.node_xyz = None
        self.node_label = None
        self.elem_conn = None
        self.elem_label = None
        self.elem_mat = None
        self.elem_type = None
        self.mat_name = None
        self.bcs_nodes = {}
        self.bcs_elems = {}
        self.sol_nodes = {}
        self.sol_elems = {}
        return

    def get_name(self):
        return self.name

    def set_node_xyz(self, _node_xyz):
        self.node_xyz = _node_xyz
        self.struct[self.name + '/Node:Coordinate'] = self.node_xyz

    def get_node_xyz(self):
        return self.node_xyz

    def set_node_label(self, _node_label):
        self.node_label = _node_label
        self.struct[self.name + '/Node:Label'] = self.node_label

    def get_node_label(self):
        return self.node_label

    def set_elem_conn(self, _elem_conn):
        self.elem_conn = _elem_conn
        self.struct[self.name + '/Element:Connectivity'] = self.elem_conn

    def get_elem_conn(self):
        return self.elem_conn

    def set_elem_mat(self, _elem_mat):
        self.elem_mat = _elem_mat
        self.struct[self.name + '/Element:Material'] = self.elem_mat

    def get_elem_mat(self):
        return self.elem_mat

    def set_elem_type(self, _type):
        if _type != 'mixed':
            if not element_type().check_type(_type):
                stdout.write(' **ERROR** dpFem.set_elem_type() Element type %s unknown!\n\n' % _type)
                stdout.write('\n E N D E D  with ERRORS \n\n')
                stdout.flush()
                exit(1)
        self.elem_type = _type
        self.struct[self.name + '/Element:Type'] = self.elem_type

    def get_elem_type(self):
        return self.elem_type

    def get_elem_mat(self):
        return self.elem_mat

    def set_elem_label(self, _elem_label):
        self.elem_label = _elem_label
        self.struct[self.name + '/Element:Label'] = self.elem_label

    def get_elem_label(self):
        return self.elem_label

    def get_struct(self):
        return self.struct
