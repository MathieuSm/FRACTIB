from odbAccess import *
from textRepr import *
from string import *
from time import *

# Define frame, step and instance
input_frame = range(0, 2)
input_step = ['0']
input_instance = ['0']

# display the reading result of odb2vtk file
print('\nBasic Information:')
print('Model:' + 'Cube.odb' + '; Mesh type:' + 'Hexahedron' + '; Number of blocks:' + str(1))
print('Convert frames: ' + str(input_frame[0]) + ' to ' + str(input_frame[-1]))
print('Step & Instance : ' + str(input_step) + ', ' + str(input_instance))

# open an ODB ( Abaqus output database )
starttime = time()
odb = openOdb('Cube.odb', readOnly=True)
print('\nODB opened')

# access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
rootassembly = odb.rootAssembly
instance = rootassembly.instances
# access attribute information
step = odb.steps
# get instance & step information : Quantity and all names
allinstancestr = str(instance)
autoins = allinstancestr.split("'")
inslen = len(autoins) / 4
instance_N = range(0, inslen)
allstepstr = str(step)
autostep = allstepstr.split("'")
steplen = len(autostep) / 4
step_N = range(0, steplen)

for i in input_step:
    if (steplen < int(i)):
        print('\nInput step exceeds the range of steps')
        os._exit(0)
for i in input_instance:
    if (inslen < int(i)):
        print('\nInput instance exceeds the range of instances')
        os._exit(0)

# step cycle
for step_i in input_step:
    n = int(step_i) * 4 + 1
    stepname = autostep[n]
    print('\nStep: ' + stepname)
    # instance cycle
    for ins_i in input_instance:
        n = int(ins_i) * 4 + 1
        instancename = autoins[n]
        print('\nInstance: ' + instancename)

        # access nodes & elements
        node = instance[instancename].nodes
        element = instance[instancename].elements
        n_nodes = len(node)
        n_elements = len(element)
        # access attribute(fieldOutputs) information
        frame = step[stepname].frames

        # compute the number of element of each block
        p_elements = n_elements + 1
        lp_elements = n_elements # last block

        # match nodes' label and its order in sequence (for empty nodes in tetra mesh)
        MLN = node[n_nodes - 1].label
        TOTAL = []
        # read node in sequence, and get the largest label of node(non-empty)
        # MLN is the max label of nodeset
        for i in node:
            TOTAL.append(i.label)
            if (i.label > MLN):
                MLN = i.label
        # match (the key)
        L = []
        n = 0
        for i in range(MLN):
            L.append(0)
        for i in TOTAL:
            L[i - 1] = n
            n += 1

        # frame cycle
        for i_frame in input_frame:

            # Detect whether the input frame is out of range
            try:
                TRY = odb.steps[stepname].frames[int(i_frame)]
            except:
                print('\nInput frame exceeds the range of frames')
                os._exit(0)
                break

            # Access a frame
            N_Frame = odb.steps[stepname].frames[int(i_frame)]
            print('\nFrame:' + str(i_frame))

            # create array for store result data temporarily
            # Vector-U,A,V,RF
            L0 = []
            # Tensors-S
            L1 = []
            for i in range(MLN):
                L0.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                L1.append([0, 0])

            print('\nReading U, RF ...')
            time1 = time()
            # Access Spatial displacement
            displacement = N_Frame.fieldOutputs['U']
            fieldValues = displacement.values
            for valueX in fieldValues:
                i = valueX.nodeLabel
                L0[i - 1][0] = valueX.data[0]
                L0[i - 1][1] = valueX.data[1]
                L0[i - 1][2] = valueX.data[2]

            # Access Reaction force
            Reaction_force = N_Frame.fieldOutputs['RF']
            fieldValues = Reaction_force.values
            for valueX in fieldValues:
                i = valueX.nodeLabel
                L0[i - 1][9] = valueX.data[0]
                L0[i - 1][10] = valueX.data[1]
                L0[i - 1][11] = valueX.data[2]
            print('Time elapsed: %.3f s' % (time() - time1))
            print('\nReading Stress ...')
            time1 = time()
            # access Stress components
            Stress = N_Frame.fieldOutputs['S']
            node_Stress = Stress.getSubset(position=ELEMENT_NODAL)
            fieldValues = node_Stress.values
            for valueX in fieldValues:
                L1[valueX.nodeLabel - 1][0] += 1
                L1[valueX.nodeLabel - 1][1] += valueX.mises
            # can first ave
            print('Time elapsed: %.3f s' % (time() - time1))

            '''============================================================'''

            print('\nPartitionning model and writing vtk files ...')
            time1 = time()
            print('Frame:' + str(i_frame) + '; Block:' + str(0))
            
            # Reorganization
            # Control&Storage
            # estimate whether the node has already existed
            stg_p = []
            
            # store the reorganized node for element
            stg_e = []
            
            # store the reorganized node for node
            stg_n = []
            for i in range(MLN):
                stg_p.append(-1)
            nodecount = 0
            
            # reorganize the node and element (reconstruct the mesh)
            for i in range(n_elements):
                for j in range(8):
                    k = element[i].connectivity[j] - 1
                    if (stg_p[k] < 0):
                        stg_p[k] = nodecount
                        stg_n.append(L[k])
                        stg_e.append(nodecount)
                        nodecount += 1
                    else:
                        stg_e.append(stg_p[k])
                        
            # compute point quantity
            n_reop = len(stg_n)
            reop_N = range(0, len(stg_n))

            # create and open a VTK(.vtu) files
            outfile = open('Cube.odb'[:-4] + '_' + stepname + '_' + instancename + 'f%03d' % int(i_frame) + '.vtu', 'w')

            # <VTKFile>, including the type of mesh, version, and byte_order
            outfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">' + '\n')
            
            # <UnstructuredGrid>
            outfile.write('<UnstructuredGrid>' + '\n')
            
            # <Piece>, including the number of points and cells
            outfile.write('<Piece NumberOfPoints="' + str(n_reop) + '"' + ' ' + 'NumberOfCells="' + str(
                lp_elements) + '">' + '\n')

            print('Writing Nodes ...')
            # <Points> Write nodes into vtk files
            displacement = N_Frame.fieldOutputs['U']
            fieldValues = displacement.values
            outfile.write('<Points>' + '\n')
            outfile.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">' + '\n')
            for i in reop_N:
                nt = stg_n[i]
                k = node[stg_n[i]].label - 1
                X, Y, Z = node[nt].coordinates[0] + L0[k][0], node[nt].coordinates[1] + L0[k][1], \
                            node[nt].coordinates[2] + L0[k][2]
                outfile.write(' ' + '%11.8e' % X + '  ' + '%11.8e' % Y + '  ' + '%11.8e' % Z + '\n')
            outfile.write('</DataArray>' + '\n')
            outfile.write('</Points>' + '\n')
            # </Points>

            print('Writing Results data ...')
            # <PointData> Write results data into vtk files
            outfile.write("<" + "PointData" + " " + "Vectors=" + '"' + "Spatial_displacement,Reaction_force" + '"' \
                            + " " + "Scalars=" + '"' + "Stress_Mises" + '"' + ">" + '\n')

            # Spatial displacement, <DataArray>
            outfile.write(
                "<" + "DataArray" + " " + "type=" + '"' + "Float32" + '"' + " " + "Name=" + '"' + "Spatial_displacement" + '"' + " " + "NumberOfComponents=" + '"' + "3" + '"' + " " + "format=" + '"' + "ascii" + '"' + ">" + '\n')
            for i in reop_N:
                k = node[stg_n[i]].label - 1
                X, Y, Z = L0[k][0], L0[k][1], L0[k][2]
                outfile.write('%11.8e' % X + ' ' + '%11.8e' % Y + ' ' + '%11.8e' % Z + '\n')
            outfile.write("</DataArray>" + '\n')
            # </DataArray>

            # Reaction force
            outfile.write(
                "<" + "DataArray" + " " + "type=" + '"' + "Float32" + '"' + " " + "Name=" + '"' + "Reaction_force" + '"' + " " + "NumberOfComponents=" + '"' + "3" + '"' + " " + "format=" + '"' + "ascii" + '"' + ">" + '\n')
            for i in reop_N:
                k = node[stg_n[i]].label - 1
                X, Y, Z = L0[k][9], L0[k][10], L0[k][11]
                outfile.write('%11.8e' % X + ' ' + '%11.8e' % Y + ' ' + '%11.8e' % Z + '\n')
            outfile.write("</DataArray>" + '\n')
            # </DataArray>

            # Stress Mises, <DataArray>
            outfile.write(
                "<" + "DataArray" + " " + "type=" + '"' + "Float32" + '"' + " " + "Name=" + '"' + "Stress_Mises" + '"' + " " + "format=" + '"' + "ascii" + '"' + ">" + '\n')
            for i in reop_N:
                k = node[stg_n[i]].label - 1
                X = L1[k][1] / L1[k][0]
                outfile.write('%11.8e' % X + '\n')
            outfile.write('</DataArray>' + '\n')
            # </DataArray>

            outfile.write("</PointData>" + '\n')
            # </PointData>

            print('Writing Cells ...')
            # <Cells> Write cells into vtk files
            outfile.write('<Cells>' + '\n')
            # Connectivity
            outfile.write('<DataArray type="Int32" Name="connectivity" format="ascii">' + '\n')
            for i in range(len(stg_e) / 8):
                outfile.write(str(stg_e[i * 8]) + ' ' + str(stg_e[i * 8 + 1]) + ' ' + str(
                    stg_e[i * 8 + 2]) + ' ' + str(stg_e[i * 8 + 3]) + ' ' + str(
                    stg_e[i * 8 + 4]) + ' ' + str(stg_e[i * 8 + 5]) + ' ' + str(
                    stg_e[i * 8 + 6]) + ' ' + str(stg_e[i * 8 + 7]) + '\n')
            outfile.write('</DataArray>' + '\n')
            
            # Offsets
            outfile.write('<DataArray type="Int32" Name="offsets" format="ascii">' + '\n')
            for i in range(len(stg_e) / 8):
                outfile.write(str(i * 8 + 8) + '\n')
            outfile.write('</DataArray>' + '\n')
            
            # Type
            outfile.write('<DataArray type="UInt8" Name="types" format="ascii">' + '\n')
            for i in range(len(stg_e) / 8):
                outfile.write(str(12) + '\n')
            outfile.write('</DataArray>' + '\n')
            outfile.write('</Cells>' + '\n')
            # </Cells>

            # </Piece>
            outfile.write('</Piece>' + '\n')
            # </UnstructuredGrid>
            outfile.write('</UnstructuredGrid>' + '\n')
            # </VTKFile>
            outfile.write('</VTKFile>' + '\n')

            outfile.close()
            print('Time elapsed: %.3f s' % (time() - time1))

            '''====================================================================='''
            print('\nCreating .pvtu file for frame ' + str(i_frame) + ' ...')

        odb.close()

print('Total time elapsed: %.3f s' % (time() - starttime))
    
