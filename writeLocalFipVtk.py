import math
import datetime

import numpy as np

# Input and ouput file names
InputName='fullInclusionXCPTransCycle'

#print 'Output Variable = ' + OutName4
outputfilename1=InputName+'-localFip.vtk'

## open txt file to write to
out1 = open(outputfilename1,'w')

fipFile = 'fullInclusionXCPTransCycle-FIP.txt'
localFipMat = np.loadtxt(fipFile,delimiter=',')

numSteps = len(localFipMat[0,:])-1
numElem  = len(localFipMat[:,0])

# get last fip
print numSteps
lastLocalFip = localFipMat[:,numSteps]
fipElemInd = localFipMat[:,0] - 1
fipElemInd = np.int_(fipElemInd)


if fipElemInd[0] == fipElemInd[1]:
        print "need to check average over gauss points"

        # get number of Gauss points
        for i in range(len(fipElemInd)-1):
                if fipElemInd[i] != fipElemInd[i+1]:
                        numGp = i + 1
                        break

        numAveFips = len(fipElemInd)/numGp                
        fipTemp = np.zeros((numAveFips,))
        fipElemIndTemp = np.zeros((numAveFips,))
        # average over guass points
        for i in range(numAveFips):
                indStart = i*numGp
                indEnd = i*numGp + numGp - 1
                fipElemIndTemp[i] = fipElemInd[startInd]
                fipTemp[i] = np.mean(lastLocalFip[startInd:endInd])

        fipElemInd = fipElemIndTemp
        lastLocalFip = fipTemp


# cell type num *(8 node hex)
cellTypeNum = 12

# data name in Visit
dataName = 'fip'

# node file name (from getNlNodesElems.py )
nodeFile = 'nlNodes.inp'
# element file name
elementFile = 'nlElements.inp'

## Open input files
nodes = np.loadtxt(nodeFile,delimiter=',')
elements = np.loadtxt(elementFile,delimiter=',',dtype=int)

numNode = len(nodes[:,0])
numElem = len(elements[:,0])

# fip for all elements
fip = np.zeros((numElem,))
fip[fipElemInd] = lastLocalFip


# file identifier
out1.write('# vtk DataFile Version 1.0 \n')

# Header
out1.write(InputName + ' ' + dataName + '\n')

# File formate
out1.write('ASCII \n')
#############################
##### Dataset structure #####
#############################

# points
out1.write('\n')
out1.write('DATASET UNSTRUCTURED_GRID \n')
out1.write('POINTS ' + str(numNode) + ' float \n')

for i in range(numNode):
        out1.write(str(nodes[i,1]) + '\t' + str(nodes[i,2]) + '\t' + str(nodes[i,3]) + '\n' )

# Cells
nodesPerElem = len(elements[0,:])-1

size = int((nodesPerElem + 1)*numElem)
out1.write('\n')
out1.write('CELLS ' + str(numElem) + ' ' + str(size) + '\n')
for i in range(numElem):
        out1.write(str(nodesPerElem) + '\t')
        for j in range(nodesPerElem):
                if j == nodesPerElem - 1:
                        out1.write(str(elements[i,j+1]-1) + '\n')
                else:
                        out1.write(str(elements[i,j+1]-1) + '\t')

# cell type
out1.write('\n')
out1.write('CELL_TYPES ' + str(numElem) + '\n')
for i in range(numElem):
        out1.write(str(cellTypeNum) + '\n')

#############################
##### Dataset dataset attributes #####
#############################
out1.write('\n')
out1.write('CELL_DATA ' + str(numElem) + '\n')
out1.write('SCALARS ' + dataName + ' float \n')
out1.write('LOOKUP_TABLE default \n')
for i in range(numElem):
        val = fip[i]
        if abs(val) < 1e-20:
                val = 1e-20
        out1.write(str(val) + '\n')

# close input and output files
out1.close()



