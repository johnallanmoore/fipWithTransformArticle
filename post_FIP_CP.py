import math
import datetime
from odbAccess import *
from abaqusConstants import *

import numpy as np

# This script one .odb file for the analysis and extracts the max
# FIP on each slip system and saves that to a file.

def getTrMat(EULER):
    TR    = np.zeros((3,3))
    # storing the sin and cos of the euler angles
    SP = math.sin(math.radians(EULER[0]))
    CP = math.cos(math.radians(EULER[0]))
    ST = math.sin(math.radians(EULER[1]))
    CT = math.cos(math.radians(EULER[1]))
    SO = math.sin(math.radians(EULER[2]))
    CO = math.cos(math.radians(EULER[2]))

    # Rotating euler angle to produce transformation matrix
    TR[0,0] =  CO*CP-SO*SP*CT
    TR[0,1] =  CO*SP+SO*CT*CP   
    TR[0,2] =  SO*ST   
    TR[1,0] = -SO*CP-SP*CO*CT 
    TR[1,1] = -SO*SP+CT*CO*CP
    TR[1,2] =  CO*ST
    TR[2,0] =  SP*ST       
    TR[2,1] = -ST*CP       
    TR[2,2] =  CT  

    return TR

def getSlipNormDir(TR):
    # making slip system data
    CMSLIP = np.zeros((3,12))
    CBSLIP = np.zeros((3,12))
    # CMSLIP Data
    CMSLIP[0,0]  =  0
    CMSLIP[1,0]  =  1
    CMSLIP[2,0]  =  0
    CMSLIP[0,1]  =  0
    CMSLIP[1,1]  =  0
    CMSLIP[2,1]  =  1
    CMSLIP[0,2]  =  1
    CMSLIP[1,2]  =  0
    CMSLIP[2,2]  =  0
    CMSLIP[0,3]  =  0
    CMSLIP[1,3]  =  0
    CMSLIP[2,3]  =  1
    CMSLIP[0,4]  =  1
    CMSLIP[1,4]  =  0
    CMSLIP[2,4]  =  0
    CMSLIP[0,5]  =  0
    CMSLIP[1,5]  =  1
    CMSLIP[2,5]  =  0
    CMSLIP[0,6]  =  0
    CMSLIP[1,6]  =  1
    CMSLIP[2,6]  =  1
    CMSLIP[0,7]  =  0
    CMSLIP[1,7]  =  1
    CMSLIP[2,7]  = -1
    CMSLIP[0,8]  =  1
    CMSLIP[1,8]  =  0
    CMSLIP[2,8]  =  1
    CMSLIP[0,9] =  1
    CMSLIP[1,9] =  0
    CMSLIP[2,9] = -1
    CMSLIP[0,10] =  1
    CMSLIP[1,10] =  1
    CMSLIP[2,10] =  0
    CMSLIP[0,11] = -1
    CMSLIP[1,11] =  1
    CMSLIP[2,11] =  0
    # CBSLIP Data
    CBSLIP[0,0]  =  1
    CBSLIP[1,0]  =  0
    CBSLIP[2,0]  =  0
    CBSLIP[0,1]  =  1
    CBSLIP[1,1]  =  0
    CBSLIP[2,1]  =  0
    CBSLIP[0,2]  =  0
    CBSLIP[1,2]  =  1
    CBSLIP[2,2]  =  0
    CBSLIP[0,3]  =  0
    CBSLIP[1,3]  =  1
    CBSLIP[2,3]  =  0
    CBSLIP[0,4]  =  0
    CBSLIP[1,4]  =  0
    CBSLIP[2,4]  =  1
    CBSLIP[0,5]  =  0
    CBSLIP[1,5]  =  0
    CBSLIP[2,5]  =  1
    CBSLIP[0,6]  =  1
    CBSLIP[1,6]  =  0
    CBSLIP[2,6]  =  0
    CBSLIP[0,7]  =  1
    CBSLIP[1,7]  =  0
    CBSLIP[2,7]  =  0
    CBSLIP[0,8]  =  0
    CBSLIP[1,8]  =  1
    CBSLIP[2,8]  =  0
    CBSLIP[0,9] =  0
    CBSLIP[1,9] =  1
    CBSLIP[2,9] =  0
    CBSLIP[0,10] =  0
    CBSLIP[1,10] =  0
    CBSLIP[2,10] =  1
    CBSLIP[0,11] =  0
    CBSLIP[1,11] =  0
    CBSLIP[2,11] =  1

    ISYS=0
    i = 0
    for ISYS in range(12):
        CMMAG=math.sqrt(CMSLIP[0,ISYS]**2 + CMSLIP[1,ISYS]**2 + CMSLIP[2,ISYS]**2)
        CBMAG=math.sqrt(CBSLIP[0,ISYS]**2 + CBSLIP[1,ISYS]**2 + CBSLIP[2,ISYS]**2)
        for i in range(3):
            CMSLIP[i,ISYS]=CMSLIP[i,ISYS]/CMMAG
            CBSLIP[i,ISYS]=CBSLIP[i,ISYS]/CBMAG

    # set indexing variables back to zero
    ISYS = 0
    i    = 0
    j    = 0
    GMSLIP = np.zeros((3,12))
    GBSLIP = np.zeros((3,12))
    for ISYS in range(12):
        for i in range(3):
            for j in range(3):
                GBSLIP[i,ISYS] = GBSLIP[i,ISYS]+TR[i,j]*CBSLIP[j,ISYS]
                GMSLIP[i,ISYS] = GMSLIP[i,ISYS]+TR[i,j]*CMSLIP[j,ISYS]

    # Rotated normal vectors for each slip system
    n1  = np.array([GMSLIP[0,0],GMSLIP[1,0],GMSLIP[2,0]])
    n2  = np.array([GMSLIP[0,1],GMSLIP[1,1],GMSLIP[2,1]])
    n3  = np.array([GMSLIP[0,2],GMSLIP[1,2],GMSLIP[2,2]])
    n4  = np.array([GMSLIP[0,3],GMSLIP[1,3],GMSLIP[2,3]])
    n5  = np.array([GMSLIP[0,4],GMSLIP[1,4],GMSLIP[2,4]])
    n6  = np.array([GMSLIP[0,5],GMSLIP[1,5],GMSLIP[2,5]])
    n7  = np.array([GMSLIP[0,6],GMSLIP[1,6],GMSLIP[2,6]])
    n8  = np.array([GMSLIP[0,7],GMSLIP[1,7],GMSLIP[2,7]])
    n9  = np.array([GMSLIP[0,8],GMSLIP[1,8],GMSLIP[2,8]])
    n10 = np.array([GMSLIP[0,9],GMSLIP[1,9],GMSLIP[2,9]])
    n11 = np.array([GMSLIP[0,10],GMSLIP[1,10],GMSLIP[2,10]])
    n12 = np.array([GMSLIP[0,11],GMSLIP[1,11],GMSLIP[2,11]])
    # Rotated direction vectors for each slip system
    t1  = np.array([GBSLIP[0,0],GBSLIP[1,0],GBSLIP[2,0]])
    t2  = np.array([GBSLIP[0,1],GBSLIP[1,1],GBSLIP[2,1]])
    t3  = np.array([GBSLIP[0,2],GBSLIP[1,2],GBSLIP[2,2]])
    t4  = np.array([GBSLIP[0,3],GBSLIP[1,3],GBSLIP[2,3]])
    t5  = np.array([GBSLIP[0,4],GBSLIP[1,4],GBSLIP[2,4]])
    t6  = np.array([GBSLIP[0,5],GBSLIP[1,5],GBSLIP[2,5]])
    t7  = np.array([GBSLIP[0,6],GBSLIP[1,6],GBSLIP[2,6]])
    t8  = np.array([GBSLIP[0,7],GBSLIP[1,7],GBSLIP[2,7]])
    t9  = np.array([GBSLIP[0,8],GBSLIP[1,8],GBSLIP[2,8]])
    t10 = np.array([GBSLIP[0,9],GBSLIP[1,9],GBSLIP[2,9]])
    t11 = np.array([GBSLIP[0,10],GBSLIP[1,10],GBSLIP[2,10]])
    t12 = np.array([GBSLIP[0,11],GBSLIP[1,11],GBSLIP[2,11]])

    n = [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12]
    t = [t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12]
    return n,t

# extract FIP only from a specific elset region
inElset = True

# elset to extract FIP in
# This is an assembly elset JR
#elsetName = 'SET-ALL'
elsetName = 'MATRIX'

# fatimi socie scaling parameter
kappa = 0.55

# yield strength
sy = 750.

pi=math.pi

# initialize the variables
EULER = np.zeros((3))

# Saving pi as a variable
PI=math.pi

# column for outputting to text file
column = 0

###############################################################################

#InputName = rndCpTransElem64kCycle
InputName = 'fullInclusionXCPTransCycle'
outputfilename1=InputName+'-FIP.txt'

 ## open txt file to write to
out1 = open(outputfilename1,'w')

outputfilename2=InputName+'-MISES-AVE.txt'
out2 = open(outputfilename2,'w')

outputfilename3=InputName+'-PE11-1-AVE.txt'
out3 = open(outputfilename3,'w')

outputfilename4=InputName+'-PE11-2-AVE.txt'
out4 = open(outputfilename4,'w')

# let you know what ODB is open
print('ODB = ' + InputName)

# opend .dob file to read from
odb = openOdb(InputName+'.odb')

# number of steps
steps = odb.steps.keys()
numSteps = len(steps)

# get element set
assembly = odb.rootAssembly.instances['PART-1-1']
elemset = assembly.elementSets[elsetName]
elemsetAve = assembly.elementSets['POLY37']


# number of elements
frame = odb.steps[steps[0]].frames[-1]
stressVar=frame.fieldOutputs['S']
if inElset == True:
    stressVar = stressVar.getSubset(region=elemset)
    stressVarAve = stressVar.getSubset(region=elemsetAve)
    numelemAve = len(stressVarAve.values)
numelem = len(stressVar.values)

storeStressN1 = np.zeros((numelem,numSteps))
storeStressN2 = np.zeros((numelem,numSteps))
storeStressN3 = np.zeros((numelem,numSteps))
storeStressN4 = np.zeros((numelem,numSteps))
storeStressN5 = np.zeros((numelem,numSteps))
storeStressN6 = np.zeros((numelem,numSteps))
storeStressN7 = np.zeros((numelem,numSteps))
storeStressN8 = np.zeros((numelem,numSteps))
storeStressN9 = np.zeros((numelem,numSteps))
storeStressN10 = np.zeros((numelem,numSteps))
storeStressN11 = np.zeros((numelem,numSteps))
storeStressN12 = np.zeros((numelem,numSteps))
storeGammaP1 = np.zeros((numelem,numSteps))
storeGammaP2 = np.zeros((numelem,numSteps))
storeGammaP3 = np.zeros((numelem,numSteps))
storeGammaP4 = np.zeros((numelem,numSteps))
storeGammaP5 = np.zeros((numelem,numSteps))
storeGammaP6 = np.zeros((numelem,numSteps))
storeGammaP7 = np.zeros((numelem,numSteps))
storeGammaP8 = np.zeros((numelem,numSteps))
storeGammaP9 = np.zeros((numelem,numSteps))
storeGammaP10 = np.zeros((numelem,numSteps))
storeGammaP11 = np.zeros((numelem,numSteps))
storeGammaP12 = np.zeros((numelem,numSteps))
storeElem = np.zeros((numelem,),dtype=int)

for sNum, s in enumerate(steps):

    # number of frames in step
    numberFrame = len(odb.steps[s].frames)

    # extract last frame 
    frame = odb.steps[s].frames[-1]

    # Stress
    stressVar  = frame.fieldOutputs['S']
    # Equivalent plastic strain in an element
    strainpVar1 = frame.fieldOutputs['SDV202']
    strainpVar2 = frame.fieldOutputs['SDV203']
    strainpVar3 = frame.fieldOutputs['SDV204']
    strainpVar4 = frame.fieldOutputs['SDV205']
    strainpVar5 = frame.fieldOutputs['SDV206']
    strainpVar6 = frame.fieldOutputs['SDV207']
    strainpVar7 = frame.fieldOutputs['SDV208']
    strainpVar8 = frame.fieldOutputs['SDV209']
    strainpVar9 = frame.fieldOutputs['SDV210']

    euler1 = frame.fieldOutputs['SDV79']
    euler2 = frame.fieldOutputs['SDV80']
    euler3 = frame.fieldOutputs['SDV81']

    if sNum == len(steps) - 1:
        stressVarAve = stressVar.getSubset(region=elemsetAve)
        strainpVar1Ave2 = strainpVar1.getSubset(region=elemsetAve)

    if sNum == len(steps) - 2:
        strainpVar1Ave1 = strainpVar1.getSubset(region=elemsetAve)

    if inElset == True:
        stressVar = stressVar.getSubset(region=elemset)
        strainpVar1 = strainpVar1.getSubset(region=elemset)
        strainpVar2 = strainpVar2.getSubset(region=elemset)
        strainpVar3 = strainpVar3.getSubset(region=elemset)
        strainpVar4 = strainpVar4.getSubset(region=elemset)
        strainpVar5 = strainpVar5.getSubset(region=elemset)
        strainpVar6 = strainpVar6.getSubset(region=elemset)
        strainpVar7 = strainpVar7.getSubset(region=elemset)
        strainpVar8 = strainpVar8.getSubset(region=elemset)
        strainpVar9 = strainpVar9.getSubset(region=elemset)

        euler1 = euler1.getSubset(region=elemset)
        euler2 = euler2.getSubset(region=elemset)
        euler3 = euler3.getSubset(region=elemset)

    numelem = len(stressVar.values)
    numelemAve = len(stressVarAve.values)

    # these just let you know its running and how long it takes
    print(s + ' ' +  str(datetime.datetime.now()))


    if sNum == len(steps) - 1:
        sumvar1 = 0.0
        sumvar2 = 0.0;
        for e in range (0,numelemAve):
            var1 = stressVarAve.values[e].mises
            sumvar1 += var1
            var2 = strainpVar1Ave2.values[e].data
            sumvar2 += var2
        avevar1 = sumvar1/numelemAve
        avevar2 = sumvar2/numelemAve

    if sNum == len(steps) - 2:
        sumvar3 = 0.0
        for e in range (0,numelemAve):
            var3 = strainpVar1Ave1.values[e].data
            sumvar3 += var3
        avevar3 = sumvar3/numelemAve

    for e in range (0,numelem):
        if sNum == 0:
            storeElem[e] = stressVar.values[e].elementLabel
        stress = np.zeros((3,3))
        #11,22,33,12,13,23,
        stress[0,0] = stressVar.values[e].data[0]
        stress[1,1] = stressVar.values[e].data[1]
        stress[2,2] = stressVar.values[e].data[2]
        stress[0,1] = stressVar.values[e].data[3]
        stress[0,2] = stressVar.values[e].data[4]
        stress[1,2] = stressVar.values[e].data[5]

        stress[1,0] = stress[0,1]
        stress[2,0] = stress[0,2]
        stress[2,1] = stress[1,2]

        strainp = np.zeros((3,3))
        #11,22,33,12,13,23,
        strainp[0,0] = strainpVar1.values[e].data
        strainp[1,1] = strainpVar5.values[e].data
        strainp[2,2] = strainpVar9.values[e].data
        strainp[0,1] = strainpVar2.values[e].data
        strainp[0,2] = strainpVar3.values[e].data
        strainp[1,2] = strainpVar6.values[e].data

        strainp[1,0] = strainp[0,1]
        strainp[2,0] = strainp[0,2]
        strainp[2,1] = strainp[1,2]

        EULER[0] = euler1.values[e].data
        EULER[1] = euler2.values[e].data
        EULER[2] = euler3.values[e].data

        TR = getTrMat(EULER)
        n,t = getSlipNormDir(TR)

        n1 = n[0]
        n2 = n[1]
        n3 = n[2]
        n4 = n[3]
        n5 = n[4]
        n6 = n[5]
        n7 = n[6]
        n8 = n[7]
        n9 = n[8]
        n10 = n[9]
        n11 = n[10]
        n12 = n[11]

        t1 = t[0]
        t2 = t[1]
        t3 = t[2]
        t4 = t[3]
        t5 = t[4]
        t6 = t[5]
        t7 = t[6]
        t8 = t[7]
        t9 = t[8]
        t10 = t[9]
        t11 = t[10]
        t12 = n[11]

        # projected plastic shear strain
        gammaP1 = 2.0*n1.dot(strainp.dot(t1))
        gammaP2 = 2.0*n2.dot(strainp.dot(t2))
        gammaP3 = 2.0*n3.dot(strainp.dot(t3))
        gammaP4 = 2.0*n4.dot(strainp.dot(t4))
        gammaP5 = 2.0*n5.dot(strainp.dot(t5))
        gammaP6 = 2.0*n6.dot(strainp.dot(t6))
        gammaP7 = 2.0*n7.dot(strainp.dot(t7))
        gammaP8 = 2.0*n8.dot(strainp.dot(t8))
        gammaP9 = 2.0*n9.dot(strainp.dot(t9))
        gammaP10 = 2.0*n10.dot(strainp.dot(t10))
        gammaP11 = 2.0*n11.dot(strainp.dot(t11))
        gammaP12 = 2.0*n12.dot(strainp.dot(t12))

        # projected normal stress
        stressN1 = n1.dot(stress.dot(n1))
        stressN2 = n2.dot(stress.dot(n2))
        stressN3 = n3.dot(stress.dot(n3))
        stressN4 = n4.dot(stress.dot(n4))
        stressN5 = n5.dot(stress.dot(n5))
        stressN6 = n6.dot(stress.dot(n6))
        stressN7 = n7.dot(stress.dot(n7))
        stressN8 = n8.dot(stress.dot(n8))
        stressN9 = n9.dot(stress.dot(n9))
        stressN10 = n10.dot(stress.dot(n10))
        stressN11 = n11.dot(stress.dot(n11))
        stressN12 = n12.dot(stress.dot(n12))

        storeGammaP1[e,sNum] = gammaP1
        storeGammaP2[e,sNum] = gammaP2
        storeGammaP3[e,sNum] = gammaP3
        storeGammaP4[e,sNum] = gammaP4
        storeGammaP5[e,sNum] = gammaP5
        storeGammaP6[e,sNum] = gammaP6
        storeGammaP7[e,sNum] = gammaP7
        storeGammaP8[e,sNum] = gammaP8
        storeGammaP9[e,sNum] = gammaP9
        storeGammaP10[e,sNum] = gammaP10
        storeGammaP11[e,sNum] = gammaP11
        storeGammaP12[e,sNum] = gammaP12

        storeStressN1[e,sNum] = stressN1
        storeStressN2[e,sNum] = stressN2
        storeStressN3[e,sNum] = stressN3
        storeStressN4[e,sNum] = stressN4
        storeStressN5[e,sNum] = stressN5
        storeStressN6[e,sNum] = stressN6
        storeStressN7[e,sNum] = stressN7
        storeStressN8[e,sNum] = stressN8
        storeStressN9[e,sNum] = stressN9
        storeStressN10[e,sNum] = stressN10
        storeStressN11[e,sNum] = stressN11
        storeStressN12[e,sNum] = stressN12
        #print stress[0,0]

if numSteps % 2 == 0:
    numFips = (numSteps - 2)/2
else:
    numFips = (numSteps - 1)/2

# saving all the fip calcs for each of the 12 slip systems
#fip = np.zeros((numelm,12))
FIP_For_Global = np.zeros(numelem)
FIP_Local_Max=np.zeros((numelem,numFips))
for e in range (0,numelem):
    #FIP_Local_Max=np.zeros((numelem,numFips))
    fip1 = 0
    fip2 = 0
    fip3 = 0
    fip4 = 0
    fip5 = 0
    fip6 = 0
    fip7 = 0
    fip8 = 0
    fip9 = 0
    fip10 = 0
    fip11 = 0
    fip12 = 0
    deltaGammaP1 = 0
    deltaGammaP2 = 0
    deltaGammaP3 = 0
    deltaGammaP4 = 0
    deltaGammaP5 = 0
    deltaGammaP6 = 0
    deltaGammaP7 = 0
    deltaGammaP8 = 0
    deltaGammaP9 = 0
    deltaGammaP10 =0
    deltaGammaP11 =0
    deltaGammaP12 =0
    #out1.write(str(storeElem[e]))
    for f in range(numFips):
        deltaGammaP1 = abs(storeGammaP1[e,2*f+2] - storeGammaP1[e,2*f+1] )
        deltaGammaP2 = abs(storeGammaP2[e,2*f+2] - storeGammaP2[e,2*f+1] )
        deltaGammaP3 = abs(storeGammaP3[e,2*f+2] - storeGammaP3[e,2*f+1] )
        deltaGammaP4 = abs(storeGammaP4[e,2*f+2] - storeGammaP4[e,2*f+1] )
        deltaGammaP5 = abs(storeGammaP5[e,2*f+2] - storeGammaP5[e,2*f+1] )
        deltaGammaP6 = abs(storeGammaP6[e,2*f+2] - storeGammaP6[e,2*f+1] )
        deltaGammaP7 = abs(storeGammaP7[e,2*f+2] - storeGammaP7[e,2*f+1] )
        deltaGammaP8 = abs(storeGammaP8[e,2*f+2] - storeGammaP8[e,2*f+1] )
        deltaGammaP9 = abs(storeGammaP9[e,2*f+2] - storeGammaP9[e,2*f+1] )
        deltaGammaP10 = abs(storeGammaP10[e,2*f+2] - storeGammaP10[e,2*f+1] )
        deltaGammaP11 = abs(storeGammaP11[e,2*f+2] - storeGammaP11[e,2*f+1] )
        deltaGammaP12 = abs(storeGammaP12[e,2*f+2] - storeGammaP12[e,2*f+1] )

        fip1 = abs(0.5*deltaGammaP1*(1 + kappa*storeStressN1[e,2*f+2]/sy))
        fip2 = abs(0.5*deltaGammaP2*(1 + kappa*storeStressN2[e,2*f+2]/sy))
        fip3 = abs(0.5*deltaGammaP3*(1 + kappa*storeStressN3[e,2*f+2]/sy))
        fip4 = abs(0.5*deltaGammaP4*(1 + kappa*storeStressN4[e,2*f+2]/sy))
        fip5 = abs(0.5*deltaGammaP5*(1 + kappa*storeStressN5[e,2*f+2]/sy))
        fip6 = abs(0.5*deltaGammaP6*(1 + kappa*storeStressN6[e,2*f+2]/sy))
        fip7 = abs(0.5*deltaGammaP7*(1 + kappa*storeStressN7[e,2*f+2]/sy))
        fip8 = abs(0.5*deltaGammaP8*(1 + kappa*storeStressN8[e,2*f+2]/sy))
        fip9 = abs(0.5*deltaGammaP9*(1 + kappa*storeStressN9[e,2*f+2]/sy))
        fip10 = abs(0.5*deltaGammaP10*(1 + kappa*storeStressN10[e,2*f+2]/sy))
        fip11 = abs(0.5*deltaGammaP11*(1 + kappa*storeStressN11[e,2*f+2]/sy))
        fip12 = abs(0.5*deltaGammaP12*(1 + kappa*storeStressN12[e,2*f+2]/sy))
        # only care about the max FIP value on each slip system
        FIP_Local_Max[e,f] = max(np.array([fip1,fip2,fip3,fip4,fip5,fip6,fip7,fip8,fip9,fip10,fip11,fip12]))
        
        #out1.write(',' + str(FIP_Local_Max[e,f]) + '\n')


## write data
for e in range(numelem):
        out1.write(str(storeElem[e]))
        for f in range(numFips):
                out1.write(', ' + str("{:e}".format(FIP_Local_Max[e,f])))
        out1.write('\n')

# write stress
out2.write(str(avevar1))
# write first strain
out3.write(str(avevar3))
# write second strain
out4.write(str(avevar2))

print 'Mises Stress in Last Step              : ' + str(avevar1)
print 'Plastic Strain X in Second to Last Step: ' + str(avevar3)
print 'Plastic Strain X in Last Step          : ' + str(avevar2)

#close odb file file
odb.close()

# close output file
out1.close()
out2.close()
out3.close()
out4.close()
        
        
