import os
import numpy as np

print('--------------------------------------')
print('Python script for extrapolating element gauss-point stress to equivalent nodal stress.')
print('Created by:  Chris McGann, Pedro Arduino: University of Washington')
print('Re-Written in python by: Abolfazl Najafi: Iran Univ. of Science and Tech.')
print('Mail: abolfazlbox@gmail.com - Webpage: najafice.github.io')
print('--------------------------------------\n')


fid = "EPPR.out"

nodes = np.loadtxt(os.path.join(os.getcwd(), 'data', 'nodesInfo.dat'))
elems = np.loadtxt(os.path.join(os.getcwd(), 'data', 'elementInfo.dat'))

#----------------------PORE PRESSURES----------------------------------

porePress = np.loadtxt(os.path.join(os.getcwd(), 'data', 'porePressure.out'))

numNodes,numSteps = porePress.shape    # 4650*2 matrix > P0 AND Pn
exPWP = np.zeros((numNodes,numSteps))  # 4650*2 matrix

iniPress = porePress[:,0] # initial pore pressure at t=0 (first column)

for i in range(numSteps):
    u = porePress[:,i]
    exPWP[:,i] = np.abs(u - iniPress)

print(">>> Done with porepressure ...")

#----------------------PORE PRESSURE RATIO----------------------------------

stress = np.loadtxt(os.path.join(os.getcwd(), 'data', 'Gstress.out'))
nElem = stress.shape[0]

sVnot = stress[:,1]

exPWP = A = np.delete(exPWP, 0, 1)  # delete second row of A   # now a 4650*1 matrix
iniStress = np.zeros((numNodes,1))  # 4650*1 matrix

# local coordinate vectors
xi  =  0.125*np.array([-1, 1, 1,-1,-1, 1, 1,-1])
et  =  0.125*np.array([-1,-1, 1, 1,-1,-1, 1, 1])
ze  =  0.125*np.array([-1,-1,-1,-1, 1, 1, 1, 1])
hst =  0.125*np.array([ 1,-1, 1,-1, 1,-1, 1,-1])
hut =  0.125*np.array([ 1, 1,-1,-1,-1,-1, 1, 1])
hus =  0.125*np.array([ 1,-1,-1, 1,-1, 1, 1,-1])
hstu = 0.125*np.array([-1, 1,-1, 1, 1,-1, 1,-1])

for i in range(nElem):   # 1 to 3720 
    x = np.zeros(8)
    y = np.zeros(8)
    z = np.zeros(8)

    for j in range(8):
        x[j] = nodes[int(elems[i][j+1]),1]
        y[j] = nodes[int(elems[i][j+1]),2]
        z[j] = nodes[int(elems[i][j+1]),3]

    # define coefficents
    a1 = np.dot(x,xi)
    a2 = np.dot(x,et)
    a3 = np.dot(x,ze)
    a4 = np.dot(x,hut)
    a5 = np.dot(x,hus)
    a6 = np.dot(x,hst)
    a7 = np.dot(x,hstu)

    b1 = np.dot(y,xi)
    b2 = np.dot(y,et)
    b3 = np.dot(y,ze)
    b4 = np.dot(y,hut)
    b5 = np.dot(y,hus)
    b6 = np.dot(y,hst)
    b7 = np.dot(y,hstu)

    c1 = np.dot(z,xi)
    c2 = np.dot(z,et)
    c3 = np.dot(z,ze)
    c4 = np.dot(z,hut)
    c5 = np.dot(z,hus)
    c6 = np.dot(z,hst)
    c7 = np.dot(z,hstu)

    # define coefficent vectors
    e1 = np.array([a1,b1,c1])
    e2 = np.array([a2,b2,c2])
    e3 = np.array([a3,b3,c3])
    e4 = np.array([a4,b4,c4])
    e5 = np.array([a5,b5,c5])
    e6 = np.array([a6,b6,c6])
    e7 = np.array([a7,b7,c7])

    # jacobian terms
    J0 = np.dot(e1,np.cross(e2,e3))
    J1 = np.dot(e1,np.cross(e2,e5)) + np.dot(e1,np.cross(e6,e3))
    J2 = np.dot(e1,np.cross(e2,e4)) + np.dot(e6,np.cross(e2,e3))
    J3 = np.dot(e5,np.cross(e2,e3)) + np.dot(e1,np.cross(e4,e3))
    J4 = np.dot(e7,np.cross(e2,e3)) + np.dot(e4,np.cross(e5,e2)) + np.dot(e4,np.cross(e3,e6))
    J5 = np.dot(e1,np.cross(e7,e3)) + np.dot(e4,np.cross(e5,e1)) + np.dot(e3,np.cross(e5,e6))
    J6 = np.dot(e1,np.cross(e2,e7)) + np.dot(e4,np.cross(e1,e6)) + np.dot(e2,np.cross(e5,e6))
    J7 = -np.dot(e1,np.cross(e5,e6))
    J8 = -np.dot(e4,np.cross(e2,e6))
    J9 = -np.dot(e4,np.cross(e5,e3))
    J10 = np.dot(e2,np.cross(e4,e7))
    J11 = -np.dot(e3,np.cross(e4,e7))
    J12 = np.dot(e3,np.cross(e5,e7))
    J13 = -np.dot(e1,np.cross(e5,e7))
    J14 = np.dot(e1,np.cross(e6,e7))
    J15 = -np.dot(e2,np.cross(e6,e7))
    J16 = 2*np.dot(e4,np.cross(e5,e6))
    J17 = np.dot(e7,np.cross(e5,e6))
    J18 = np.dot(e4,np.cross(e7,e6))
    J19 = np.dot(e4,np.cross(e5,e7))

    nodeStress = sVnot[i]*(J0 + (J1*xi + J2*et + J3*ze + J7 + J8 + J9)/3.0
                + (J4*hut + J5*hus + J6*hst + J10*ze + J11*et + J12*xi + J13*ze + J14*et + J15*xi)/9.0
                + (J16*hstu + J17*hut + J18*hus + J19*hst)/27.0)

    for j in range(8):
        iniStress[int(elems[i][j+1])] += nodeStress[j]

iniStress *= 2

ppRatio = abs(exPWP/iniStress)
# print(ppRatio.shape)

with open(fid, 'w') as f:
    f.write('Pore Water Pressure Ratio\n')
    f.write('---------------------------------\n')
    
    for j in range(numNodes):
        f.write("{} \t {:<12.8e}\n".format(j+1, ppRatio.item(j)))
    f.close()
print('>>> Done with porepressure ratio ...')
#----------------------END OF RESULT FILE ----------------------------------
