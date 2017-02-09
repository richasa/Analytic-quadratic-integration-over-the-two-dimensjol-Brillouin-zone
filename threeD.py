""" 
Analytic qudtratic interpolation over the three-dimensional Brillouin zone "G wiesenekker,G received 1990" 
Redistributions of source code must retain the above copyright notice. 
"""
import numpy as np
import math
import twoD
from scipy import integrate

def energyf (k):
  x = k[0]
  y = k[1]
  z = k[2]
  E = x**2+y**2+z**2
  return E;

def energy2Df (k, q):
  x = k[0]
  y = k[1]
  E = q[0] + q[1]*k[0] + q[2]*k[1] + q[3]*k[2] + q[4]*k[0]**2 + q[5]*k[0]*k[1] + q[6]* k[0]*k[2] + q[7]* k[1]**2 + q[8]*k[1]*k[2] + q[9]* k[2]**2
  return E;

FI3 = [[1,0,0,0,0,0,0,0,0,0],
       [-3,-1,0,0,0,4,0,0,0,0],
       [-3,0,-1,0,0,4,0,0,0,0],
       [-3,0,0,-1,0,0,4,0,0,0],
       [2,2,0,0,-4,0,0,0,0,0],
       [4,0,0,0,-4,-4,0,4,0,0],
       [4,0,0,0,-4,0,-4,0,0,4],
       [2,0,2,0,0,-4,0,0,0,0],
       [4,0,0,0,0,-4,-4,0,4,0],
       [2,0,0,2,0,0,-4,0,0,0]]

def getConstantsEi(E,k):
  matriseA = [[1, k[0][0], k[0][1], k[0][2], k[0][0]**2, k[0][0]*k[0][1], k[0][0]*k[0][2], k[0][1]**2, k[0][1]*k[0][2], k[0][2]**2 ],
              [1, k[1][0], k[1][1], k[1][2], k[1][0]**2, k[1][0]*k[1][1], k[1][0]*k[1][2], k[1][1]**2, k[1][1]*k[1][2], k[1][2]**2 ],
              [1, k[2][0], k[2][1], k[2][2], k[2][0]**2, k[2][0]*k[2][1], k[2][0]*k[2][2], k[2][1]**2, k[2][1]*k[2][2], k[2][2]**2 ],
              [1, k[3][0], k[3][1], k[3][2], k[3][0]**2, k[3][0]*k[3][1], k[3][0]*k[3][2], k[3][1]**2, k[3][1]*k[3][2], k[3][2]**2 ],
              [1, k[4][0], k[4][1], k[4][2], k[4][0]**2, k[4][0]*k[4][1], k[4][0]*k[4][2], k[4][1]**2, k[4][1]*k[4][2], k[4][2]**2 ],
              [1, k[5][0], k[5][1], k[5][2], k[5][0]**2, k[5][0]*k[5][1], k[5][0]*k[5][2], k[5][1]**2, k[5][1]*k[5][2], k[5][2]**2 ],
              [1, k[6][0], k[6][1], k[6][2], k[6][0]**2, k[6][0]*k[6][1], k[6][0]*k[6][2], k[6][1]**2, k[6][1]*k[6][2], k[6][2]**2 ],
              [1, k[7][0], k[7][1], k[7][2], k[7][0]**2, k[7][0]*k[7][1], k[7][0]*k[7][2], k[7][1]**2, k[7][1]*k[7][2], k[7][2]**2 ],
              [1, k[8][0], k[8][1], k[8][2], k[8][0]**2, k[8][0]*k[8][1], k[8][0]*k[8][2], k[8][1]**2, k[8][1]*k[8][2], k[8][2]**2 ],
              [1, k[9][0], k[9][1], k[9][2], k[9][0]**2, k[9][0]*k[9][1], k[9][0]*k[9][2], k[9][1]**2, k[9][1]*k[9][2], k[9][2]**2 ]]

  matriseB =[E[0],
             E[1],
             E[2],
             E[3],
             E[4],
             E[5],
             E[6],
             E[7],
             E[8],
             E[9]]
  q = np.linalg.solve(matriseA, matriseB)
  return q;

def surfaces(e,k,E):
  """
                  (x2,y2)
                  / \ 
                 /   \ 
                /     \ 
               /       \   
      (x0,y0) /_________\(x1,y1)

  """
  k = [[0,0,0],
      [1,0,0],
      [0,1,0],
      [0,0,1],
      [0.5,0,0],
      [0,0.5,0],
      [0,0,0.5],
      [0.5,0.5,0],
      [0,0.5,0.5],
      [0.5,0,0.5]]
  q = np.dot(FI,E)
  q = np.around(q, 10)
  return 0;

def make2Dk(z, q):
  k2 =[[0, 0 , z],
      [1.0-z, 0, z],
      [0, 1.0-z, z],
      [(1.0-z)/2, 0, z],
      [(1.0-z)/2, (1.0-z)/2, z],
      [0, (1.0-z)/2, z]];
  E2 = [0,0,0,0,0,0]
  for i in range(0, 6):
    E2[i] = energy2Df(k2[i],q)

  k2 =[[0, 0 ],
      [1.0-z, 0],
      [0, 1.0-z],
      [(1.0-z)/2, 0],
      [(1.0-z)/2, (1.0-z)/2],
      [0, (1.0-z)/2]];
  return k2, E2;

def triangelIntegral(e, E, k, f):

  v = surfaces(e,k,E)
  sum = 0
  for i in range(0, 10):
    sum += f[i]*v[i]
  return sum;

#creating the test case in the artikel.
E = [0,0,0,0,0,0,0,0,0,0]
k = [[0,0,0],
      [1,0,0],
      [0,1,0],
      [0,0,1],
      [0.5,0,0],
      [0,0.5,0],
      [0,0,0.5],
      [0.5,0.5,0],
      [0.5,0,0.5],
      [0,0.5,0.5]]

for i in range(0, 10):
  E[i] = energyf (k[i])
q = getConstantsEi(E,k)

#prdouce the z point and call the twoD function
#print  triangelIntegral()
e = 0.1
f = [1,1,1,1,1,1]
def triangelIntegralInz(z):
  k2,E2 = make2Dk(z,q)
  #k2 = np.around(k2, 10)
  #E2 = np.around(E2, 10)

  return twoD.triangelIntegral(e,E2,k2,f);

##just testing ................
#k2,E2 = make2Dk( 0.14917311663,q)
#k2 = np.around(k2, 10)
#E2 = np.around(E2, 10)
#print twoD.triangelIntegral(e,E2,k2,f)
#print twoD.triangelIntegral(e,E2,k2,f)

function2d = lambda z: triangelIntegralInz(z)
print integrate.quad(function2d, 0, 0.3)
