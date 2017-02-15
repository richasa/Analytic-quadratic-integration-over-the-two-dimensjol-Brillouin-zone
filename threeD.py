""" 
Analytic qudtratic interpolation over the three-dimensional Brillouin zone "G wiesenekker,G received 1990" 
Redistributions of source code must retain the above copyright notice. 
"""
import numpy as np
import math
import twoD
from scipy import integrate


miniX = 0
maxX = 1
miniY = 0
maxY = 1
miniZ = 0
maxZ = 1
Deltak = 1

def creatkGrid():

  width = maxX - miniX
  height = maxY - miniY
  length = maxZ - miniZ 
  #rotation grid was a nessery step in testing. but not any more. 
  #rotation = 0.0
  nx = int(float(width/Deltak)) 
  ny = int(float(height/Deltak)) 
  nz = int(float(height/Deltak)) 

  kGrid = [[[0 for i in xrange(ny +1)] for i in xrange(nx +1)]for i in xrange(nz +1)]
  EGrid = [[[0 for i in xrange(ny +1)] for i in xrange(nx +1)]for i in xrange(nz +1)]
  fGrid = [[[0 for i in xrange(ny +1)] for i in xrange(nx +1)]for i in xrange(nz +1)]

  for ix in range(0,nx + 1):
    for iy in range(0,ny + 1):
      for iz in range(0,nz + 1):
  #    x = (miniX + Deltak *ix) * np.cos(rotation) + -(miniY + Deltak*iy) * np.sin(rotation)  
  #    y = (miniX + Deltak *ix) * np.sin(rotation) + (miniY + Deltak*iy) * np.cos(rotation)  
        kGrid[ix][iy][iz] = ix , iy , iz
        EGrid[ix][iy][iz] = energyf (kGrid[ix][iy][iz])
        fGrid[ix][iy][iz] = functionf(kGrid[ix][iy][iz])

  return kGrid, EGrid, fGrid;
def findmid(a,b):
  point = [0,0,0]
  point[0] = (b[0] + a[0])/2
  point[1] = (b[1] + a[1])/2
  point[2] = (b[2] + a[2])/2
  return point;
def cubeToTetrahedral (c):
  #print c
  c01 = findmid(c[0],c[1])
  c02 = findmid(c[0],c[2])
  c03 = findmid(c[0],c[3])
  c04 = findmid(c[0],c[4])
  c05 = findmid(c[0],c[5])
  c06 = findmid(c[0],c[6])
  c07 = findmid(c[0],c[7])
  c12 = findmid(c[1],c[2])
  c15 = findmid(c[1],c[5])
  c16 = findmid(c[1],c[6])
  c23 = findmid(c[2],c[3])
  c26 = findmid(c[2],c[6])
  c36 = findmid(c[3],c[6])
  c37 = findmid(c[3],c[6])
  c45 = findmid(c[4],c[5])
  c46 = findmid(c[4],c[6])
  c47 = findmid(c[4],c[7])
  c56 = findmid(c[5],c[6])
  c67 = findmid(c[6],c[7])

  m = [[c[0],c[1],c[5],c[6], c01, c05, c06, c15, c16, c56],
       [c[0],c[4],c[5],c[6], c04, c05, c06, c45, c46, c56],
       [c[0],c[1],c[2],c[6], c01, c02, c06, c12, c16, c26],
       [c[0],c[4],c[6],c[7], c04, c06, c07, c46, c47, c67],
       [c[0],c[2],c[3],c[6], c02, c03, c06, c23, c26, c36],
       [c[0],c[3],c[6],c[7], c03, c06, c07, c36, c37, c67]]
  return m;

def energyf (k):
  x = k[0]
  y = k[1]
  z = k[2]
  E =  x**2+y**2+z**2 
  return E;

def functionf (k

  ):
  x = k[0]
  y = k[1]
  z = k[2]
  f =  1
  return f;

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
ks = [[0,0,0],
      [1,0,0],
      [0,1,0],
      [0,0,1],
      [0.5,0,0],
      [0,0.5,0],
      [0,0,0.5],
      [0.5,0.5,0],
      [0,0.5,0.5],
      [0.5,0,0.5]]

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

def surfaces(E):
  """
                  (x2,y2)
                  / \ 
                 /   \ 
                /     \ 
               /       \   
      (x0,y0) /_________\(x1,y1)

  """

  q = np.dot(FI3,E)
  print q
  #q = np.around(q, 10)
  return q;

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
k = [[1,3,2],
      [-3,2,1],
      [-2,-3,0],
      [2,3,-1],
      [-1,2.5,1.5],
      [1.5,3,0.5],
      [-0,5,0,1],
      [-2.5,-0.5,0.5],
      [-0.5,2.5,0],
      [0,0,-0.5]]
k = ks
for i in range(0, 10):
  E[i] = energyf (k[i])
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



#q = np.around(q, 10)

#prdouce the z point and call the twoD function
#print  triangelIntegral()
e = 0.02
f = [1,1,1,1,1,1]
q = np.dot(FI3,E)

def triangelIntegralInz(z):
  k2,E2 = make2Dk(z,q)
  return twoD.triangelIntegral(e,E2,k2,f);

##just testing ................
#k2,E2 = make2Dk( 0.14917311663,q)
#k2 = np.around(k2, 10)
#E2 = np.around(E2, 10)
#print twoD.triangelIntegral(e,E2,k2,f)
#print twoD.triangelIntegral(e,E2,k2,f)

#function2d = lambda z: triangelIntegralInz(z)
#print  integrate.quad(function2d, 0, 1)
c =creatkGrid()
#print c[0]


def totalInegral():
  totalInegralSum = 0;

  kGrid ,EGrid, fGrid = creatkGrid()
  sum = 0.0

  #integrating over the  triangels
  for ix in range(0, len(kGrid) - 1):
    for iy in range (0, len(kGrid[ix]) - 1):
      for iz in range (0, len(kGrid[ix][iy]) - 1):
        k = [kGrid[ix][iy][iz], kGrid[ix + 1][iy][iz], kGrid[ix + 1][iy + 1][iz], kGrid[ix][iy + 1][iz], kGrid[ix][iy][iz +1 ], kGrid[ix + 1][iy][iz + 1], kGrid[ix + 1][iy+1][iz+1], kGrid[ix][iy + 1][iz + 1]]
        Mat = cubeToTetrahedral(k)
        for i in range(0, 6):
          E = [energyf(Mat[i][0]), energyf(Mat[i][1]), energyf(Mat[i][2]), energyf(Mat[i][3]), energyf(Mat[i][4]), energyf(Mat[i][5]), energyf(Mat[i][6]), energyf(Mat[i][7]), energyf(Mat[i][8]), energyf(Mat[i][9])]
          q = np.dot(FI3,E)
          print q
          function2d = lambda z: triangelIntegralInz(z)
          sum += integrate.quad(function2d, 0, 1)[0]
          #print sum
  return sum

print totalInegral()