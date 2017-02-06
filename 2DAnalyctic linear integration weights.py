"""
Analytic linear interpolation "G wiesenekker,G received 1987" 2D test
dimensionaless k vekort
the answer given in (1/2pi)^2 *(2m / hbar^2)^n where  n = s +1. fn(k^(2s))  
"""
import numpy as np

# testng the function
e = 1
miniX = -2
maxX = 2
miniY = -2
maxY = 2
Deltak = 1
Det = Deltak**2
FI = [[1, 0, 0 ],
      [-1, 1, 0],
      [-1, 0, 1]]
F = [[0,0],
     [1,0],
     [0,1]]
cosT = np.cos(1)
sinT = np.sin(1)

R = [[cosT,-sinT],
     [sinT,cosT]]
F = np.dot(F,R)
FF = [[1,0,0],
       [1,F[1][0],F[1][1]],
       [1,F[2][0],F[2][1]]] 
print F

FI  =  np.linalg.inv(FF)
#print nF
#F =  [[0.0, 0.0],
#     [0.5, np.sqrt(3)/2],
#     [-np.sqrt(3)/2, 5]]

#Fx = [[1,0.0, 0.0],
#     [1,0.5, np.sqrt(3)/2],
#     [1,-np.sqrt(3)/2, 5]]
#FI = np.linalg.inv(Fx)
#Det = Deltak**2* np.linalg.det(Fx)
#enter the function f here 
def functionf (x,y):
  #f = (x**2) + y**2;
  #f = 1+x+y+y**2 + y*x+x**2
  #f = y**2
  #f = 1
  #f = x*y
  #f = y #look tom
  #f =1+ y + 2*x +0.5*x**2 +0.5*y**2+ x*y
  f =  1+x**2+y**2+x+y+x*y
  #if x >0:
  #   f =   x#**2+y**2+x+y+x*y
  #else:
  #  f= y
  return f;

#enter the function E here 
def energyf (x,y):
  #E = x**2 + y**2;
  #E = 2+ x*y
  #E = 2*y + x**2
  #E =  x**2
  #E =  1 + 2*x
  #E = 0.9 + x**2
  #E = x**3
  #E = x +y
  #E = 2*x + y +0.5*x**2-0.5*y**2+ x*y
  #E = 0.5*x**2+0.7*y**2+ x*y 
  #E =  x*y + 0.5*x**2 + 0.5*y**2 + x
  #E = x*y +x**2+y**2
  #E = 2*x*y +x**2 +y**2
  #E = 2*x*y +x**2 -y**2
  #E  = -2*x*y -x**2 - y**2 + 2
  #E = 1 + 4*x + x**2
  #E = 1 + 4*x + x**2 + y**2
  #E =  0.3+ 2*x+ x**2 #look tom
  #E = 0.5*y**2 + 0.5*y
  #E = 0.5*x + 2*x**2 + -5*y**2
  #E =0.2 + -2*y  + -2*x+2*y**2
  #E =  -0.25 +x*y
  #E = 2*x+0.2*y+ y**2
  #E =  x+y
  #E = x
  E = x+2*x**2+2*y**2 +2*y
  return E;


# function creat k grid. triangel
"""  ___________
    |/\/\/\/\/\/|
    |-----------|
    |/\/\/\/\/\/|
    |-----------|
    |/\/\/\/\/\/|
    |-----------|
    |/\/\/\/\/\/|
    -------------    
"""
def creatkGrid():

  width = maxX - miniX
  height = maxY - miniY
  #only to avoid triangel dx and dy = 0! this depends on the grid we chose 
  rotation = 0.1
 
  nx = int(float(width/Deltak)) 
  ny = int(float(height/Deltak)) 


  kGrid = [[0 for i in xrange(ny +1)] for i in xrange(nx +1)]
  EGrid = [[0 for i in xrange(ny +1)] for i in xrange(nx +1)]
  fGrid = [[0 for i in xrange(ny +1)] for i in xrange(nx +1)]

  for ix in range(0,nx + 1):
    for iy in range(0,ny + 1):
      x = (miniX + Deltak *ix) * np.cos(rotation) + -(miniY + Deltak*iy) * np.sin(rotation)  
      y = (miniX + Deltak *ix) * np.sin(rotation) + (miniY + Deltak*iy) * np.cos(rotation)  
      kGrid[ix][iy] = x , y
      EGrid[ix][iy] = energyf (x,y)
      fGrid[ix][iy] = functionf(x, y)
  return kGrid, EGrid, fGrid;
#transform to an afineTriangel ref my paper()
def ABv(A,B,v):
  vn = [0,0,0]
  vn[0] = v[0]
  vn[1] = A[0][0]*v[1] + A[0][1]*v[2] + B[0]*v[0]
  vn[2] = A[1][0]*v[1] + A[1][1]*v[2] + B[1]*v[0]
  return vn

# Function sortE return E1<E2<E3
def sortE(E, f, k):
  for i in range(1, len(E)):
        j = i
        while j > 0 and E[j] < E[j-1]:
            
            f[j], f[j-1] = f[j-1], f[j]
            k[j], k[j-1] = k[j-1], k[j]
            E[j], E[j-1] = E[j-1], E[j]
            j -= 1
  return E, f, k;
def getAfineConstat(q,p,k,A,B):


	return

def getConstantsPi(E,k,f):
  # p0 + p1 x0 + p2 y0 = f0
  # p0 + p1 x1 + p2 y1 = f1
  # p0 + p1 x2 + p2 y2 = f1

  matriseA =[ [1, k[0][0], k[0][1]],
              [1, k[1][0], k[1][1]],
              [1, k[2][0], k[2][1]]]
  matriseB =[f[0],
             f[1],
             f[2]]

  p = np.linalg.solve(matriseA, matriseB)
  return p;
#ref -(6b)
def getConstantsEi(E,k):
  # q0 + q1 x0 + q2 y0 = E0
  # q0 + q1 x1 + q2 y1 = E1
  # q0 + q1 x2 + q2 y2 = E1

  matriseA =[ [1, k[0][0], k[0][1]],
              [1, k[1][0], k[1][1]],
              [1, k[2][0], k[2][1]]]
  matriseB =[E[0],
             E[1],
             E[2]]
  q = np.linalg.solve(matriseA, matriseB)
  return q;

# Function returns the integral of I(E)for E1<e<E2 ref(14)
def IEIntegral1(e, E, k):
  #kt = k[0] + (e - E[0])/(E[2]-E[0]) * (k[2] - k[0]) 
  kt = np.add(k[0], np.multiply((e - E[0])/(E[2]-E[0]), (np.subtract(k[2],k[0])))) 
  #ku = (e - E[0])/(E[1]-E[0]) * (k[1]-k[0]). - (e - E[0]) / (E[2]-E[0]) * (k[2]- k[0]) 
  ku = np.subtract(np.multiply((e-E[0])/(E[1]-E[0]),(np.subtract(k[1],k[0]))), np.multiply((e-E[0])/(E[2]-E[0]) ,np.subtract(k[2],k[0])) )
  # Function returns the integral of I(E) for E2<e<E3
  return kt, ku;
#ref (16)
def IEIntegral2(e, E, k):
  #kt = k[2] + (e - E[2])/(E[2]-E[0]) * (k[2] - k[0])
  kt = np.add(k[2], np.multiply((e - E[2])/(E[2]-E[0]), (np.subtract(k[2],k[0])))) 

  #ku = (e - E[2])/(E[2]-E[1]) * (k[2]-k[1]). - (e - E[2]) / (E[2]-E[0]) * (k[2]- k[0]) 
  ku = np.subtract(np.multiply((e-E[2])/(E[2]-E[1]),(np.subtract(k[2],k[1]))), np.multiply((e-E[2])/(E[2]-E[0]) ,np.subtract(k[2],k[0])) )
  # Function returns the integral of I(E) for E2<e<E3
  return kt, ku;

#ref(14)
def determinantJacobian1(e, E, k):
  # the Jacobian (e - E[0])/((E[1]-E[0])*(E[2]-E[0])) * [[k20 - k00]  , [k10 - k00]]
  #                                                     [[k21 - k01]  , [k11 - k01]] 
  j =  [[k[2][0] - k[0][0], k[1][0] - k[0][0]],
       [k[2][1] - k[0][1], k[1][1] - k[0][1]]]
  j = np.linalg.det(j)
  if j < 0 :
    j = j * -1
  return (e - E[0])/((E[1]-E[0])*(E[2]-E[0])) * j;

#ref(17)
def determinantJacobian2(e, E, k):
  # the Jacobian (E[2] - e)/((E[2]-E[1])*(E[2]-E[0])) * [[k20 - k10]  , [k20 - k00]]
  #                                                     [[k21 - k11]  , [k21 - k01]] 
  j = [[k[2][0] - k[1][0], k[2][0] - k[0][0]],
      [k[2][1] - k[1][1], k[2][1] - k[0][1]]]
  j = np.linalg.det(j)
  if j < 0 :
    j = j * -1
  return (E[2] - e)/((E[2]-E[1])*(E[2]-E[0])) * j;

#ref(19)
def triangelIntegral(e,E,k,f):
  D = [[1,k[0][0],k[0][1]],
       [1,k[1][0],k[1][1]],
       [1,k[2][0],k[2][1]]]
  E, f, k = sortE(E, f, k);
  k = F
  det =  Det

  if E[0] <= e and e < E[1]:
    
    kt, ku = IEIntegral1(e, E, k);
    Jacobian = determinantJacobian1(e,E,k)
  elif E[1] <= e and e < E[2]:

    kt, ku = IEIntegral2(e, E, k);
    Jacobian = determinantJacobian2(e,E,k)   
  else:
    kt = 0,0
    ku = 0,0
    Jacobian = 0

  v = [0,0,0]
  v[0] =  Jacobian 
  v[1] =  Jacobian*(kt[0]+ 0.5*ku[0])
  v[2] =  Jacobian * (kt[1] + 0.5 * ku[1]) 
  vn = np.dot(v,FI)

  sum = f[0] * vn[0] + f[1]*vn[1]+ f[2] *vn[2]
  sum = sum * det
  return sum;

def totalInegral():
  totalInegralSum = 0;

  #creating the k grid and calc the energy 
  kGrid ,EGrid, fGrid = creatkGrid()
  numrows = len(kGrid)
  numcols = len(kGrid)
  sum = 0.0
  #integrating over the  triangels
  for ix in range(0, len(kGrid) - 1):
    for iy in range (0, len(kGrid[ix]) - 1): # q and p must be here 

      #summing over the odd triangel
      E = [EGrid[ix][iy], EGrid[ix + 1][iy], EGrid[ix][iy + 1]]
      k = [kGrid[ix][iy], kGrid[ix + 1][iy], kGrid[ix][iy + 1]]
      f = [fGrid[ix][iy], fGrid[ix + 1][iy], fGrid[ix][iy + 1]]
      sum += triangelIntegral(e,E,k,f)
      #summing over the partall triangel
      E = [EGrid[ix + 1][iy + 1], EGrid[ix + 1][iy], EGrid[ix ][iy + 1]]
      k = [kGrid[ix + 1][iy + 1], kGrid[ix + 1][iy], kGrid[ix ][iy + 1]]
      f = [fGrid[ix + 1][iy + 1], fGrid[ix + 1][iy], fGrid[ix ][iy + 1]]
      sum += triangelIntegral(e,E,k,f)
  return sum

print totalInegral()
