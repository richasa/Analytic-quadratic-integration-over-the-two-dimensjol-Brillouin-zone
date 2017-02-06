""" 
Analytic qudtratic interpolation "G wiesenekker,G received 1987" 2D test
dimensionaless k vekort
the answer given in (1/2pi)^2 *(2m / hbar^2)^n where  n = s +1. fn(k^(2s))  
Redistributions of source code must retain the above copyright notice, 
"""

import numpy as np
import math
# testng the function
#initializing the integrating grid, e must be an array, but here we used 1 point test. Deltak is the triangel width and height.  
e = 1
miniX = -4
maxX = 4
miniY = -4
maxY = 4
Deltak =1
Det = (2*Deltak)**2
pi = 3.14159265359
epsilion = 1*10**-9
FI =[[1,0,0,0,0,0],
     [-3,-1,0,4,0,0],
     [-3,0,-1,0,0,4],
     [2,2,0,-4,0,0],
     [4,0,0,-4,4,-4],
     [2,0,2,0,0,-4]]

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

#just apply any desired funksjon f.
def functionf (x,y):
  f =  1+ x+3*x-5*y+7*x**2+11*y**2
  return f;

#also any desired funksjon E.
def energyf (x,y,z):
  #E = x**2 + y**2;
  E =  4*x**2
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

  #rotation grid was a nessery step in testing. but not any more. 
  rotation = 0.00

  kGrid = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
  EGrid = [[energyf (KGrid[0])],[energyf (KGrid[1])],[energyf (KGrid[2])],[energyf (KGrid[3])]] 
  fGrid = [[functionf (KGrid[0])],[functionf(KGrid[1])],[functionf(KGrid[2])],[functionf (KGrid[3])]]

  return kGrid, EGrid, fGrid;

# Function sortu return u0<u1<u2<u3<u4<u5
def sortU(u):
  for i in range(1, len(u)):
        j = i
        while j > 0 and u[j] < u[j-1]:
            
            u[j], u[j-1] = u[j-1], u[j]
            j -= 1
  return u;
#ref (46) integration over the affine transformation. with respect to  A and B matrixes. 
def ABv(A,B,v):
  vn = [0,0,0,0,0,0]
  vn[0] = v[0]
  vn[1] = A[0][0]*v[1] + A[0][1]*v[2] + B[0]*v[0]
  vn[2] = A[1][0]*v[1] + A[1][1]*v[2] + B[1]*v[0]
  vn[3] = (A[0][0]**2)*v[3] + (A[0][1]**2)*v[5] + (B[0]**2)*v[0] + 2*(A[0][0]*A[0][1]*v[4] + A[0][0]*B[0]*v[1] + A[0][1]*B[0]*v[2])
  vn[4] = A[0][0]*A[1][0]*v[3] + A[0][1]*A[1][1]*v[5] + B[0]*B[1]*v[0] +(A[0][0]*A[1][1] + A[1][0]*A[0][1])*v[4] + (A[0][0]*B[1] + A[1][0]*B[0])*v[1] + (A[0][1]*B[1] + A[1][1]*B[0])*v[2]
  vn[5] = (A[1][0]**2)*v[3] + (A[1][1]**2)*v[5] + (B[1]**2)*v[0] + 2*(A[1][0]*A[1][1]*v[4] + A[1][0]*B[1]*v[1] + A[1][1]*B[1]*v[2])
  return vn


#ref (20 a) this is not necessary and not used function, since we used weight method for integration I = f(k)*w 
#this was for function was for testing where intergration  I = pi*vi  
def getConstantsAB(k,corners):
  matriseA =[ [k[0][0],k[0][1],1],
              [k[1][0],k[1][1],1],
              [k[2][0],k[2][1],1]]

  matriseB1 =[0,1,0]
  matriseB2 =[0,0,1]

  A1 = np.linalg.solve(matriseA, matriseB1)
  A2 = np.linalg.solve(matriseA, matriseB2)
  #print A1,A2
  A = [[A1[0],A1[1]],
       [A2[0],A2[1]]]
  B = [-A1[2],-A2[2]] #---------------------------maybe minus
  AI = np.linalg.inv(A)
  #for i in range(0, 6):
  #  k[i] = k[i][0]*A[0][0] + k[i][1]*A[0][1] -B[0], k[i][0]*A[1][0] + k[i][1]*A[1][1] -B[1]
  #print k
  return AI,B;

#ref -(20b) here is where we find which energy surface we are dealing with. all the six point is required to solv the 6 unkown problem [q1,q2,q3,q4,q5,q6].
def getConstantsEi(E,k):
  # q0 + q1 x0 + q2 y0 + q3 X0^2 + q4 x0y0 + q5 y0^2 = E0
  # q1 + q1 x1 + q2 y1 + q3 X1^2 + q4 x1y1 + q5 y1^2 = E1
  # q2 + q1 x2 + q2 y2 + q3 X2^2 + q4 x2y2 + q5 y2^2 = E2
  # q3 + q1 x3 + q2 y3 + q3 X3^2 + q4 x3y3 + q5 y3^2 = E3
  # q4 + q1 x4 + q2 y4 + q3 X4^2 + q4 x4y4 + q5 y4^2 = E4
  # q5 + q1 x5 + q2 y5 + q3 X5^2 + q4 x5y5 + q5 y5^2 = E5
  matriseA = [[1, k[0][0], k[0][1], k[0][0]**2, k[0][0]*k[0][1], k[0][1]**2 ],
              [1, k[1][0], k[1][1], k[1][0]**2, k[1][0]*k[1][1], k[1][1]**2 ],
              [1, k[2][0], k[2][1], k[2][0]**2, k[2][0]*k[2][1], k[2][1]**2 ],
              [1, k[3][0], k[3][1], k[3][0]**2, k[3][0]*k[3][1], k[3][1]**2 ],
              [1, k[4][0], k[4][1], k[4][0]**2, k[4][0]*k[4][1], k[4][1]**2 ],
              [1, k[5][0], k[5][1], k[5][0]**2, k[5][0]*k[5][1], k[5][1]**2 ]]

  matriseB =[E[0],
             E[1],
             E[2],
             E[3],
             E[4],
             E[5]]
  q = np.linalg.solve(matriseA, matriseB)
  return q;
#ref -(not createdyet) here is where we run all the Affine transformationes
def getAfineConstat(q,k,A,B):
  #ref A1
  if (q[4] != 0):
    s = [[q[3] , 0.5*q[4]],
        [0.5*q[4], q[5]]]

    #A here the eigen vector which diagonal matrix s / we multyplay the second col with -1 to make nA matrix is also the inverse matix itself. 
    eigvals, nA = np.linalg.eig(s)
    nA = [[ nA[0][0] , -nA[0][1]],
          [nA[1,0], -nA[1][1]]]

    q1 = q[1]*nA[0][0] + q[2]*nA[0][1]
    q2 = q[1]*nA[1][0] + q[2]*nA[1][1]
    q3 = q[3]*nA[0][0]**2 + q[4]*nA[0][0]*nA[1][0] + q[5]*nA[1][0]**2
    q5 = q[3]*nA[0][1]**2 + q[4]*nA[0][1]*nA[1][1] + q[5]*nA[1][1]**2

    q = [q[0],q1,q2,q3,0,q5] 
    q = np.around(q, 14)
    #update k after the first Affine transformation.
    for i in range(0, 6):
      k[i] = k[i][0]*nA[0][0] + k[i][1]*nA[0][1] , k[i][0]*nA[1][0] + k[i][1]*nA[1][1]
    A = nA

    #ref (A.19)
  if(q[5] != 0 and q[2] != 0):

    dB = 0 , - q[2]/(2*q[5])
    B = B[0] +dB[0]*A[0][0]+dB[1]*A[0][1], B[1]+dB[0]*A[1][0]+dB[1]*A[1][1]
    #B = np.dot(B,nA)
    for i in range(0, 6):
      k[i] = k[i][0] - dB[0] , k[i][1] -dB[1]  
    q0 = q[0] - (q[2]**2 /(4*q[5]))
    q = [q0,q[1],0,q[3],q[4],q[5]]
  if (q[3] == 0 and q[5] != 0):
    #"x = y' ; y = x' 
    nA = [[0,1],
          [1,0]]
    for i in range(0, 6):
      k[i] = k[i][0]*nA[0][0] + k[i][1]*nA[0][1] , k[i][0]*nA[1][0] + k[i][1]*nA[1][1]
    #B = B[0]*nA[0][0] + B[1]*nA[0][1] , B[0]*nA[1][0] + B[1]*nA[1][1]
    q = [q[0],q[2],q[1],q[5],q[4],0]
    #update A
    A = np.dot(A,nA)
  #ref (A.17)
  if(q[3] !=0 and q[1] != 0 ):
    # "x = x' -q2/2q4   y=y"
    #update B,q
    dB = - q[1]/(2*q[3]) , 0 
    B = B[0] +dB[0]*A[0][0]+dB[1]*A[0][1], B[1]+dB[0]*A[1][0]+dB[1]*A[1][1]
    for i in range(0, 6):
      k[i] = k[i][0] -dB[0] , k[i][1] - dB[1]

    q0 = q[0] - (q[1]**2 /(4*q[3]))
    q = [q0,0,q[2],q[3],q[4],q[5]]

 
  #ref (A1.10)
  #Here we have only 3 forms
  # q0 + q3 + q5 (q3 != 0  and q5 !=0 )
  # q0 + q2 + a3 (q3 != 0, q2 may be 0) 
  # q0 +  q1 + q2 (q1 and (or) q2 may be 0)

  #ref A1.11-14
  #start returning surfaces.

  if (q[3] != 0 and q[5] != 0):
    if (q[3]*q[5] > 0):
      #"return elipse + return"
      return q, k, A, B

    else:

      #"return hyperbola A1.12 + return"
      nA = [[1/np.sqrt(2) , np.sqrt(abs(q[5]) /abs(2*q[3])) ],
           [ np.sqrt(abs(q[3]) /abs(2*q[5])) , -1/np.sqrt(2)]]
      nAI = np.linalg.inv(nA)

      #update k,B,A,q ref (A1.12)
      for i in range(0, 6):
        k[i] = (k[i][0])*nA[0][0] + k[i][1]*nA[0][1] , k[i][0]*nA[1][0] + k[i][1]*nA[1][1] 
      A = np.dot(A,nA)

      q4 = q[3]*np.sqrt(abs(q[5]/abs(q[3]))) - q[5]*np.sqrt(abs(q[3]/abs(q[5])))
      q = [q[0],0,0,0,q4,0]
      return q, k, A, B

  if (q[3] != 0):
    #print "paralel lines+ return" , q[3]
    return q, k, A, B
  if (q[1] != 0 or q[2] != 0):
    #"stright line q2' != 0 and q3 = 0 A1.13 and 14"
    cosT = q[1]/(q[1]**2 + q[2]**2)**0.5
    sinT = q[2]/(q[1]**2 + q[2]**2)**0.5
    nA = [[cosT,sinT],
        [sinT,-cosT]]
    #print A
    #update k,A,B,q
    for i in range(0, 6):
      k[i] = k[i][0]*nA[0][0] + k[i][1]*nA[0][1] , k[i][0]*nA[1][0] + k[i][1]*nA[1][1]
    B = B[0]*nA[0][0] + B[1]*nA[0][1] , B[0]*nA[1][0] + B[1]*nA[1][1]
    #print "reach"
    q1 = cosT*q[1] + sinT*q[2]
    q = [q[0],q1,0,0,0,0]
    A = np.dot(A,nA)
    return q, k, A, B
  return q, k, A, B


"""------------------------------------------Surface integration functions--------------------------------------------------------"""
def elipse (e, q, corners,dx,dy,nx,ny):
  print "reach"
  c = [dx[0] * corners[0][1] - dy[0]*corners[0][0],
       dx[1] * corners[1][1] - dy[1]*corners[1][0],
       dx[2] * corners[2][1] - dy[2]*corners[2][0]]
  #print c
  constantA = ((e - q[0]) / q[3]) **0.5
  constantB =  ((e - q[0]) / q[5]) **0.5
  #ref 25 my papper.   
  b = [nx[0]*constantA,
       nx[1]*constantA,
       nx[2]*constantA]

  a = [ny[0]*constantB,
       ny[1]*constantB,
       ny[2]*constantB]

  #ref 33 my papper
  t = [0,0,0,0,0,0]
  w = [0,0,0,0,0,0]
  u = [0,0,0,0,0,0]
  tStatus = [False, False, False, False, False, False]
  z = 0; 
  for i in range(0, 3):
    if (4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i])) >= 0 :
      if abs(b[i]+c[i] ) != 0:
        u[z] = 2*np.arctan((2*a[i] + np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) / (2*(b[i] + c [i])))
        u[z+1] =2*np.arctan( (2*a[i] - np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) / (2*(b[i] + c [i])))
        tStatus[z] = True
        tStatus[z+1] = True
      else:
        if a[i] != 0:
          if abs(2*a[i] + np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) > 0:
            u[z+1] = 2*np.arctan((c[i]-b[i])/(2*a[i]))
            tStatus[z+1] = True
          elif abs(2*a[i] - np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) > 0:
            u[z] = 2*np.arctan((c[i]-b[i])/(2*a[i]))
            tStatus[z] = True
        if  (2*a[i] + np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i])))  != 0:
          if (2*a[i] + np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) > 0:
            u[z] = pi
          else:
            u[z] = -pi
          tStatus[z] = True
        if (2*a[i] - np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i])))  != 0:
          if  (2*a[i] - np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) > 0:
            u[z+1] = pi
          else :
            u[z+1] = -pi

          tStatus[z+1] = True

    z += 2
  z = 0
  n = 0
  for i in range (0,6):
    if tStatus[i] == True:
      point = constantA*np.cos(u[i]),constantB*np.sin(u[i])
      if(inTriangel(corners,point)):
        w[n] = u[i]
        n += 1
  uf = []
  w1 = [w[i] for i in xrange(n)]
  w1 = sortU(w1)
  z = 0
  #print w1
  for i in range(0, n):
    if (i + 1 < n ) and (w1[i] != w1[i + 1]):
      point = constantA*np.cos( w1[i]+0.0001), constantB*np.sin(w1[i]+0.0001)
      if (inTriangel(corners,point)):
        uf.append( w1[i])
        uf.append( w1[i+1])
        z += 2 
    elif(i < n) and (w1[i] != w1[0]):
      point = constantA*np.cos(w1[0]-0.0001), constantB*np.sin(w1[0]-0.0001)
      if (inTriangel(corners,point)):
        uf.append( w1[n - 1] ) 
        uf.append( w1[0] + 2*pi)
        z += 2
  if z == 0:
    return [0,0,0,0,0,0]
  #ref (27)
  Jacobian = ( 1.0 / (2*(q[3]*q[5])**0.5))
  sum1 = 0
  v = [0,0,0,0,0,0]

  for i in range(0, z-1, 2):
    #ref (28) 
    v[0] +=  Jacobian * (uf[i+1] - uf[i])
    v[1] +=  Jacobian * ((e - q[0]) / q[3]) ** 0.5 * (np.sin(uf[i+1])-(np.sin(uf[i])))
    v[2] +=  Jacobian * ((e - q[0]) / q[5]) ** 0.5 * (np.cos(uf[i]) - np.cos(uf[i+1]))
    v[3] +=  Jacobian * abs((e - q[0]) / q[3]) * (((0.5 * uf[i+1] + 0.25 * np.sin(2*uf[i+1]))-(0.5 * uf[i] + 0.25 * np.sin(2*uf[i]))))
    v[4] +=  Jacobian * ((e - q[0]) / q[3]) **0.5 * ((e - q[0]) / q[5]) **0.5 *(( -0.25 * np.cos(2*uf[i+1]))-(-0.25 * np.cos(2*uf[i])))
    v[5] +=  Jacobian * abs((e - q[0]) / q[5]) * ((0.5 * uf[i+1] - 0.25 * np.sin(2*uf[i+1]))-(0.5 * uf[i] - 0.25 * np.sin(2*uf[i])))
  return v

"""--------------------------------------------------------------------------------------------------------------------------"""
def hyperbola (e, q, corners,dx,dy,nx,ny):
  #print A, B,q
  #ref 36-37 my papper    
  a = [nx[0],
       nx[1],
       nx[2]]
  
  b = [dy[0]*corners[0][0] - dx[0] * corners[0][1],
       dy[1]*corners[1][0] - dx[1] * corners[1][1],
       dy[2]*corners[2][0] - dx[2] * corners[2][1]]

  constantA = ((e - q[0]) / q[4])
  c = [ny[0]*constantA,
       ny[1]*constantA,
       ny[2]*constantA]

  #ref(38) my paper
  u = [0,0,0,0,0,0]
  w = [0,0,0,0,0,0]
  uStatus = [False, False, False, False, False, False]
  z = 0;
  for i in range(0, 3):
    if (b[i]**2 - 4*a[i]*c[i]) > 0 :
      if a[i] != 0:
        u[z] = ((-b[i] + np.sqrt(b[i]**2 - 4*a[i]*c[i]))/ (2*a[i]))
        u[z+1] = ((-b[i] - np.sqrt(b[i]**2 - 4*a[i]*c[i])) / (2*a[i]))
        uStatus[z] = True
        uStatus[z+1] = True 
      if ( b[i] != 0 and a[i] == 0):
        u[z] = -c[i]/b[i]
        uStatus[z] = True
       
    z += 2
  n = 0
  for i in range (0,6):
    if uStatus[i] == True:
      if u[i] != 0:
        point = u[i], constantA/(u[i])
        if(inTriangel(corners,point)):
          w[n] = u[i]   
          n += 1
  uf = []
  w1 = [w[i] for i in xrange(n)]
  w1 = sortU(w1)

  z = 0
  for i in range(0, n):
    if (i + 1 < n ):
      point = ( w1[i]+0.000001), constantA/(w1[i]+0.000001)
      if (inTriangel(corners,point)):
        uf.append( w1[i])
        uf.append( w1[i+1])
        z += 2 
    elif(i < n):
      point = ( w1[0]-0.000001), constantA/(w1[0]-0.000001)
      if (inTriangel(corners,point)):
        uf.append( w1[n - 1])
        uf.append( w1[0])
        z += 2
  if z == 0:
    return [0,0,0,0,0,0];
  #ref (29)
  sum1 = 0
  #Jacobian = ( 1.0 / abs(q[4]*uf))

  v = [0,0,0,0,0,0]
  for i in range(0, z-1, 2):
    if (uf[i] <= 0 and uf[i+1] <= 0):
      v[0] -= (1/abs(q[4])) * (np.log(abs(uf[i+1])) - np.log(abs(uf[i])))
      v[1] -= (1/abs(q[4])) * (uf[i+1] - uf[i])  
      v[2] -= (1/abs(q[4])) * constantA*(-1/uf[i+1] + 1/uf[i]) 
      v[3] -= (1/abs(q[4])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2) 
      v[4] -= (1/abs(q[4])) * constantA*(np.log(abs(uf[i+1])) - np.log(abs(uf[i]))) 
      v[5] -= (1/abs(q[4])) * ((e-q[0])**2 / q[4]**2) * (-1/(2*uf[i+1]**2) + 1/(2*uf[i]**2))

      #ref (28)
    else:
      v[0] += (1/abs(q[4])) * (np.log(abs(uf[i+1])) - np.log(abs(uf[i])))
      v[1] += (1/abs(q[4])) * (uf[i+1] - uf[i])  
      v[2] += (1/abs(q[4])) * constantA*(-1/uf[i+1] + 1/uf[i]) 
      v[3] += (1/abs(q[4])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2) 
      v[4] += (1/abs(q[4])) * constantA*(np.log(abs(uf[i+1])) - np.log(abs(uf[i]))) 
      v[5] += (1/abs(q[4])) * ((e-q[0])**2 / q[4]**2) * (-1/(2*uf[i+1]**2) + 1/(2*uf[i]**2))

      #ref (28)
  #ref(46)
  return v;
"""--------------------------------------------------------------------------------------------------------------------------"""
def parabola(e, q, corners,dx,dy,nx,ny):
  #ref 36-37 my papper    
  a = [-ny[0]*(q[3]/q[2]),
       -ny[1]*(q[3]/q[2]),
       -ny[2]*(q[3]/q[2])]
  
  b = [nx[0],
       nx[1],
       nx[2]]

  #dxj*yj - dyj*xj     
  cT = [dx[0] * corners[0][1] - dy[0]*corners[0][0],
       dx[1] * corners[1][1] - dy[1]*corners[1][0],
       dx[2] * corners[2][1] - dy[2]*corners[2][0]]

  constantA = ((e - q[0]) / q[2])
  c = [ny[0]*constantA - cT[0],
       ny[1]*constantA - cT[1],
       ny[2]*constantA - cT[2]]

  #ref(42) my paper
  u = [0,0,0,0,0,0]
  w = [0,0,0,0,0,0]
  uStatus = [False, False, False, False, False, False]
  z = 0;
  for i in range(0, 3):
    if (b[i]**2 - 4*a[i]*c[i]) > 0 :
      if a[i] != 0:
        u[z] = ((-b[i] + np.sqrt(b[i]**2 - 4*a[i]*c[i]))/ (2*a[i]))
        u[z+1] = ((-b[i] - np.sqrt(b[i]**2 - 4*a[i]*c[i])) / (2*a[i]))
        uStatus[z] = True
        uStatus[z+1] = True
      if ( b[i] != 0 and a[i] == 0):
        u[z] = -c[i]/b[i]
        uStatus[z] = True
        #print "rech", u[z]
    z += 2
  #chicking if the intersection point is actuly inside the triangle.
  n = 0
  for i in range (0,6):
    if uStatus[i] == True:
      point = u[i], constantA - (q[3]/q[2])*u[i]**2
      if(inTriangel(corners,point)):
        #print point
        w[n] = u[i]   
        n += 1
  uf = []
  w1 = [w[i] for i in xrange(n)]
  w1 = sortU(w1)
  z = 0
  for i in range(0, n):
    if (i + 1 < n ):
      point = ( w1[i]+0.0001),constantA - (q[3]/q[2])*(w1[i]+0.0001)**2
      if (inTriangel(corners,point)):
        #print point
        uf.append( w1[i])
        uf.append( w1[i+1])
        z += 2 
    elif(i < n):
      point = ( w1[0]-0.0001), constantA - (q[3]/q[2])*(w1[0]-0.0001)**2
      if (inTriangel(corners,point)):
        uf.append( w1[n - 1])
        uf.append( w1[0])
        z += 2
  if z == 0:
    return [0,0,0,0,0,0];
  #ref (32)
  sum1 = 0
  v =[0,0,0,0,0,0]
  #print p
  for i in range(0, z-1, 2):
    #ref (28)
    v[0] += (1/abs(q[2])) * (uf[i+1] - uf[i])
    v[1] += (1/abs(q[2])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)  
    v[2] += (1/abs(q[2])) * (constantA*uf[i+1]- (1.0/3)*(q[3]/q[2]*uf[i+1]**3) - (constantA*uf[i]- (1.0/3)*(q[3]/q[2]*uf[i]**3)))
    v[3] += (1/abs(q[2])) * ((1.0/3)*uf[i+1]**3 - (1.0/3)*uf[i]**3 )
    v[4] += (1/abs(q[2])) * (((0.5)*constantA * uf[i+1]**2 - (1.0/4)*(q[3]/q[2])*uf[i+1]**4) - ((0.5)*constantA * uf[i]**2 - (1.0/4)*(q[3]/q[2])*uf[i]**4))
    v[5] += (1/abs(q[2])) * ( ( ((e-q[0])**2 /(q[2]**2))*uf[i+1]-(2.0/3)*(e-q[0])/q[2]*(q[3]/q[2])*uf[i+1]**3+(1.0/5)*(q[3]**2/q[2]**2)*uf[i+1]**5 )-( ((e-q[0])**2 /(q[2]**2))*uf[i]-(2.0/3)*(e-q[0])/q[2]*(q[3]/q[2])*uf[i]**3+(1.0/5)*(q[3]**2/q[2]**2)*uf[i]**5 ))
  #ref(46)
  #vn = ABv(A,B,v)
  return v;
"""-------------------------------------------------------------------------------------------------------------------------------"""
def straightLine(e, q, corners,dx,dy,nx,ny):
  #shift if e is on the side of the triangel
  b = [ny[0],
       ny[1],
       ny[2]]

  #dxj*yj - dyj*xj     
  cT = [dx[0] * corners[0][1] - dy[0]*corners[0][0],
       dx[1] * corners[1][1] - dy[1]*corners[1][0],
       dx[2] * corners[2][1] - dy[2]*corners[2][0]]

  constantA = ((e - q[0]) / q[1])
  c = [nx[0]*constantA - cT[0],
       nx[1]*constantA - cT[1],
       nx[2]*constantA - cT[2]]

  #ref(42) my paper
  u = [0,0,0]
  w = [0,0,0]
  uStatus = [False, False, False]
  z = 0;
  for i in range(0, 3):
    #print b[i]
    if b[i] != 0:
      #print b[i]
      u[i] = -c[i]/b[i]
      uStatus[i] = True

  #chicking if the intersection point is actuly inside the triangle.
  n = 0
  for i in range (0,3):
    if uStatus[i] == True:
      point = constantA, u[i]
      if(inTriangel(corners,point)):
        w[n] = u[i]   
        n += 1
  uf = []
  w1 = [w[i] for i in xrange(n)]
  w1 = sortU(w1)
  z = 0
  for i in range(0, n):
    if (i + 1 < n ):
      point = constantA, (w1[i]+0.0001)
      if (inTriangel(corners,point)):
        #print point
        uf.append( w1[i])
        uf.append( w1[i+1])
        z += 2 
    elif(i < n):
      point = constantA, (w1[0]-0.0001)
      if (inTriangel(corners,point)):
        uf.append( w1[n - 1])
        uf.append( w1[0])
        z += 2
  if z == 0:
    return [0,0,0,0,0,0];
  #ref (35)
  sum1 = 0
  v = [0,0,0,0,0,0]

  for i in range(0, z-1, 2):
    shift = inTriangelStraight(corners,(constantA,(uf[i]+0.0001)))
    v[0] += shift*(1/abs(q[1])) * (uf[i+1] - uf[i])
    v[1] += shift*(1/abs(q[1])) * constantA * (uf[i+1] - uf[i])
    v[2] += shift*(1/abs(q[1])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)
    v[3] += shift*(1/abs(q[1])) * constantA**2 * (uf[i+1] - uf[i])
    v[4] += shift*(1/abs(q[1])) * constantA * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)
    v[5] += shift*(1/abs(q[1])) * (1.0/3 * uf[i+1]**3 - 1.0/3 * uf[i]**3)
  #ref(46)
  return v;

"""-------------------------------------------------------------------------------------------------------------------------------"""
def degenerat(e, q, corners,dx,dy,nx,ny):
  shift = 1
  b = [ny[0],
       ny[1],
       ny[2]]

  #dxj*yj - dyj*xj     
  cT = [dx[0] * corners[0][1] - dy[0]*corners[0][0],
       dx[1] * corners[1][1] - dy[1]*corners[1][0],
       dx[2] * corners[2][1] - dy[2]*corners[2][0]]

  constantAP = ((e - q[0]) / (q[3]*(e-q[0]))**0.5)
  constantAM = -constantAP 
  cP = [nx[0]*constantAP - cT[0],
       nx[1]*constantAP - cT[1],
       nx[2]*constantAP - cT[2]]
  cM = [nx[0]*constantAM - cT[0],
       nx[1]*constantAM - cT[1],
       nx[2]*constantAM - cT[2]]

  #ref(42) my paper
  uP = [0,0,0]
  wP = [0,0,0]
  #ref(42) my paper
  uM = [0,0,0]
  wM = [0,0,0]
  uStatusP = [False, False, False]
  uStatusM = [False,False, False]
  z = 0;
  for i in range(0, 3):
    if b[i] != 0:
      uP[i] = -cP[i]/b[i]
      uM[i] = -cM[i]/b[i]
      uStatusP[i] = True
      uStatusM[i] = True
  #chicking if the intersection point is actuly inside the triangle.
  nP = 0
  for i in range (0,3):
    if uStatusP[i] == True:
      point = constantAP, uP[i]
      if(inTriangel(corners,point)):
        #print point
        wP[nP] = uP[i]   
        nP += 1
  nM = 0
  for i in range (0,3):
    if uStatusM[i] == True:
      point = constantAM, uM[i]
      if(inTriangel(corners,point)):
        #print point
        wM[nM] = uM[i]   
        nM += 1

  ufP = []
  w1P = [wP[i] for i in xrange(nP)]
  w1P = sortU(w1P)
  ufM = []
  w1M = [wM[i] for i in xrange(nM)]
  w1M = sortU(w1M)
  zP = 0
  for i in range(0, nP):
    if (i + 1 < nP ):
      point = constantAP, (w1P[i]+0.0001)
      if (inTriangel(corners,point)):
        #print point
        ufP.append( w1P[i])
        ufP.append( w1P[i+1])
        zP += 2 
    elif(i < nP):
      point = constantAP, (w1P[0]-0.0001)
      if (inTriangel(corners,point)):
        ufP.append( w1P[nP - 1])
        ufP.append( w1P[0])
        zP += 2
  zM = 0
  for i in range(0, nM):
    if (i + 1 < nM ):
      point = constantAM, (w1M[i]+0.0001)
      if (inTriangel(corners,point)):
        #print point
        ufM.append( w1M[i])
        ufM.append( w1M[i+1])
        zM += 2 
    elif(i < nM):
      point = constantAM, (w1M[0]-0.0001)
      if (inTriangel(corners,point)):
        ufM.append( w1M[nM - 1])
        ufM.append( w1M[0])
        zM += 2
  if zP ==  0 and zM == 0:
    return [0,0,0,0,0,0];
  #ref (40)/ my paper (47- 60)
  sum1P = 0
  j = 0.5/(q[3]*(e-q[0]))**0.5
  vp = [0,0,0,0,0,0]
  for i in range(0, zP-1, 2):
    shift = inTriangelStraight(corners,(constantAP,(ufP[i]+0.0001)))
    vp[0] += shift*j * (ufP[i+1] - ufP[i])
    vp[1] += shift*j * constantAP * (ufP[i+1] - ufP[i])
    vp[2] += shift*j * (0.5*ufP[i+1]**2 - 0.5*ufP[i]**2)
    vp[3] += shift*j * constantAP**2 * (ufP[i+1] - ufP[i])
    vp[4] += shift*j * constantAP * (0.5*ufP[i+1]**2 - 0.5*ufP[i]**2)
    vp[5] += shift*j * (1.0/3 * ufP[i+1]**3 - 1.0/3 * ufP[i]**3)

  sum1M = 0
  j = 0.5/(q[3]*(e-q[0]))**0.5
  
  vm = [0,0,0,0,0,0]
  for i in range(0, zM-1, 2):  
    shift = inTriangelStraight(corners,(constantAM,(ufM[i]+0.0001)))
    vm[0] += shift*j * (ufM[i+1] - ufM[i])
    vm[1] += shift*j * constantAM * (ufM[i+1] - ufM[i])
    vm[2] += shift*j * (0.5*ufM[i+1]**2 - 0.5*ufM[i]**2)
    vm[3] += shift*j * constantAM**2 * (ufM[i+1] - ufM[i])
    vm[4] += shift*j* constantAM * (0.5*ufM[i+1]**2 - 0.5*ufM[i]**2)
    vm[5] += shift*j * (1.0/3 * ufM[i+1]**3 - 1.0/3 * ufM[i]**3)
  
  v = [0,0,0,0,0,0]
  for i in range(0, 6):
    v[i] = (vp[i]+vm[i])
  return v;

"""-------------------------------------------------(linear test)-------------------------------------------------------------------------------"""
# Function check if a point p inside the triangel! 
#p = p0 + (c1 - c0) * s + (c2 - c0) * t
#The point p is inside the triangle if 0 <= s <= 1 and 0 <= t <= 1 and s + t <= 1.
def inTriangel(c,p):

  A = [[c[1][0] - c[0][0] , c[2][0] - c[0][0]],
      [c[1][1] - c[0][1] , c[2][1] - c[0][1]]]

  B = [p[0] - c[0][0],
     p[1] - c[0][1] ]
  #print A,B,"look here"
  st = np.linalg.solve(A, B)
  st = np.around(st, 10)

  if  ( 0 <= st[0] <=  1) and  ( 0 <= st[1] <=  1):
    if st[0] + st[1] <= 1:
      #print "status True : ",st
      return True


  #print "status false : ",st
  return False;

def inTriangelStraight(c,p):

  A = [[c[1][0] - c[0][0] , c[2][0] - c[0][0]],
      [c[1][1] - c[0][1] , c[2][1] - c[0][1]]]

  B = [p[0] - c[0][0],
     p[1] - c[0][1] ]
  #print A,B,"look here"
  st = np.linalg.solve(A, B)
  st = np.around(st, 10)
  if  ( 0 <= st[0] <=  1) and  ( 0 <= st[1] <=  1):
    if st[0]  == 0 or st[1] == 0:
      #print "status True : ",st
      return 0.5
    else :
      return 1

"""-------------------------------------------------------------------------------------------------------------------"""

# Function returns the integral of I(E)for a surface " Ellipse, Hyperbola, Parabola, Stright line, Degenerate"
def surfaces(e, corners,k,E):
  """
                  (x2,y2)
                  / \ 
                 /   \ 
                /     \ 
               /       \   
      (x0,y0) /_________\(x1,y1)

  """
  A1,B1 = getConstantsAB(k,corners)
  k = [[0,0],
      [1,0],
      [0,1],
      [0.5,0],
      [0.5,0.5],
      [0,0.5]]
  q = np.dot(FI,E)
  q = np.around(q, 10)
  #print ("q",q)

  A = [[1,0],
       [0,1]]
  B = [0,0]
 
  q, k, A, B = getAfineConstat(q,k,A,B)


  k = np.around(k, 10)
  A = np.around(A, 10)
  B = np.around(B, 10)
  q = np.around(q, 10)

  corners = [k[0], k[1], k[2]]

  dx = [corners[1][0] - corners[0][0],
          corners[2][0] - corners[1][0],
          corners[0][0] - corners[2][0]]

  dy = [corners[1][1] - corners[0][1],
          corners[2][1] - corners[1][1],
          corners[0][1] - corners[2][1]]

  nx = [-dy[0],
        -dy[1],
          -dy[2]]

  ny = [dx[0],
          dx[1],
          dx[2]]
  #ref (25) various forms for the qudratis surface.
  #elipse,hyperbola,parabola,straight line, degenerat1, degenerate2
  vn = [0,0,0,0,0,0]
  if q[3]*q[5] > 0 and q[2] == 0 and q[4] == 0 and q[1] == 0:
    vn = elipse(e, q, corners,dx,dy,nx,ny)
  elif q[4] != 0 and (q[2],q[3],q[5] == 0,0,0 ):
    vn = hyperbola(e, q, corners,dx,dy,nx,ny)
  elif q[2] != 0 and q[3] != 0 and (q[1], q[4], q[5] == 0,0,0):
    vn = parabola(e, q, corners,dx,dy,nx,ny)
  elif q[1] != 0 and (q[2],q[3], q[4], q[5] == 0,0,0,0):
    vn = straightLine(e, q, corners,dx,dy,nx,ny)
  elif q[3] != 0 and (q[1],q[2], q[4], q[5] == 0,0,0,0):
    vn = degenerat(e, q, corners,dx,dy,nx,ny)
  else :
    return [0,0,0,0,0,0];
  sum = 0.0
  Det = np.linalg.det(A1)

  vn = ABv(A,B,vn)
  vn = np.dot(vn,FI) * Det
  #print vn
  return vn

def triangelIntegral(e, E, k, f, corners):

  v = surfaces(e, corners,k,E)
  sum = 0
  for i in range(0, 6):
    sum += f[i]*v[i]
  return sum;


def totalInegral():
  totalInegralSum = 0;

  kGrid ,EGrid, fGrid = creatkGrid()
  #summing over the odd triangel
  E = EGrid
  k = kGrid
  f = FGrid
  corners = [k[0], k[1], k[2], k[3]]
  sum = triangelIntegral(e,E,k,f, corners)
  return sum

print totalInegral()
