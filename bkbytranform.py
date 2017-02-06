
""" found the error modife the W must be minus not plus abs ! 
Analytic qudtratic interpolation "G wiesenekker,G received 1987" 2D test
dimensionaless k vekort
the answer given in (1/2pi)^2 *(2m / hbar^2)^n where  n = s +1. fn(k^(2s))  
"""
import numpy as np
import math
# testng the function
e = 1
miniX = -2
maxX = 2
miniY = -2
maxY = 2
Deltak = 0.1
pi = 3.14159265359
#A2.8
FI = [[1,0,0,0,0,0],
      [-3,-1,0,4,0,0],
      [-3,0,-1,0,0,4],
      [2,2,0,-4,0,0],
      [4,0,0,-4,4,-4],
      [2,0,2,0,0,-4]];

#enter the function f here 
def functionf (x,y):
  #f = x**2 + y**2;
  #f = 1+x+y+y**2 + y*x+x**2
  #f = x+y
  #f = x
  f = x + y**2
  return f;

#enter the function E here 
def energyf (x,y):
  #E = x**2 + y**2;
  #E = 2+x*y;
  E = x + y
  #E = 1 + 2*x
  #E = 0.9 + x**2
  #E = x**3
  #E = 1+x+y+5*y**2 + 2*y*x+x**2
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
  rotation = 0.5
 
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

# Function sortu return u0<u1<u2<u3<u4<u5
def sortU(u):
  for i in range(1, len(u)):
        j = i
        while j > 0 and u[j] < u[j-1]:
            
            u[j], u[j-1] = u[j-1], u[j]
            j -= 1
  return u;


#ref (20 a)
def getConstantsPi(E,k,f):
  matriseA =[ [1, k[0][0], k[0][1], k[0][0]**2, k[0][0]*k[0][1], k[0][1]**2 ],
              [1, k[1][0], k[1][1], k[1][0]**2, k[1][0]*k[1][1], k[1][1]**2 ],
              [1, k[2][0], k[2][1], k[2][0]**2, k[2][0]*k[2][1], k[2][1]**2 ],
              [1, k[3][0], k[3][1], k[3][0]**2, k[3][0]*k[3][1], k[3][1]**2 ],
              [1, k[4][0], k[4][1], k[4][0]**2, k[4][0]*k[4][1], k[4][1]**2 ],
              [1, k[5][0], k[5][1], k[5][0]**2, k[5][0]*k[5][1], k[5][1]**2 ]]

  matriseB =[f[0],
             f[1],
             f[2],
             f[3],
             f[4],
             f[5]]
  #print matriseB,"f"
  p = np.linalg.solve(matriseA, matriseB)
  return p;

#ref -(20b)
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

#ref -(not createdyet)
def getAfineConstat(q,p,k):
  #ref A1
  if (q[4] != 0):
    print"no equl 0"
    print q
    s = [[q[3] , 0.5*q[4]],
        [0.5*q[4], q[5]]]
    eigvals, A = np.linalg.eig(s)

    AT = np.transpose(A)
    q1, q2 = np.dot([q[0],q[1]],AT)

    m = np.dot(A,np.dot(s,AT))
    q3 = m[0][0]
    q5 = m[1][1]
    q = [q[0],q1,q2,q3,0,q5] 
    q = np.around(q, 6)
  if (q[3] == 0 and q[5] != 0):
    print "x = y' ; y = x' our tartgeer"
  elif(q[3] !=0 and q[1] != 0 ):
    print "x = x' -q2/2q4   Y=y'"
    #new -E is A1.1
  elif(q[5] != 0 and q[2] != 0):
    print "x=x' ; y = y'-q3/2q6"
  #print q , "q here"

  #ref A1.11-14
  #start returning surfaces.
  if (q[3] != 0 and q[5] != 0 and q[2] == 0 and q[1] == 0):
    if (q[3]*q[5] > 0):
      print "return elipse + return"

    else:
      print "rerutn hyperbola A1.12 + return" 
  if (q[3] != 0):
    print "paralel lines+ return" , q[3]
  if (q[1] != 0 or q[2] != 0):
    #"stright line q2' != 0 and q3 = 0 A1.13 and 14"
    cosT = q[1]/(q[1]**2 + q[2]**2)**0.5
    sinT = q[2]/(q[1]**2 + q[2]**2)**0.5
    for i in range(0, 6):
      k[i] = k[i][0]*cosT + k[i][1]*sinT , k[i][0]*sinT - k[i][1]*cosT
      #here must find the A and program it .
      #print "reach"
    A = [[cosT,sinT],
        [sinT,-cosT]]
    q1 = cosT*q[1] + sinT*q[2]
    q = [q[0],q1,0,0,0,0]
    return q,p,k,A

  # a11*x0 + a12*y0 + b0 = 0
  # a11*x1 + a12*y1 + b0 = 1
  # a11*x2 + a12*y2 + b0 = 0

  # a21*x0 + a22*y0 + b1 = 0
  # a21*x1 + a22*y1 + b1 = 0
  # a21*x2 + a22*y2 + b1 = 1
  """
  matriseA1 = [[corners[0][0], corners[0][1], 1],
              [corners[1][0], corners[1][1], 1 ],
              [corners[2][0], corners[2][1], 1 ]]

  matriseB1 =[0,1,0]

  matriseA2 = [[corners[0][0], corners[0][1], 1],
              [corners[1][0], corners[1][1], 1 ],
              [corners[2][0], corners[2][1], 1 ]]

  matriseB2 =[0,0,1]

  a1 = np.linalg.solve(matriseA1, matriseB1)
  a2 = np.linalg.solve(matriseA2, matriseB2)
  a1 = np.around(a1, 6)
  a2 = np.around(a2, 6)
  print a1, a2, "lokkkkkkkkkkkkkkkkkk"
  print matriseA1.diagonalize()
  return a1,a2;
 """
"""------------------------------------------Surface integration functions--------------------------------------------------------"""
def elipse (e, q, p, corners,dx,dy,nx,ny,fk,k,E):
 	#dxj*yj - dyj*xj     
  c = [dx[0] * corners[0][1] - dy[0]*corners[0][0],
  	   dx[1] * corners[1][1] - dy[1]*corners[1][0],
  	   dx[2] * corners[2][1] - dy[2]*corners[2][0]]
  	
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
  tStatus = [False, False, False, False, False, False]
  z = 0; 
  for i in range(0, 3):
    if (4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i])) > 0 :
      t[z] = (2*a[i] + np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) / (2*(b[i] + c [i]))
      t[z+1] = (2*a[i] - np.sqrt(4*a[i]**2 - 4*(c[i] - b[i])*(b[i] + c[i]))) / (2*(b[i] + c [i]))
      tStatus[z] = True
      tStatus[z+1] = True
    z += 2

  u = 2*np.arctan(t)
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
    return 0;
  #ref (27)
  #print uf 
  #check for new qn and make sure u have the correct E value this can make a diffrent. 
  Jacobian = ( 1.0 / (2*(q[3]*q[5])**0.5))
  sum1 = 0
  v = [0,0,0,0,0,0]
  #print p ,"p is here"
  qn = np.dot(FI,E)
  print qn
  for i in range(0, z-1, 2):
    #ref (28) 
    sum1 += p[0] * Jacobian * (uf[i+1] - uf[i])
    sum1 += p[1] * Jacobian * ((e - q[0]) / q[3]) ** 0.5 * (np.sin(uf[i+1])-(np.sin(uf[i])))
    sum1 += p[2] * Jacobian * ((e - q[0]) / q[5]) ** 0.5 * (np.cos(uf[i]) - np.cos(uf[i+1]))
    sum1 += p[3] * Jacobian * abs((e - q[0]) / q[3]) * (((0.5 * uf[i+1] + 0.25 * np.sin(2*uf[i+1]))-(0.5 * uf[i] + 0.25 * np.sin(2*uf[i]))))
    sum1 += p[4] * Jacobian * ((e - q[0]) / q[3]) **0.5 * ((e - q[0]) / q[5]) **0.5 *(( -0.25 * np.cos(2*uf[i+1]))-(-0.25 * np.cos(2*uf[i])))
    sum1 += p[5] * Jacobian * abs((e - q[0]) / q[5]) * ((0.5 * uf[i+1] - 0.25 * np.sin(2*uf[i+1]))-(0.5 * uf[i] - 0.25 * np.sin(2*uf[i])))

    v[0] +=  Jacobian * (uf[i+1] - uf[i])
    v[1] +=  Jacobian * ((e - qn[0]) / qn[3]) ** 0.5 * (np.sin(uf[i+1])-(np.sin(uf[i])))
    v[2] +=  Jacobian * ((e - qn[0]) / qn[5]) ** 0.5 * (np.cos(uf[i]) - np.cos(uf[i+1]))
    v[3] +=  Jacobian * abs((e - qn[0]) / qn[3]) * (((0.5 * uf[i+1] + 0.25 * np.sin(2*uf[i+1]))-(0.5 * uf[i] + 0.25 * np.sin(2*uf[i]))))
    v[4] +=  Jacobian * ((e - qn[0]) / qn[3]) **0.5 * ((e - qn[0]) / qn[5]) **0.5 *(( -0.25 * np.cos(2*uf[i+1]))-(-0.25 * np.cos(2*uf[i])))
    v[5] +=  Jacobian * abs((e - qn[0]) / qn[5]) * ((0.5 * uf[i+1] - 0.25 * np.sin(2*uf[i+1]))-(0.5 * uf[i] - 0.25 * np.sin(2*uf[i])))
    #print "sum",sum1 , sum
    #fix the order of f,



    #print k , "k"
    #print fk,"f"
    #getAfineConstat(corners)
    print v ,"v"
    print np.dot(((np.dot(v,FI))),fk),"sim old"
    # print FI
    #wf = 0 
    #for i in range(0, 6):
      #wf += np.dot(fk[i],((np.dot(v,FI[i]))))
      #print i," : ", ((np.dot(v,FI[i]))), fk[i]#k is missing here maybe it make the diffrent!!!
    #print wf
    print (sum1)
    print "--------------------------------"
  return (sum1);
"""--------------------------------------------------------------------------------------------------------------------------"""
def hyperbola (e, q, p, corners,dx,dy,nx,ny):

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
      u[z] = ((-b[i] + np.sqrt(b[i]**2 - 4*a[i]*c[i]))/ (2*a[i]))
      u[z+1] = ((-b[i] - np.sqrt(b[i]**2 - 4*a[i]*c[i])) / (2*a[i]))
      uStatus[z] = True
      uStatus[z+1] = True
    z += 2
  #chicking if the intersection point is actuly inside the triangle.
  n = 0
  for i in range (0,6):
    if uStatus[i] == True:
      point = u[i], constantA/u[i]
      if(inTriangel(corners,point)):
        w[n] = u[i]   
        n += 1
  uf = []
  w1 = [w[i] for i in xrange(n)]
  w1 = sortU(w1)
  z = 0
  for i in range(0, n):
    if (i + 1 < n ):
      point = ( w1[i]+0.0001), constantA/(w1[i]+0.0001)
      if (inTriangel(corners,point)):
        uf.append( w1[i])
        uf.append( w1[i+1])
        z += 2 
    elif(i < n):
      point = ( w1[0]-0.0001), constantA/(w1[0]-0.0001)
      if (inTriangel(corners,point)):
        uf.append( w1[n - 1])
        uf.append( w1[0])
        z += 2
  if z == 0:
    return 0;
  #ref (29)
  sum1 = 0
  #Jacobian = ( 1.0 / abs(q[4]*uf))
  for i in range(0, z-1, 2):
    #ref (28)
    sum1 += p[0] * (1/abs(q[4])) * (np.log(abs(uf[i+1])) - np.log(abs(uf[i])))
    sum1 += p[1] * (1/abs(q[4])) * (uf[i+1] - uf[i])  
    sum1 += p[2] * (1/abs(q[4])) * constantA*(-1/uf[i+1] + 1/uf[i]) 
    sum1 += p[3] * (1/abs(q[4])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2) 
    sum1 += p[4] * (1/abs(q[4])) * constantA*(np.log(abs(uf[i+1])) - np.log(abs(uf[i]))) 
    sum1 += p[5] * (1/abs(q[4])) * ((e-q[0])**2 / q[4]**2) * (-1/(2*uf[i+1]**2) + 1/(2*uf[i]**2))
  return abs(sum1);
"""--------------------------------------------------------------------------------------------------------------------------"""
def parabola(e, q, p, corners,dx,dy,nx,ny):
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
      u[z] = ((-b[i] + np.sqrt(b[i]**2 - 4*a[i]*c[i]))/ (2*a[i]))
      u[z+1] = ((-b[i] - np.sqrt(b[i]**2 - 4*a[i]*c[i])) / (2*a[i]))
      uStatus[z] = True
      uStatus[z+1] = True
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
  #print w1
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
    return 0;
  #ref (32)
  sum1 = 0
  #print p
  for i in range(0, z-1, 2):
    #ref (28)
    sum1 += p[0] * (1/abs(q[2])) * (uf[i+1] - uf[i])
    sum1 += p[1] * (1/abs(q[2])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)  
    sum1 += p[2] * (1/abs(q[2])) * (constantA*uf[i+1]- (1.0/3)*(q[3]/q[2]*uf[i+1]**3) - (constantA*uf[i]- (1.0/3)*(q[3]/q[2]*uf[i]**3)))
    sum1 += p[3] * (1/abs(q[2])) * ((1.0/3)*uf[i+1]**3 - (1.0/3)*uf[i]**3 )
    sum1 += p[4] * (1/abs(q[2])) * (((0.5)*constantA * uf[i+1]**2 - (1.0/4)*(q[3]/q[2])*uf[i+1]**4) - ((0.5)*constantA * uf[i]**2 - (1.0/4)*(q[3]/q[2])*uf[i]**4))
    sum1 += p[5] * (1/abs(q[2])) * ( ( ((e-q[0])**2 /(q[2]**2))*uf[i+1]-(2.0/3)*(e-q[0])/q[2]*(q[3]/q[2])*uf[i+1]**3+(1.0/5)*(q[3]**2/q[2]**2)*uf[i+1]**5 )-( ((e-q[0])**2 /(q[2]**2))*uf[i]-(2.0/3)*(e-q[0])/q[2]*(q[3]/q[2])*uf[i]**3+(1.0/5)*(q[3]**2/q[2]**2)*uf[i]**5 ))
   
  return abs(sum1);
"""-------------------------------------------------------------------------------------------------------------------------------"""
def straightLine(e, q, p, corners,dx,dy,nx,ny,A):  
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
  uStatus = [True, True, True]
  z = 0;
  for i in range(0, 3):
    u[i] = -c[i]/b[i]
  #chicking if the intersection point is actuly inside the triangle.
  n = 0
  for i in range (0,3):
    if uStatus[i] == True:
      point = constantA, u[i]
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
    return 0;
  #ref (35)
  sum1 = 0
  v = [0,0,0,0,0,0]
  print A
  for i in range(0, z-1, 2):
    v[0] = (1/abs(q[1])) * (uf[i+1] - uf[i])
    v[1] = (1/abs(q[1])) * constantA * (uf[i+1] - uf[i])
    v[2] = (1/abs(q[1])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)
    v[3] = (1/abs(q[1])) * constantA**2 * (uf[i+1] - uf[i])
    v[4] = (1/abs(q[1])) * constantA * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)
    v[5] = (1/abs(q[1])) * (1.0/3 * uf[i+1]**3 - 1.0/3 * uf[i]**3)

    sum1 += p[0] * (1/abs(q[1])) * (uf[i+1] - uf[i])
    sum1 += p[1] * (1/abs(q[1])) * constantA * (uf[i+1] - uf[i])
    sum1 += p[2] * (1/abs(q[1])) * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)
    sum1 += p[3] * (1/abs(q[1])) * constantA**2 * (uf[i+1] - uf[i])
    sum1 += p[4] * (1/abs(q[1])) * constantA * (0.5*uf[i+1]**2 - 0.5*uf[i]**2)
    sum1 += p[5] * (1/abs(q[1])) * (1.0/3 * uf[i+1]**3 - 1.0/3 * uf[i]**3)
  #ref(46)
  vn = [0,0,0,0,0,0]
  vn[0] = v[0]
  vn[1] = A[0][0]*v[1] + A[0][1]*v[2]
  vn[2] = A[1][0]*v[1] + A[1][1]*v[2]
  vn[3] = (A[0][0]**2)*v[3] + (A[0][1]**2)*v[5] + 2*(A[0][0]*A[0][1]*v[4] )
  vn[4] = A[0][0]*A[1][0]*v[3] + A[0][1]*A[1][1]*v[5] + (A[0][0]*A[1][1] + A[1][0]*A[0][1])*v[4]
  vn[5] = (A[1][0]**2)*v[3] + (A[1][1]**2)*v[5] + 2*(A[1][0]*A[1][1]*v[4])
  sum1 = 0
  for i in range(0, 6):
    sum1 += p[i]*vn[i]
  return (sum1);

"""-------------------------------------------------------------------------------------------------------------------------------"""
def degenerat(e, q, p, corners,dx,dy,nx,ny):  
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
  uStatusP = [True, True, True]
  uStatusM = [True, True, True]
  z = 0;
  for i in range(0, 3):
    uP[i] = -cP[i]/b[i]
    uM[i] = -cM[i]/b[i]
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
    return 0;
  #ref (40)/ my paper (47- 60)
  sum1P = 0
  j = 0.5/(q[3]*(e-q[0]))**0.5
  for i in range(0, zP-1, 2):
    sum1P += p[0] * j * (ufP[i+1] - ufP[i])
    sum1P += p[1] * j * constantAP * (ufP[i+1] - ufP[i])
    sum1P += p[2] * j * (0.5*ufP[i+1]**2 - 0.5*ufP[i]**2)
    sum1P += p[3] * j * constantAP**2 * (ufP[i+1] - ufP[i])
    sum1P += p[4] * j* constantAP * (0.5*ufP[i+1]**2 - 0.5*ufP[i]**2)
    sum1P += p[5] * j * (1.0/3 * ufP[i+1]**3 - 1.0/3 * ufP[i]**3)
  sum1M = 0
  j = 0.5/(q[3]*(e-q[0]))**0.5

  for i in range(0, zM-1, 2):
    sum1M += p[0] * j * (ufM[i+1] - ufM[i])
    sum1M += p[1] * j * constantAM * (ufM[i+1] - ufM[i])
    sum1M += p[2] * j * (0.5*ufM[i+1]**2 - 0.5*ufM[i]**2)
    sum1M += p[3] * j * constantAM**2 * (ufM[i+1] - ufM[i])
    sum1M += p[4] * j* constantAM * (0.5*ufM[i+1]**2 - 0.5*ufM[i]**2)
    sum1M += p[5] * j * (1.0/3 * ufM[i+1]**3 - 1.0/3 * ufM[i]**3)
   
  return abs(sum1M + sum1P);
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

"""-------------------------------------------------------------------------------------------------------------------"""

# Function returns the integral of I(E)for a surface " Ellipse, Hyperbola, Parabola, Stright line, Degenerate"
def surfaces(e, q, p, corners,f,k,E):
  """
                  (x2,y2)
                  / \ 
                 /   \ 
                /     \ 
               /       \   
      (x0,y0) /_________\(x1,y1)
  """ 
  
  #need to find u for the intersection of E and the triganel sides.
  #ref 41 & 42
  #ref 24-35 my paper
  #corners = [(0,0),(1,1),(3,1)]
  #corners = [(0.0, 0.0), (1.0, 3.0), (0.0, 0.0)]
  print q ," q before"
  #getAfineConstat(q,p,k)
  q, p, k, A = getAfineConstat(q,p,k)
  #importan becasue we change the x and the
  corners = [k[0], k[1], k[2]]
  print q , "q after"

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
  if q[3]*q[5] > 0 and q[2] == 0 and q[4] == 0:
	  sum = elipse(e, q, p, corners,dx,dy,nx,ny,f,k,E)
  elif q[4] != 0 and (q[2],q[3],q[5] == 0,0,0 ):
    sum = hyperbola(e, q, p, corners, dx, dy, nx, ny)
  elif q[2] != 0 and q[3] != 0 and (q[1], q[4], q[5] == 0,0,0):
    sum = parabola(e, q, p, corners, dx, dy, nx, ny)
  elif q[1] != 0 and (q[2],q[3], q[4], q[5] == 0,0,0,0):
    sum = straightLine(e, q, p, corners, dx, dy, nx, ny,A)
  elif q[3] != 0 and (q[1],q[2], q[4], q[5] == 0,0,0,0):
    sum = degenerat(e, q, p, corners, dx, dy, nx, ny)
  else :
  	return 0;

  return sum;

def triangelIntegral(e, E, k, f, corners):

  #E, f, k = sortE(E, f, k);
  p = getConstantsPi(E,k,f)
  q = getConstantsEi(E,k)
  #rounding to 10^-10 decimal. 
  p = np.around(p, 10)
  q = np.around(q, 10)

  sum = surfaces(e, q, p, corners,f,k,E)
  return sum;


def totalInegral():
  totalInegralSum = 0;

  #creating the k grid and calc the energy 
  kGrid ,EGrid, fGrid = creatkGrid()

  numrows = len(kGrid)
  numcols = len(kGrid)
  sum = 0.0
  #integrating over the  triangels
  i = 0
  #i think the second raw must be fixed flip the triangel number 2, it  give a m inus evedence.
  for ix in range(0, len(kGrid) - 1, 2):
    for iy in range (0, len(kGrid[ix]) - 1, 2): # q and p must be here 
      i = i+1 
      #summing over the odd triangel
      E = [EGrid[ix][iy], EGrid[ix + 2][iy], EGrid[ix ][iy + 2], EGrid[ix + 1][iy], EGrid[ix + 1][iy + 1], EGrid[ix][iy + 1]]
      k = [kGrid[ix][iy], kGrid[ix + 2][iy], kGrid[ix ][iy + 2], kGrid[ix + 1][iy], kGrid[ix + 1][iy + 1], kGrid[ix][iy + 1]]
      f = [fGrid[ix][iy], fGrid[ix + 2][iy], fGrid[ix ][iy + 2], fGrid[ix + 1][iy], fGrid[ix + 1][iy + 1], fGrid[ix][iy + 1]]
      corners = [k[0], k[1], k[2]]
      sum += triangelIntegral(e,E,k,f, corners)

      #summing over the partall triangel
      E = [EGrid[ix + 2][iy +2], EGrid[ix][iy + 2],EGrid[ix + 2][iy], EGrid[ix + 1][iy + 2], EGrid[ix + 1][iy+1], EGrid[ix + 2][iy + 1]]
      k = [kGrid[ix + 2][iy +2], kGrid[ix][iy + 2],kGrid[ix + 2][iy], kGrid[ix + 1][iy + 2], kGrid[ix + 1][iy+1], kGrid[ix + 2][iy + 1]]
      f = [fGrid[ix + 2][iy +2], fGrid[ix][iy + 2],fGrid[ix + 2][iy], fGrid[ix + 1][iy + 2], fGrid[ix + 1][iy+1], fGrid[ix + 2][iy + 1]]
      corners = [k[0], k[1], k[2]]
      sum += triangelIntegral(e, E, k, f, corners)
  #print "i", i
  return sum

print totalInegral()
matriseA1 = [[2, 1],
            [7, 3]]

eigvals, eigvecs = np.linalg.eig(matriseA1)
print eigvals, eigvecs 
invp = np.linalg.inv(eigvecs)

print np.dot(eigvecs,np.dot(matriseA1,invp)) ,"first"

# test 1 triangel
"""
kGrid ,EGrid, fGrid = creatkGrid()
E = [0,1,1,0.25,0.5,0.25]
k = [[0,0],[0.0000000001,1.00000000001],[1.000000000001,0.000000001],[0.50000000000002,0.0000000000005],[0.5000000000002,0.50000000004],[0,0.500000000006]]
f = [1,1,1,1,1,1]
corners = [k[0], k[1], k[2]]
triangelIntegral(e,E,k,f, corners)
"""