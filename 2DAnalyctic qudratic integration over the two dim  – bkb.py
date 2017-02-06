
""" found the error modife the W must be minus not plus abs ! 
Analytic qudtratic interpolation "G wiesenekker,G received 1987" 2D test
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
Deltak = 0.1
pi = 3.14159265359
#enter the function f here 
def functionf (x,y):
  f = x**2 + y**2;
  #f = 1
  return f;

#enter the function E here 
def energyf (x,y):
  E = x**2 + y**2;
  #E = y**2
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

# Function sortE return E0<E1<E2<E3<E4<E5
def sortE(E, f, k):
  for i in range(1, len(E)):
        j = i
        while j > 0 and E[j] < E[j-1]:
            f[j], f[j-1] = f[j-1], f[j]
            k[j], k[j-1] = k[j-1], k[j]
            E[j], E[j-1] = E[j-1], E[j]
            j -= 1
  return E, f, k;
def sortE2(E, f, k):
  for i in range(1, len(E)):
        j = i
        while j > 0 and E[j] < E[j-1]:
            #f[j], f[j-1] = f[j-1], f[j] # no need for linear 
            k[j], k[j-1] = k[j-1], k[j]
            E[j], E[j-1] = E[j-1], E[j]
            j -= 1
  return E, f, k;
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
  # p0 + p1 x0 + p2 y0 + p3 X0^2 + p4 x0y0 + p5 y0^2 = f0
  # p1 + p1 x1 + p2 y1 + p3 X1^2 + p4 x1y1 + p5 y1^2 = f1
  # p2 + p1 x2 + p2 y2 + p3 X2^2 + p4 x2y2 + p5 y2^2 = f2
  # p3 + p1 x3 + p2 y3 + p3 X3^2 + p4 x3y3 + p5 y3^2 = f3
  # p4 + p1 x4 + p2 y4 + p3 X4^2 + p4 x4y4 + p5 y4^2 = f4
  # p5 + p1 x5 + p2 y5 + p3 X5^2 + p4 x5y5 + p5 y5^2 = f5
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

"""------------------------------------------Surface integration functions--------------------------------------------------------"""
def elipse (e, q, p, corners):
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
 	#u = sortU(u) wrong! becasue it change the order of lines in u.
 	#print u
 	#check if the point on the side of the triangel.

 	n = 0
 	sum = 0	
 	for i in range(0, 6):

 		if tStatus[i] == True :

 			if i/2 < 2 :
 				z = 1
 			else:
 				z = -2
 			#print ".............................."
 			#print corners
 			#print i
 			#print i/2
 			#print u[i]
 			#print a[i/2] ,"here"
 			#print b[i/2] 
 			#print corners
 			#print corners[(i/2)+z], corners[(i/2)] 
 			#print (constantA*np.cos(u[i]), constantB*np.sin(u[i]))
 			#print (constantA*np.cos(u[i]) - corners[i/2][0])/dx[i/2], (constantB*np.sin(u[i]) - corners[i/2][1])/dy[i/2]

 			s = abs((constantA*np.cos(u[i]) - corners[i/2][0])/dx[i/2])
 			limt = ((constantA*np.cos(u[i]) - corners[i/2][0])*(corners[i/2 + z][0]- corners[i/2][0]))
 			s = np.around ( s, 10)
 			limt = np.around (limt, 10)

 			if  ( s <=  1) and  0 <=  limt:
 				w[n] = u[i]
				n = n + 1
  				E1 = [corners[0][0]**2+corners[0][1]**2, corners[1][0]**2+corners[1][1]**2, corners[2][0]**2+corners[2][1]**2]
  				corners1 = [(corners[0][0], corners[0][1]),(corners[1][0], corners[1][1]),(corners[2][0], corners[2][1])]
  				sum = triangelIntegral1(e,E1, corners1,E1)
  				#print "............................."
  				#print corners , (constantA*np.cos(u[i]), constantB*np.sin(u[i]))
				#print sum
		else :
			if i/2 < 2 :
 				z = 1
 			else:
 				z = -2

		#print "case t ", t, "txT ", ((constantA*np.cos(u[i]) - corners[i/2][0])*(corners[i/2 + z][0]- corners[i/2][0]))
		#print i

	#calculating the value of u form the intersiction of the eclipse and the triangel.
	#fuck this shit just ma a loop and end this !!! 
	uf = 0
	#print u
	w1 = [w[i] for i in xrange(n)]
	w1 = sortU(w1)
	#print "wq afte sorting",w1
	for i in range(0, n):
		if (i + 1 < n ) and (w1[i] != w1[i + 1]) :
			#print "------------------up "
			point = constantA*np.cos( w1[i]+0.01), constantB*np.sin(w1[i]+0.01)
			#point = constantA*np.cos((w[i]+w[i+1])/2), constantB*np.sin((w[i]+w[i+1])/2)
			if (inTriangel(corners[0],corners[1],corners[2],point)):
				if abs((w1[i +1]) - (w1[i])) < pi :
					uf += abs((w1[i +1]) - (w1[i]))
				else :
					uf += abs((w1[i +1]) - (w1[i])) -2*pi
			#print point
			#print "plothing ", corners, point
		elif (i + 1 >= n) and (w1[i] != w1[0]):
			#print "------------------down " "this sutiaont must be fixed! "
			point = constantA*np.cos(w1[0]-0.01), constantB*np.sin(w1[0]-0.01)
			if (inTriangel(corners[0],corners[1],corners[2],point)):
				if abs((w1[0]) - ((w1[i]))) < pi :
					uf += abs((w1[0]) - ((w1[i])))
				else :
					uf += abs((w1[0]) - ((w1[i]))) - 2*pi
			#print point,"-------"
			#print "plothing ", corners, point
		#print w
		#print "..............",uf
		
	#uf = abs((w[0]) - (w[n-1]))
	uf = abs(uf)
	
	#print w
	#print 
	#print n
	#print " uf = ",uf


  	if n == 0:
  		return 0;

  	#ref (27)
  	Jacobian = ( 1.0 / (2*(q[3]*q[5])**0.5))

  	#ref (28) 
  	sum1 = p[0] * Jacobian * uf
  	sum1 += p[1] * Jacobian * ((e - q[0]) / q[3]) ** 0.5 * np.sin(uf)
  	sum1 += p[2] * Jacobian * ((e - q[0]) / q[5]) ** 0.5 * (-np.cos(uf))
  	sum1 += p[3] * Jacobian * abs((e - q[0]) / q[3]) * (0.5 * uf + 0.25 * np.sin(2*uf))
  	sum1 += p[4] * Jacobian * ((e - q[0]) / q[3]) **0.5 * ((e - q[0]) / q[5]) **0.5 *( -0.25 * np.cos(2*uf))
  	sum1 += p[5] * Jacobian * abs((e - q[0]) / q[5]) * (0.5 * uf - 0.25 * np.sin(2*uf))
  	#print "sum",sum1 , sum
  	return sum1;
def  hyperbola (e, q, p, corners):
  """
                  (x2,y2)
                  / \ 
                 /   \ 
                /     \ 
               /       \   
      (x0,y0) /_________\(x1,y1)
  """ 
  
  #need to find u for the intersection of E and the triganel sides.
  #ref 41 & 43
  #ref 36-40 my paper
    #corners = [(0,0),(1,1),(3,1)]
    #corners = [(0.0, 0.0), (1.0, 3.0), (0.0, 0.0)]
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
  #u = sortU(u) wrong! becasue it change the order of lines in u.
  #print u
  #check if the point on the side of the triangel.

  n = 0
  sum = 0 
  for i in range(0, 6):

    if tStatus[i] == True :

      if i/2 < 2 :
        z = 1
      else:
        z = -2
      #print ".............................."
      #print corners
      #print i
      #print i/2
      #print u[i]
      #print a[i/2] ,"here"
      #print b[i/2] 
      #print corners
      #print corners[(i/2)+z], corners[(i/2)] 
      #print (constantA*np.cos(u[i]), constantB*np.sin(u[i]))
      #print (constantA*np.cos(u[i]) - corners[i/2][0])/dx[i/2], (constantB*np.sin(u[i]) - corners[i/2][1])/dy[i/2]

      s = abs((constantA*np.cos(u[i]) - corners[i/2][0])/dx[i/2])
      limt = ((constantA*np.cos(u[i]) - corners[i/2][0])*(corners[i/2 + z][0]- corners[i/2][0]))
      s = np.around ( s, 10)
      limt = np.around (limt, 10)

      if  ( s <=  1) and  0 <=  limt:
        w[n] = u[i]
        n = n + 1
          E1 = [corners[0][0]**2+corners[0][1]**2, corners[1][0]**2+corners[1][1]**2, corners[2][0]**2+corners[2][1]**2]
          corners1 = [(corners[0][0], corners[0][1]),(corners[1][0], corners[1][1]),(corners[2][0], corners[2][1])]
          sum = triangelIntegral1(e,E1, corners1,E1)
          #print "............................."
          #print corners , (constantA*np.cos(u[i]), constantB*np.sin(u[i]))
        #print sum
    else :
      if i/2 < 2 :
        z = 1
      else:
        z = -2

    #print "case t ", t, "txT ", ((constantA*np.cos(u[i]) - corners[i/2][0])*(corners[i/2 + z][0]- corners[i/2][0]))
    #print i

  #calculating the value of u form the intersiction of the eclipse and the triangel.
  #fuck this shit just ma a loop and end this !!! 
  uf = 0
  #print u
  w1 = [w[i] for i in xrange(n)]
  w1 = sortU(w1)
  #print "wq afte sorting",w1
  for i in range(0, n):
    if (i + 1 < n ) and (w1[i] != w1[i + 1]) :
      #print "------------------up "
      point = constantA*np.cos( w1[i]+0.01), constantB*np.sin(w1[i]+0.01)
      #point = constantA*np.cos((w[i]+w[i+1])/2), constantB*np.sin((w[i]+w[i+1])/2)
      if (inTriangel(corners[0],corners[1],corners[2],point)):
        if abs((w1[i +1]) - (w1[i])) < pi :
          uf += abs((w1[i +1]) - (w1[i]))
        else :
          uf += abs((w1[i +1]) - (w1[i])) -2*pi
      #print point
      #print "plothing ", corners, point
    elif (i + 1 >= n) and (w1[i] != w1[0]):
      #print "------------------down " "this sutiaont must be fixed! "
      point = constantA*np.cos(w1[0]-0.01), constantB*np.sin(w1[0]-0.01)
      if (inTriangel(corners[0],corners[1],corners[2],point)):
        if abs((w1[0]) - ((w1[i]))) < pi :
          uf += abs((w1[0]) - ((w1[i])))
        else :
          uf += abs((w1[0]) - ((w1[i]))) - 2*pi
      #print point,"-------"
      #print "plothing ", corners, point
    #print w
    #print "..............",uf
    
  #uf = abs((w[0]) - (w[n-1]))
  uf = abs(uf)
  
  #print w
  #print 
  #print n
  #print " uf = ",uf


    if n == 0:
      return 0;

    #ref (27)
    Jacobian = ( 1.0 / (2*(q[3]*q[5])**0.5))

    #ref (28) 
    sum1 = p[0] * Jacobian * uf
    sum1 += p[1] * Jacobian * ((e - q[0]) / q[3]) ** 0.5 * np.sin(uf)
    sum1 += p[2] * Jacobian * ((e - q[0]) / q[5]) ** 0.5 * (-np.cos(uf))
    sum1 += p[3] * Jacobian * abs((e - q[0]) / q[3]) * (0.5 * uf + 0.25 * np.sin(2*uf))
    sum1 += p[4] * Jacobian * ((e - q[0]) / q[3]) **0.5 * ((e - q[0]) / q[5]) **0.5 *( -0.25 * np.cos(2*uf))
    sum1 += p[5] * Jacobian * abs((e - q[0]) / q[5]) * (0.5 * uf - 0.25 * np.sin(2*uf))
    #print "sum",sum1 , sum
    return sum1;

"""-------------------------------------------------(linear test)-------------------------------------------------------------------------------"""
# Function check if a point p inside the triangel! 
#p = p0 + (c1 - c0) * s + (c2 - c0) * t
#The point p is inside the triangle if 0 <= s <= 1 and 0 <= t <= 1 and s + t <= 1.
def inTriangel(c0,c1,c2,p):

	A = [[c1[0] - c0[0] , c2[0] - c0[0]],
    	[c1[1] - c0[1] , c2[1] - c0[1]]]

	B = [p[0] - c0[0],
		 p[1] - c0[1] ]
	#print A,B,"look here"
	st = np.linalg.solve(A, B)
	st = np.around(st, 10)
	if  ( 0 <= st[0] <=  1) and  ( 0 <= st[1] <=  1):
		if st[0] + st[1] <= 1:
			#print "status True : ",st
			return True

	#print "status false : ",st
	return False;

def getConstantsPi2(E,k,f):
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
def getConstantsEi2(E,k):
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
def triangelIntegral1(e,E,k,f):
  E, f, k = sortE2(E, f, k);
  p = getConstantsPi2(E,k,f)
  #print E, e
  #q = getConstantsEi(E,k) #is not needed for linear intorpolation. 

  if E[0] <= e and e <= E[1]:
    
    kt, ku = IEIntegral1(e, E, k);
    Jacobian = determinantJacobian1(e,E,k)
  elif E[1] <= e and e <= E[2]:

    kt, ku = IEIntegral2(e, E, k);
    Jacobian = determinantJacobian2(e,E,k)   
  else:
    kt = 0,0
    ku = 0,0
    Jacobian = 0
  
  #sum = p0*v0 + p1 *v1+ p2*v2  
  sum = p[0] * Jacobian + p[1]*Jacobian*(kt[0]+ 0.5*ku[0])+ p[2] *Jacobian * (kt[1] + 0.5 * ku[1]) ;
  return sum;
"""-------------------------------------------------------------------------------------------------------------------"""

# Function returns the integral of I(E)for a surface " Ellipse, Hyperbola, Parabola, Stright line, Degenerate"
def surfaces(e, q, p, corners):

  if q[3]*q[5] > 0 and q[2] == 0 and q[4] == 0:
	sum = elipse(e, q, p, corners)

  else :
  	return 0;

  return sum;

def triangelIntegral(e, E, k, f, corners):

  E, f, k = sortE(E, f, k);
  p = getConstantsPi(E,k,f)
  q = getConstantsEi(E,k)
  #rounding to 10^-10 decimal. 
  p = np.around(p, 10)
  q = np.around(q, 10)

  sum = surfaces(e, q, p, corners)
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
  for ix in range(0, len(kGrid) - 1, 2):
    for iy in range (0, len(kGrid[ix]) - 1, 2): # q and p must be here 
      i = i+1 
      #summing over the odd triangel
      E = [EGrid[ix][iy], EGrid[ix + 1][iy], EGrid[ix + 2][iy], EGrid[ix][iy + 1], EGrid[ix][iy +2], EGrid[ix + 1][iy + 1]]
      k = [kGrid[ix][iy], kGrid[ix + 1][iy], kGrid[ix + 2][iy], kGrid[ix][iy + 1], kGrid[ix][iy +2], kGrid[ix + 1][iy + 1]]
      f = [fGrid[ix][iy], fGrid[ix + 1][iy], fGrid[ix + 2][iy], fGrid[ix][iy + 1], fGrid[ix][iy +2], fGrid[ix + 1][iy + 1]]
      corners = [k[0], k[2], k[4]]
      sum += triangelIntegral(e,E,k,f, corners)

      #summing over the partall triangel
      E = [EGrid[ix][iy + 2], EGrid[ix + 1][iy + 2], EGrid[ix + 2][iy + 2], EGrid[ix + 2][iy + 1], EGrid[ix + 2][iy], EGrid[ix + 1][iy + 1]]
      k = [kGrid[ix][iy + 2], kGrid[ix + 1][iy + 2], kGrid[ix + 2][iy + 2], kGrid[ix + 2][iy + 1], kGrid[ix + 2][iy], kGrid[ix + 1][iy + 1]]
      f = [fGrid[ix][iy + 2], fGrid[ix + 1][iy + 2], fGrid[ix + 2][iy + 2], fGrid[ix + 2][iy + 1], fGrid[ix + 2][iy], fGrid[ix + 1][iy + 1]]
      corners = [k[0], k[2], k[4]]
      sum += triangelIntegral(e, E, k, f, corners)
  #print "i", i
  return sum

print totalInegral()
