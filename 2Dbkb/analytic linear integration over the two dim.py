"""
Analytic linear interpolation "G wiesenekker,G received 1987"
"""
import numpy as np

# testing the function
e = 3.1
E = [1.0,12.0,6.0]
k = [[1.4,2.0],[3.0,4.2],[5.3,6.2]]
f = [1.0,2.0,3.0]

# Function sortE return E1<E2<E3
def sortE(E, f, k):
  for i in range(1, len(E)):
        j = i
        while j > 0 and E[j] < E[j-1]:
            E[j], E[j-1] = E[j-1], E[j]
            f[j], f[j-1] = f[j-1], f[j]
            k[j], k[j-1] = k[j-1], k[j]
            j -= 1
  return E, f, k;

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

# Function returns the integral of I(E)for E1<e<E2
def IEIntegral1(e, E, k):
  #kt = k[0] + (e - E[0])/(E[2]-E[0]) * (k[2] - k[0]) 
  kt = np.add(k[0], np.multiply((e - E[0])/(E[2]-E[0]), (np.subtract(k[2],k[0])))) 

  #ku = (e - E[0])/(E[1]-E[0]) * (k[1]-k[0]). - (e - E[0]) / (E[2]-E[0]) * (k[2]- k[0]) 
  ku = np.subtract(np.multiply((e-E[0])/(E[1]-E[0]),(np.subtract(k[1],k[0]))), np.multiply((e-E[0])/(E[2]-E[0]) ,np.subtract(k[2],k[0])) )
  # Function returns the integral of I(E) for E2<e<E3
  return kt, ku;

def IEIntegral2(e, E, k):
  #kt = k[2] + (e - E[2])/(E[2]-E[0]) * (k[2] - k[0]) 
  kt = np.add(k[2], np.multiply((e - E[2])/(E[2]-E[0]), (np.subtract(k[2],k[0])))) 

  #ku = (e - E[2])/(E[2]-E[1]) * (k[2]-k[1]). - (e - E[2]) / (E[2]-E[0]) * (k[2]- k[0]) 
  ku = np.subtract(np.multiply((e-E[2])/(E[2]-E[1]),(np.subtract(k[2],k[1]))), np.multiply((e-E[2])/(E[2]-E[0]) ,np.subtract(k[2],k[0])) )
  # Function returns the integral of I(E) for E2<e<E3
  return kt, ku;

def determinantJacobian1(e, E, k):
  # the Jacobian (e - E[0])/((E[1]-E[0])*(E[2]-E[0])) * [[k20 - k00]  , [k10 - k00]]
  #                                                     [[k21 - k01]  , [k11 - k01]] 
  j = np.multiply((e - E[0])/((E[1]-E[0])*(E[2]-E[0])) , [[k[2][0] - k[0][0], k[1][0] - k[0][0]],[k[2][1] - k[0][1], k[1][1] - k[0][1]]])
  return np.linalg.det(j);

def determinantJacobian2(e, E, k):
  # the Jacobian (E[2] - e)/((E[2]-E[1])*(E[2]-E[0])) * [[k20 - k10]  , [k20 - k00]]
  #                                                     [[k21 - k11]  , [k21 - k01]] 
  j = np.multiply((E[2] - e)/((E[2]-E[1])*(E[2]-E[0])) , [[k[2][0] - k[1][0], k[2][0] - k[0][0]],[k[2][1] - k[1][1], k[2][1] - k[0][1]]])
  return np.linalg.det(j);

def integralSum(e,E,k,f):
  sortE(E, f, k);
  p = getConstantsPi(E,k,f)

  if E[0] < e and e < E[1]:

    kt, ku = IEIntegral1(e, E, k);
    Jacobian = determinantJacobian1(e,E,k)
  elif E[2] < e and e < E[2]:

    kt, ku = IEIntegral2(e, E, k);
    Jacobian = determinantJacobian2(e,E,k)   
  else:
    kt = 0,0
    ku = 0,0
    Jacobian = 0
  #sum = p0*v0 + p1 *v1+ p2*v2  
  sum = p[0] * Jacobian + p[1]*Jacobian*(kt[0]+ 0.5*ku[0])+ Jacobian * (kt[1] + 0.5 * ku[1]) ;
  return sum;


sum = integralSum(e,E,k,f)
print sum
