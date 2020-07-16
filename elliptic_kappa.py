import numpy as np
import scipy.special as sp

N     = 18*80
N1    = N-1
nfold = 3
A = 8.075921715351875e+01
B = 2.591801747383097e+02

al1 = A + np.sqrt(A**2 + B**2)
al2 = 0;
al3 = -A + np.sqrt(A**2 + B**2)

p = np.sqrt((al3-al2)/(al3+al1))
q = np.sqrt(1-al2/al3)
r = 1/2*np.sqrt(al3+al1)
K = sp.ellipk(p)
s = np.linspace(0,2*nfold*K/r,N+1)
tau = 8.094090000000000e+00*2*nfold*K/r  #8.093946633549898e+00


sn = np.zeros(N+1)
cn = np.zeros(N+1)
dn = np.zeros(N+1)

for i in range(0,N+1):
    sn[i],cn[i], dn[i],temp = sp.ellipj(r*s[i],p)

kappa = np.sqrt(al3*(1-q**2*sn**2))
