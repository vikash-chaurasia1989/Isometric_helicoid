import numpy as np
import mathutils
import scipy.special as sp
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

global N,fac,wd,th,thickness

th = np.pi/2
fac = 5
wd  = fac/17




N     = 18*8
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
h = 1/N*2*nfold*K/r;

rx = np.zeros(N+1)
ry = np.zeros(N+1)
rz = np.zeros(N+1)
tx = np.zeros(N+1)
ty = np.zeros(N+1)
tz = np.zeros(N+1)
nx = np.zeros(N+1)
ny = np.zeros(N+1)
nz = np.zeros(N+1)
bx = np.zeros(N+1)
by = np.zeros(N+1)
bz = np.zeros(N+1)

sn = np.zeros(N+1)
cn = np.zeros(N+1)
dn = np.zeros(N+1)

for i in range(0,N+1):
    sn[i],cn[i], dn[i],temp = sp.ellipj(r*s[i],p)
    #print(p**2*sn[i]**2 + dn**2)

kappa = np.sqrt(al3*(1-q**2*sn**2))

#=== correcting sign of kappa in the relevant intervals of s

#=== correcting sign of kappa in the relevant intervals of s

for m in range(0,nfold-1,2):
    kappa[int((2*m+1)*N/(2*nfold)+1):int((2*m+1)*N/(2*nfold) + N/(nfold))] = -kappa[int((2*m+1)*N/(2*nfold)+1):int((2*m+1)*N/(2*nfold) + N/(nfold))]

m = nfold-1
kappa[int((2*m+1)*N/(2*nfold)+1):N+1] = -kappa[int((2*m+1)*N/(2*nfold)+1):N+1]


print(nfold)

'''
nfold = nfold*2
temp = np.arange(int(N/nfold)-1,int(3*N/nfold))
for i in range(0,int(((nfold/2-1)/2))):
     ind = int(4*(i-1)*N/nfold)+1 + temp
     kappa[ind] = -kappa[ind]

ind = np.arange(int((nfold-1)*N/nfold),(N+1))
kappa[ind] = -kappa[ind]
'''




# === frenet frame integration
initial = np.array([0, -0.1276, 0, 0.9280, 0, -0.3726, 0.3726, 0, 0.9280, 0, 1, 0],'d')
initial = np.array([0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 1],'d')

rx[0] = initial[0]
ry[0] = initial[1]
rz[0] = initial[2]

tx[0] = initial[3]
ty[0] = initial[4]
tz[0] = initial[5]

nx[0] = initial[6]
ny[0] = initial[7]
nz[0] = initial[8]

bx[0] = initial[9]
by[0] = initial[10]
bz[0] = initial[11]

i = 0

tx[i+1] = tx[i] + kappa[i]*h*nx[i]
ty[i+1] = ty[i] + kappa[i]*h*ny[i]
tz[i+1] = tz[i] + kappa[i]*h*nz[i]

nx[i+1] = nx[i] + (-kappa[i]*tx[i] + tau*bx[i])*h
ny[i+1] = ny[i] + (-kappa[i]*ty[i] + tau*by[i])*h
nz[i+1] = nz[i] + (-kappa[i]*tz[i] + tau*bz[i])*h

bx[i+1] = bx[i] - tau*nx[i]*h
by[i+1] = by[i] - tau*ny[i]*h
bz[i+1] = bz[i] - tau*nz[i]*h

#== central difference integration for t,n, and b

for i in range(1,N):

    tx[i+1] = tx[i-1] + kappa[i]*2*h*nx[i]
    ty[i+1] = ty[i-1] + kappa[i]*2*h*ny[i]
    tz[i+1] = tz[i-1] + kappa[i]*2*h*nz[i]


    nx[i+1] = nx[i-1] + (-kappa[i]*tx[i] + tau*bx[i])*2*h
    ny[i+1] = ny[i-1] + (-kappa[i]*ty[i] + tau*by[i])*2*h
    nz[i+1] = nz[i-1] + (-kappa[i]*tz[i] + tau*bz[i])*2*h

    bx[i+1] = bx[i-1] - tau*nx[i]*2*h
    by[i+1] = by[i-1] - tau*ny[i]*2*h
    bz[i+1] = bz[i-1] - tau*nz[i]*2*h

for i in range(1,N+1):
    rx[i] = rx[i-1]+h*tx[i-1]
    ry[i] = ry[i-1]+h*ty[i-1]
    rz[i] = rz[i-1]+h*tz[i-1]

#===== Integration complete ===

#==== Constructing midline and edges
# Transforming the midline
rx = rx - np.mean(rx)
ry = ry - np.mean(ry)
rz = rz - np.mean(rz)


#== Half width of the strip

l = sum(((rx[0:N-1]-rx[1:N])**2 + (ry[0:N-1]-ry[1:N])**2 +(rz[0:N-1]-rz[1:N])**2)**.5)

rx,ry,rz = fac*rx/l,fac*ry/l,fac*rz/l

points = []
points.append((float(rx[0]),float(ry[0]),float(rz[0])))
for i in range(N):
     points.append((float(rx[i+1]),float(ry[i+1]),float(rz[i+1])))

#== Half width of the strip
x1 = rx-wd*bx
y1 = ry-wd*by
z1 = rz-wd*bz

x2 = rx+wd*bx
y2 = ry+wd*by
z2 = rz+wd*bz
