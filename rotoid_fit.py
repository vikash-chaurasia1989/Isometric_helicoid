#=== Main solver
# importing library
import numpy as np
import math
import matplotlib
import scipy.optimize as optimize

# shortcuts
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import pi as PI
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import curve_fit

from mayavi.mlab import *
from plotplytest import *


#==== toroid curvefit

def fun_rotoid(t,a,b,c,p):

    data[:,0] = (a - b*np.cos(p*t))*np.cos(t) - c*np.sin(t)
    data[:,1] = (a - b*np.cos(p*t))*np.sin(t) + c*np.cos(t)
    data[:,2] =  b*np.sin(p*t)
    return data

def plotplytest(x,y,z):
    fig = go.Figure(data=[go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode= 'lines+markers',
    marker=dict(
        size=2,
        color=z,                # set color to an array/list of desired values
        colorscale= 'Viridis',   # choose a colorscale
        opacity=1

        ),

        #colorbar=dict(thickness=10, tickvals=[np.min(z), np.max(z)], ticktext=['Low', 'High']),
        )])

    # tight layout
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    #fig.add_trace(colorbar_trace)
    fig.show()
    return


'''
=========================================================================
 Name and path of the data file
=========================================================================
'''

branch = 1

if branch==1:

    nfold = 3
    str1 =  '3fold_N'
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/'

    strtau = path+'tau_branch_1.txt'
    tau1 = np.loadtxt(strtau)
    N = 72
    N1=N-1
    h = 1


elif branch==2:

    nfold = 5
    str1 = '5pi_knot_N'
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/'

    strtau = path + 'tau_branch_2.txt'
    tau1 = np.loadtxt(strtau)
    N = 120
    N1= N-1
    h = 1/N
else:

    nfold = 7
    str1 = '7pi_knot_N'
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/'

    strtau = path + 'tau_branch_3.txt'
    tau1 = np.loadtxt(strtau)

#
i = 7
tau  =  int(np.around(tau1[i+10]*(10**10),decimals=0))
#tau  =   np.around(tau1[i+10]*(10**10),decimals=1)

strb =  path + 'b_' + str1 + str(N)+ '_tau_'+ str(tau)+ '.txt'

#if i==0:
    #strb = path + 'b_5pi_knot_N120_tau_119842333032.1598.txt'
ln1 = np.size(tau1)


points = np.loadtxt(strb)

bx = points[:,0]
by = points[:,1]
bz = points[:,2]


bxp = np.transpose(np.zeros(N+1))
byp = np.transpose(np.zeros(N+1))
bzp = np.transpose(np.zeros(N+1))


#========================================================
#               b'
#========================================================
h = 1/N
bxp[0]   = (bx[0]+ bx[N-1])/h
byp[0]   = (by[0]+ by[N-1])/h
bzp[0]   = (bz[0]+ bz[N-1])/h


ig = np.arange(1,N+1)

bxp[ig] = (bx[ig]-bx[ig-1])/h
byp[ig] = (by[ig]-by[ig-1])/h
bzp[ig] = (bz[ig]-bz[ig-1])/h


#========================================================
#               Constructing tangent  t= bxb'
#========================================================
tx = np.transpose(np.zeros((N+1)))
ty = np.transpose(np.zeros((N+1)))
tz = np.transpose(np.zeros((N+1)))

#== tangent ==
tx[0:N] = np.multiply(by[0:N],bzp[0:N]) - np.multiply(bz[0:N],byp[0:N])
ty[0:N] = np.multiply(bz[0:N],bxp[0:N]) - np.multiply(bx[0:N],bzp[0:N])
tz[0:N] = np.multiply(bx[0:N],byp[0:N]) - np.multiply(by[0:N],bxp[0:N])

#========================================================

#========================================================
#           Constructing midline
#========================================================
#== Midline
rx = np.transpose(np.zeros((N+1)))
ry = np.transpose(np.zeros((N+1)))
rz = np.transpose(np.zeros((N+1)))



for i in range(N):
    rx[i+1] = rx[i] + tx[i]
    ry[i+1] = ry[i] + ty[i]
    rz[i+1] = rz[i] + tz[i]


# Transforming the midline
rx = rx - np.mean(rx)
ry = ry - np.mean(ry)
rz = rz - np.mean(rz)


#== Half width of the strip

l = sum(((rx[0:N-1]-rx[1:N])**2 + (ry[0:N-1]-ry[1:N])**2 +(rz[0:N-1]-rz[1:N])**2)**.5)

rx,ry,rz = rx/l,ry/l,rz/l

#== transforming the midline such that it is parallel to x by plane

ind = np.array([0,np.int(N/nfold)-1,np.int(2*N/nfold)-1])    # index of three symmetry points

v1x  = rx[ind[0]] - rx[ind[1]]
v1y  = ry[ind[0]] - ry[ind[1]]
v1z  = rz[ind[0]] - rz[ind[1]]

v2x  = rx[ind[0]] - rx[ind[2]]
v2y  = ry[ind[0]] - ry[ind[2]]
v2z  = rz[ind[0]] - rz[ind[2]]


norx = v1y*v2z - v1z*v2y
nory = v1z*v2x - v1x*v2z
norz = v1x*v2y - v2x*v1y

a = norx/(norx**2 + nory**2 + norz**2)**.5
b = nory/(norx**2 + nory**2 + norz**2)**.5
c = norz/(norx**2 + nory**2 + norz**2)**.5

u1,u2 = b,-a
th = math.acos(c)

rot = np.array([[np.cos(th) + u1**2*(1-np.cos(th)), u1*u2*(1-np.cos(th)), u2*np.sin(th)],\
[u1*u2*(1-np.cos(th)), np.cos(th) + u2**2*(1-np.cos(th)),-u1*np.sin(th) ],\
[-u2*np.sin(th), u1*np.sin(th), np.cos(th)]])


for i in range(0,N+1):
    temp = np.dot(rot,np.array([[rx[i]],[ry[i]],[rz[i]]]))

    rx[i]= temp[0]
    ry[i]= temp[1]
    rz[i]= temp[2]


#plotplytest(rx,ry,rz)

fig1 = plt.figure()
ax = fig1.add_subplot(111,projection='3d')
ax.plot(rx,ry,rz)
#plt.show()

#=== calling the curve fitting routine

data = np.zeros((N+1,3))
data[:,0] = rx
data[:,1] = ry
data[:,2] = rz



t = np.linspace(0,2*PI,N+1)
# initial guess
a = np.mean((rx**2 + ry**2 + rz**2)**.5)
b = .5*a
c = .5*a
p = 3

temp = fun_rotoid(t,a,b,c,p)
optimize.curve_fit(fun_rotoid,t,data,a,b,c,p)
