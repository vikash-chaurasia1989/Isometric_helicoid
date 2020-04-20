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
from mayavi.mlab import *
#from mayavi import mlab

#==== Importing my subroutine for mobius surface and midline

#from plot_mobius_data import *

#==== Global variables

global N,N1,u,v,w,f_b,h,tau,rho,lm
global bx,by,bz,bxp,byp,bzp,bx2p,by2p,bz2p,bx4p,by4p,bz4p,bext,f_b,id,tau,jac,ind,lmp,path,p2




#=== Initialization
#N = 240
#tau1 = np.array([16,16])#np.arange(16,16,.2)
tau1 = np.arange(15.8,20,.1)
tau1 = np.arange(19.9,25,.1)
tau1 = np.arange(22.165,25,.05)
#tau1 = np.arange(23.665,25,.05)
#tau1 = np.arange(23.065,22.165,-.05)
#tau1 = np.array([22.15,22.15])
tau1 = np.arange(22.165,20.2,-.05)

temp1 = np.arange(15.3,22.1,.1)
temp2 = np.arange(22.165,25,.05)

tau1  = np.concatenate((temp1,temp2))

tau1 = np.arange(8.1,21,.1)

#tau1 = np.arange(12.7,25,1)

err  =np.zeros(len(tau1))


#for p2 in range(45,len(tau1)):
for p2 in range(3,len(tau1)):
#for p2 in range(2,1,-1):
    tau = tau1[p2]
    #===============================================================================
    #======================     Initial guess ======================================
    #===============================================================================
    def initial_guess():

        global N,N1,u,v,w,f_b,h,tau,rho,lm
        global bx,by,bz,bxp,byp,bzp,bx2p,by2p,bz2p,bx4p,by4p,bz4p,bext,f_b,id,tau,jac,path, p2,tau1

        #N  = 100

        # reading data file
        tau2 = int(np.around(tau1[p2-1]*(10**10),decimals=0))
        path ='/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/'
        #path ='/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_discrete/data_3pi/'
        #str0 = '3fold_N72_tau_' + str(tau2) + '.txt'
        #str0 = '7pi_knot_49_tau_'+ str(tau2) + '.txt'
        #str0  = '7pi_knot_N49_tau_160000000000.txt'

        if(p2==0):
            str0 = '7pi_knot_N84_tau_153000000000.txt'#'7pi_knot_N84_tau_237000000000.txt'
            str0 = '7pi_knot_discrete_N_84_tau_154500000000.txt'
        else:
            str0 = '7pi_knot_discrete_N_84_tau_' + str(tau2) + '.txt'

        str0 = '3fold_N72_tau_' + str(tau2) + '.txt'
        #str0 = '3fold_discrete_N72_tau_' + str(tau2) + '.txt'
        #print(str0)
        #str0 = '7pi_knot_N84_tau_' + str(tau2) + '.txt'
        #str0 = '7pi_knot_84_tau_154000000000.0.txt'
        #str0 = '7pi_knot_discrete_N_84_tau_154500000000.txt'
        strb    = path + 'b_' + str0
        strlm   = path + 'lm_' + str0
        strrho  = path + 'rho_' + str0
        struvw  = path + 'uvw_' + str0
        strbext = path + 'bext_' + str0


        points = np.loadtxt(strb)

        N = len(points)-1
        N1 = N-1


        #   Initialization
        bx = np.transpose(np.zeros(N+1))
        by = np.transpose(np.zeros(N+1))
        bz = np.transpose(np.zeros(N+1))

        bxp = np.transpose(np.zeros(N+1))
        byp = np.transpose(np.zeros(N+1))
        bzp = np.transpose(np.zeros(N+1))


        bx2p = np.transpose(np.zeros(N+1))
        by2p = np.transpose(np.zeros(N+1))
        bz2p = np.transpose(np.zeros(N+1))

        bx4p = np.transpose(np.zeros(N-1))
        by4p = np.transpose(np.zeros(N-1))
        bz4p = np.transpose(np.zeros(N-1))

        bext = np.transpose(np.zeros(6))


        bx = points[:,0]
        by = points[:,1]
        bz = points[:,2]
        '''
        s = np.linspace(0,1,N+1)

        th = PI/2 + .1*np.cos(7*PI*s)
        ph = PI*s

        bx = np.multiply(np.sin(th),np.cos(th))
        by = np.multiply(np.sin(th),np.sin(th))
        bz = np.cos(th)
        '''

        N  = len(bx)-1
        N1 = N-1
        h  = 1/N
        #=== parameters u v and w for closure constraint
        u =  -5.1566319025341907e-01
        v =  -8.00693449731905190e+00
        w =   3.9967030985498e-01

        points = np.loadtxt(struvw)
        u = points[0]
        v = points[1]
        w = 0*points[2]

        '''
        temp = .001/(500*h**3)

        u = u*temp
        v = v*temp
        w = 0*temp
        '''
        #== lagrange multiplier rho for unity constraint
        rho = np.loadtxt(strrho)
        #rho = points[:,0]
        #rho = 0*np.ones(N1)
        lm = np.loadtxt(strlm)

        #lm = 0.0*np.ones(N1+2)

        #lm = 0*np.loadtxt(strlm)
        #c = np.array([temp[0]])
        #lm = np.concatenate((temp,c))

        #lm  = points[:,0]
        #

        #==== Smoothening of the lagrange multipliers
        '''
        lm[0:N1] = .5*(lm[0:N1]+ lm[1:N1+1])
        lm[N1+1] = lm[0]
        #lm[N1]   = lm[1]
        lm[0:N1] = .5*(lm[0:N1]+ lm[1:N1+1])
        lm[N1+1] = lm[0]
        lm[N1]   = lm[1]

        rho[0:N1-2] = .5*(rho[0:N1-2]+ rho[1:N1-1])
        rho[N1-1]   =  rho[0]
        #rho[N1-2]   =  rho[1]
        rho[0:N1-2] = .5*(rho[0:N1-2]+ rho[1:N1-1])
        rho[N1-1]   =  rho[0]
        '''

        #=== polynomial fit for the boundary points ==
        nfit = 4
        #-- At i = 1
        px = np.polyfit(np.linspace(1,nfit+1,nfit+1),bx[0:nfit+1],nfit)
        bx1 = px[-1]

        py = np.polyfit(np.linspace(1,nfit+1,nfit+1),by[0:nfit+1],nfit)
        by1 = py[-1]

        pz = np.polyfit(np.linspace(1,nfit+1,nfit+1),bz[0:nfit+1],nfit)
        bz1 = pz[-1]

        #-- At i = N+1
        px = np.polyfit(np.linspace(1,nfit+1,nfit+1),bx[N:N-nfit-1:-1],nfit)
        bxN2 = px[-1]

        py = np.polyfit(np.linspace(1,nfit+1,nfit+1),by[N:N-nfit-1:-1],nfit)
        byN2 = py[-1]

        pz = np.polyfit(np.linspace(1,nfit+1,nfit+1),bz[N:N-nfit-1:-1],nfit)
        bzN2 = pz[-1]

        bext = np.zeros(6)
        bext[0] = bx1
        bext[1] = by1
        bext[2] = bz1

        bext[3] = bxN2
        bext[4] = byN2
        bext[5] = bzN2

        bext = np.loadtxt(strbext)
        #bext= points[:,0]

        #lm,rho = fun_lambda(bx,by,bz,bext)

        b = np.zeros(3*N1)

        b[0:3*N1-2:3]   = bx[1:N]
        b[1:3*N1-1:3]   = by[1:N]
        b[2:3*N1:3]     = bz[1:N]

        #== arranging all the variables in a column vector

        f = np.zeros(5*N1+11)

        f[0:3*N1]           = b
        f[3*N1:4*N1]        = rho
        f[4*N1:5*N1+2]      = lm
        f[5*N1+2:5*N1+5]    = [u,v,w]
        f[5*N1+5:5*N1+11]   = bext


        #=== reassigning boundary values to bx by and bz

        bx[0] = 1.000000000000000
        by[0] = 0.000000000000000
        bz[0] = 0.000000000000000

        bx[N] = -1.000000000000000
        by[N] =  0.000000000000000
        bz[N] =  0.000000000000000

        bz[1]   = 0.000000000000000
        bz[N-1] = 0.000000000000000

        #=== Intialize arrays for function and jacobian evaluation

        f_b = np.zeros(5*N1+11)
        jac = np.zeros((5*N1+11,5*N1+11))
        id  = np.array([(1,0,0),(0,1,0),(0,0,1)])

       # deleting the entries corresponding to bz[1] and bz[N-1]
        f = np.delete(f,2)
        f = np.delete(f,3*N1-2)


        return f


    #===============================================================================
    #===============================================================================
    #===============================================================================
    #======================    fun_curve ======================================
    #===============================================================================



    def fun_curve(x):


        global N,N1,u,v,w,f_b,h,tau,rho,lm
        global bx,by,bz,bxp,byp,bzp,bx2p,by2p,bz2p,bx4p,by4p,bz4p,bext,f_b,tau,jac,ind,lmp


        #=== Initialize arrays ===



        bx[1] = x[0]
        by[1] = x[1]

        bx[2:N]     = x[2:3*N1-3:3]
        by[2:N]     = x[3:3*N1-2:3]
        bz[2:N-1]   = x[4:3*N1-4:3]


        rho         = x[3*N1-2:4*N1-2]
        lm          = x[4*N1-2:5*N1]

        u           = x[5*N1]
        v           = x[5*N1+1]
        w           = x[5*N1+2]

        bext        = x[5*N1+3:5*N1+9]


        #=== calculation of derivatives ===

        ig = np.arange(1,N)

        #===========================================================================
        # ========= b' ==== O(h**2) ===========
        #===========================================================================
        bxp[0]  = (bx[1]-bext[0])/(2*h)
        byp[0]  = (by[1]-bext[1])/(2*h)
        bzp[0]  = (bz[1]-bext[2])/(2*h)

        bxp[ig] = (bx[ig+1]-bx[ig-1])/(2*h)
        byp[ig] = (by[ig+1]-by[ig-1])/(2*h)
        bzp[ig] = (bz[ig+1]-bz[ig-1])/(2*h)

        bxp[N]  = (bext[3]-bx[N-1])/(2*h)
        byp[N]  = (bext[4]-by[N-1])/(2*h)
        bzp[N]  = (bext[5]-bz[N-1])/(2*h)

        lmp = (lm[ig+1]-lm[ig-1])/(2*h)
        #===========================================================================
        #===========================================================================


        bx2p[0] = (bx[1] + bext[0] -2*bx[0])/h**2
        by2p[0] = (by[1] + bext[1] -2*by[0])/h**2
        bz2p[0] = (bz[1] + bext[2] -2*bz[0])/h**2

        bx2p[ig] = (bx[ig+1] + bx[ig-1] - 2*bx[ig])/h**2
        by2p[ig] = (by[ig+1] + by[ig-1] - 2*by[ig])/h**2
        bz2p[ig] = (bz[ig+1] + bz[ig-1] - 2*bz[ig])/h**2

        bx2p[N] = (bext[3] + bx[N-1]  - 2*bx[N])/h**2
        by2p[N] = (bext[4] + by[N-1]  - 2*by[N])/h**2
        bz2p[N] = (bext[5] + bz[N-1]  - 2*bz[N])/h**2

        #===========================================================================
        #===========================================================================

        bx4p[0] = (bx[3] - 4*bx[2] + 6*bx[1] - 4*bx[0] + bext[0])/h**4
        by4p[0] = (by[3] - 4*by[2] + 6*by[1] - 4*by[0] + bext[1])/h**4
        bz4p[0] = (bz[3] - 4*bz[2] + 6*bz[1] - 4*bz[0] + bext[2])/h**4


        ig = np.arange(2,N-1)

        bx4p[ig-1] = (bx[ig+2] - 4*bx[ig+1] + 6*bx[ig] - 4*bx[ig-1] + bx[ig-2])/h**4
        by4p[ig-1] = (by[ig+2] - 4*by[ig+1] + 6*by[ig] - 4*by[ig-1] + by[ig-2])/h**4
        bz4p[ig-1] = (bz[ig+2] - 4*bz[ig+1] + 6*bz[ig] - 4*bz[ig-1] + bz[ig-2])/h**4

        bx4p[N-2]  = (bext[3]  - 4*bx[N]  + 6*bx[N-1]  - 4*bx[N-2] + bx[N-3] )/h**4
        by4p[N-2]  = (bext[4]  - 4*by[N]  + 6*by[N-1]  - 4*by[N-2] + by[N-3] )/h**4
        bz4p[N-2]  = (bext[5]  - 4*bz[N]  + 6*bz[N-1]  - 4*bz[N-2] + bz[N-3] )/h**4



       #===============================================================================
       #===============================================================================
       #===============================================================================

        # array for closure constraint
        f_c = np.array([0.0,0.0,0.0])
        om = np.array([ (0, -w,  v),
                        (w,  0,  -u),
                        (-v,    u,    0 )  ])

       #===  p =  1  ===================

        p = 1


        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

      #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)


        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

        jac[3*p-3:3*p,5*N1+5:5*N1+8] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]             = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,5*N1+5:5*N1+8]         = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-3:3*p]             =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])

        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]/2-bz[p+1]), (by[p-1]/2-by[p+1])),\
                                                    ((bz[p-1]/2-bz[p+1]),       0,           -(bx[p-1]/2-bx[p+1])),\
                                                    (-(by[p-1]/2-by[p+1]),       (bx[p-1]/2-bx[p+1]),      0) ])
       #======================================================================


       #==== p = 2 =================================================

        p = 2


        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

        #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)


        jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

        #jac[3*p-3:3*p,5*N1+6:5*N1+9] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])

        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                    ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                    (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])



        for p in range(3,N-2):

             gmx = v*bzp[p] - w*byp[p]
             gmy = w*bxp[p] - u*bzp[p]
             gmz = u*byp[p] - v*bxp[p]


             f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
             f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
             f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


             f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

             f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


             f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                                bz[p]*bxp[p] - bx[p]*bzp[p],\
                                                bx[p]*byp[p] - by[p]*bxp[p]])
            #==================== Jacobian entry ================================
             jac[3*p-3:3*p,3*p-9:3*p-6] = id*1*.001/h**4
             jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
             jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
             jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
             jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

             #--- Constraint gradient ---

             #--- unit vector constraint --
             jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
             #--- unit speed constraint --
             jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
             jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

             #-- differentiation with rho
             jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

             #-- differentiation with lambda
             jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                             (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                             (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
             #-- differentiation with u v w
             jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                             ( bzp[p],      0,       -bxp[p]),\
                                                             (-byp[p],     bxp[p],     0    )])

             #---- closure ----

             jac[5*N1+2:5*N1+5,3*p-3:3*p]  =    np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                         ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                         (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])
            #=================================================================================================================================
            #=====================================  For loop ends here =======================================================================
            ##=================================================================================================================================

        #=== Entry at p = N-2
        p = N-2

        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

        #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)

        jac[3*p-3:3*p,3*p-9:3*p-6] = id*1*.001/h**4
        jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        #jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

        #jac[3*p-3:3*p,5*N1+6:5*N1+9] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])

        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                    ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                    (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])
        #--- CLOSURE CONSTRAINT WITH RESPECT TO bext(1:3) and bext(4:6);


        jac[5*N1+2:5*N1+5,5*N1+8:5*N1+11] =   0.5*np.array([(0,     -bz[N],   by[N]),\
                                                            (bz[N],   0,     -bx[N]),\
                                                            (-by[N],  bx[N],   0) ])
        #===============================================================================================================
        #==== Entry at p = N-2 ends here =========
        #===============================================================================================================

        #===============================================================================================================
        #==== Entry at p = N-1 begins here =========
        #===============================================================================================================
        p = N-1

        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

        #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)

        jac[3*p-3:3*p,3*p-9:3*p-6] = id*1*.001/h**4
        jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        #jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        #jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4
        jac[3*p-3:3*p,5*N1+8:5*N1+11] = id*1*.001/h**4
        #jac[3*p-3:3*p,5*N1+6:5*N1+9] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])
        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                    ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                    (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])


        #===============================================================================================================
        #==== Entry at p = N-1 ends here =========
        #===============================================================================================================

        #===============================================================================================================
        #==== Entry at p = N begins here =========
        #===============================================================================================================
        p = N
        jac[4*N1+p-1,3*p-9:3*p-6]   = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])   #differentiation with b(N-1)

        f_b[5*N1+2:5*N1+5] =  h*f_c +   h*np.array([by[p]*bzp[p]  - bz[p]*byp[p],\
                                                    bz[p]*bxp[p] - bx[p]*bzp[p],\
                                                    bx[p]*byp[p] - by[p]*bxp[p] ] )

        #===============================================================================================================
        #==== Entry at p = N ends here =========
        #===============================================================================================================
        # differentiation with bext(4:6)
        p = N+1
        jac[4*N1+p-1,5*N1+8:5*N1+11]    =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])



        #-- |b'(0)| = tau
        f_b[4*N1]   = 1/2*(bxp[0]**2 + byp[0]**2 + bzp[0]**2 - tau**2)

         # |b'(1)| = tau
        f_b[5*N1+1] = 1/2*(bxp[N]**2 + byp[N]**2 + bzp[N]**2 - tau**2)

         # ---- b'(1) + b'(N+1) = 0 ----

        f_b[5*N1+5] = bxp[0] + bxp[N]    # bxp1**2  + byp1**2  + bzp1**2  - tau**2 ;
        f_b[5*N1+6] = byp[0] + byp[N]    #bxpN1**2 + bypN1**2 + bzpN1**2 - tau**2 ;
        f_b[5*N1+7] = bzp[0] + bzp[N]    #

        jac[5*N1+5:5*N1+8,0:3]           =  1/(2*h)*id
        jac[5*N1+5:5*N1+8,5*N1+5:5*N1+8] = -1/(2*h)*id

        jac[5*N1+5:5*N1+8,3*N1-3:3*N1]   = -1/(2*h)*id
        jac[5*N1+5:5*N1+8,5*N1+8:5*N1+11] = 1/(2*h)*id
        # ---- b''(1) + b''(N+1) = 0 ----

        f_b[5*N1+8]   = bx2p[0] + bx2p[N]
        f_b[5*N1+9]   = by2p[0] + by2p[N]
        f_b[5*N1+10]  = bz2p[0] + bz2p[N]

        jac[5*N1+8:5*N1+11,0:3]            = 1/h**2*id
        jac[5*N1+8:5*N1+11,5*N1+5:5*N1+8]  = 1/h**2*id

        jac[5*N1+8:5*N1+11,3*N1-3:3*N1]     = 1/h**2*id
        jac[5*N1+8:5*N1+11,5*N1+8:5*N1+11] = 1/h**2*id



        J = jac.copy()
        F = f_b.copy()

        #=== deleting the entries corresponding to bz[1] and bz[N-1]

        F = np.delete(F,2)
        F = np.delete(F,3*N1-2)

        J = np.delete(J,2,0)
        J = np.delete(J,2,1)
        J = np.delete(J,3*N1-2,0)
        J = np.delete(J,3*N1-2,1)



        print(sum(F**2))


        return F

    #===============================================================================
    #======================    fun_jac ======================================
    #===============================================================================


    def fun_jac(x):


        global N,N1,u,v,w,f_b,h,tau,rho,lm
        global bx,by,bz,bxp,byp,bzp,bx2p,by2p,bz2p,bx4p,by4p,bz4p,bext,f_b,tau,jac,ind,lmp


        #=== Initialize arrays ===



        bx[1] = x[0]
        by[1] = x[1]

        bx[2:N]     = x[2:3*N1-3:3]
        by[2:N]     = x[3:3*N1-2:3]
        bz[2:N-1]   = x[4:3*N1-4:3]


        rho         = x[3*N1-2:4*N1-2]
        lm          = x[4*N1-2:5*N1]

        u           = x[5*N1]
        v           = x[5*N1+1]
        w           = x[5*N1+2]

        bext        = x[5*N1+3:5*N1+9]


        #=== calculation of derivatives ===

        ig = np.arange(1,N)

        #===========================================================================
        # ========= b' ==== O(h**2) ===========
        #===========================================================================
        bxp[0]  = (bx[1]-bext[0])/(2*h)
        byp[0]  = (by[1]-bext[1])/(2*h)
        bzp[0]  = (bz[1]-bext[2])/(2*h)

        bxp[ig] = (bx[ig+1]-bx[ig-1])/(2*h)
        byp[ig] = (by[ig+1]-by[ig-1])/(2*h)
        bzp[ig] = (bz[ig+1]-bz[ig-1])/(2*h)

        bxp[N]  = (bext[3]-bx[N-1])/(2*h)
        byp[N]  = (bext[4]-by[N-1])/(2*h)
        bzp[N]  = (bext[5]-bz[N-1])/(2*h)

        lmp = (lm[ig+1]-lm[ig-1])/(2*h)
        #===========================================================================
        #===========================================================================


        bx2p[0] = (bx[1] + bext[0] -2*bx[0])/h**2
        by2p[0] = (by[1] + bext[1] -2*by[0])/h**2
        bz2p[0] = (bz[1] + bext[2] -2*bz[0])/h**2

        bx2p[ig] = (bx[ig+1] + bx[ig-1] - 2*bx[ig])/h**2
        by2p[ig] = (by[ig+1] + by[ig-1] - 2*by[ig])/h**2
        bz2p[ig] = (bz[ig+1] + bz[ig-1] - 2*bz[ig])/h**2

        bx2p[N] = (bext[3] + bx[N-1]  - 2*bx[N])/h**2
        by2p[N] = (bext[4] + by[N-1]  - 2*by[N])/h**2
        bz2p[N] = (bext[5] + bz[N-1]  - 2*bz[N])/h**2

        #===========================================================================
        #===========================================================================

        bx4p[0] = (bx[3] - 4*bx[2] + 6*bx[1] - 4*bx[0] + bext[0])/h**4
        by4p[0] = (by[3] - 4*by[2] + 6*by[1] - 4*by[0] + bext[1])/h**4
        bz4p[0] = (bz[3] - 4*bz[2] + 6*bz[1] - 4*bz[0] + bext[2])/h**4


        ig = np.arange(2,N-1)

        bx4p[ig-1] = (bx[ig+2] - 4*bx[ig+1] + 6*bx[ig] - 4*bx[ig-1] + bx[ig-2])/h**4
        by4p[ig-1] = (by[ig+2] - 4*by[ig+1] + 6*by[ig] - 4*by[ig-1] + by[ig-2])/h**4
        bz4p[ig-1] = (bz[ig+2] - 4*bz[ig+1] + 6*bz[ig] - 4*bz[ig-1] + bz[ig-2])/h**4

        bx4p[N-2]  = (bext[3]  - 4*bx[N]  + 6*bx[N-1]  - 4*bx[N-2] + bx[N-3] )/h**4
        by4p[N-2]  = (bext[4]  - 4*by[N]  + 6*by[N-1]  - 4*by[N-2] + by[N-3] )/h**4
        bz4p[N-2]  = (bext[5]  - 4*bz[N]  + 6*bz[N-1]  - 4*bz[N-2] + bz[N-3] )/h**4



       #===============================================================================
       #===============================================================================
       #===============================================================================

        # array for closure constraint
        f_c = np.array([0.0,0.0,0.0])
        om = np.array([ (0, -w,  v),
                        (w,  0,  -u),
                        (-v,    u,    0 )  ])

       #===  p =  1  ===================

        p = 1


        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

      #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)


        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

        jac[3*p-3:3*p,5*N1+5:5*N1+8] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]             = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,5*N1+5:5*N1+8]         = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-3:3*p]             =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])

        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]/2-bz[p+1]), (by[p-1]/2-by[p+1])),\
                                                    ((bz[p-1]/2-bz[p+1]),       0,           -(bx[p-1]/2-bx[p+1])),\
                                                    (-(by[p-1]/2-by[p+1]),       (bx[p-1]/2-bx[p+1]),      0) ])
       #======================================================================


       #==== p = 2 =================================================

        p = 2


        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

        #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)


        jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

        #jac[3*p-3:3*p,5*N1+6:5*N1+9] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])

        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                    ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                    (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])



        for p in range(3,N-2):

             gmx = v*bzp[p] - w*byp[p]
             gmy = w*bxp[p] - u*bzp[p]
             gmz = u*byp[p] - v*bxp[p]


             f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
             f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
             f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


             f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

             f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


             f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                                bz[p]*bxp[p] - bx[p]*bzp[p],\
                                                bx[p]*byp[p] - by[p]*bxp[p]])
            #==================== Jacobian entry ================================
             jac[3*p-3:3*p,3*p-9:3*p-6] = id*1*.001/h**4
             jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
             jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
             jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
             jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

             #--- Constraint gradient ---

             #--- unit vector constraint --
             jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
             #--- unit speed constraint --
             jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
             jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

             #-- differentiation with rho
             jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

             #-- differentiation with lambda
             jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                             (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                             (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
             #-- differentiation with u v w
             jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                             ( bzp[p],      0,       -bxp[p]),\
                                                             (-byp[p],     bxp[p],     0    )])

             #---- closure ----

             jac[5*N1+2:5*N1+5,3*p-3:3*p]  =    np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                         ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                         (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])
            #=================================================================================================================================
            #=====================================  For loop ends here =======================================================================
            ##=================================================================================================================================

        #=== Entry at p = N-2
        p = N-2

        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

        #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)

        jac[3*p-3:3*p,3*p-9:3*p-6] = id*1*.001/h**4
        jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        #jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4

        #jac[3*p-3:3*p,5*N1+6:5*N1+9] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])

        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                    ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                    (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])
        #--- CLOSURE CONSTRAINT WITH RESPECT TO bext(1:3) and bext(4:6);


        jac[5*N1+2:5*N1+5,5*N1+8:5*N1+11] =   0.5*np.array([(0,     -bz[N],   by[N]),\
                                                            (bz[N],   0,     -bx[N]),\
                                                            (-by[N],  bx[N],   0) ])
        #===============================================================================================================
        #==== Entry at p = N-2 ends here =========
        #===============================================================================================================

        #===============================================================================================================
        #==== Entry at p = N-1 begins here =========
        #===============================================================================================================
        p = N-1

        gmx = v*bzp[p] - w*byp[p]
        gmy = w*bxp[p] - u*bzp[p]
        gmz = u*byp[p] - v*bxp[p]


        f_b[3*p-3] = .001*bx4p[p-1] + lmp[p-1]*bxp[p] + lm[p]*bx2p[p] - rho[p-1]*bx[p] + gmx
        f_b[3*p-2] = .001*by4p[p-1] + lmp[p-1]*byp[p] + lm[p]*by2p[p] - rho[p-1]*by[p] + gmy
        f_b[3*p-1] = .001*bz4p[p-1] + lmp[p-1]*bzp[p] + lm[p]*bz2p[p] - rho[p-1]*bz[p] + gmz


        f_b[3*N1+p-1] = 1/2*(bx[p]**2 + by[p]**2 + bz[p]**2 -1)

        f_b[4*N1+p]   = 1/2*(bxp[p]**2 + byp[p]**2 + bzp[p]**2 - tau**2)


        f_c  = f_c     +         np.array([by[p]*bzp[p] - bz[p]*byp[p],\
                                           bz[p]*bxp[p] - bx[p]*bzp[p],\
                                           bx[p]*byp[p] - by[p]*bxp[p]])

        #==================== Jacobian entry ================================

        ind = [3*p-3,3*p-2,3*p-1]#np.arange(3*p-3,3*p)

        jac[3*p-3:3*p,3*p-9:3*p-6] = id*1*.001/h**4
        jac[3*p-3:3*p,3*p-6:3*p-3] = id*(-4*.001/h**4 -lmp[p-1]/(2*h) +lm[p]/h**2) - om/(2*h)
        jac[3*p-3:3*p,3*p-3:3*p]   = id*(6*.001/h**4-2*lm[p]/h**2-rho[p-1])
        #jac[3*p-3:3*p,3*p:3*p+3]   = id*(-4*.001/h**4 +lmp[p-1]/(2*h) +lm[p]/h**2) + om/(2*h)
        #jac[3*p-3:3*p,3*p+3:3*p+6] = id*1*.001/h**4
        jac[3*p-3:3*p,5*N1+8:5*N1+11] = id*1*.001/h**4
        #jac[3*p-3:3*p,5*N1+6:5*N1+9] = id*1*.001/h**4

        #--- Constraint gradient ---

        #--- unit vector constraint --
        jac[3*N1+p-1,3*p-3:3*p]         = np.array([bx[p], by[p], bz[p]])
        #--- unit speed constraint --
        jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-3:3*p]         =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])

        #-- differentiation with rho

        jac[3*p-3:3*p,3*N1+p-1] = -1*np.array([bx[p],by[p],bz[p]])

        #-- differentiation with lambda
        jac[3*p-3:3*p,4*N1+p-1:4*N1+p+2]   = np.array([ (-1/(2*h)*bxp[p],  bx2p[p], 1/(2*h)*bxp[p]),\
                                                        (-1/(2*h)*byp[p],  by2p[p], 1/(2*h)*byp[p]),\
                                                        (-1/(2*h)*bzp[p],  bz2p[p], 1/(2*h)*bzp[p]) ])
        #-- differentiation with u v w
        jac[3*p-3:3*p,5*N1+2:5*N1+5]       = -1*np.array([ (0,          -bzp[p],   byp[p]), \
                                                        ( bzp[p],      0,       -bxp[p]),\
                                                        (-byp[p],     bxp[p],     0    )])
        #---- closure ----

        jac[5*N1+2:5*N1+5,3*p-3:3*p]  =   np.array([(0,                 -(bz[p-1]-bz[p+1]), (by[p-1]-by[p+1])),\
                                                    ((bz[p-1]-bz[p+1]),       0,           -(bx[p-1]-bx[p+1])),\
                                                    (-(by[p-1]-by[p+1]),       (bx[p-1]-bx[p+1]),      0) ])


        #===============================================================================================================
        #==== Entry at p = N-1 ends here =========
        #===============================================================================================================

        #===============================================================================================================
        #==== Entry at p = N begins here =========
        #===============================================================================================================
        p = N
        jac[4*N1+p-1,3*p-9:3*p-6]   = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])   #differentiation with b(N-1)

        f_b[5*N1+2:5*N1+5] =  h*f_c +   h*np.array([by[p]*bzp[p]  - bz[p]*byp[p],\
                                                    bz[p]*bxp[p] - bx[p]*bzp[p],\
                                                    bx[p]*byp[p] - by[p]*bxp[p] ] )

        #===============================================================================================================
        #==== Entry at p = N ends here =========
        #===============================================================================================================
        # differentiation with bext(4:6)
        p = N+1
        jac[4*N1+p-1,5*N1+8:5*N1+11]    =  1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])
        jac[4*N1+p-1,3*p-9:3*p-6]       = -1/(2*h)*np.array([bxp[p-1], byp[p-1], bzp[p-1]])



        #-- |b'(0)| = tau
        f_b[4*N1]   = 1/2*(bxp[0]**2 + byp[0]**2 + bzp[0]**2 - tau**2)

         # |b'(1)| = tau
        f_b[5*N1+1] = 1/2*(bxp[N]**2 + byp[N]**2 + bzp[N]**2 - tau**2)

         # ---- b'(1) + b'(N+1) = 0 ----

        f_b[5*N1+5] = bxp[0] + bxp[N]    # bxp1**2  + byp1**2  + bzp1**2  - tau**2 ;
        f_b[5*N1+6] = byp[0] + byp[N]    #bxpN1**2 + bypN1**2 + bzpN1**2 - tau**2 ;
        f_b[5*N1+7] = bzp[0] + bzp[N]    #

        jac[5*N1+5:5*N1+8,0:3]           =  1/(2*h)*id
        jac[5*N1+5:5*N1+8,5*N1+5:5*N1+8] = -1/(2*h)*id

        jac[5*N1+5:5*N1+8,3*N1-3:3*N1]   = -1/(2*h)*id
        jac[5*N1+5:5*N1+8,5*N1+8:5*N1+11] = 1/(2*h)*id
        # ---- b''(1) + b''(N+1) = 0 ----

        f_b[5*N1+8]   = bx2p[0] + bx2p[N]
        f_b[5*N1+9]   = by2p[0] + by2p[N]
        f_b[5*N1+10]  = bz2p[0] + bz2p[N]

        jac[5*N1+8:5*N1+11,0:3]            = 1/h**2*id
        jac[5*N1+8:5*N1+11,5*N1+5:5*N1+8]  = 1/h**2*id

        jac[5*N1+8:5*N1+11,3*N1-3:3*N1]     = 1/h**2*id
        jac[5*N1+8:5*N1+11,5*N1+8:5*N1+11] = 1/h**2*id



        J = jac.copy()
        F = f_b.copy()

        #=== deleting the entries corresponding to bz[1] and bz[N-1]

        F = np.delete(F,2)
        F = np.delete(F,3*N1-2)

        J = np.delete(J,2,0)
        J = np.delete(J,2,1)
        J = np.delete(J,3*N1-2,0)
        J = np.delete(J,3*N1-2,1)



        #print(sum(F**2))


        return J




    #===============================================================================
    #======================    fun_lambda ======================================
    #===============================================================================

    def fun_lambda(bx,by,bz,bext):
        #=== calculation of derivatives ===
        global N,N1,h,u,v,w
        global bxp,byp,bzp,bx2p,by2p,bz2p,bx4p,by4p,bz4p,tau

        ig = np.arange(1,N)

        #===========================================================================
        # ========= b' ==== O(h**2) ===========
        #===========================================================================
        bxp[0]  = (bx[1]-bext[0])/(2*h)
        byp[0]  = (by[1]-bext[1])/(2*h)
        bzp[0]  = (bz[1]-bext[2])/(2*h)

        bxp[ig] = (bx[ig+1]-bx[ig-1])/(2*h)
        byp[ig] = (by[ig+1]-by[ig-1])/(2*h)
        bzp[ig] = (bz[ig+1]-bz[ig-1])/(2*h)

        bxp[N]  = (bext[3]-bx[N-1])/(2*h)
        byp[N]  = (bext[4]-by[N-1])/(2*h)
        bzp[N]  = (bext[5]-bz[N-1])/(2*h)


        #===========================================================================
        #===========================================================================


        bx2p[0] = (bx[1] + bext[0] -2*bx[0])/h**2
        by2p[0] = (by[1] + bext[1] -2*by[0])/h**2
        bz2p[0] = (bz[1] + bext[2] -2*bz[0])/h**2

        bx2p[ig] = (bx[ig+1] + bx[ig-1] - 2*bx[ig])/h**2
        by2p[ig] = (by[ig+1] + by[ig-1] - 2*by[ig])/h**2
        bz2p[ig] = (bz[ig+1] + bz[ig-1] - 2*bz[ig])/h**2

        bx2p[N] = (bext[3] + bx[N-1]  - 2*bx[N])/h**2
        by2p[N] = (bext[4] + by[N-1]  - 2*by[N])/h**2
        bz2p[N] = (bext[5] + bz[N-1]  - 2*bz[N])/h**2

        #===========================================================================
        #===========================================================================

        bx4p[0] = (bx[3] - 4*bx[2] + 6*bx[1] - 4*bx[0] + bext[0])/h**4
        by4p[0] = (by[3] - 4*by[2] + 6*by[1] - 4*by[0] + bext[1])/h**4
        bz4p[0] = (bz[3] - 4*bz[2] + 6*bz[1] - 4*bz[0] + bext[2])/h**4


        ig = np.arange(2,N-1)

        bx4p[ig-1] = (bx[ig+2] - 4*bx[ig+1] + 6*bx[ig] - 4*bx[ig-1] + bx[ig-2])/h**4
        by4p[ig-1] = (by[ig+2] - 4*by[ig+1] + 6*by[ig] - 4*by[ig-1] + by[ig-2])/h**4
        bz4p[ig-1] = (bz[ig+2] - 4*bz[ig+1] + 6*bz[ig] - 4*bz[ig-1] + bz[ig-2])/h**4

        bx4p[N-2]  = (bext[3]  - 4*bx[N]  + 6*bx[N-1]  - 4*bx[N-2] + bx[N-3] )/h**4
        by4p[N-2]  = (bext[4]  - 4*by[N]  + 6*by[N-1]  - 4*by[N-2] + by[N-3] )/h**4
        bz4p[N-2]  = (bext[5]  - 4*bz[N]  + 6*bz[N-1]  - 4*bz[N-2] + bz[N-3] )/h**4

        lm = 1.5*.001*(bx2p*bx2p + by2p*by2p + bz2p*bz2p)/tau**2 - tau

        rho = .001*(bx4p*bx[1:N] + by4p*by[1:N] + bz4p*bz[1:N]) - (tau**2)* lm[1:N]\
              +(v*bz[1:N]-w*by[1:N])*bx[1:N] + (w*bx[1:N]-u*bz[1:N])*by[1:N]\
              +(u*by[1:N]-v*bx[1:N])*bz[1:N]

        return lm,rho



    #===============================================================================
    #====================  Function ends here
    #===============================================================================




    #===============================================================================

    #===============================================================================


    #===============================================================================
    #======================    Solving the equations ===============================
    #===============================================================================


    x0 = initial_guess()

    #x2 = fsolve(fun_curve,x0,fprime = fun_jac,xtol=10**(-7),maxfev=1000 )
    #x2 = optimize.leastsq(fun_curve, x0 )
    #x2 = optimize.root(fun_curve,x0,method = 'lm',jac=fun_jac, options={'xatol': 1e-8, 'disp': True})
    x2 = optimize.least_squares(fun_curve,x0,method = 'lm', jac = fun_jac, ftol=1.49012e-25)
    #===============================================================================
    #======================    Post processing ======================================
    #===============================================================================

    '''
    fig1 = plt.figure()
    ax = fig1.add_subplot(111,projection='3d')
    ax.plot(bx,by,bz)
    plt.show()
    '''

    #plot3d(bx,by,bz,np.sin(bz), tube_radius=0.025, colormap='Spectral')
    #plot3d(bx,by,bz, tube_radius=0.025)
    #show()
    #===============================================================================
    #======================     main solver   ======================================
    #===============================================================================

    #===============================================================================
    #======================    fun_curve ======================================
    #===============================================================================

    #==============================================================================
    #======================= saving data file ====================================

    path ='/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/'

    tau2    = int(np.around(tau*(10**10),decimals=0))
    #str0    = '7pi_knot_N' +str(len(bx)-1)+ '_tau_'+ str(tau2) + '.txt'
    str0 = '3fold_N72_tau_' + str(tau2) + '.txt'
    strb    = path + 'b_' + str0
    strlm   = path + 'lm_' + str0
    strrho  = path + 'rho_' + str0
    struvw  = path + 'uvw_' + str0
    strbext = path + 'bext_' + str0


    #=== saving b ==
    with open(strb, 'w') as fh:
        for i in range(0,len(bx)):
            #fh.write('%25s  %25s %s\n' %(bx[i],  by[i],bz[i],fmt="%.8g" ))
            fh.write('%30.16E %30.16E %30.16E \r\n' % (bx[i],by[i],bz[i]))
    fh.close()
    #=== saving lm ==
    with open(strlm, 'w') as fh:
        for i in range(0,len(lm)):
            #fh.write('%25s  %25s %s\n' %(bx[i],  by[i],bz[i],fmt="%.8g" ))
            fh.write('%30.16E  \r\n' % (lm[i]))
    fh.close()
    #=== saving rho ==
    with open(strrho, 'w') as fh:
        for i in range(0,len(rho)):
            #fh.write('%25s  %25s %s\n' %(bx[i],  by[i],bz[i],fmt="%.8g" ))
            fh.write('%30.16E \r\n' % (rho[i]))
    fh.close()
    #=== saving uvw ==
    with open(struvw, 'w') as fh:
            fh.write('%30.16E \r\n' % (u))
            fh.write('%30.16E \r\n' % (v))
            fh.write('%30.16E \r\n' % (w))
    fh.close()
    #=== saving bext ==
    with open(strbext, 'w') as fh:
        for i in range(0,6):
            #fh.write('%25s  %25s %s\n' %(bx[i],  by[i],bz[i],fmt="%.8g" ))
            fh.write('%30.16E \r\n' % (bext[i]))
    fh.close()


    F = x2.fun
    err[p2] = sum(F**2)

    print(p2)
#createplots(strb)
#=== saving error ====


with open(path+'err.txt', 'w') as fh:
         for i in range(1,len(tau1)):
             #fh.write('%25s  %25s %s\n' %(bx[i],  by[i],bz[i],fmt="%.8g" ))
             fh.write('%30.16E \r\n' % (err[i]))
fh.close()
