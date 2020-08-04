function [] = parameters5()

global N N1 h tau sig fac branch ht bt gamma t1 nstep strpath path bx1 by1 bz1 bx2 by2 bz2 d21


 strpath = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
 %strpath  = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';

 
 path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent6/'    ;

 % path = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent6/';


N = 105;
N1 = N-1;
h = 1/N;
tau = 16;

%== prefactor for aspect ratio
sig = 0.01;
n = tau/2/pi;
fac = (asinh(n*pi*sig))/(4*pi^3*n^3);

branch =2;
%==========================================================================
ht = 0.001;

bt    =  -5*10^-5;
gamma =  -1*10^-5;

t1 = 0:ht:1;

nstep = 4000;%2500;


branch = 1;
path1 = [strpath num2str(branch) '/'];

str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];

x = load([path1 str0]);

bx1(2:N,1)     = x(1:3:3*N1-2,1)    ;
by1(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz1(2:N,1)     = x(3:3:3*N1,1)      ;

bx1(1,1) = 1;
by1(1,1) = 0;
bz1(1,1) = 0;

bx1(N+1,1) = -1;
by1(N+1,1) =  0;
bz1(N+1,1) =  0;

branch = 2;
path1 = [strpath num2str(branch) '/'];

str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];

x = load([path1 str0]);

bx2(2:N,1)     = x(1:3:3*N1-2,1)    ;
by2(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz2(2:N,1)     = x(3:3:3*N1,1)      ;

bx2(1,1) = 1;
by2(1,1) = 0;
bz2(1,1) = 0;

bx2(N+1,1) = -1;
by2(N+1,1) =  0;
bz2(N+1,1) =  0;

%== Distance between the equiilbrium points on branch 1 and 2 ===
d21 =  sqrt(sum((bx1-bx2).^2 + (by1-by2).^2 + (bz1-bz2).^2)) ;



end