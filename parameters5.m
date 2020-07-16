function [] = parameters5()

global N N1 h tau sig fac branch ht bt gamma t1 

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

bt = 5*10^-5;
gamma = 1*10^-5;

t1 = 0:ht:1;

end