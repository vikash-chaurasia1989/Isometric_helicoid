function [] = parameters()

global N N1 h tau c ht tspan nnt


N = 105;
N1 = N-1;
h = 1/N;
tau = 16;
c =  10^-5;





ht = 0.002;  % time step size
nnt =100;
tspan = 0:ht:nnt*ht;  % timespan for pseudodynamics



end