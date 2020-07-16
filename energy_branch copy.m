%=== calculate energy of equilibrium solution for given tau 
clear all
clc

global N tau sig N1

N = 105;
N1 = N-1;
h = 1/N;
tau = 16;
sig = 0.01;


branch =2;

path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
strtau = ['tau_branch_' num2str(branch) '.txt'];
tau1 = load([path strtau]);
[val,ind] = min(abs(tau1-tau));
tau2 = tau1(ind);               % nearest value to input tau for which we have the data
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];
x = load([path str0]);  % initial point for pseudodynamics


    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;

    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
   
  bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
  by(2:N,1)     = x(2:3:3*N1-1,1)    ;
  bz(2:N,1)     = x(3:3:3*N1,1)      ;
  
  
  

E(branch) = energy_b(bx,by,bz);




%==== branch 1 ===

branch =1;

path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
strtau = ['tau_branch_' num2str(branch) '.txt'];
tau1 = load([path strtau]);
[val,ind] = min(abs(tau1-tau));
tau2 = tau1(ind);               % nearest value to input tau for which we have the data
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];
x = load([path str0]);  % initial point for pseudodynamics


    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;

    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
   
  bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
  by(2:N,1)     = x(2:3:3*N1-1,1)    ;
  bz(2:N,1)     = x(3:3:3*N1,1)      ;
  
  

E(branch) = energy_b(bx,by,bz);

  