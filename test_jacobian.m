%=================

clear all
clc

 global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

 global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext tau id sig E

  format longE
  
% tau = 8.1;
% N = 77;
% N1 = N-1;

id = eye(3,3);

%x = initial_guess_loop8();

% 
% path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';
% 
% %N = 72;
% tau = 16;
% 
% N = 105;
% h = 1 / N
% N1 = N - 1;
% sig = 0.01;
% 
% tau2 = tau;
% 
% str1 = ['5pi_knot_N' num2str(N) '_tau_' num2str(10^10 * tau2) '.txt'];
% str2 = ['3fold_N' num2str(N) '_tau_' num2str(10^10 * tau2) '.txt'];
% str3 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10 * tau2) '.txt'];
% strb1 = [path 'b_' str1];
% strb2 = [path 'b_' str2];
% strb3 = [path 'b_' str3];
% %==== load b ===========
% temp = load(strb1);
% 
% bx  = temp(:, 1);
% by  = temp(:, 2);
% bz  = temp(:, 3);
% 
% 
% E = energy_b(bx,by,bz);
% 
% 
% x(1:3:3 * N1 - 2, 1) = bx(2:N, 1);
% x(2:3:3 * N1 - 1, 1) = by(2:N, 1);
% x(3:3:3 * N1, 1)     = bz(2:N, 1);
% 
% %=== analytical ===
% [F,J] = fun_evolution(x);
% 

tau = 11;

%=== numerical ===

 x = initial_guess_coil();

 [F,J] = fun_orientable(x);
 
%f = fun_evolution_c(x);

%x = initial_guess_loop8();

j = jacobian(@fun_orientableC,x);


ind = 1:3*N1;

plot(J(ind,1)-j(ind,1));
