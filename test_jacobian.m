%=================

clear all
clc

 global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

 global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext tau id

  format longE
  
% tau = 8.1;
% N = 77;
% N1 = N-1;

id = eye(3,3);

x = initial_guess_loop8();


%=== analytical ===
[F,J] = fun_jacobian8(x);

%=== numerical ===

x = initial_guess_loop8();

f = fun_curve8(x);

x = initial_guess_loop8();

j = jacobian(@fun_curve8,x);


ind = 1:3*N1;

plot(J(ind,1)-j(ind,1));
