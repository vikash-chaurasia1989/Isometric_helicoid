%=================

clear all
clc

 global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 bxm1 bym1 bzm1 d id bx0 by0 bz0

 global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext tau id sig E path branch 
 
 

  format longE
  
  
  parameters5();
  
  %=== calculate the Frechet distance between bx0 and the solution at the
%previous step
%=== calculate the Frechet distance between bx0 and the solution at the
%previous step
p1= 1;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'] ;

%== previous steps ==
x0 =  load([path str0]);

bxm1(2:N,1)     = x0(1:3:3*N1-2,1)    ;
bym1(2:N,1)     = x0(2:3:3*N1-1,1)    ;
bzm1(2:N,1)     = x0(3:3:3*N1,1)      ;

bxm1(1,1) = 1;
bym1(1,1) = 0;
bzm1(1,1) = 0;

bxm1(N+1,1) = -1;
bym1(N+1,1) =  0;
bzm1(N+1,1) =  0;

bx0 = bxm1;
by0 = bym1;
bz0 = bzm1;

id = eye(3,3);
%===== initial guess for the current step ==

 

str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1)   '.txt'] ;


%== previous steps ==
x0 =  load([path str0]);

 
bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;
% initial guess for x2 =====

bx(2:N,1)     = x0(1:3:3*N1-2,1)    ;
by(2:N,1)     = x0(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x0(3:3:3*N1,1)      ;


d =  4*sum((bxm1-bx).^2 + (bym1-by).^2 + (bzm1-bz).^2) ;


var_in = zeros(5*N1+4,1);

var_in(1:3*N1,1) = x0(1:3*N1,1); 



[F,J] = fun_pullback2(var_in);
%
% %== Constraint gradient
%
A = J(3*N1+1:5*N1+4,1:3*N1);
var_in(3*N1+1:5*N1+4,1) = pinv(A')*F(1:3*N1,1);

var_in(5*N1+4,1) = 100;

[F,J] = fun_pullback2(var_in);

j = jacobian(@curve_pullback4,var_in) ;

F1 = curve_pullback4(var_in);

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

% tau = 11;
% 
% %=== numerical ===
% 
%  x = initial_guess_coil();
% 
%  [F,J] = fun_orientable(x);
%  
% %f = fun_evolution_c(x);
% 
% %x = initial_guess_loop8();
% 
% j = jacobian(@fun_orientableC,x);
% 
% 
% ind = 1:3*N1;
% 
% plot(J(ind,1)-j(ind,1));
