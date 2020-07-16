clear all
clc

global    N1  tau   h tx ty tz  s nx ny nz k0 al w

global bx by bz   N     rx ry rz   bext or kappa A B tau1
set(0,'defaultlinelinewidth',2)


or =    -1;

branch = 1;
N = 72;
N1 = N-1;
h = 1/N;
tau =8.3; 

% strpath = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
% %strpath = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
 path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
% 
  strtau = ['tau_branch_' num2str(branch) '.txt'];
% 
 tau2 = load([path strtau]);
% 
 [val,ind] = min(abs(tau2-tau));
% 
 tau3 = tau2(ind);
  
 tau = tau3;

  
 
str0 = ['b_3fold_N72_tau_' num2str(tau*10^10) '.txt'];

x = load([path str0]);

%----- reading input and constructing bx1 by1 bz1 ----


% 
% bx1(2:N,1)     = x(1:3:3*N1-2,1)    ;
% by1(2:N,1)     = x(2:3:3*N1-1,1)    ;
% bz1(2:N,1)     = x(3:3:3*N1,1)      ;


bx1 = x(:,1);
by1 = x(:,2);
bz1 = x(:,3);

% N2 = 120;
% s2 = linspace(0,1,N+1)';
% s1 = linspace(0,1,N2+1)';
% 
% bx1 = interp1q(s2,bx2,s1);
% by1 = interp1q(s2,by2,s1);
% bz1 = interp1q(s2,bz2,s1);

%N = N2;
N1 = N-1;
h = 1/N;
bx1(1,1) = 1;
by1(1,1) = 0;
bz1(1,1) = 0;

bx1(N+1,1) = -1;
by1(N+1,1) =  0;
bz1(N+1,1) =  0;


bext(1:3,1) = -1*[bx1(N,1);by1(N,1);bz1(N,1)];
bext(4:6,1) = -1*[bx1(2,1);by1(2,1);bz1(2,1)];

%---- derivative calculation for b2 --- bN  ---

ig = 2:N;

%=== O(h)====
bx1p(1,1) = (bx1(1,1)-bext(1,1))/(h);
by1p(1,1) = (by1(1,1)-bext(2,1))/(h);
bz1p(1,1) = (bz1(1,1)-bext(3,1))/(h);

bx1p(ig,1) = (bx1(ig,1)-bx1(ig-1,1))/(h);
by1p(ig,1) = (by1(ig,1)-by1(ig-1,1))/(h);
bz1p(ig,1) = (bz1(ig,1)-bz1(ig-1,1))/(h);


bx1p(N+1,1) = (bx1(N+1,1)-bx1(N,1))/(h);
by1p(N+1,1) = (by1(N+1,1)-by1(N,1))/(h);
bz1p(N+1,1) = (bz1(N+1,1)-bz1(N,1))/(h);

  %--  b''  (O(h^2))----
   
   bx12p(1,1) = (bx1(2,1) + bext(1,1) -2*bx1(1,1))/h^2 ;
   by12p(1,1) = (by1(2,1) + bext(2,1) -2*by1(1,1))/h^2 ;
   bz12p(1,1) = (bz1(2,1) + bext(3,1) -2*bz1(1,1))/h^2 ;
   
   bx12p(ig,1) = (bx1(ig+1,1)  + bx1(ig-1,1) - 2*bx1(ig,1))/h^2 ;
   by12p(ig,1) = (by1(ig+1,1)  + by1(ig-1,1) - 2*by1(ig,1))/h^2 ;
   bz12p(ig,1) = (bz1(ig+1,1)  + bz1(ig-1,1) - 2*bz1(ig,1))/h^2 ;
   
   bx12p(N+1,1) = (bext(4,1) + bx1(N,1)   - 2*bx1(N+1,1))/h^2 ;
   by12p(N+1,1) = (bext(5,1) + by1(N,1)   - 2*by1(N+1,1))/h^2 ;
   bz12p(N+1,1) = (bext(6,1) + bz1(N,1)   - 2*bz1(N+1,1))/h^2 ;
%      

kappa2 = sqrt((bx12p.*bx12p + by12p.*by12p + bz12p.*bz12p -tau^4)/tau^2);

i = 1:N+1   ;

%
tx1(i,1) = by1(i,1).*bz1p(i,1) - bz1(i).*by1p(i,1);
ty1(i,1) = bz1(i,1).*bx1p(i,1) - bx1(i).*bz1p(i,1);
tz1(i,1) = bx1(i,1).*by1p(i,1) - by1(i).*bx1p(i,1);
tx1 = tx1/tau;
ty1 = ty1/tau;
tz1 = tz1/tau;

%--------------------

% %--- position vector using integration of tangent
% % initializ1ation
%h=1;
rx1(1,1) = 0; ry1(1,1) = 0; rz1(1,1) = 0;



for i = 1:N
    
    rx1(i+1,1) = rx1(i,1) + h*tx1(i+1,1);
    ry1(i+1,1) = ry1(i,1) + h*ty1(i+1,1);
    rz1(i+1,1) = rz1(i,1) + h*tz1(i+1,1);
    
end


%==========================================================================
% 
rx = zeros(N+1,1);
ry = zeros(N+1,1);
rz = zeros(N+1,1);
tx = zeros(N+1,1);
ty = zeros(N+1,1);
tz = zeros(N+1,1);
nx = zeros(N+1,1);
ny = zeros(N+1,1);
nz = zeros(N+1,1);
bx = zeros(N+1,1);
by = zeros(N+1,1);
bz = zeros(N+1,1);


%== initial point same as the saved data 

tx(1,1) = tx1(1,1);
ty(1,1) = ty1(1,1);
tz(1,1) = tz1(1,1);

bx(1,1) = bx1(1,1);
by(1,1) = by1(1,1);
bz(1,1) = bz1(1,1);

nx(1,1) = by(1,1)*tz(1,1) - bz(1,1)*ty(1,1);
ny(1,1) = bz(1,1)*tx(1,1) - bx(1,1)*tz(1,1);
nz(1,1) = bx(1,1)*ty(1,1) - by(1,1)*tx(1,1);

% rx = rx';
% ry = ry';
% rz = rz';
% 
% tx = tx';
% ty = ty';
% tz = tz';
% 
% nx = nx';
% ny = ny';
% nz = nz';
% 
% bx = bx';
% by = by';
% bz = bz';
 
s1 = linspace(0,1,N+1)';
N2 = 120;
s2 = linspace(0,1,N2+1)';

kappa3 = interp1q(s1,kappa2,s2);

kappa3 = kappa3-min(kappa3);

kappa3(N2/6+1:(N2/6+N2/3)) = - kappa3(N2/6+1:(N2/6+N2/3));
kappa3(N2-N2/6+1:N2+1) = -kappa3(N2-N2/6+1:N2+1);

 N = N2;
N1 = N-1;
h = 1/N;
 
tau1 = tau;
var_in = [3,14^2 ];

%var_in = [1.324912993438339e+01, 4.811607621617519e+01 , 1.302531929059680e+00];

%--------------------------- Solver ---------------------------------------


algo = 1

stralgo = {'levenberg-marquardt' ,'trust-region-reflective' ,'trust-region-dogleg'};

options             = optimset('Display','iter', 'Algorithm',stralgo{algo},'Jacobian','off', 'TOlFun',10^(-15),'TOlX',10^(-15),'MaxFunEvals',695000  ) ;
options.MaxIter     = 50000  ;
%     [x,fval,exitflag,output,qd1] =  fsolve(@fun_curve2 ,var_initial,options)          ;
  [x,fval,exitflag,output,qd1] =  fsolve(@fun_frenet ,var_in,options)          ;
%[x,fval,exitflag,output,qd1] =  fsolve(@fun_p ,.1,options)          ;

 


