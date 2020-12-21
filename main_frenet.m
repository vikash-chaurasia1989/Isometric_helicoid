clear all
clc

global    N1  tau   h tx ty tz  s nx ny nz k0 al w r

global bx by bz   N     rx ry rz   bext or kappa A B tau1 nfold
set(0,'defaultlinelinewidth',2)


or =    -1;

branch =1;
N = 120;
N1 = N-1;
h = 1/N;
tau =10.1% 8.0940900000;

nfold = 3;

strpath = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
% strpath = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch';
path = [strpath num2str(branch) '/'];

% strtau = ['tau_branch_' num2str(branch) '.txt'];
%
% tau1 = load([path strtau]);
%
% [val,ind] = min(abs(tau1-tau));
%
% tau2 = tau1(ind);
%
% tau = tau2;

% str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];

  str0 =  'branch_1_N120_tau_80940990000.txt'; 

%str0  = 'branch_2_N120_tau_119787000480_symmetry.txt';

%  path = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
%  str0 = 'b_5pi_knot_N120_tau_119842333032.1598.txt';
x = load([path str0]);

%----- reading input and constructing bx1 by1 bz1 ----

bx1(1,1) = 1;
by1(1,1) = 0;
bz1(1,1) = 0;

bx1(N+1,1) = -1;
by1(N+1,1) =  0;
bz1(N+1,1) =  0;



bx1(2:N,1)     = x(1:3:3*N1-2,1)    ;
by1(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz1(2:N,1)     = x(3:3:3*N1,1)      ;

% bx1 = x(:,1);
% by1 = x(:,2);
% bz1 = x(:,3);

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

bx12p(ig,1) = (bx1(ig+1,1) + bx1(ig-1,1) - 2*bx1(ig,1))/h^2 ;
by12p(ig,1) = (by1(ig+1,1) + by1(ig-1,1) - 2*by1(ig,1))/h^2 ;
bz12p(ig,1) = (bz1(ig+1,1) + bz1(ig-1,1) - 2*bz1(ig,1))/h^2 ;

bx12p(N+1,1) = (bext(4,1) + bx1(N,1) - 2*bx1(N+1,1))/h^2 ;
by12p(N+1,1) = (bext(5,1) + by1(N,1) - 2*by1(N+1,1))/h^2 ;
bz12p(N+1,1) = (bext(6,1) + bz1(N,1) - 2*bz1(N+1,1))/h^2 ;
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





k0 = 1;

s1 = linspace(0,1,N+1)';
N2 = 120;
s2 = linspace(0,1,N2+1)';

kappa3 = interp1q(s1,kappa2,s2);

kappa3 = kappa3-min(kappa3);

for m = 1:2:nfold-1
    
    kappa3((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold))) = -kappa3((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold)));
    
end
m = nfold;
kappa3((2*m-1)*N/(2*nfold)+1:N+1) = -kappa3((2*m-1)*N/(2*nfold)+1:N+1);

N = N2;

s = linspace(0,1 ,N+1);

tau1 = tau;

var_in = [ 4.037960852507271e+01,  1.679359081062158e+04] ;   % for 3pi solution

var_in = [ 4.037960852507271e+02,  1.679359081062158e+05] ;   % for 3pi solution

%var_in = [2.768311267779191e+01,8.564919953409128e+02];

tau  = 15.1;
N = 120*2;
h = 1/N;
%--------------------------- Solver ---------------------------------------
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

nfold = 10;
algo = 1;

stralgo = {'levenberg-marquardt' ,'trust-region-reflective' ,'trust-region-dogleg'};

options             = optimset('Display','iter', 'Algorithm',stralgo{algo},'Jacobian','off', 'TOlFun',10^(-32),'TOlX',10^(-32),'MaxFunEvals',695000  ) ;
options.MaxIter     = 50000  ;
%     [x,fval,exitflag,output,qd1] =  fsolve(@fun_curve2 ,var_initial,options)          ;
[x,fval,exitflag,output,qd1] =  fsolve(@fun_frenet5  ,var_in,options)          ;
%[x,fval,exitflag,output,qd1] =  fsolve(@fun_p ,.1,options)          ;



% %=== kappa = +- ==== .. Correcting sign of kappa in suitable intervals
%
%  kappa(N/6+1:(N/6+N/3)) = - kappa(N/6+1:(N/6+N/3));
%  kappa(N-N/6+1:N+1) = -kappa(N-N/6+1:N+1);

s1 = linspace(0,1,length(kappa));
s2 = linspace(0,1,length(kappa3));


plot(s1,kappa)
hold on
plot(s2,kappa3)



%=== for 5pi knot ====



kappa2 = [kappa2(end  ,1) kappa2(1:end-1,1)']';

N = 120;

kappa2 = [kappa2(12:end,1)' kappa2(1:11,1)']';
kappa2 = kappa2 - min(kappa2);

for m = 1:2:nfold-1
    
    kappa2((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold))) = -kappa2((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold)));
    
end
m = nfold;
kappa2((2*m-1)*N/(2*nfold)+1:N+1) = -kappa2((2*m-1)*N/(2*nfold)+1:N+1);



%==== constructing band from frenet frame data =====

N = length(kappa)-1;
h = 1/N;

% bz1 = x(:,3);

bext(1:3,1) = -1*[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -1*[bx(2,1);by(2,1);bz(2,1)];

%---- derivative calculation for b2 --- bN  ---

ig = 2:N;

%=== O(h)====
bxp(1,1) = (bx(1,1)-bext(1,1))/(h);
byp(1,1) = (by(1,1)-bext(2,1))/(h);
bzp(1,1) = (bz(1,1)-bext(3,1))/(h);

bxp(ig,1) = (bx(ig,1)-bx(ig-1,1))/(h);
byp(ig,1) = (by(ig,1)-by(ig-1,1))/(h);
bzp(ig,1) = (bz(ig,1)-bz(ig-1,1))/(h);


bxp(N+1,1) = (bx(N+1,1)-bx(N,1))/(h);
byp(N+1,1) = (by(N+1,1)-by(N,1))/(h);
bzp(N+1,1) = (bz(N+1,1)-bz(N,1))/(h);


% 
% 
% i = 1:N+1   ;
% 
% %
% tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
% ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
% tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
% tx = tx/tau;
% ty = ty/tau;
% tz = tz/tau;

%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
% rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;
% 
% 
% 
% for i = 1:N
%     
%     rx(i+1,1) = rx(i,1) + h*tx(i,1);
%     ry(i+1,1) = ry(i,1) + h*ty(i,1);
%     rz(i+1,1) = rz(i,1) + h*tz(i,1);
%     
% end




s1 = linspace(0,1,length(kappa));
s2 = linspace(0,1,length(kappa2));


plot(s1,kappa);
hold on
plot(s2,kappa2);
hold on

grid on
set(gca,'Fontsize',30)
set(gca,'XMinorTick','on','YMinorTick','on')
xlhand = get(gca,'xlabel');
set(xlhand,'string','$s$ ','Interpreter','LaTex','fontsize',30)
ylhand = get(gca,'ylabel');
set(ylhand,'string','$\kappa$','Interpreter','LaTex','fontsize',30)

title(' ')



% %============ Perturbation about elliptic function ===========
% 
% al1 = 2*A + 2*sqrt(A^2+B);
% al2 = 0;
% al3 = -2*A +2*sqrt(A^2+B);
%  
% p = sqrt((al3-al2)/(al3+al1));
% q = sqrt(1-al2/al3);
% r = 1/2*sqrt(al3+al1);
% 
% 
% [K,E] = ellipke(p);
%  s = linspace(0,2*nfold*K/r,N+1);
% %s = linspace(0,1,N+1);
%   tau = tau1*2*nfold*K/r;
%  
%  h = 1/N*2*nfold*K/r;
% 
% [sn,cn,dn] = ellipj(r*s,p,eps);
%  %kappa =  sqrt(al3*(1-q^2*sn.^2));
% 
%  dn = sqrt(1-p^2*sn.^2);
%  
% %== satisfying the equilibrium equation =====
% kp = -al3*q^2*r*sn.*cn.*dn./kappa'   ;
% err = kp.^2 +kappa'.^4/4 + A*kappa'.^2 -B; 
%  
% % 
% u = al3*(1-q^2*sn.^2);
%  up = -2*q^2*al3*r*sn.*cn.*dn;%4*q^4*al3^2*r^2*sn.^2.*cn.^2.*dn.^2;
% 
% err = up.^2+ (u+al1).*(u-al2).*(u-al3);
% 
% % 
% %1.380983773188140e+01
