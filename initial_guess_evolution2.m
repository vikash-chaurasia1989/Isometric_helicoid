function f =   initial_guess_evolution2()


global   N1  N  qd   id   h   rho lm  tau tau1 p1 fac tau2

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold sig branch
format longE


%
% if(p1>0)
%
%     [val,ind] = min(abs(tau2-tau1(p1)));
%
%      tau3  = tau2(ind(1));
%
%
%     if(branch==1)
%         N = 72;
%         str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau3) '.txt'];
%         %  str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '_sig_' num2str(sig) '.txt'];
%
%         path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
%     elseif(branch==2)
%         N = 120;
%         str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau3) '.txt'];
%         path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
%         %
%     else
%         N = 84;
%         str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau3) '.txt'];
%         path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
%     end
%
%
% end
%
% path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
%
% str0 =   '5pi_knot_N180_tau_342000000000.txt';
%
% path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Constant torsion/Full_discrete_formulation/data_9pi/';
%
% str0 = '9pi_knot_N_432_tau_80.9361_z0.txt';
%
%
% str0 = '9pi_knot_N_432_tau_184.7472_z0.txt';
%
% path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch1/';

% strb =  'branch_1_N105_tau_80940900000.txt'% 'b_branch_2_N240_symmetric.txt'%  'b_branch_3_N336_symmetric.txt';

% str0 =   '9pi_knot_N162_tau_185050200000.txt';
%
%
% strb =  [path 'b_' str0] ;
% strlm = [path 'lm_' str0];
% strrho = [path 'rho_' str0];
% struvw = [path 'uvw_' str0];
%
%
% %==== load b ===========
% temp = load(strb);
% N   = length(temp)-1;
% N1 = N-1;
%
% h = 1/N;
% % N = N/2;
% bx1 = temp(1:N+1,1);
% by1 = temp(1:N+1,2);
% bz1 = temp(1:N+1,3);
%
%
%
% temp =     load(struvw);
% %
% u  =    temp(1,1);
% v  =    temp(2,1);
% w =     temp(3,1);
%
% %===============================
%
% %=== load rho ==========
% rho1 =      load(strrho);
%
%
% % %=== load lm ====
% lm1 =    load(strlm);
%
%
% %================  Interpolation of the initial guess ==============
%
% %======== b =========
% % %
% N2 = 153;        % new number of points 1
%
%
%  s1 = (linspace(0,1,N+1))';
%  s2 = (linspace(0,1,N2+1))';
%
%  bx = interp1q(s1,bx1,s2)  ;
%  by = interp1q(s1,by1,s2)  ;
%  bz = interp1q(s1,bz1,s2)  ;
%
%  s1 = (linspace(0,1,N))';
%  s2 = (linspace(0,1,N2))';
%
%
% lm  = (interp1q(s1,lm1,s2));
%
%  s1 = (linspace(0,1,N-1))';
%  s2 = (linspace(0,1,N2-1))';
%
%
%  rho = (interp1q(s1,rho1,s2));
%
%
%
% N = N2;
% N1 = N-1;
%
% h = 1/N;
%
% %
% %   s = linspace(0,1,N);
% % %
% % lm = 0.01*sin(18*pi*s);
% % lm = lm';
% % rho = lm(1:N-1,1);
% %
% % rho(1,1) = 0;
% % rho(end,1) = 0;
%
% mod = sqrt(bx.^2+by.^2+bz.^2);
%
% % bx = bx./mod;
% % by = by./mod;
% % bz = bz./mod;
%
% b(1:3:3*N1-2) = bx(2:N,1)  ;
% b(2:3:3*N1-1) = by(2:N,1)  ;
% b(3:3:3*N1)   = bz(2:N,1)  ;
% %====================
%
% %f = b';
%
%  % lag = fun_lagrangian6(b');
%
%    f =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;         % for fun_jacobian6
%
% %  f =   [b(1:3*N1)  lag'  ]'  ;

% dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2+(bz(1:end-1,1)-bz(2:end,1)).^2);
% tau = sum(dl);

% % % % % %
N = 180;
N1 = N-1;
h = 1/N;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

 %path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

  
%str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau1(p1))) '.txt'];

str0 =     'branch_1_N180_tau_80939466335.txt';
;%'branch_1_N120_tau_80940990000.txt';% 'branch_1_N120_tau_80940900000.txt';% 'branch_1_N105_tau_80940900000.txt'; % 'branch_1_N240_tau_80940900000.txt';%

x = load([path str0]);

 f = x;

% 
% %
% %==========================================================================
% %==========================================================================
%             %   Interpolation of the saved data 
% %==========================================================================
% %==========================================================================

bx1(1,1) = 1;
by1(1,1) = 0;
bz1(1,1) = 0;

bx1(N+1,1) = -1;
by1(N+1,1) =  0;
bz1(N+1,1) =  0;


bx1(2:N,1)     = x(1:3:3*N1-2,1)    ;
by1(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz1(2:N,1)     = x(3:3:3*N1,1)      ;


rho1           = x(3*N1+1:4*N1,1)   ;
lm1            = x(4*N1+1:5*N1+1,1) ;

u              = x(5*N1+2 ,1)       ;
v              = x(5*N1+3,1)        ;
w              = x(5*N1+4,1)        ;


%==== interpolation of the saved data ==
N2 = 180;        % new number of points 1
N1 = N2-1;

 s1 = (linspace(0,1,N+1))';
 s2 = (linspace(0,1,N2+1))';

 bx = interp1q(s1,bx1,s2)  ;
 by = interp1q(s1,by1,s2)  ;
 bz = interp1q(s1,bz1,s2)  ;

 s1 = (linspace(0,1,N))';
 s2 = (linspace(0,1,N2))';


lm  = (interp1q(s1,lm1,s2));

 s1 = (linspace(0,1,N-1))';
 s2 = (linspace(0,1,N2-1))';


 rho = (interp1q(s1,rho1,s2));
 
 
 mod = sqrt(bx.^2+by.^2+bz.^2);

% bx = bx./mod;
% by = by./mod;
% bz = bz./mod;

N = N2;
N1 = N-1;
h = 1/N;

b(1:3:3*N1-2) = bx(2:N,1)  ;
b(2:3:3*N1-1) = by(2:N,1)  ;
b(3:3:3*N1)   = bz(2:N,1)  ;

 f =  [b(1:3*N1)  rho' lm'  u  v   w  ]'  ;         % for fun_jacobian6







%f(1:5*N1+4,1) = temp(1:5*N1+4,1);















% f(3*N1+1,1) = 0;
%  f(4*N1 ,1) = 0;
%  f(5*N1+4,1) = 0;
%
%--- boundary points --

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


%------------ Array initialization ----------------

bxp = zeros(N+1,1);
byp = zeros(N+1,1);
bzp = zeros(N+1,1);

bx2p = zeros(N+1,1);
by2p = zeros(N+1,1);
bz2p = zeros(N+1,1);

bx4p = zeros(N-1,1);
by4p = zeros(N-1,1);
bz4p = zeros(N-1,1);



%-- hessian
id = [1 0 0;0 1 0;0 0 1];


figure(1)

plotbrowser on
title('Initial guess')
% plot3(bx,by,bz,'LineWidth',4)
plot(lm,'-o')
set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
box on
grid on

end