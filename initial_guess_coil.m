function x =   initial_guess_coil()


global   N1  N  qd   id   h   rho lm  tau tau1 p1 fac tau2

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold sig branch
format longE

% 
% % %
%  path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch4/';
% % % 
%  strb =   'b_coil_N108_tau_116788578830.5211.txt';% 'b_branch_2_N240_symmetric.txt';%  
% % 
%  N = 108;
% temp = load([path strb]);
% 
% bx   = temp(:,1);
% by   = temp(:,2);
% bz   = temp(:,3);
% 
% bx(1,1) = 1;
% by(1,1) = 0;
% bz(1,1) = 0;
% 
% bx(N+1,1) = -1;
% by(N+1,1) = 0;
% bz(N+1,1) = 0;
N= 100 ;
h = 1/N;
t1 = linspace(0,2*pi,N+1);
t1 = t1';

m = 2;
n = -3;
al = pi/3;
bt = 1/2*acos(-5/7);

S_al = [cos(al) -sin(al) 0; sin(al) cos(al) 0; 0 0 1];
S_bt = [cos(bt) -sin(bt) 0; sin(bt) cos(bt) 0; 0 0 1];


u = 1/3^.5;
v =1/6^.5;
w = v/u;

C = [u -v  w ; 
     u -v -w ;
     u 2*v 0];




U = [ 0 -u  u ;
      u  0 -u ;
     -u  u  0 ];

id = eye(3,3);
for i = 1:N+1
    
    t = t1(i);
    
    Rnt = [1  0 0 ; 0 cos(n*t) -sin(n*t);0 sin(n*t) cos(n*t)];
        
    Qmt = id + sin(m*t)*U + (1-cos(m*t))*U*U;
    
    
    
    temp = Qmt*C*S_al*Rnt*S_bt;
    
    bx(i,1) = temp(1,1);
    by(i,1) = temp(2,1);
    bz(i,1) = temp(3,1);
    
    
    
end



lm = zeros(N,1);
rho= zeros(N-1,1);



 s1 = (linspace(0,1,N))';
 
 lm = .1*sin(9*pi*s1);
 
 rho = lm(1:N-1,1);
 
u = 0;
v = 0;
w = 0;

% 
% % % %=== lagrange multipliers from saved data --
% % % 
% % % % % % % % %
% N = 105;
% N1=N-1;
% tau2 = 15.425;
% branch = 3;
% 
% path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
% str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];
% 
% x = load([path str0]);
% 
% rho1           = x(3*N1+1:4*N1,1)   ;
% lm1            = x(4*N1+1:5*N1+1,1) ;
% 
% u             = x(5*N1+2 ,1)       ;
% v             = x(5*N1+3,1)        ;
% w             = x(5*N1+4,1)        ;
% 
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
%  %mod = sqrt(bx.^2+by.^2+bz.^2);
% 
% % bx = bx./mod;
% % by = by./mod;
% % bz = bz./mod;
% 
% N = 336 ;
% N1=N-1;
% h = 1/N;
% 
% 
% 
% b(1:3:3*N1-2) = bx(2:N,1)  ;
% b(2:3:3*N1-1) = by(2:N,1)  ;
% b(3:3:3*N1)   = bz(2:N,1)  ;
% 
% 
% %====================
% 
% %f = b';
% 
%  %lag = fun_lagrangian6(b');
% 
%   x =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;         % for fun_jacobian6

 
 
%== for saving data at intermediate steps 
% b(1:3:3*N1-2) = bx(2:N,1)  ;
% b(2:3:3*N1-1) = by(2:N,1)  ;
% b(3:3:3*N1)   = bz(2:N,1)  ;
% 
% x  =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;
%  path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%        % str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
%         
%         str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_symmetry.txt'];
%         
%         fileID = fopen([path str0],'w');
%         fprintf(fileID,'%30.16E   \r\n',x );
%         fclose(fileID);
%  
%  
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

N1 = N-1;

b(1:3:3*N1-2) = bx(2:N,1)  ;
b(2:3:3*N1-1) = by(2:N,1)  ;
b(3:3:3*N1)   = bz(2:N,1)  ;
% %====================
%
% %f = b';
%
%  % lag = fun_lagrangian6(b');
%
     x =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;         % for fun_jacobian6
%
% %  f =   [b(1:3*N1)  lag'  ]'  ;




 
%
% % % f(3*N1+1,1) = 0;
% % %  f(4*N1 ,1) = 0;
% % %f(5*N1+2,1) = 0;
 %x(5*N1+3,1) = 0;
 %    x(5*N1+4,1) = 0;
%
%--- boundary points --

% %==== resizing the saved data ===
% % 
% %  
%   bx1(2:N,1)     = x(1:3:3*N1-2,1)    ;
%   by1(2:N,1)     = x(2:3:3*N1-1,1)    ;
%   bz1(2:N,1)     = x(3:3:3*N1,1)      ;
%   
% bx1(1,1) = 1;
% by1(1,1) = 0;
% bz1(1,1) = 0;
% 
% bx1(N+1,1) = -1;
% by1(N+1,1) = 0;
% bz1(N+1,1) = 0;

%   
%   
%   
%   rho1          = x(3*N1+1:4*N1,1)   ;
%   lm1           = x(4*N1+1:5*N1+1,1) ;
%   
%   
%   
%   u             = x(5*N1+2 ,1)       ;
%   v             = x(5*N1+3,1)        ;
%   w             = x(5*N1+4,1)        ;
%   
%   
%   
%   N2 = 120;        % new number of points 1
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
% b(1:3:3*N1-2) = bx(2:N,1)  ;
% b(2:3:3*N1-1) = by(2:N,1)  ;
% b(3:3:3*N1)   = bz(2:N,1)  ;
% %====================
% 
% %f = b';
% 
%  % lag = fun_lagrangian6(b');
% 
%   x =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;         % for fun_jacobian6
%   
%   
%   
%   
  
  
  
  
% 
% bx(1,1) = 1;
% by(1,1) = 0;
% bz(1,1) = 0;
% 
% bx(N+1,1) = -1;
% by(N+1,1) =  0;
% bz(N+1,1) =  0;
% 

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