function f =   initial_guess_evolution2()


global   N1  N  qd   id   h   rho lm  tau tau1 p1 fac tau2 id

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold sig branch algo
format longE

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
% % %====================
% %
% % %f = b';
% %
% %  % lag = fun_lagrangian6(b');
% %
%     x =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;         % for fun_jacobian6
%
% %  f =   [b(1:3*N1)  lag'  ]'  ;

% dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2+(bz(1:end-1,1)-bz(2:end,1)).^2);
% tau = sum(dl);

% % % % % %
if (branch==1||branch==2||branch==3)
    N=105;
end

if(branch==4||branch==9)
    N = 180%120;
end
if(branch==5)
    N = 140;
end

if(branch== 6)
    N = 90;
end

if(branch== 8)
    N = 150;
end

 %N = 405 ;
N= 600;
N1 = N-1;
h = 1/N;

%--- boundary points --

id = eye(3,3);


path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
 % path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
 
  
% str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau1(p1-1))) '.txt'];

 %path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent6/';
%  str0 = 'branch_3_N105_tau_160000000000_step_1500.txt';
%  str0 =  'branch_1_N105_tau_160000000000_unstable.txt';
%   str0 = 'b_5pi_unknot_discrete_N_105_tau_150177413261.5587.txt';%'b_5pi_unknot_N_120_tau_150.1774_z0.txt';
%  %str0 =  'branch_unknot_5pi_N120_tau_151250000000.txt' ;
% %str0 =     'branch_1_N180_tau_80939466335.txt';
   %'branch_1_N120_tau_80940990000.txt';% 'branch_1_N120_tau_80940900000.txt';% 'branch_1_N105_tau_80940900000.txt'; % 'branch_1_N240_tau_80940900000.txt';%
% 
%  str0 =  'branch_unknot_7pi_N105_tau_200000000000.txt';%   'branch_unknot_5pi_N105_tau_162500000000.txt'
%  str0 = 'branch_unknot_5pi_N120_tau_210000000000.txt';
%  'branch_unknot_5pi_N105_tau_214805796283.txt';
% %  str0 = ['branch_unknot_5pi_N' num2str(N) '_tau_' num2str(round(10^10*tau1(p1))) '.txt'];
 

% % 
% if(algo==3)
%        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau1(p1-1))) '.txt'];
%        x = load([path str0]);
%          x(3*N1+1:end,1)  = 0;%rho; 
%        
% else
%        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau1(p1))) '.txt'];
%         x = load([path str0]);
%       
% end
%        f=x;
   tau1 = 78.2  ;
  str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau1)) '.txt'];
  str0 =   'branch_15_N600_tau_783585500000_cos.txt';
% %tau = 78.25;
% %    str0 =   'branch_4_N150_tau_150000000000.txt';
% %  
% %  %  'branch_1_N240_tau_80936000000.txt';
% % % 
  %  str0 =  'branch_15_N400_tau_782000000000_knot2.txt';
     x = load([path str0]);
%      % x(3*N1+1:5*N1+4,1)  = 0;%rho;
%      y = x(1:3*N1,1);
%    
       f = x;
% % %      
     
bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


%  
%     N = 1000;
%     N1 = N-1;
%      strb = 'b_25_N1000_tau_783300000000_cos.txt';% 'b_25_knot_N300_tau_782000000000_knot2.txt';%'branch_15_N250_tau_782000000000_unknot1.txt';% 'b_25_knot_N1000_tau_782000000000_knot1.txt';
% %     
%     temp = load([path strb]);
%     
%     bx = temp(:,1);
%     by = temp(:,2);
%     bz = temp(:,3);
%     
%     mod = sqrt(bx.^2+by.^2+bz.^2);
%     bx = bx./mod;
%     by = by./mod;
%     bz = bz./mod;
%     
 
% 
% %=== generating initial guess using spherical angle 
% 
% 
% N = 500;
% N1  = N-1;
   nfold = 25;
   s = linspace(0,1,N+1);
   s = s';
% % % 
% ph1 = -sin(2*nfold*pi*s);
% al = -atan((nfold-2)*pi/N);
% 
%  ph  = -(1:N+1)'*sin(al) + ph1*cos(al);
% 
% %ph  = (1:N+1)'*13*pi/(N+1) + ph1*cos(al);
% 
% 
% 
% %ph = 3*pi*s+ sin(9*pi*s);
% th = pi/2- sin(nfold*pi*s);
% 
% %th = pi/2 + .2*sin(5*s*pi);
% 
% 
% bx = sin(th).*cos(ph);
% by = sin(th).*sin(ph);
% bz = cos(th)         ;
% % 
% x = zeros(5*N1+4,1);
% x(1:3:3*N1-2,1) = bx(2:N,1);
% x(2:3:3*N1-1,1) = by(2:N,1);
% x(3:3:3*N1,1)   = bz(2:N,1);
% % % % 
% % % % %  x(1:3:3*N1-2) = x1(2:end-1,1)  ;
% % % % %  x(2:3:3*N1-1) = x1(2:end-1,2)  ;
% % % % %  x(3:3:3*N1)   = x1(2:end-1,3)  ;
% % % % % 
% % % % %  x = x';
% % % % % % 
% % % % % %  %  f = x;
% % % % % % %   %== calculating lagrange multiplier 
% % % % %   
%   [fval0,jac] = fun_jacobian6([x(1:3*N1,1)' zeros(1,2*N1+4)]');
% % % 
%   gradG = jac(3*N1+1:5*N1+4,1:3*N1);
% % %   
% % % %== lagrange multipliers corresponding to bx0, by0, and bz0 === 
% % % 
%   lag = pinv(gradG')*fval0(1:3*N1,1); %(gradG*gradG')\(gradG*fval0(1:3*N1,1));  % lagrange multiplier corresponding to x0
% % % 
%   f = [x(1:3*N1,1)'  lag']';
% % % 
% % % % 
% % % 
%   x(3*N1+1:4*N1,1) =  sin(2*nfold*pi*s(1:end-2,1))       ;
%   x(4*N1+1:5*N1+1,1) = 10*sin(2*nfold*pi*s(1:end-1,1))   ;
% % % %    
 %   f = [x(1:3*N1,1)'  zeros(1,2*N1+4)]';
% % % 
%     y = x(1:3*N1,1);
%  %  
    % f = x;
%    


bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

   %  
% % 
% % %
% %==========================================================================
% %==========================================================================
%             %   Interpolation of the saved data 
% %==========================================================================
% % %==========================================================================
% % 
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
N2 = 600;        % new number of points 1
 

 s1 = (linspace(0,1,N+1))';
 s2 = (linspace(0,1,N2+1))';

 
 bx = zeros(N2+1,1);
 by = zeros(N2+1,1);
 bz = zeros(N2+1,1);
 
 bx = interp1q(s1,bx1,s2)  ;
 by = interp1q(s1,by1,s2)  ;
 bz = interp1q(s1,bz1,s2)  ;
 
%  bx = cos(25*pi*s2);
%  by = sin(25*pi*s2);
%  bz = sin(25*pi*s2);

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
% 
         %   f = [b(1:3*N1)   zeros(1,2*N1+4)]';




% % 
% % %f(1:5*N1+4,1) = temp(1:5*N1+4,1);
% % 
% % 
% % 
% % 
% 
% 
% 
% 
% 
% 
% 
 
% f(3*N1+1,1) = 0;
%  f(4*N1 ,1) = 0;
%  f(5*N1+4,1) = 0;
%

%------------ Array initialization ----------------


bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

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


% figure(1)
% 
% plotbrowser on
% title('Initial guess')
% % plot3(bx,by,bz,'LineWidth',4)
% plot(lm,'-o')
% set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
% box on
% grid on

end