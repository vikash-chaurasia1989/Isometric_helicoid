   function f =   initial_guess_loop2()
 
  
  global   N1  N  qd   id   h   rho lm  tau tau1 p1 fac
 
  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold sig branch
  format longE
  %--------Array initialization -------
  path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_2/';
 % path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_discrete/data_7pi/';
 % path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Constant torsion/Full_discrete_formulation/data_5pi/';

   %  
   
% str0 =   '3fold_N72_tau_80995000000.txt' ;
 
%  
%  str0 = '3fold_discrete_N72_tau_81000000000.txt' ;
% %  str0 = '3fold_N72_tau_81000000000.txt' ;
%  
%  %str0 = '5pi_knot_discrete_N120_tau_119850000000.txt';
%  str0 =   '5pi_knot_N120_tau_119842313030.txt';
%  
%  str0 = '7pi_knot_discrete_N_84_tau_153000000000.txt';
% %  
% 
% % branch =1;
 
 
 
%  if(branch==1)
%      tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%      tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%      tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
%      
%      
%      
%      tau1 =      [tau_3 tau_5 tau_7];
%      
%      
%      N = 72;
%      N1 = N-1;
%      h = 1/N;
%      nfold=3;
%      
%  elseif(branch==2)
%      %================ Branch 2 =====================
%      
%      tau1 = [ 11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
%      
%      
%      N = 120;
%      N1 = N-1;
%      h = 1/N;
%      nfold=5;
%      
%     % =============== Branch 3 =====================
%      
%  else
%      tau1 = [15.45:.1:25 25.1:.1:40];
%      
%      N = 84;
%      N1 = N-1;
%      h = 1/N;
%      
%      nfold=7;
%  end
 % 
 % %  
% tau = tau1(p1);
 
if(p1>0)
    tau2 =  tau1(p1);
    
    if(nfold==3)
        str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
      %  str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '_sig_' num2str(sig) '.txt'];

        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    elseif(nfold==5)
        str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
        %
    else
        str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
    end
    
    
end

 % path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_discrete/data_3pi/';
 % path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Constant torsion/Full_discrete_formulation/data_5pi/';

   %  
   
% str0 =   '3fold_N72_tau_80995000000.txt' ;
 
 
 %str0 = '3fold_discrete_N72_tau_81000000000.txt' ;
 

   
strb =  [path 'b_' str0] ;
strlm = [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];
 
 
%==== load b ===========
temp = load(strb);
N   = length(temp)-1;
N1 = N-1;

h = 1/N;
% N = N/2;
bx = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);

 

b(1:3:3*N1-2) = bx(2:N,1)  ;
b(2:3:3*N1-1) = by(2:N,1)  ;
b(3:3:3*N1)    = bz(2:N,1)  ;
%====================
  

 temp =     load(struvw);
% 
u  =    temp(1,1);
v  =    temp(2,1);
w =     temp(3,1);

  %===============================

%=== load rho ==========   
   rho =      load(strrho);
  
 %     rho(1,1) = 2*rho(2,1)-rho(3,1);
%     rho(N-1,1) = rho(1,1);
% %=== load lm ====  
      lm =    load(strlm);
     s = linspace(0,1,N );
%       %  lm(1:N ,1) = .3*sin(10*pi*s);
%       lm(2,1) = 2*lm(3,1)-lm(4,1);
        %  lm(1,1) = 2*lm(2,1)-lm(3,1);
% %       
% %       lm(N-1,1) = 2*lm(N-2,1) - lm(N-3,1);
        % lm(N,1) = 2*lm(N-1,1) - lm(N-2,1);
        p3= 15;
  %  lm(N:-1:N-p3,1) = lm(1:p3+1,1);
%    lm(N-1,1) = lm(2,1);
%    lm(N-3,1) = lm(3,1);
%    lm(N-4,1) = lm(4,1);

%==== for studying effect of aspect ratio on the lagrange multipliers ====
  
  %fac= asinh(tau*pi*sig)/(4*pi^3*tau^3);

  lm = lm/.001*fac;
  rho = rho/.001*fac;
  u = u/.001*fac;
  v =  v/.001*fac;
  w =  w/.001*fac;
 % w =   2.650573192256430e+04




    f =  [b(1:3*N1)  rho' lm' u v    1*w  ]'  ;         % for fun_jacobian6

  
  %--- fix bz(2) = 0; bz(N) = 0 ---
%   
%   f(3,:) = [];
%   f(3*N1-1,:) = [];
%   
  
  %--- boundary points --
  
  bx(1,1) = 1;
  by(1,1) = 0;
  bz(1,1) = 0;
  
  bx(N+1,1) = -1;
  by(N+1,1) =  0;
  bz(N+1,1) =  0;
  
 
%--- for fun_curve9 ---
% f = [b(1:3*N1)  rho' lm' u v w  bx1 by1 bz1]'  ;
 
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