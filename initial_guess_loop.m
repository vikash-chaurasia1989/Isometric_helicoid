   function f =   initial_guess_loop()
 
  
  global   N1  N  qd   id   h   rho lm  tau tau1 p1
 
  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p
  
  %--------Array initialization -------
     path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_discrete/data_3pi/';
 % path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/';

   %  
   
% str0 =   '3fold_N72_tau_80995000000.txt' ;
 
 
 str0 = '3fold_discrete_N72_tau_81000000000.txt' ;
%  str0 = '3fold_N72_tau_81000000000.txt' ;


%  
   if(p1>1)
        tau2 = tau1(p1-1);
  str0 =['3fold_N72_tau_' num2str(10^10*tau2) '.txt'];
 
  end
 %   str0 =  '3fold_N72_tau_210_1.txt';
%   
%  if(p1>32)
% str0 = ['3fold_N72_tau_' num2str(10*tau1(p1-1)) '_1.txt'];
%  end
 str0 = '5pi_knot_N120_tau_330000000000.txt';
 path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_2/';


strb =  [path 'b_' str0] ;
strlm = [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];
strbext = [path 'bext_' str0];

%path2 =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_discrete/data_3pi/';

%strb = [path2 'b_3fold_discrete_N72_tau_81000000000.txt'];%

 
%==== load b ===========
temp = load(strb);
N   = length(temp)-1;
N1 = N-1;

h = 1/N;
% N = N/2;
bx = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);

% mod = sqrt(bx.^2+by.^2+bz.^2);
% 
% bx = bx./mod;
% by = by./mod;
% bz = bz./mod;

b(1:3:3*N1-2) = bx(2:N,1)  ;
b(2:3:3*N1-1) = by(2:N,1)  ;
b(3:3:3*N1)    = bz(2:N,1)  ;
%====================
%     rho(1:N-1,1) = .01*ones(N1,1)*0;
%
     u = 0;
     v = 0;
     w = .1;
   
% 
%  rho = [(rho(1:end-1,1)+rho(2:end,1))'/2 rho(1,1)]';
%   rho = (rho(1:end-1,1)+rho(2:end,1))/2;
%             rho(N-1,1) = rho(1,1);
%              rho(N-2,1) = rho(2,1);
% %   

s = linspace(0,1,N+1);

%    lm(1:N,1) = -1/(h*tau*(N))*ones(N,1)*0;
%    
%    lm(1:N+1,1) = .3*sin(6*pi*s);
%    rho = -lm(1:N-1,1);
% %=== load uvw ========= 
% 

%--- higher order interpolation ---
nfit = 4;

%-- At i = 1;
px = polyfit(linspace(1,nfit+1,nfit+1),bx(1:nfit+1,1)',nfit);
bx1 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),by(1:nfit+1,1)',nfit); 
by1 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),bz(1:nfit+1,1)',nfit); 
bz1 = px(end);

%-- At i = N+1
px = polyfit(linspace(1,nfit+1,nfit+1),bx(N+1:-1:N+1-nfit,1)',nfit);
bxN2 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),by(N+1:-1:N+1-nfit,1)',nfit);
byN2 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),bz(N+1:-1:N+1-nfit,1)',nfit);
bzN2 = px(end);


  temp =     load(struvw);
% 
u  =    temp(1,1);
v  =    temp(2,1);
w =     temp(3,1);
 %===============================

%=== load rho ==========   
   rho =      load(strrho);
  
%   rho(3,1) = 2*(rho(4,1)- rho(5,1));
  %  rho(1,1) = rho(end);
%   rho(2,1) = .5*(rho(1,1)+rho(3,1));
%   rho(1:5,1) = rho(end:-1:end-4,1);
  %  rho(end,1) = 2*rho(end-1,1)-rho(end-2,1);
   
%=== load lm ====  
      lm =   load(strlm);
   % lm(end) = [];
   % lm(1:5,1) = lm(end:-1:end-4,1);
   %  lm(1,1) = 2*lm(2,1)-lm(3,1);
  % lm(3,1) = 2*lm(4,1)-lm(5,1);
   %lm(2,1) = 2*lm(3,1)-lm(4,1);
   %lm(1,1) = (2*lm(2,1) - lm(3,1));
   
    %     lm(N+1,1) = 2*lm(N,1)-lm(N-1,1);
   %  lm(N+1,1) = lm(1,1);%2*lm(N,1)-lm(N-1,1);
      % lm(N,1) = .5*(lm(N+1,1) + lm(N,1));
%    % lm = 1.547869306128353e+00*ones(N,1);
%-- extrapolated points from saved data file 

%strbext = 'bext_3fold_N72_tau_81000000000.txt';
  
  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
 %bext = load(strbext);

bx1 = bext(1);
by1 = bext(2);
bz1 = bext(3);

bxN2 = bext(4);
byN2 = bext(5);
bzN2 = bext(6);
%==============

 
 % f = [b(1:3*N1)  rho' lm' u v w  bx1 by1 bz1 bxN2 byN2 bzN2 ]'  ;  % for fun_jacobian
   f = [b(1:3*N1)  rho' lm' u v 0*w  ]'  ;         % for fun_jacobian6

  
  %--- fix bz(2) = 0; bz(N) = 0 ---
  
  f(3,:) = [];
  f(3*N1-1,:) = [];
  
  
  %--- boundary points --
  
  bx(1,1) = 1;
  by(1,1) = 0;
  bz(1,1) = 0;
  
  bx(N+1,1) = -1;
  by(N+1,1) =  0;
  bz(N+1,1) =  0;
  
  bz(2,1) = 0;
  bz(N,1) = 0;
%  
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
       plot3(bx,by,bz,'LineWidth',4)
       set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
       box on
       grid on
 
  end