       function f =   initial_guess_refine()
 
  
  global   N1  N  qd   id   h   rho lm  tau
 
  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p
  
  %--------Array initialization -------
 %--- This file only uses saved data for initial guess 
   %  
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_2/';
    % path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/';

   %  
   
 str0 =  '3fold_N36_tau_81000000000.txt' ;
 
% ==========   b ==================  
%tau1 = 8.1;
   % str0 = ['3fold_N72_tau_' num2str(10*tau1) '.txt'];
    
    %-- 5pi --
  %  str0 = ['5pi_knot_N120_tau_' num2str(10*tau1(p1)) '.txt'];
  
  N = 36;
 % tau = 21.5;
 %str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
 
% str0 = '7pi_knot_N49_tau_215000000000.txt';
 
 strb = [path 'b_' str0];
 %strb = [path 'b_' '7pi_knot_N49_tau_228000000000_1.txt']
%strr = ['r_' str0];

strlm =  [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];
strbext = [path 'bext_' str0];



       temp = load(strb);
 
           bx = temp(:,1);
           by = temp(:,2);
           bz = temp(:,3);
 
            N = length(bx)-1;
  
%=========== rho =================

rho = load(strrho);
 %  v = w;       
        temp = load(struvw);
       
       u  = temp(1,1);
       v  = temp(2,1);
       w = temp(3,1);

      
%=========== \lambda ==============

lm = load(strlm);   
    
 
%============ Refinement of the uploaded data =====

%   %    %--- refinement --  

%======== b =========
% %
 s1 = (linspace(0,1,N+1))';
 s2 = (linspace(0,1,2*N+1))';
 
 bx = interp1q(s1,bx,s2)  ;
 by = interp1q(s1,by,s2)  ;
 bz = interp1q(s1,bz,s2)  ;

  s1 = (linspace(0,1,N))';
 s2 = (linspace(0,1,2*N))';
 
 
lm  = (interp1q(s1,lm,s2));
 
 s1 = (linspace(0,1,N-1))';
 s2 = (linspace(0,1,2*N-1))';
 
 
 rho = (interp1q(s1,rho,s2));
 
  

 N = length(bx)-1;

 N1 = N-1;
%  
%       
           

%---- coarsing the data ---
% 
%        N2 = N/2;
%        
%        bx1 = bx(1:2:N+1);
%        by1 = by(1:2:N+1);
%        bz1 = bz(1:2:N+1);
%        
%        bx = bx1;
%        by = by1;
%        bz = bz1;
%        
%        N = length(bx)-1;
%        N1 = N-1;
%        
       %--- plotting initial guess 
       figure(1)
       
       plotbrowser on 
       title('Initial guess')
       plot3(bx,by,bz,'LineWidth',4)
       set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
       box on
       grid on
  
   b(1:3:3*N1-2) = bx(2:N,1)  ;
   b(2:3:3*N1-1) = by(2:N,1)  ;
   b(3:3:3*N1)    = bz(2:N,1)  ;
   
 
   h = 1/N ;
 
 
%--- linear interpolation
bx1 = 2*bx(1,1)-bx(2,1);
by1 = 2*by(1,1)-by(2,1);
bz1 = 2*bz(1,1)-bz(2,1);

bxN2 = 2*bx(N+1,1)-bx(N,1);
byN2 = 2*by(N+1,1)-by(N,1);
bzN2 = 2*bz(N+1,1)-bz(N,1);


%--- higher order interpolation ---
nfit = 4;

%-- At i = 1;
px = polyfit(linspace(1,nfit+1,nfit+1),bx(1:nfit+1,1)',nfit);
bx1 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),by(1:nfit+1,1)',nfit);;
by1 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),bz(1:nfit+1,1)',nfit);;
bz1 = px(end);

%-- At i = N+1
px = polyfit(linspace(1,nfit+1,nfit+1),bx(N+1:-1:N+1-nfit,1)',nfit);
bxN2 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),by(N+1:-1:N+1-nfit,1)',nfit);
byN2 = px(end);

px = polyfit(linspace(1,nfit+1,nfit+1),bz(N+1:-1:N+1-nfit,1)',nfit);
bzN2 = px(end);


%-- extrapolated points from saved data file 
bext = load(strbext);

bx1 = bext(1);
by1 = bext(2);
bz1 = bext(3);

bxN2 = bext(4);
byN2 = bext(5);
bzN2 = bext(6);

%--- for fun_curve7 ---
  f = [b(1:3*N1)  rho' lm' u v w  bx1 by1 bz1 bxN2 byN2 bzN2 ]'  ;
  
  
  %--- fix bz(2) = 0; bz(N) = 0 ---
  
  f(3,:) = [];
  f(3*N1-1,:) = [];
  
  
  %--- boundary points --
  
  bx(1,1) = 1;
  by(1,1) = 0;
  bz(1,1) = 0;
  
  bx(N+1,1) = -1;
  by(N+1,1) =   0;
  bz(N+1,1) =   0;
  
  bz(2,1) = 0;
  bz(N,1) = 0;
 
%--- for fun_curve81 ---
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
 
   

%   %-- hessian
    qd  = zeros(3*N1+3,3*N1+3);
    id = [1 0 0;0 1 0;0 0 1];
 
  end