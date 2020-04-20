%--- In this file, we read saved data and compute bending energy for each
%given \tau 

%==========================================================================
%==========================================================================

%====== In this file, we calculate bending energy from the saved
%equilibrium solutions and check their stability as well. Stability check
%is the additional feature in this program if compared to bending_energy

clear all
clc

   
global  N  N1  u v w    f_b     id   h tau rho lm   
global bx by bz   bx2p by2p bz2p   bext

 addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_2')
   %==== data of number of twist and tau
%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
  col6 = [1, 0.929, 0.278];
  col = {col1,col2,col3,col4,col5,col6} ;
%--- identity matrix 

id = [1 0 0;0 1 0;0 0 1];
%================  3pi twist ==========
 


%===== tau for 3 fold ====
 
 
tau_3 = [ 8.0957344694550208E+00  8.09574 8.095745 8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.05:8.5 8.6:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
    13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
    13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];

tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];

    
tau1 =      [tau_3 tau_5 tau_7];
N = 72;
h = 1/N;

N1 = N-1;
 
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';

for i = 1:110%length(tau1)
    
     %==== Number of twist =====
     
     if(tau1(i)<=tau_3(end))
         num_t(i) = 3;
     elseif(tau_5(1)<=tau1(i) && tau1(i)<=tau_5(end))
         num_t(i) = 5;
     else
         num_t(i) = 7;
         
     end
  
     
     
    tau = tau1(i);
    
    
     str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];

    
    strb = [path 'b_' str0];
    strr = [path  'r_' str0];
    
    strlm = [path  'lm_' str0];
    strrho = [path 'rho_' str0];
    struvw = [path 'uvw_' str0];
    strbext = [path 'bext_' str0];
    

    
     %--- load b
    temp = load(strb);

  bx = temp(1:N+1,1);
  by = temp(1:N+1,2);
  bz = temp(1:N+1,3);
 
   b(1:3:3*N1-2) = bx(2:N,1)  ;
   b(2:3:3*N1-1) = by(2:N,1)  ;
   b(3:3:3*N1)    = bz(2:N,1)  ;
   
   
  %--- load u, v, and w -- 
   temp = load(struvw);
       
       u  = temp(1,1);
       v  = temp(2,1);
       w =  temp(3,1);

%=====================================================

 %---- load \lambda and \rho ----
 
   rho = load(strrho);
   lm  = load(strlm);

%-- extrapolated points from saved data file 
%bext = load(strbext);

  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
bx1 = bext(1);
by1 = bext(2);
bz1 = bext(3);

bxN2 = bext(4);
byN2 = bext(5);
bzN2 = bext(6);


%===== Using extrapolation instead of equlibrium solutions at the ends 


%==========================================================================
%------- Array for input of the function fun_jacobian ---------------------
%==========================================================================

  f = [b(1:3*N1)  rho' lm' u v w  bx1 by1 bz1 bxN2 byN2 bzN2 ]'  ;
  
  
  %--- fix bz(2) = 0; bz(N) = 0 ---
  
  f(3,:)      = [];
  f(3*N1-1,:) = [];
  
  %--- fixing the end points ----
  
  bx(1,1) = 1;
  by(1,1) = 0;
  bz(1,1) = 0;
  
  bx(N+1,1) = -1;
  by(N+1,1) =  0;
  bz(N+1,1) =  0;
  
  bz(2,1)   =  0;
  bz(N,1)   =  0;
  
  
 
    
%==========================================================================
%====================  Stability Check ====================================

bx2p = zeros(N+1,1);
by2p = zeros(N+1,1);
bz2p = zeros(N+1,1);

 
 st_3pi(i) =1;% fun_stability(f);

%==========================================================================


%==========================================================================
%=======================  Bending energy calculation ======================

%--- 

%--- higher order interpolation ---
nfit = 5;

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
%bext = load(strbext);
% 
 bext(1,1) = bx1;
 bext(2,1)    = by1;
 bext(3,1)    = bz1;

bext(4,1)     = bxN2;
bext(5,1)     = byN2; 
bext(6,1)     = bzN2; 
%--- 
%==========================================================================
  ig = 2:N;
   
   
   %--  b'   (O(h^2))----
   
   %-- Central difference
    
     
   %--  b''  (O(h^2))----
   
   bx2p(1,1) = (bx(2,1) + bext(1,1) -2*bx(1,1))/h^2 ;
   by2p(1,1) = (by(2,1) + bext(2,1) -2*by(1,1))/h^2 ;
   bz2p(1,1) = (bz(2,1) + bext(3,1) -2*bz(1,1))/h^2 ;
   
   bx2p(ig,1) = (bx(ig+1,1) + bx(ig-1,1) - 2*bx(ig,1))/h^2 ;
   by2p(ig,1) = (by(ig+1,1) + by(ig-1,1) - 2*by(ig,1))/h^2 ;
   bz2p(ig,1) = (bz(ig+1,1) + bz(ig-1,1) - 2*bz(ig,1))/h^2 ;
   
   bx2p(N+1,1) = (bext(4,1) + bx(N,1)   - 2*bx(N+1,1))/h^2 ;
   by2p(N+1,1) = (bext(5,1) + by(N,1)   - 2*by(N+1,1))/h^2 ;
   bz2p(N+1,1) = (bext(6,1) + bz(N,1)   - 2*bz(N+1,1))/h^2 ;
%    




  bpp = bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p;
  
  
   n = tau/2/pi;
%  E_3pi(i) = h*sum(bpp(1:N))/tau^2 - tau^2;
  
    E_3pi(i) = h*sum(bpp(5:N-5))/(16*pi^4*n^2) - n^2;
  
 % E_3pi(i) = h*sum(bpp(1:N));
  
  
  figure(1)
  
  if(st_3pi(i)>0)
     
      if(num_t(i)==3)
          colr = col{1};
      elseif(num_t(i)==5)
          colr = col{2};
      else
          colr = col{3};
      end
      
     
     % plot((2*pi)/tau1(i),E_3pi(i),'o','Color',colr,'MarkerSize',8)
%       plot(tau1(i)/(2*pi),E_3pi(i),'o','Color',colr,'MarkerSize',8)
      hold on
  else
      % plot((2*pi)/tau1(i),E_3pi(i),'or','MarkerSize',8)
  %     plot(tau1(i)/(2*pi),E_3pi(i),'or','MarkerSize',8)
       hold on 
  end
  
%==========================================================================  
 
 % th_3fold(i,:) = 180/pi*acos(bx(1,1)*bx +  by(1,1)*by + bz(1,1)*bz) ;
%   
%   figure(i*100)
% plotbrowser on 
% plot(th_3fold(i,:))
end
hold on 
plot(tau1/(2*pi),E_3pi)
title('Bending energy')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 
%     
% figure(3)
% plotbrowser on
% plot(tau1/(2*pi),E_3pi,'LineWidth',4)
% title('Bending energy')
% set(gca,'FontSize',25,'LineWidth',1)
% box on
% grid on
% hold on 

 