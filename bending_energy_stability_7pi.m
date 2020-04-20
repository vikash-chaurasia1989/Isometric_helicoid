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

 addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
 addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi')
 addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')
  %==== data of number of twist and tau
%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6

  col = {col1,col2,col3,col4,col5} ;
%--- identity matrix 

id = [1 0 0;0 1 0;0 0 1];
%================  3pi twist ==========
 
 
%==== 7 pi twist ====

tau0 = 15.45:.1:25;


 tau_7 = tau0(1:71);
 tau_9 = tau0(72:73);
 tau_11 = tau0(74:82);
 tau_13 = tau0(83);
 tau_15 = tau0(84:end);
  
%  tau0 = [15.3:.1:22.1 22.165:.05:25];
%  tau1 = tau0;
%   
%   N = 84;
%   N1 = N-1;
%   h = 1/N;  
%    tau1 = tau1;
   
 tau1 = [15.45:.1:25 25.1:.1:40];

   
  N = 84;
  N1 = N-1;
  h = 1/N;  
  
  path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_2/';

for i = 1:length(tau1)
    p1 = i;
    
    
    if(tau1(i)<=tau_7(end))
        num_t(i) = 7;
    elseif(tau_9(1)<=tau1(i) && tau1(i)<=tau_9(end))
        num_t(i) = 9;
    elseif(tau_11(1)<=tau1(i) && tau1(i)<=tau_11(end))
        num_t(i) =11;
    elseif(tau_13(1)<=tau1(i) && tau1(i)<=tau_13(end))
        num_t(i) =13;
    else
        num_t(i) = 15;
        
    end
     
        
    
    tau = tau1(i);
    
    
     str0 =['7pi_knot_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];

    
   
   
strb =  [path 'b_' str0] ;
strlm = [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];
 
     
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

  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

bx1 = bext(1);
by1 = bext(2);
bz1 = bext(3);

bxN2 = bext(4);
byN2 = bext(5);
bzN2 = bext(6);


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
  
  
 bx2p = zeros(N+1,1);
 by2p = zeros(N+1,1);
 bz2p = zeros(N+1,1);

    
%==========================================================================
%====================  Stability Check ====================================
 
 st_7pi(i) = 1;%fun_stability(f);

%==========================================================================


%==========================================================================
%=======================  Bending energy calculation ======================
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

  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
  
  
 bext(1,1) = bx1;
 bext(2,1)    = by1;
 bext(3,1)    = bz1;

bext(4,1)     = bxN2;
bext(5,1)     = byN2; 
bext(6,1)     = bzN2; 
%--- 


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
  
    E_7pi(i) = h*sum(bpp(1:N))/(16*pi^4*n^2) - n^2;  
  
  
  figure(1)
  
  if(st_7pi(i)==1)
      
       if(num_t(p1)==7)
          colr = col{1};
       elseif(num_t(p1)==9)
          colr = col{2};
       elseif(num_t(p1)==11)
           colr = col{3};
       elseif(num_t(p1) ==15)
           colr = col{4};
       else 
           colr = [0,0,0];
       end
      
     % plot(tau1(i)/(2*pi),E_7pi(i),'square','color',colr,'MarkerSize',8)
      plot(p1,E_7pi(i),'square','color',colr,'MarkerSize',8)

      hold on
  else
      % plot(tau1(i)/(2*pi),E_7pi(i),'square','color','r','MarkerSize',8)
       hold on 
  end
   

end
hold on 
%plot(tau1/2/pi,E_7pi,':k')
  title('Bending energy')
  set(gca,'FontSize',25,'LineWidth',1)
  box on
  grid on
% hold on 
% ind = length(tau1);
% figure(3)
% plotbrowser on
% plot(tau1(1:ind)/(2*pi),E_5pi(1:ind),'r','LineWidth',4)
% hold on 
% plot(tau1(ind:end)/(2*pi),E_5pi(ind:end),'--r','LineWidth',4)
% hold on 
% plot(tau1(ind)/(2*pi),E_5pi(ind),'.r','LineWidth',4)
% 
% title('Bending energy')
% set(gca,'FontSize',25,'LineWidth',1)
% box on
% grid on



