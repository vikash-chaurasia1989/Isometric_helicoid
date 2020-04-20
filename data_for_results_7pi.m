clear all
clc

 global   N1  N  qd   id   h   rho lm  tau bx by bz h id
 global bxp byp bzp bx2p by2p bz2p 

addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')
%==== data of number of twist and tau
%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6

  col = {col1,col2,col3,col4,col5} ;%==== data of number of twist and tau

%=========================== 7 fold =======================================

 
   
 tau_7 = [15.46 15.5:.1:16.1];
 
% tau_9 = [17.47837  17.479 17.48 17.4844];
tau_9 = 16.2:.1:24.9;

 tau_11 = [ 17.482926 17.4844 17.485 17.5:.1:19  19.5:.5:25];
 
 % % 
  tau1 = [tau_7 tau_9];% [tau_5 tau_7 tau_9 tau_11];
  N = 49;
  N1 = N-1;
  h = 1/N;
  
  
    %== check stability in loop
 %--- 3x3 identity matrix --

id = [1 0 0;0 1 0;0 0 1];

 %===== number of twist ====
 
 for i = 1:length(tau1)
     
     if(tau1(i)<=tau_7(end))
         num_t(i) = 7;
     elseif(tau_9(1)<=tau1(i) && tau1(i)<=tau_9(end))
         num_t(i) = 9;
     elseif(tau_9(1)<=tau1(i) && tau1(i)<=tau_9(end))
         num_t(i) =9;
     else
         num_t(i) = 11;
        
     end
     
     p1 = i;
    tau = tau1(p1);
    
    
 % % %     

    str0 =['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
    
   % str0  = '3fold_N72_tau_170.txt';
    
    %-- 5pi --
    %  str0 = ['5pi_knot_N120_tau_' num2str(10*tau1(p1)) '.txt'];
    
    
    strb = ['b_' str0];
    strr = ['r_' str0];
    
    strlm = ['lm_' str0];
    strrho = ['rho_' str0];
    struvw = ['uvw_' str0];
    strbext = ['bext_' str0];
    
    
    %--- load b
    temp = load(strb);
    
    N = length(temp)-1;
    N1 = N-1;
    h = 1/N;

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
bext = load(strbext);

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
  
  
  stbl(p1) = 1;%fun_stability(f);
  
  figure(1)
  hold on 
  plotbrowser on 
  
  if(stbl(p1)==1)
       if(num_t(p1)==7)
          colr = col{5};
      elseif(num_t(p1)==9)
          colr = col{2};
      elseif(num_t(p1)==8)
          colr = col{2};
      else
          colr = col{3};
      end
      plot(tau/(2*pi),num_t(p1),'square','Color',colr)
 else
      plot(tau/(2*pi),num_t(p1),'or')
  end
  hold on 
      
  end
  hold on 
  plot(tau1/(2*pi),num_t,'--k')
    set(gca,'FontSize',25,'LineWidth',1)
         box on
         grid on
%=== index of unstable points ==

% p1 = 4     6     8     9    10    11    12    18    21    25    26    59    62    65    66    73    82    89

