clear all
clc

 global   N1  N  qd   id   h   rho lm  tau bx by bz h id
 global bxp byp bzp bx2p by2p bz2p 

  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
%==== data of number of twist and tau
%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6

  col = {col1,col2,col3,col4,col5} ;

%=========================== 5 fold =======================================

%===== tau for 3 fold ====
 
tau_3 = [8.1:.1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
             13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
         
tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
             13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
         
tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21];

%tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.1:.1:13.9 14:.5:21];
tau1 =      [tau_3 tau_5 tau_7];
N = 72;
h = 1/N;

N1 = N-1;
  %-- index of unstable points --
  %indx =  11    12    30    38    51    54    57    66    69    70    75    76    79    80   100
   
 %===== number of twist ====
 id = [1 0 0;0 1 0;0 0 1];

 for i = 1:length(tau1)
     
     
     %==== Number of twist =====
     
     if(tau1(i)<=tau_3(end))
         num_t(i) = 3;
     elseif(tau_5(1)<=tau1(i) && tau1(i)<=tau_5(end))
         num_t(i) = 5;
     else
         num_t(i) = 7;
         
     end
  
     p1 = i;
       tau = tau1(p1);

       %===== Checking stability =====
       
     str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
 
    
    
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
  
  
  stbl(p1) = fun_stability(f);
  
  figure(1)
  hold on 
  plotbrowser on 
  
  if(stbl(p1)==1)
      if(num_t(p1)==3)
          colr = col{1};
      elseif(num_t(p1)==5)
          colr = col{2};
      else
          colr = col{3};
      end
      
      plot(tau/(2*pi),num_t(p1),'o','Color',colr)
  else
      plot(tau/(2*pi),num_t(p1),'or')
  end
  hold on 
      
  
  
  
 end
 
   
%=========================================================================  
   hold on 
  plot(tau1/(2*pi),num_t,'--k')
    set(gca,'FontSize',25,'LineWidth',1)
         box on
         grid on
%=== index of unstable points ==

% p1 = 4     6     8     9    10    11    12    18    21    25    26    59    62    65    66    73    82    89

