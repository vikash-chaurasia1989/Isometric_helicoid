%==== This file checks stability 


clear all
clc

  global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext

  format longE
  tic

 % addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
%******************************** Main file for two chain interaction *********

for branch = 3:3 
 %===== branch 1   
  
 if(branch==1)
     tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
         13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
     tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
         13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
     
     tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
     
     
     
     tau1 =      [tau_3 tau_5 tau_7];
     
     
     N = 72;
     N1 = N-1;
     h = 1/N;
     nfold=3;
     
 elseif(branch==2)
     %================ Branch 2 =====================
     
     tau1 = [ 11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
     
     
     N = 120;
     N1 = N-1;
     h = 1/N;
     nfold=5;
     
    % =============== Branch 3 =====================
     
 else
     tau1 = [15.45:.1:25 25.1:.1:40];
     
     N = 84;
     N1 = N-1;
     h = 1/N;
     
     nfold=7;
 end
%==========================================================================  
if(nfold==3)
    str1 =  '3fold_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_2/';
elseif(nfold==5)
    str1 = '5pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_2/';
%  
else
    str1 = '7pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_2/';
end


  st = zeros(1,length(tau1));
%==========================================================================  
for p1  = 1:length(tau1)
    sv = 1;
    tau =   tau1(p1);
    
    
    str0 = [str1  num2str(N) '_tau_' num2str((10^10*tau)) '.txt'];
    
    if(nfold==7)
     str0 = [str1  num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    end
    
    strb =  [path 'b_' str0] ;
    strlm = [path 'lm_' str0];
    strrho = [path 'rho_' str0];
    struvw = [path 'uvw_' str0];
    
    
    %==== load b ===========
    temp = load(strb);
    
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
      %=== load lm ====
    lm =   load(strlm);
      
    f = [b(1:3*N1)  rho' lm' u v w  ]'  ;         % for fun_jacobian6
    
    
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
    
    
    %======= calling fun_stability2 ====
    
    st(p1) = fun_stability2(f);
    
    
    
end
plot(st,'-o','LineWidth',.5)
hold on 

end