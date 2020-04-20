clear all
clc

  global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig

  format longE
  tic

 branch =  1;   
    
 if(branch==1)   % 8.0957344694550208E+00   
%          tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%      tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%      tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
         
       tau_3 = [ 8.0957344694550208E+00  8.09574 8.095745 8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.05:8.5 8.6:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2...
                  13.3 13.31 13.32  13.33 13.335 13.3375 13.33755];
       tau_5 =   [ 13.337575 13.3376 13.33775];     
     
       tau_7 = [ 13.338 13.3385 13.34 13.35  13.36 13.37 13.38 13.39 13.395  13.4:.1:13.9 14.5:.5:18 18.1:.1:22.9 23:.5:27 27.01:.01:27.05 27.0525 27.0535 27.05355 27.053575];  
       tau_9 = [27.05375 27.054 27.0545  27.055 27.06:.01:27.13 27.1325 27.132575 27.133];
       tau_11 = [27.1335 27.135 27.14 27.15 27.2:.1:27.4 27.5:.5:33.5 33.6:.1:33.9 33.95 33.96 33.965 33.966 33.9665 33.96675 33.96676 33.96677 ];
%      
      
     
     tau1 =      [tau_3 tau_5 tau_7 tau_9 tau_11];
     
     
     N = 72;
     N1 = N-1;
     h = 1/N;
     nfold=3;
     
 elseif(branch==2) % 1.1984230303215979E+01
     %================ Branch 2 =====================
     
    % tau1 = [1.1984233303215979E+01   11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
     tau1 = [1.1984233303215979E+01   11.984235 11.98425  11.9845 11.985 12:.1:16.5 16.55 16.56 16.561 16.5625 16.565 16.57 16.575...
              16.6: 16.7 16.8 16.81 16.815 16.8175 16.8176 16.81775 16.818 16.82 16.825 16.85 16.86 16.865 16.8675 16.8685 16.8695 16.869575 16.8696 16.86975 16.8699 16.87 16.875 16.9:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];

     
     N = 120;
     N1 = N-1;
     h = 1/N;
     nfold=5;
     
    % =============== Branch 3 =====================
     
 else
     tau1 = [15.45:.1:20.95 20.96:.01:21.25  21.35:.1:25 25.1:.1:40];
     
     N = 84;
     N1 = N-1;
     h = 1/N;
     
     nfold=7;
 end
%==========================================================================  
p1 = 7;

tau = tau1(p1);

n = tau/2/pi

sig1 = [.01 .02 .03 .04 .05 .06];

str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '_sig_' num2str(.01) '.txt'];

path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';

strb = ['b_' str0];
strr = ['r_' str0];

strlm = ['lm_' str0];
strrho = ['rho_' str0];
struvw = ['uvw_' str0];

lm = load([path strlm]);
rho = load([path strrho]);
 

figure(1)
plotbrowser on
plot(lm,'LineWidth',1)
title('\lambda')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 


figure(2)
plotbrowser on
plot(rho,'LineWidth',1)
title('\rho')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 


for p1  =1:length(sig1)-1
    
    
    sig = sig1(p1);
    
    if(nfold==3)
        %  str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
        str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '_sig_' num2str(sig) '.txt'];
        
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    elseif(nfold==5)
        str0 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
        %
    else
        str0 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
        path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
    end
    
    
    
    strb = ['b_' str0];
    strr = ['r_' str0];
    
    strlm = ['lm_' str0];
    strrho = ['rho_' str0];
    struvw = ['uvw_' str0];
    
    
       fac = asinh(n*pi*sig)/(4*pi^3*n^3);
    
    lm = load([path strlm]);
    rho = load([path strrho]);
    
    lm = lm/fac;
    rho = rho/fac;
    
    figure(1)
    plotbrowser on
    plot(lm,'LineWidth',1)
    title('\lambda')
    set(gca,'FontSize',25,'LineWidth',1)
    box on
    grid on
    
    figure(2)
    plotbrowser on
    plot(rho,'LineWidth',1)
    title('\rho')
    set(gca,'FontSize',25,'LineWidth',1)
    box on
    grid on


    
end