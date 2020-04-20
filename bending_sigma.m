%--- 
clear all
clc
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')


 global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext

  format longE
  
 
%-- This is the final file for calculating Linking number
%-- Here, I have used codes developed by Zin Arai and modified it for my
%purpose. The accuracy is impressive. Way better than the double integrals
%I was using to calculate Writhe and then the linking number using the
%Caligerano theorem  Lk = 2*(\tau/(2\pi) + Wr);

%---- read r and b file and create ruling data ----
  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
  col6 = [0.341, 0.505, 0.819]/2  ;         % -- color for mode n = 6
  col7 = [0.831, 0.545, 0.247]/2  ;      % -- color for mode n = 3

  
  col = {col1,col2,col3,col4,col5,col6,col7} ;

%temp = load('05_fold.txt');


for branch =1:1 
 %===== branch 1   
  
 if(branch==1)
%           tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%               13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%           tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%               13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%           tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
%      
%      
%      
%           tau1 =      [tau_3 tau_5 tau_7];
%      
   
%         
%           tau_3 = [ 8.0957344694550208E+00  8.09574 8.095745 8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.05:8.5 8.6:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%        tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%        tau_7 = [13.224  13.225 13.23 13.3 13.31 13.32  13.33 13.3 13.34 13.35  13.36 13.37 13.38 13.39 13.395  13.4:.1:13.9 14.5:.5:21 21.5:.5:33.5 33.6:.1:33.9 33.95 33.96 33.965 33.966 33.9665 33.96675 33.96676 33.96677 ];
% %      
%     
%      tau1 =      [tau_3 tau_5 tau_7];
%      
      


     N = 72;
     N1 = N-1;
     h = 1/N;
     nfold=3;
     
    str1 =  '3fold_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    
    strtau = 'tau_branch_1.txt';

    tau1 = load([path strtau]);

 elseif(branch==2)
     %================ Branch 2 =====================
     
     %tau1 = [ 11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
     
     
     N = 120;
     N1 = N-1;
     h = 1/N;
     nfold=5;
     
    str1 = '5pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
    
     strtau = 'tau_branch_2.txt';

    tau1 = load([path strtau]);
    
     % =============== Branch 3 =====================
     
 else
     tau1 = [15.45:.1:25 25.1:.1:40];
     
     N = 84;
     N1 = N-1;
     h = 1/N;
     
     nfold=7;
    str1 = '7pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
      strtau = 'tau_branch_3.txt';

    tau1 = load([path strtau]);
 
 end
%==========================================================================  
 
%==================

h = 1/N      ;

 
%strLk = [path 'Linking_branch' num2str(branch) '.txt'];
strLk = ['Linking_branch_' num2str(branch) '.txt'];

  
 Lk = load([path strLk]);
 st = zeros(1,length(tau1));
 E = zeros(1,length(tau1));

 
 
 rx = zeros(N+1,1);
 ry = zeros(N+1,1);
 rz = zeros(N+1,1);
 
 
 
 
 
for p1 =  1:length(tau1)
    
    tau =  tau1(p1);
    
    
    
    %-- 3 fold --
    
    str0 = [str1  num2str(N) '_tau_' num2str((10^10*tau)) '.txt'];
    
    if(nfold==7)
    %    str0 = [str1  num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    end
    
    
    
    strb = [path 'b_' str0];
    strr = [path  'r_' str0];
    
    strlm = [path  'lm_' str0];
    strrho = [path 'rho_' str0];
    struvw = [path 'uvw_' str0];
     
    
    
    %--- load b
    temp = load(strb);
    
    bx  = temp(1:N+1,1);
    by = temp(1:N+1,2);
    bz = temp(1:N+1,3);
    
    
    bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
    bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
    
    
    %  dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);
    %
    % tau = sum(dl);
    %
    
    i = 1:N;
    
    tx = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
    ty = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
    tz = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
    
    
    
    tx(N+1,1) = tx(1,1);
    ty(N+1,1) = ty(1,1);
    tz(N+1,1) = tz(1,1);
    
    tx = tx/tau/h;
    ty = ty/tau/h;
    tz = tz/tau/h;
    
    
    rx(1,1) = 0;
    ry(1,1) = 0;
    rz(1,1) = 0;
    
    for i=1:N
        
        rx(i+1,1) = rx(i,1) + h*tx(i,1);
        ry(i+1,1) = ry(i,1) + h*ty(i,1);
        rz(i+1,1) = rz(i,1) + h*tz(i,1);
    end
    
    
    
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
    
    E(p1) = h*sum(bpp(1:N))/2 - 8*pi^4*n^4;
    
    sig = .1;
    
    
    
    
    % E_3pi(i) = h*sum(bpp(1:N));
    
    %===Loading the linking number
    
    
    
    figure(1)
    plotbrowser on
    
    
    
    icol =  round((Lk(p1)-1)/2);
     plot(tau1(p1)/(2*pi),E(p1),'-ok','LineWidth',.5,'color',col{icol})
    
    % plot(tau1(p1)/(2*pi),Lk(p1),'-ok','LineWidth',.5)
    % else
    %   plot(tau1(p1)/(2*pi),Lk(p1),'-or','LineWidth',.5)
    
     hold on
     set(gca,'FontSize',25,'LineWidth',.5)
    %hold on
    
end

plot(tau1/(2*pi),E,'--k','LineWidth',.5)
hold on 

% sig = .1;
% n = tau1/2/pi;
% fac = asinh(n*pi*sig)./(n*pi*sig);
% plot(tau1/(2*pi),fac.*E,'--k','LineWidth',.5)
% hold on 
% sig = .5;
% n = tau1/2/pi;
% fac = asinh(n*pi*sig)./(n*pi*sig);
% 
% plot(tau1/(2*pi),fac.*E,'--b','LineWidth',.5)

title('Linking number')

box on
grid on
hold on

end
