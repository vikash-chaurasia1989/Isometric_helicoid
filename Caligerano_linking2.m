clear all
clc
 
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


for branch =3:3 
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
 
 
%==================
 
%==================

h = 1/N      ;

 
  
 Lk = zeros(1,length(tau1));
 st = zeros(1,length(tau1));
 
 
 
 rx = zeros(N+1,1);
 ry = zeros(N+1,1);
 rz = zeros(N+1,1);
 
 
 
 
 
for p1 =  1:length(tau1)
    
    
    if(branch==1)
        N=72;
        N1 = N-1;
    elseif(branch==2)
        N=120;
        N1=N-1;
    else
        N=84;
        N1 = N-1;
    end
    
    tau =  tau1(p1);
      %-- 3 fold --
    
      str0 = [str1  num2str(N) '_tau_' num2str((10^10*tau)) '.txt'];
    
    if(nfold==7)
     str0 = [str1  num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    end
    
%     %-- 5 fold ----    
%     str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
%  
%     if(tau>13.9)
%             str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
%     end
%     
   % str0  = '3fold_N72_tau_170.txt';
    
    %-- 5pi --
    %  str0 = ['5pi_knot_N120_tau_' num2str(10*tau1(p1)) '.txt'];
    
    
    %  7pi 
    
    % str0 =['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
    
      
 
    
    strb = [path 'b_' str0];
    strr = [path  'r_' str0];
    
    strlm = [path  'lm_' str0];
    strrho = [path 'rho_' str0];
    struvw = [path 'uvw_' str0];
    strbext = [path 'bext_' str0];
    
    
    
 %==========================================================================
%==========================================================================



    bx = zeros(N+1,1);
    by = zeros(N+1,1);
    bz = zeros(N+1,1);
   
    rx = zeros(N+1,1);
    ry = zeros(N+1,1);
    rz = zeros(N+1,1);
    
    tx = zeros(N+1,1);
    ty = zeros(N+1,1);
    tz = zeros(N+1,1);
    
    rx0 = zeros(N+1,1);
    ry0 = zeros(N+1,1);
    rz0 = zeros(N+1,1);
    
    rx1 = zeros(N+1,1);
    ry1 = zeros(N+1,1);
    rz1 = zeros(N+1,1);
    
    
    %================================================
  
    %--- load b
    temp = load(strb);
    
    bx = temp(:,1);
    by = temp(:,2);
    bz = temp(:,3);
    
    
    %============== Refining data for more accuarate results =====
    %-- load r --
    
    
 s1 = (linspace(0,1,N+1))';
 %
 N = 500;
 N1 = N-1;
 h = 1/N;
 
 s2 = (linspace(0,1,N+1))';
 
 bx = interp1q(s1,bx,s2)  ;
 by = interp1q(s1,by,s2)  ;
 bz = interp1q(s1,bz,s2)  ;
 
 dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);
 tau2 = sum(dl);
 %========  creating tangent from binormal ----
    
      %-- Central difference
   %---O(h^2)---
 bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

   
   bx1 = bext(1);
   by1 = bext(2);
   bz1 = bext(3);
   
   bxN2 = bext(4);
   byN2 = bext(5);
   bzN2 = bext(6);
   
   
     bxp(1,1) = (bx(2,1)-bext(1,1))/(2*h);
     byp(1,1) = (by(2,1)-bext(2,1))/(2*h);
     bzp(1,1) = (bz(2,1)-bext(3,1))/(2*h);
     
     ig = 2:N;
 
     bxp(ig,1) = (bx(ig+1,1)-bx(ig,1))/(h);
     byp(ig,1) = (by(ig+1,1)-by(ig,1))/(h);
     bzp(ig,1) = (bz(ig+1,1)-bz(ig,1))/(h);
     
 
    bxp(N+1,1) = (bext(4,1)-bx(N,1))/(2*h);
    byp(N+1,1) = (bext(5,1)-by(N,1))/(2*h);
    bzp(N+1,1) = (bext(6,1)-bz(N,1))/(2*h);
    
   %----- t = bxb'/tau
   
   %----------------- Post processing -----
%---- Tangent ti = bi \times bi+1

i = 1:N+1;


% tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
% ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
% tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
%
tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);

%
%   tx(N+1,1) = tx(1,1);
%   ty(N+1,1) = ty(1,1);
%   tz(N+1,1) = tz(1,1);
% %
% %
 tx = tx/tau;
 ty = ty/tau;
 tz = tz/tau;
 
 
 
 % %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;

% i = 1 ;
%
%     rx(i+1) = rx(i) + h*tx(i);
%     ry(i+1) = ry(i) + h*ty(i);
%     rz(i+1) = rz(i) + h*tz(i);
%
% for i = 2:N-1
%     rx(i+1) = rx(i-1) + 2*h*tx(i);
%     ry(i+1) = ry(i-1) + 2*h*ty(i);
%     rz(i+1) = rz(i-1) + 2*h*tz(i);
% end
%


for i = 1:N

     rx(i+1,1) = rx(i,1) + h*tx(i,1);
     ry(i+1,1) = ry(i,1) + h*ty(i,1);
     rz(i+1,1) = rz(i,1) + h*tz(i,1);

end
    
    
%     plot3(rx,ry,rz)
%     hold on 
    %---- create edge of the binormal scroll using b and r
    wd = .01;
    
  %  wd = .1;
    v = linspace(-wd,wd,10);
    
 %--- edge   ----
 
    
    rx0 = rx - v(1)*bx;
    ry0 = ry - v(1)*by;
    rz0 = rz - v(1)*bz;
    
    tx0 = tx - v(1)*bxp;
    ty0 = ty - v(1)*byp;
    tz0 = tz - v(1)*bzp;
   
    
%==========================================================================

%-------- self linking number ---
 
 
 
for i = 1:N     
     % 
     j = 1:N+1     ;             % j array index
     j = j(j~=i);   % indexes not containing i
      
         %--------------------- 1. Curve r(s)  -------------------------------
         
      rij  =  ((rx(j,1)-rx(i,1)).^2 + (ry(j,1)-ry(i,1)).^2 +(rz(j,1)-rz(i,1)).^2).^.5 ;
       
      %-- for self linking number
      tjix = ty(j,1).*tz(i,1) - tz(j,1).*ty(i,1);
      tjiy = tz(j,1).*tx(i,1) - tx(j,1).*tz(i,1);
      tjiz = tx(j,1).*ty(i,1) - ty(j,1).*tx(i,1);
         
      
      temp =   (tjix.*(rx(j,1)-rx(i,1)) + tjiy.*(ry(j,1)-ry(i,1)) + tjiz.*(rz(j,1)-rz(i,1)))./rij.^3 ;
        
      ker(i,1) = trapz(temp);
      
      
      
end

Wr  =  h^2*trapz(ker)/(4*pi);
Lk(p1)  =  (2*(tau/(2*pi) + Wr ))  ;


 
    
 icol =  round((Lk(p1)-1)/2);  
% plot(tau1(p1)/(2*pi),Lk(p1),'-ok','LineWidth',.5,'color',col{icol})
    
  plot(tau1(p1)/(2*pi),Lk(p1),'-ok','LineWidth',.5)
% else
%   plot(tau1(p1)/(2*pi),Lk(p1),'-or','LineWidth',.5)
 
hold on 
set(gca,'FontSize',25,'LineWidth',.5)
hold on 

end

plot(tau1/(2*pi),Lk,'--k','LineWidth',.5)
  
title('Linking number')

box on
grid on
hold on

end
