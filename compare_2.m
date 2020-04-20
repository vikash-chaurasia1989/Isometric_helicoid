  clear all
  clc
  
  
  global   N1  N  qd   id   h   rho lm  tau tau1 p1
 
  global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p nfold
  
  %--------Array initialization -------
    %  
   
% str0 =   '3fold_N72_tau_80995000000.txt' ;
  %===== branch 1   
  
tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
             13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;          
tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
             13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
         
tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];



%tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.1:.1:13.9 14:.5:21];
tau1 =      [tau_3 tau_5 tau_7];


N = 72;
  N1 = N-1;
  h = 1/N;
  nfold=3;
  
  p1 = 7;
  
%  
if(p1>0)
    tau2 = tau1(p1);
    
    if(nfold==3)
        str1 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path1 = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi/';
   
        str2 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path2 = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_2/';
  
    
    elseif(nfold==5)
        str1 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path1 = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi/';
       
        str2 =['5pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path2 = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_2/';
       
        %
    else
        str1 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path1 = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi/';
    
        str2 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau2) '.txt'];
        path2 = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_2/';

    
    end
    
    
end

 
   
strb =  [path1 'b_' str1] ;
strlm = [path1 'lm_' str1];
strrho = [path1 'rho_' str1];
struvw = [path1 'uvw_' str1];
 
 
%==== load b ===========
temp = load(strb);
N   = length(temp)-1;
N1 = N-1;

h = 1/N;
% N = N/2;
bx = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);


%=== load rho ==========   
   rho =      load(strrho);
  
    
%=== load lm ====  
   lm =   load(strlm);
      

  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
  
  %---- derivative calculation for b2 --- bN  ---
 
   ig = 2:N;
   
 
   
 
  %=== O(h)====
     bxp(1,1) = (bx(1,1)-bext(1,1))/(h);
     byp(1,1) = (by(1,1)-bext(2,1))/(h);
     bzp(1,1) = (bz(1,1)-bext(3,1))/(h);
 
     bxp(ig,1) = (bx(ig,1)-bx(ig-1,1))/(h);
     byp(ig,1) = (by(ig,1)-by(ig-1,1))/(h);
     bzp(ig,1) = (bz(ig,1)-bz(ig-1,1))/(h);
     
 
    bxp(N+1,1) = (bx(N+1,1)-bx(N,1))/(h);
    byp(N+1,1) = (by(N+1,1)-by(N,1))/(h);
    bzp(N+1,1) = (bz(N+1,1)-bz(N,1))/(h);
    
  %     
   %   lmp = (lm(ig+1,1)-lm(ig-1))/(2*h);
     lmp = (lm(ig,1)-lm(ig-1,1))/(h);
   
          
   %--  b''  (O(h^2))----
   
   bx2p(1,1) = (bx(2,1) + bext(1,1) -2*bx(1,1))/h^2 ;
   by2p(1,1) = (by(2,1) + bext(2,1) -2*by(1,1))/h^2 ;
   bz2p(1,1) = (bz(2,1) + bext(3,1) -2*bz(1,1))/h^2 ;
   
   bx2p(ig,1) = (bx(ig+1,1)  + bx(ig-1,1) - 2*bx(ig,1))/h^2 ;
   by2p(ig,1) = (by(ig+1,1)  + by(ig-1,1)  - 2*by(ig,1))/h^2 ;
   bz2p(ig,1) = (bz(ig+1,1)  + bz(ig-1,1)  - 2*bz(ig,1))/h^2 ;
   
   bx2p(N+1,1) = (bext(4,1) + bx(N,1)   - 2*bx(N+1,1))/h^2 ;
   by2p(N+1,1) = (bext(5,1) + by(N,1)   - 2*by(N+1,1))/h^2 ;
   bz2p(N+1,1) = (bext(6,1) + bz(N,1)   - 2*bz(N+1,1))/h^2 ;
%    
     
%    %=====================================================
   




   bx4p(1,1) = (bx(4,1) - 4*bx(3,1) + 6*bx(2,1)-4*bx(1,1)+bext(1,1))/h^4;
   by4p(1,1) = (by(4,1) - 4*by(3,1) + 6*by(2,1)-4*by(1,1)+bext(2,1))/h^4;
   bz4p(1,1) = (bz(4,1) - 4*bz(3,1) + 6*bz(2,1)-4*bz(1,1)+bext(3,1))/h^4;
   
   ig = 3:N-1;
   
   bx4p(ig-1,1) = (bx(ig+2,1) - 4*bx(ig+1,1) + 6*bx(ig,1) - 4*bx(ig-1,1) + bx(ig-2,1))/(h^4);
   by4p(ig-1,1) = (by(ig+2,1) - 4*by(ig+1,1) + 6*by(ig,1) - 4*by(ig-1,1) + by(ig-2,1))/(h^4);
   bz4p(ig-1,1) = (bz(ig+2,1) - 4*bz(ig+1,1) + 6*bz(ig,1) - 4*bz(ig-1,1) + bz(ig-2,1))/(h^4);
 
   bx4p(N-1,1)   = (bext(4,1)  - 4*bx(N+1,1)  + 6*bx(N,1)  - 4*bx(N-1,1) + bx(N-2,1) )/h^4;
   by4p(N-1,1)   = (bext(5,1)  - 4*by(N+1,1)  + 6*by(N,1)  - 4*by(N-1,1) + by(N-2,1) )/h^4;
   bz4p(N-1,1)   = (bext(6,1)  - 4*bz(N+1,1)  + 6*bz(N,1)  - 4*bz(N-1,1) + bz(N-2,1) )/h^4;
        
        
figure(1)
plotbrowser on
plot3(bx,by,bz,'--ok','LineWidth',1)
title('binormal')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 

  

figure(2)
plotbrowser on
plot(bxp,'--ok','LineWidth',1)
title('bxp')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 

figure(3)
plotbrowser on
plot(bx2p,'--ok','LineWidth',1)
title('bx2p')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on

figure(4)
plotbrowser on
plot(bx4p,'--ok','LineWidth',1)
title('bx4p')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on
   
figure(5)
plotbrowser on
plot(lm,'-ok','LineWidth',1)
title('\lambda')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on


figure(6)
plotbrowser on
plot(rho,'-ok','LineWidth',1)
title('\rho')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 
   
   
   
   
   
 %========== data for the correct discretization scheme =========
 
 
 
    
strb =  [path2 'b_' str2] ;
strlm = [path2 'lm_' str2];
strrho = [path2 'rho_' str2];
struvw = [path2 'uvw_' str2];
 
 
%==== load b ===========
temp = load(strb);
N   = length(temp)-1;
N1 = N-1;

h = 1/N;
% N = N/2;
bx = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);

%=== load rho ==========   
   rho =      load(strrho);
  
    
%=== load lm ====  
   lm =   load(strlm);
      
  bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
  bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
  
   %---- derivative calculation for b2 --- bN  ---
 
   ig = 2:N;
   
 
   
 
  %=== O(h)====
     bxp(1,1) = (bx(1,1)-bext(1,1))/(h);
     byp(1,1) = (by(1,1)-bext(2,1))/(h);
     bzp(1,1) = (bz(1,1)-bext(3,1))/(h);
 
     bxp(ig,1) = (bx(ig,1)-bx(ig-1,1))/(h);
     byp(ig,1) = (by(ig,1)-by(ig-1,1))/(h);
     bzp(ig,1) = (bz(ig,1)-bz(ig-1,1))/(h);
     
 
    bxp(N+1,1) = (bx(N+1,1)-bx(N,1))/(h);
    byp(N+1,1) = (by(N+1,1)-by(N,1))/(h);
    bzp(N+1,1) = (bz(N+1,1)-bz(N,1))/(h);
    
  %     
   %   lmp = (lm(ig+1,1)-lm(ig-1))/(2*h);
     lmp = (lm(ig,1)-lm(ig-1,1))/(h);
   
          
   %--  b''  (O(h^2))----
   
   bx2p(1,1) = (bx(2,1) + bext(1,1) -2*bx(1,1))/h^2 ;
   by2p(1,1) = (by(2,1) + bext(2,1) -2*by(1,1))/h^2 ;
   bz2p(1,1) = (bz(2,1) + bext(3,1) -2*bz(1,1))/h^2 ;
   
   bx2p(ig,1) = (bx(ig+1,1)  + bx(ig-1,1) - 2*bx(ig,1))/h^2 ;
   by2p(ig,1) = (by(ig+1,1)  + by(ig-1,1)  - 2*by(ig,1))/h^2 ;
   bz2p(ig,1) = (bz(ig+1,1)  + bz(ig-1,1)  - 2*bz(ig,1))/h^2 ;
   
   bx2p(N+1,1) = (bext(4,1) + bx(N,1)   - 2*bx(N+1,1))/h^2 ;
   by2p(N+1,1) = (bext(5,1) + by(N,1)   - 2*by(N+1,1))/h^2 ;
   bz2p(N+1,1) = (bext(6,1) + bz(N,1)   - 2*bz(N+1,1))/h^2 ;
%    
     
%    %=====================================================
   




   bx4p(1,1) = (bx(4,1) - 4*bx(3,1) + 6*bx(2,1)-4*bx(1,1)+bext(1,1))/h^4;
   by4p(1,1) = (by(4,1) - 4*by(3,1) + 6*by(2,1)-4*by(1,1)+bext(2,1))/h^4;
   bz4p(1,1) = (bz(4,1) - 4*bz(3,1) + 6*bz(2,1)-4*bz(1,1)+bext(3,1))/h^4;
   
   ig = 3:N-1;
   
   bx4p(ig-1,1) = (bx(ig+2,1) - 4*bx(ig+1,1) + 6*bx(ig,1) - 4*bx(ig-1,1) + bx(ig-2,1))/(h^4);
   by4p(ig-1,1) = (by(ig+2,1) - 4*by(ig+1,1) + 6*by(ig,1) - 4*by(ig-1,1) + by(ig-2,1))/(h^4);
   bz4p(ig-1,1) = (bz(ig+2,1) - 4*bz(ig+1,1) + 6*bz(ig,1) - 4*bz(ig-1,1) + bz(ig-2,1))/(h^4);
 
   bx4p(N-1,1)   = (bext(4,1)  - 4*bx(N+1,1)  + 6*bx(N,1)  - 4*bx(N-1,1) + bx(N-2,1) )/h^4;
   by4p(N-1,1)   = (bext(5,1)  - 4*by(N+1,1)  + 6*by(N,1)  - 4*by(N-1,1) + by(N-2,1) )/h^4;
   bz4p(N-1,1)   = (bext(6,1)  - 4*bz(N+1,1)  + 6*bz(N,1)  - 4*bz(N-1,1) + bz(N-2,1) )/h^4;
        
        
figure(1)
plotbrowser on
plot3(bx,by,bz,'--or','LineWidth',1)
title('binormal')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 

 
figure(2)
plotbrowser on
plot(bxp,'--or','LineWidth',1)
title('bxp')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 

figure(3)
plotbrowser on
plot(bx2p,'--or','LineWidth',1)
title('bx2p')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on

figure(4)
plotbrowser on
plot(bx4p,'--or','LineWidth',1)
title('bx4p')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on
   
figure(5)
plotbrowser on
plot(lm,'-or','LineWidth',1)
title('\lambda')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on


figure(6)
plotbrowser on
plot(rho,'-or','LineWidth',1)
title('\rho')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 
   
   
   
 
 