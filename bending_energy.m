%--- In this file, we read saved data and compute bending energy for each
%given \tau 
clear all
clc


%================  3pi twist ==========
tau1 = [8.1:.1:9.5 10:1:19];

N = 72;
h = 1/N;

for i = 1:length(tau1)
    
    tau = tau1(i);
    
    
    strb = ['b_3fold_N72_tau_' num2str(10*tau) '.txt'];
    strbext = ['bext_3fold_N72_tau_' num2str(10*tau) '.txt'];
    
   strlm = ['lm_3fold_N72_tau_' num2str(10*tau) '.txt'];

    
   %==== b==== 
     temp = load(strb);
 
           bx = temp(:,1);
           by = temp(:,2);
           bz = temp(:,3);
   %=============
   
   %=== b_extrapolated =====
    bext = load(strbext);
    
    
    %=== \lambda =====
    
    lm = load(strlm);
    
 
   
       ig = 2:N;

   
    %========== b' O(h)===============
     bxp(1,1) = (bx(2,1)-bext(1,1))/(2*h);
     byp(1,1) = (by(2,1)-bext(2,1))/(2*h);
     bzp(1,1) = (bz(2,1)-bext(3,1))/(2*h);
 
     bxp(ig,1) = (bx(ig+1,1)-bx(ig-1,1))/(2*h);
     byp(ig,1) = (by(ig+1,1)-by(ig-1,1))/(2*h);
     bzp(ig,1) = (bz(ig+1,1)-bz(ig-1,1))/(2*h);
     
 
    bxp(N+1,1) = (bext(4,1)-bx(N,1))/(2*h);
    byp(N+1,1) = (bext(5,1)-by(N,1))/(2*h);
    bzp(N+1,1) = (bext(6,1)-bz(N,1))/(2*h);
    
   %==============================
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


%--- 

  bpp = bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p;
  
  
  E(i) = h*sum(bpp(1:N))/tau^2 - tau^2;
  
  
 
 % th_3fold(i,:) = 180/pi*acos(bx(1,1)*bx +  by(1,1)*by + bz(1,1)*bz) ;
%   
%   figure(i*100)
% plotbrowser on 
% plot(th_3fold(i,:))
end
    
figure(3)
plotbrowser on
plot(tau1/(2*pi),E,'LineWidth',4)
title('Bending energy')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
hold on 



%==== 5 pi twist ====


%================  3pi twist ==========
 
tau2 = [12 12.5 13 14:.5:25] ;
    tau2 = [11.9 12:.1:13.9 14:.5:25];
    
    
 tau2 = [12:.1:13.9 14:.5:25];
 tau2 = [11.97871 11.98071 11.98571 11.99071 11.99571  tau2];
%tau2 = 14:.5:25;
N = 120;
h = 1/N;

for i = 1:length(tau2)
    
    tau = tau2(i);
    
    
 %str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
   str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];

   if(tau>13.9)
          str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
   end
   
   if(tau ==11.9)
                 str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];

   end

 
strb = ['b_' str0];
strr = ['r_' str0];

strlm = ['lm_' str0];
strrho = ['rho_' str0];
struvw = ['uvw_' str0];
strbext = ['bext_' str0];
    
     
   %==== b==== 
     temp = load(strb);
 
           bx = temp(:,1);
           by = temp(:,2);
           bz = temp(:,3);
   %=============
   
   %=== b_extrapolated =====
    bext = load(strbext);
    
    
    %=== \lambda =====
    
    lm = load(strlm);
    
 
   
       ig = 2:N;

   
    %========== b' O(h)===============
     bxp(1,1) = (bx(2,1)-bext(1,1))/(2*h);
     byp(1,1) = (by(2,1)-bext(2,1))/(2*h);
     bzp(1,1) = (bz(2,1)-bext(3,1))/(2*h);
 
     bxp(ig,1) = (bx(ig+1,1)-bx(ig-1,1))/(2*h);
     byp(ig,1) = (by(ig+1,1)-by(ig-1,1))/(2*h);
     bzp(ig,1) = (bz(ig+1,1)-bz(ig-1,1))/(2*h);
     
 
    bxp(N+1,1) = (bext(4,1)-bx(N,1))/(2*h);
    byp(N+1,1) = (bext(5,1)-by(N,1))/(2*h);
    bzp(N+1,1) = (bext(6,1)-bz(N,1))/(2*h);
    
   %==============================
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


%--- 

  bpp = bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p;
  
  
  E_5pi(i) = h*sum(bpp(1:N))/tau^2 - tau^2;
  
  
th_5fold(i,:) = 180/pi*acos(bx(1,1)*bx +  by(1,1)*by + bz(1,1)*bz) ; 

  
end
    
ind = length(tau2);
figure(3)
plotbrowser on
plot(tau2(1:ind)/(2*pi),E_5pi(1:ind),'r','LineWidth',4)
hold on 
plot(tau2(ind:end)/(2*pi),E_5pi(ind:end),'--r','LineWidth',4)
hold on 
plot(tau2(ind)/(2*pi),E_5pi(ind),'.r','LineWidth',4)

title('Bending energy')
set(gca,'FontSize',25,'LineWidth',1)
box on
grid on
