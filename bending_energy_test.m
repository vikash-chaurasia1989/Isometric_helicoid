%--- In this file, we read saved data and compute bending energy for each
%given \tau 
clear all
clc

tau1 = [8.1:.1:9.5 10:1:19];

N = 72;
h = 1/N;

for i = 1:1%length(tau1)
    
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
    
    lm = (lm(1:end-1,1)+lm(2:end,1))/2;
         lm(N+1,1) = lm(1,1);
         
         lm = (lm(1:end-1,1)+lm(2:end,1))/2;
         lm(N+1,1) = lm(1,1);
    
   
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
    
%            bxp(1,1) = 1/h*(-25/12*bx(1,1) + 4*bx(2,1)-3*bx(3,1)+4/3*bx(4,1)-1/4*bx(5,1));
%        byp(1,1) = 1/h*(-25/12*by(1,1) + 4*by(2,1)-3*by(3,1)+4/3*by(4,1)-1/4*by(5,1));
%        bzp(1,1) = 1/h*(-25/12*bz(1,1) + 4*bz(2,1)-3*bz(3,1)+4/3*bz(4,1)-1/4*bz(5,1));
%        
%        bxp(2,1) = 1/h*(1/12*bext(1,1)-2/3*bx(1,1) + 2/3*bx(3,1)-1/12*bx(4,1));
%        byp(2,1) = 1/h*(1/12*bext(2,1)-2/3*by(1,1) + 2/3*by(3,1)-1/12*by(4,1));
%        bzp(2,1) = 1/h*(1/12*bext(3,1)-2/3*bz(1,1) + 2/3*bz(3,1)-1/12*bz(4,1));
%      
%        ig = 3:N-1;
%                 
%        bxp(3:N-1,1) = (-bx(ig+2,1)+8*bx(ig+1,1) -8*bx(ig-1,1)+bx(ig-2,1))/(12*h);
%        byp(3:N-1,1) = (-by(ig+2,1)+8*by(ig+1,1) -8*by(ig-1,1)+by(ig-2,1))/(12*h);
%        bzp(3:N-1,1) = (-bz(ig+2,1)+8*bz(ig+1,1) -8*bz(ig-1,1)+bz(ig-2,1))/(12*h);
%       
%        bxp(N,1) = 1/h*(1/12*bx(N-2,1)-2/3*bx(N-1,1) + 2/3*bx(N+1,1)-1/12*bext(4,1));
%        byp(N,1) = 1/h*(1/12*by(N-2,1)-2/3*by(N-1,1) + 2/3*by(N+1,1)-1/12*bext(5,1));
%        bzp(N,1) = 1/h*(1/12*bz(N-2,1)-2/3*bz(N-1,1) + 2/3*bz(N+1,1)-1/12*bext(6,1));
   
          
   %==============================
    %--  b''  (O(h^2))----
        ig = 2:N;
       
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
  
  
  %=== Validating \lambda = 3/2*\kappa^2 + c 
  
  kappa2 = bpp/tau^2-tau^2;
  
  bbpp = bx.*bx2p + by.*by2p + bz.*bz2p;
  bpbp = bxp.*bxp + byp.*byp + bzp.*bzp;
  
  figure(4)
   %plot(lm- (3/2*kappa2))
 %  plot(bbpp+tau^2);
  
   plot(bx.*bxp + by.*byp + bz.*bzp)
    
  hold on 
  
end
    
% figure(3)
% plotbrowser on
% plot(tau1/(2*pi),E,'LineWidth',4)
% title('Bending energy')
% set(gca,'FontSize',25,'LineWidth',1)
% box on
% grid on
