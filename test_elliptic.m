%=== we test how elliptic functions behave in this file ===
clear all
clc


N = 240;

 

p =  0.2981; 
[K,E] = ellipke(p);
 s = linspace(0,1,N+1);

 m = 5;
 

[sn,cn,dn] = ellipj(2*m*s*K,p);



plot(s,sn)
hold on 
plot(s,sin(s*m*pi))


 
% 
% p = x;
% 
% %fun_frenet(p);
% %plot3(rx,ry,rz)
% 
% sig = 2*tau;
% 
% 
% 
% 
% x1 = k0*s/(2*p);
% dx1 = k0/(2*p)*1/h;    % dx = k0/(2*p)*ds 
% 
% [sn,cn,dn] = ellipj(x1,p);
% 
% mu = 1/4*sqrt((p^(-2)-sig^2)^2 + 4*sig^2);
% m = 16*mu^2/(p^(-2)+sig^2)^2;
% lm1 = (sig^2 -p^(-2)+2)/4;
% 
%  %r = (1/mu)*sqrt(1/m-sn.^2);
% 
% L = sig*p*(p^(-2)-sig^2-2)/8/mu;
% M = 16*mu^2/(p^(-2)+sig^2)^2;
% N2 = 2*mu*p*sig/(p^(-2)+sig^2) + L ;
% 
% 
% 
%   dth = -p*sig/mu*(lm1 + (p^(-2)-1)*(p^(-2) -sig^2)/2/(p^(-2)+sig^2)./(1-m*sn.^2));
% % 
% % h = 1/N;
%  th(1) = 0;
% for i =1:N
%     
%     th(i+1) = th(i) +dx1*dth(i);
%     
% end    
% 
% 
% 
% 
% dth = L - N2./(1-M*sn.^2);
% dz  = (p^(-2) -sig^2 -2*sn.^2)/4/mu;
% 
%  r = sqrt(dn.^2 - sig^2/(4*mu^2))./mu;
% 
% th(1) = 0; z(1) = 0;
% 
% for i = 1:N
%     
%     th(i+1) = th(i) + h/(2*p)*dth(i);
%     z(i+1) = z(i) + h*dz(i);
%     
% end
%     
% x = r.*cos(th);
% y = r.*sin(th);
% 
% 
% plot(x,y)
% 
 
%plot(sin(pi/2+s/max(s)*3*pi)*6.112)
%  al1 = A +sqrt(A^2+B^2);
% al2 = 0;
% al3 = -A +sqrt(A^2+B^2);
%  
% p = sqrt((al3-al2)/(al3+al1));
% q = sqrt(1-al2/al3);
% r = 1/2*sqrt(al3+al1);
% 
% 
% [K,E] = ellipke(p);
% %s = linspace(0,6*K/r,N+1);
% s = linspace(0,1,N+1);
%  
% %tau = tau1;
% 
% [sn,cn,dn] = ellipj(r*s*6*K/r,p);
% 
% 
% kappa =  sqrt(al3*(1-q^2*sn.^2));

