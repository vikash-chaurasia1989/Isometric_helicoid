clear all
clc
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi')
  addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi')


format longE

%--- In this file, we read the saved data and compute self linking and edge
%linking number using the Caligariano theorem 

%Lk = 2*(tau/(2*pi) + Wr), where Wr is writhe of the curve %
%---- plot the saved data 

     %-- 3pi --
    %str0 = ['3fold_N72_tau_' num2str(10*tau) '.txt'];
    
    %--- Array initialization 
    %    N = 144;
    %   tau = 8.1;
    
    %===== 3 fold  =====
   %  tau1 = [8.1:.1:9.5 10:1:19];
   %     N = 72;
    %===== 5 fold =====   
  %===== tau for 3 fold ====
 
tau_3 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
             13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725 13.216865473  13.216865474355 ] ;
         
tau_5 = [13.216865474 13.2168655 13.21686575 13.216866 13.21687 13.2169 13.217 13.2175 13.22 13.23];
tau_7 = [13.3:.1:13.9 14.5:.5:21];

%tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.1:.1:13.9 14:.5:21];
tau_1 = [tau_3 tau_5 tau_7];
N = 72;
h = 1/N;

N1 = N-1;
 
         
         for p1 =  1:1%length(tau1)
             
         
 %========== 3pi ==================  
       
 N = 72;
 h = 1/N;
 
 N1 = N-1;
 
 
 tau =  13.21686547250001;%tau1(p1);
 
% tau=13.22;%

 
 str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
 
 
%==========================================================================
%==========================================================================


%===== 5pi ========================
%  N = 120;
%  h = 1/N;
%  
%  N1 = N-1;
 
 %     str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
%  
%     if(tau>13.9)
%             str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
%     end
    
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
 
    
    strb = ['b_' str0];
    strr = ['r_' str0];
    
    strlm = ['lm_' str0];
    strrho = ['rho_' str0];
    struvw = ['uvw_' str0];
    strbext = ['bext_' str0];
    
    
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
   bext = load(strbext);
   
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
Lk_s(p1)  =  (2*(tau/(2*pi) + Wr ))  ;

%==========================================================================


    
%==========================================================================

%-------- self linking number ---
 
 
%  
% for i = 1:N     
%      % 
%      j = 1:N+1     ;             % j array index
%    %  j = j(j~=i);   % indexes not containing i
%       
%          %--------------------- 1. Curve r(s)  -------------------------------
%          
%       rij  =  ((rx0(j,1)-rx(i,1)).^2 + (ry0(j,1)-ry(i,1)).^2 +(rz0(j,1)-rz(i,1)).^2).^.5 ;
%        
%       %-- for self linking number
%       tjix = ty0(j,1).*tz(i,1) - tz0(j,1).*ty(i,1);
%       tjiy = tz0(j,1).*tx(i,1) - tx0(j,1).*tz(i,1);
%       tjiz = tx0(j,1).*ty(i,1) - ty0(j,1).*tx(i,1);
%          
%       
%       temp =   (tjix.*(rx0(j,1)-rx(i,1)) + tjiy.*(ry0(j,1)-ry(i,1)) + tjiz.*(rz0(j,1)-rz(i,1)))./rij.^3 ;
%         
%       ker(i,1) = trapz(temp);
%       
%       
%       
% end
% 
% Wr  =   h^2*trapz(ker)/(4*pi);
% Lk_edge(p1)  = 2*(tau/(2*pi) + Wr )  ;

%==========================================================================

         end
         
%          figure(1)
%          plotbrowser on
%          plot(tau1/(2*pi),round_odd(Lk_s),'--o','LineWidth',1)
%          title('self linking')
%          set(gca,'FontSize',25,'LineWidth',1)
%          box on
%          grid on
%   

Lk_s