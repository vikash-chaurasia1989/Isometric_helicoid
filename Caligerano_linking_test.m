clear all
clc

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
    tau1 = [12 12.5 13 14:.5:25] ;
    
  
    tau1 = [12:.1:13.9 14:.5:25];
    
    tau1 = [11.97871 11.98071 11.98571 11.99071 11.99571  tau1];
       N = 120      ; 
       
       h = 1/N      ;   
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
    
%==========================================================================
%==========================================================================
         
         for p1 = 1:length(tau1)
             
             
         tau = tau1(p1);
     

   % str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];
    
   N = 120;
   h = 1/N;
    str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
 
    if(tau>13.9)
            str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
    end
    
   % str0  = '3fold_N72_tau_170.txt';
    
    %-- 5pi --
    %  str0 = ['5pi_knot_N120_tau_' num2str(10*tau1(p1)) '.txt'];
    
    
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
    
  
    %===== b'=====
    
      
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
    
    
    
    
    
    
    
    
    %========= Generating full binormal by going twice over the arc length
    %===
%     
%       bx(N+2:2*N+1,1) = -bx(2:N+1,1);
%       by(N+2:2*N+1,1) = -by(2:N+1,1);
%       bz(N+2:2*N+1,1) = -bz(2:N+1,1);
% 
%     
%       bxp(N+2:2*N+1,1) = -bxp(2:N+1,1);
%       byp(N+2:2*N+1,1) = -byp(2:N+1,1);
%       bzp(N+2:2*N+1,1) = -bzp(2:N+1,1);
% 
%           
%     
%     N = 2*N;
    %================================================
    
    
    
    %---------------------------------

    %-- load r --
    
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
   
%     plot3(rx,ry,rz)
%     hold on 
    %---- create edge of the binormal scroll using b and r
    wd = .005;
    
 
    v = linspace(-wd,wd,10);
    
 %--- edge   ----
 
    
    rx0 = rx - v(1)*bx;
    ry0 = ry - v(1)*by;
    rz0 = rz - v(1)*bz;
    
    tx0 = tx - v(1)*bxp;
    ty0 = ty - v(1)*byp;
    tz0 = tz - v(1)*bzp;
   
 
 
%==========================================================================

 
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
Lk_self(p1)  =  2*((tau/(2*pi) + Wr ))  ;
    
%==========================================================================

%--------edge linking ---
 
 
 
for i = 1:N     
     % 
     j = 1:N+1     ;             % j array index
   %  j = j(j~=i);   % indexes not containing i
      
         %--------------------- 1. Curve r(s)  -------------------------------
         
      rij  =  ((rx0(j,1)-rx(i,1)).^2 + (ry0(j,1)-ry(i,1)).^2 +(rz0(j,1)-rz(i,1)).^2).^.5 ;
       
      %-- for self linking number
      tjix = ty0(j,1).*tz(i,1) - tz0(j,1).*ty(i,1);
      tjiy = tz0(j,1).*tx(i,1) - tx0(j,1).*tz(i,1);
      tjiz = tx0(j,1).*ty(i,1) - ty0(j,1).*tx(i,1);
         
      
      temp =   (tjix.*(rx0(j,1)-rx(i,1)) + tjiy.*(ry0(j,1)-ry(i,1)) + tjiz.*(rz0(j,1)-rz(i,1)))./rij.^3 ;
        
      ker(i,1) = trapz(temp);
      
      
      
end

Wr  =   h^2*trapz(ker)/(4*pi);
Lk_edge(p1)  =   ((tau/(2*pi) + Wr ) ) ;

%==========================================================================
%===========  Linking number using the solid angle formulation ============

%-- creating array for linking number subroutine --

for i = 1:N+1
    
    gamma0(1,i) = rx0(i,1);
    gamma0(2,i) = ry0(i,1);
    gamma0(3,i) = rz0(i,1);
    
    gamma1(1,i) = rx(i,1);
    gamma1(2,i) = ry(i,1);
    gamma1(3,i) = rz(i,1);
end
n= 0;

 
 % gamma0 = gamma1;

length_0 = size(gamma0, 2);
length_1 = size(gamma1, 2);

% gamma0 = intval(gamma0);
% gamma1 = intval(gamma1);

for i = 1:length_0
    for j = 1:length_1
        a = gamma1(:, j)                    + 0     - gamma0(:, i);
        b = gamma1(:, j)                    - gamma0(:, mod(i, length_0) + 1);
        c = gamma1(:, mod(j, length_1) + 1) - gamma0(:, mod(i, length_0) + 1);
        d = gamma1(:, mod(j, length_1) + 1) + 0     - gamma0(:, i);
        n = n + solid_angle_quadrilateral(a, b, c, d);
        
       
    end
    
              
        
end

 
%===================================================
%n = n / (4 * midrad(pi, 1e-14))

n = n/(4*pi) ;

Lk_solid(p1) =  n*2 ;

         end
         
 
         figure(1)
         plotbrowser on
         plot(tau1/(2*pi),Lk_solid,'--o','LineWidth',1)
         hold on 
         plot(tau1/(2*pi),Lk_edge,'--o','LineWidth',1)

         title('solid angle')
         set(gca,'FontSize',25,'LineWidth',1)
         box on
         grid on
         hold on 
         
         
         
         figure(2)
         plotbrowser on
         plot(tau1/(2*pi),Lk_edge,'--o','LineWidth',1)
         title('edge linking')
         set(gca,'FontSize',25,'LineWidth',1)
         box on
         grid on
         hold on 
         
         
         
         figure(3)
         plotbrowser on
         plot(tau1/(2*pi),Lk_self,'--o','LineWidth',1)
         title('self linking')
         set(gca,'FontSize',25,'LineWidth',1)
         box on
         grid on
         hold on 

 
         
         
         
 