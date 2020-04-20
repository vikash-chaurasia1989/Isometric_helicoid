clear all
clc

%--- In this file, we read the saved data and plot quantities in frame of
%reference of the tangent vector of the midline 
%---- plot the saved data 

     %-- 3pi --
    %str0 = ['3fold_N72_tau_' num2str(10*tau) '.txt'];
    
    %--- Array initialization 
        N = 144;
    tau = 8.1;
    
    
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
    
    

    h =1/N;% mean(len);

    str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];
    
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
    
    %-- load r --
    
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
    
    
    
    %---- create edge of the binormal scroll using b and r
    wd = .01;
    
    wd = .1;
    v = linspace(-wd,wd,10);
    
    figure(1)
    rx0 = rx + v(1)*bx;
    ry0 = ry + v(1)*by;
    rz0 = rz + v(1)*bz;
    
    plotbrowser on
    title('Initial guess')
    plot3(rx0,ry0,rz0,'LineWidth',4)
    set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
    box on
    grid on
    
  
    
%  
%    
%  
% for i=1:N+1
%     
%     %--
%  a = tx(i)/norm(i)  ;
%  b = ty(i)/norm(i)  ;
%  c = tz(i)/norm(i)  ;
% 
%  
% 
% %-- rotating the mid-plane to align with x-y plane 
% %--- axis of rotation u = (u1, u2, u3)
% %  u = n x k  = (b,-a,0)
% u1 =  b;
% u2 = -a;
% 
% th = acos(c);
% 
% 
% %------ Rotation matrix that will transform each point 
% Ru = [(cos(th) + u1^2*(1-cos(th))) u1*u2*(1-cos(th))                u2*sin(th);
%        u1*u2*(1-cos(th))           (cos(th) + u2^2*(1-cos(th)))    -u1*sin(th);
%        -u2*sin(th)                 u1*sin(th)                       cos(th)   ] ;
%     
%     
%     
% %==================================    
%     temp = Ru*[  rx0(i,1);
%                         ry0(i,1);
%                         rz0(i,1)];
%                     
%      rx1(i,1) = temp(1);
%      ry1(i,1) = temp(2);
%      rz1(i,1) = temp(3);
%      
%      
%      %-- rotating binormal --
%      
%  %==================================    
%     temp = Ru*[  bx(i,1);
%                         by(i,1);
%                         bz(i,1)];
%                     
%      bx1(i,1) = temp(1);
%      by1(i,1) = temp(2);
%      bz1(i,1) = temp(3);    
%      
% end


rx2 = rx;
ry2 = ry;
rz2 = rz;

bx2 = bx;
by2 = by;
bz2 = bz;
   
 
for i=2:N 
    
    %--
% angle between hinge i-1, i, and i+1 

v1x = rx(i,1)-rx(i-1,1);
v1y = ry(i,1)-ry(i-1,1);
v1z = rz(i,1)-rz(i-1,1);

 v2x = rx(i+1,1)-rx(i,1);
v2y = ry(i+1,1)-ry(i,1);
v2z = rz(i+1,1)-rz(i,1);


th = acos((v1x*v2x + v1y*v2y + v1z*v2z)/sqrt((v1x^2 + v1y^2 + v1z^2)*((v2x^2 + v2y^2 + v2z^2)))) ;

th = pi-th;
%-- axis of rotation 
ux = bx(i,1) ;
uy = by(i,1) ;
uz = bz(i,1) ;

R = [(cos(th) + ux^2*(1-cos(th)))         (ux*uy*(1-cos(th))-uz*sin(th))      (ux*uz*(1-cos(th)) + uy*sin(th)) ;
       (ux*uy*(1-cos(th)) +uz*sin(th))    (cos(th) + uy^2*(1-cos(th)))          (uy*uz*(1-cos(th)) - ux*sin(th))  ;
       (ux*uz*(1-cos(th)) - uy*sin(th))    (uy*uz*(1-cos(th)) + ux*sin(th))    (cos(th) + uz^2*(1-cos(th))) ];

%-- rotating the midline ----------------------- 
   
temp = R*[rx(i+1,1) ;
                ry(i+1,1)  ;
                rz(i+1,1) ];
rx2(i+1,1) = temp(1);
ry2(i+1,1) = temp(2);
rz2(i+1,1) = temp(3);


%-- rotating the hinge or the binormal vector ---
temp = R*[bx(i+1,1) ;
                by(i+1,1)  ;
                bz(i+1,1) ];
bx2(i+1,1) = temp(1);
by2(i+1,1) = temp(2);
bz2(i+1,1) = temp(3);


th2(i-1) = 180/pi*acos((bx2(i,1)*bx2(i+1,1) + by2(i,1)*by2(i+1,1) + bz2(i,1)*bz2(i+1,1)));
 
end




figure(2)
plotbrowser on
title('Initial guess')
plot3(rx2,ry2,rz2,'LineWidth',4)
set(gca,'DefaultAxesFontSize',4,'FontSize',25,'LineWidth',2)
box on
grid on
 %---- Rotating ---- 
    