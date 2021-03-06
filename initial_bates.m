%==== initial guess for trifoil knot ====

clear all
clc

N= 100 ;
t1 = linspace(0,2*pi,N+1);
t1 = t1';

m = 2;
n = -3;
al = pi/3;
bt = 1/2*acos(-5/7);

S_al = [cos(al) -sin(al) 0; sin(al) cos(al) 0; 0 0 1];
S_bt = [cos(bt) -sin(bt) 0; sin(bt) cos(bt) 0; 0 0 1];


u = 1/3^.5;
v =1/6^.5;
w = v/u;

C = [u -v  w ; 
     u -v -w ;
     u 2*v 0];




U = [ 0 -u  u ;
      u  0 -u ;
     -u  u  0 ];

id = eye(3,3);
for i = 1:N+1
    
    t = t1(i);
    
    Rnt = [1  0 0 ; 0 cos(n*t) -sin(n*t);0 sin(n*t) cos(n*t)];
        
    Qmt = id + sin(m*t)*U + (1-cos(m*t))*U*U;
    
    
    
    temp = Qmt*C*S_al*Rnt*S_bt;
    
    bx(i,1) = temp(1,1);
    by(i,1) = temp(2,1);
    bz(i,1) = temp(3,1);
    
    
    
end


modb = sqrt(bx.^2+by.^2+bz.^2);

bx = bx./modb;
by = by./modb;
bz = bz./modb;

ig = 1:N;

tx = by(ig,1).*bz(ig+1,1) - bz(ig,1).*by(ig+1,1);
ty = bx(ig,1).*bz(ig+1,1) - bz(ig,1).*bx(ig+1,1);
tz = bx(ig,1).*by(ig+1,1) - by(ig,1).*bx(ig+1,1);


rx(1,1) = 0;
ry(1,1) = 0;
rz(1,1) = 0;


for i = 1:N 
rx(i+1,1) = rx(i,1) + tx(i,1);
ry(i+1,1) = ry(i,1) + ty(i,1);
rz(i+1,1) = rz(i,1) + tz(i,1);
end

 bxg  = [bx(N-1,1) bx(N,1) bx' bx(2,1) bx(3,1)];
 byg  = [by(N-1,1) by(N,1) by' by(2,1) by(3,1)];
 bzg  = [bz(N-1,1) bz(N,1) bz' bz(2,1) bz(3,1)];
% 


%  
% 
% rxg  = [rx(N-1,1) rx(N,1) rx' rx(2,1) rx(3,1)];
% ryg  = [ry(N-1,1) ry(N,1) ry' ry(2,1) ry(3,1)];
% rzg  = [rz(N-1,1) rz(N,1) rz' rz(2,1) rz(3,1)];
% 
% 
%  
h = 1/N;

ig = 3:N+3;




bxp= (bxg(ig+1)-bxg(ig))/(h);
byp= (byg(ig+1)-byg(ig))/(h);
bzp= (bzg(ig+1)-bzg(ig))/(h);


% 
% rx2p = (rxg(ig+1) + rxg(ig-1) -2*rxg(ig))/h^2;
% ry2p = (ryg(ig+1) + ryg(ig-1) -2*ryg(ig))/h^2;
% rz2p = (rzg(ig+1) + rzg(ig-1) -2*rzg(ig))/h^2;
% 
% rx3p = (-1/2*rxg(ig-2) + rxg(ig-1) -rxg(ig+1) +1/2*rxg(ig+2))/h^3;
% ry3p = (-1/2*ryg(ig-2) + ryg(ig-1) -ryg(ig+1) +1/2*ryg(ig+2))/h^3;
% rz3p = (-1/2*rzg(ig-2) + rzg(ig-1) -rzg(ig+1) +1/2*rzg(ig+2))/h^3;
% 
% 
% 
% %  mod2 = (ty'.*rz2p - ry2p.*tz').^2 + (rx2p.*tz' - tx'.*rz2p).^2 + (tx'.*ry2p - rx2p.*ty').^2;
% %  tau2(p1,:) = (rx3p.*(ty'.*rz2p - ry2p.*tz') + ry3p.*(rx2p.*tz' - tx'.*rz2p) + rz3p.*(tx'.*ry2p - rx2p.*ty'))./mod2;
% 
% 
% mod2 = (ryp.*rz2p - ry2p.*rzp).^2 + (rx2p.*rzp - rxp.*rz2p).^2 + (rxp.*ry2p - rx2p.*ryp).^2;
% tau  = (rx3p.*(ryp.*rz2p - ry2p.*rzp) + ry3p.*(rx2p.*rzp - rxp.*rz2p) + rz3p.*(rxp.*ry2p - rx2p.*ryp))./mod2;
% 
% 
%  
% modt = sqrt(rxp.^2 + ryp.^2 + rzp.^2);
% modn = sqrt(rx2p.^2 + ry2p.^2 + rz2p.^2);
% 
% tx = rxp./modt;
% ty = ryp./modt;
% tz = rzp./modt;
% 
% nx = rx2p./modn;
% ny = ry2p./modn;
% nz = rz2p./modn;
% 
% %====  Bates curve 
% 
% 
% 
% 
% 
% 
% bx = ty.*nz - tz.*ny;
% by = tx.*nz - tz.*nx;
% bz = tx.*ny - ty.*nx;
% 




