%==== initial guess for trifoil knot ====

clear all
clc

N= 100;
t = linspace(0,2*pi,N+1);
t = t';


rx = sin(t) + 2*sin(2*t);
ry = cos(t) - 2*cos(2*t);
rz = -sin(3*t);


rxg  = [rx(N-1,1) rx(N,1) rx' rx(2,1) rx(3,1)];
ryg  = [ry(N-1,1) ry(N,1) ry' ry(2,1) ry(3,1)];
rzg  = [rz(N-1,1) rz(N,1) rz' rz(2,1) rz(3,1)];


 
h = 1/N;

ig = 3:N+3;




rxp= (rxg(ig+1)-rxg(ig-1))/(2*h);
ryp= (ryg(ig+1)-ryg(ig-1))/(2*h);
rzp= (rzg(ig+1)-rzg(ig-1))/(2*h);



rx2p = (rxg(ig+1) + rxg(ig-1) -2*rxg(ig))/h^2;
ry2p = (ryg(ig+1) + ryg(ig-1) -2*ryg(ig))/h^2;
rz2p = (rzg(ig+1) + rzg(ig-1) -2*rzg(ig))/h^2;

rx3p = (-1/2*rxg(ig-2) + rxg(ig-1) -rxg(ig+1) +1/2*rxg(ig+2))/h^3;
ry3p = (-1/2*ryg(ig-2) + ryg(ig-1) -ryg(ig+1) +1/2*ryg(ig+2))/h^3;
rz3p = (-1/2*rzg(ig-2) + rzg(ig-1) -rzg(ig+1) +1/2*rzg(ig+2))/h^3;



%  mod2 = (ty'.*rz2p - ry2p.*tz').^2 + (rx2p.*tz' - tx'.*rz2p).^2 + (tx'.*ry2p - rx2p.*ty').^2;
%  tau2(p1,:) = (rx3p.*(ty'.*rz2p - ry2p.*tz') + ry3p.*(rx2p.*tz' - tx'.*rz2p) + rz3p.*(tx'.*ry2p - rx2p.*ty'))./mod2;


mod2 = (ryp.*rz2p - ry2p.*rzp).^2 + (rx2p.*rzp - rxp.*rz2p).^2 + (rxp.*ry2p - rx2p.*ryp).^2;
tau  = (rx3p.*(ryp.*rz2p - ry2p.*rzp) + ry3p.*(rx2p.*rzp - rxp.*rz2p) + rz3p.*(rxp.*ry2p - rx2p.*ryp))./mod2;


 
modt = sqrt(rxp.^2 + ryp.^2 + rzp.^2);
modn = sqrt(rx2p.^2 + ry2p.^2 + rz2p.^2);

tx = rxp./modt;
ty = ryp./modt;
tz = rzp./modt;

nx = rx2p./modn;
ny = ry2p./modn;
nz = rz2p./modn;

%====  Bates curve 






bx = ty.*nz - tz.*ny;
by = tx.*nz - tz.*nx;
bz = tx.*ny - ty.*nx;





