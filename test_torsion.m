%-- test torsion 
clear all
clc


temp = load('b_3fold_N36_tau_809.txt');

bx = temp(:,1);
by = temp(:,2);
bz = temp(:,3);



N = length(bx)-1;

h = 1/N;

bext = load('bext_3fold_N36_tau_809.txt');

bx1 = bext(1);
by1 = bext(2);
bz1 = bext(3);

bxN2 = bext(4);
byN2 = bext(5);
bzN2 = bext(6);


%-- augmented matrix ---

  ig = 2:N;
   
   
   %--  b'   (O(h^2))----
   %-- boundary points ---
   %- forward Euler --
     bxp(1,1) = (bx(2,1)-bext(1,1))/(2*h);
     byp(1,1) = (by(2,1)-bext(2,1))/(2*h);
     bzp(1,1) = (bz(2,1)-bext(3,1))/(2*h);
   
%      bxp(1,1) = (bx(2,1)-bx(1,1))/(h);
%      byp(1,1) = (by(2,1)-by(1,1))/(h);
%      bzp(1,1) = (bz(2,1)-bz(1,1))/(h);
   
   
%    bxp(1,1) = 1/h*(-25/12*bxg(2,1) + 4*bxg(3,1) -3*bxg(4,1)+4/3*bxg(5,1)-1/4*bxg(6,1));
%    byp(1,1) = 1/h*(-25/12*byg(2,1) + 4*byg(3,1) -3*byg(4,1)+4/3*byg(5,1)-1/4*byg(6,1));
%    bzp(1,1) = 1/h*(-25/12*bzg(2,1) + 4*bzg(3,1) -3*bzg(4,1)+4/3*bzg(5,1)-1/4*bzg(6,1));
%    %-- Central difference
   %---O(h^2)---
     bxp(ig,1) = (bx(ig+1,1)-bx(ig-1,1))/(2*h);
     byp(ig,1) = (by(ig+1,1)-by(ig-1,1))/(2*h);
     bzp(ig,1) = (bz(ig+1,1)-bz(ig-1,1))/(2*h);
     
      %--- Backward Euler --
    bxp(N+1,1) = (bext(4,1)-bx(N,1))/(2*h);
    byp(N+1,1) = (bext(5,1)-by(N,1))/(2*h);
    bzp(N+1,1) = (bext(6,1)-bz(N,1))/(2*h);
    
    
    
    %--- recovering midline ---
    
    %---- Tangent ti = bi \times bi+1

i = 1:N+1;


% tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
% ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
% tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
%

tau = 8.09;

tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);

tx = tx/tau;
ty = ty/tau;
tz = tz/tau;

%
%   tx(N+1,1) = tx(1,1);
%   ty(N+1,1) = ty(1,1);
%   tz(N+1,1) = tz(1,1);
% %
% %
% tx = tx/tau;
% ty = ty/tau;
% tz = tz/tau;

%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1) = 0; ry(1) = 0; rz(1) = 0;

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

     rx(i+1) = rx(i) + h*tx(i);
     ry(i+1) = ry(i) + h*ty(i);
     rz(i+1) = rz(i) + h*tz(i);

end


%========= Load r from data =======


temp = load('r_3fold_N36_tau_82.txt');

rx = temp(:,1);
ry = temp(:,2);
rz = temp(:,3);

rx = rx';
ry = ry';
rz = rz';
%================

rxg = [rx(N-1) rx(N) rx  rx(2) rx(3)] ;
ryg = [ry(N-1)  ry(N)  ry  ry(2) ry(3)] ;
rzg = [rz(N-1)  rz(N) rz  rz(2) rz(3)] ;



rxp =  (rxg(4:N+4) - rxg(2:N+2))/(2*h);
ryp =  (ryg(4:N+4) - ryg(2:N+2))/(2*h);
rzp =  (rzg(4:N+4) - rzg(2:N+2))/(2*h);


rx2p = (rxg(4:N+4) + rxg(2:N+2)-2*rxg(3:N+3))/h^2;
ry2p = (ryg(4:N+4) + ryg(2:N+2)-2*ryg(3:N+3))/h^2;
rz2p = (rzg(4:N+4) + rzg(2:N+2)-2*rzg(3:N+3))/h^2;


rx3p = (rxg(5:N+5) -2*rxg(4:N+4) + 2*rxg(2:N+2) - rxg(1:N+1))/(2*h^3);
ry3p = (ryg(5:N+5) -2*ryg(4:N+4) + 2*ryg(2:N+2) - ryg(1:N+1))/(2*h^3);
rz3p = (rzg(5:N+5) -2*rzg(4:N+4) + 2*rzg(2:N+2) - rzg(1:N+1))/(2*h^3);



nx  =  ryp.*rz2p - rzp.*ry2p ;
ny  =  rzp.*rx2p - rxp.*rz2p ;
nz  =  rxp.*ry2p - ryp.*rx2p ;



%------

tau2 = (nx.*rx3p + ny.*ry3p + nz.*rz3p)./(nx.^2+ny.^2+nx.^2);

    