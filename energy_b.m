function E = energy_b(bx,by,bz)

%== energy given binormal and torsion ==
global sig 

% bx = b(:,1);
% by = b(:,2);
% bz = b(:,3);

N = length(bx)-1;
h = 1/N;

tau  = sum(sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2+(bz(1:end-1,1)-bz(2:end,1)).^2));


 
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
  
  
  bpp = bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p;

  %-- curvature 
  n = tau/2/pi;
  
  % aspect ratio  
  %sig = 0.01;
  
  E  = (asinh(n*pi*sig))/(4*pi^3*n^3)*(h*sum(bpp(1:N))/2 - 8*pi^4*n^4) ;
    




end