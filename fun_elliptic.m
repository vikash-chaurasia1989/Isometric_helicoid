function fval =   fun_elliptic( var_in)

global k0 kappa  k or s  h N   A   f_t f_b f_n   tau   f_r tx ty tz nx ny nz bx by bz  rx ry rz   p   B tau1 r nfold


A = var_in(1);
B = var_in(2);

%tau1 = var_in(3);

al1 =  2*A + 2*sqrt(A^2+B);
al2 =  0                  ;
al3 = -2*A +2*sqrt(A^2+B) ;

p = sqrt((al3-al2)/(al3+al1));
%p = 0.99;
q = sqrt(1-al2/al3);
r = 1/2*sqrt(al3+al1);


[K,E] = ellipke(p);

s = linspace(0,nfold*K/r,N+1);

tau = tau1*nfold*K/r;

h = 1/N*nfold*K/r;


[sn,cn,dn] = ellipj(r*s,p);


kappa =  sqrt(al3*(1-q^2*sn.^2));


%=== kappa = +- ==== .. Correcting sign of kappa in suitable intervals

% for m = 1:2:nfold-1
%
%     kappa((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold))) = -kappa((2*m-1)*N/(2*nfold)+1:((2*m-1)*N/(2*nfold) + N/(nfold)));
%
% end
%     m = nfold;
%     kappa((2*m-1)*N/(2*nfold)+1:N+1) = -kappa((2*m-1)*N/(2*nfold)+1:N+1);
%  kappa = kappa';
%

% 
%nfold = 3;
temp = N/nfold:3*N/nfold;
for i = 1:((nfold/2-1)/2)
    
    ind = 4*(i-1)*N/nfold+1 + temp;%((2*i-1)*N/nfold + 1):(((2*i+1)*N/nfold));
    
    kappa(ind) = -kappa(ind);
    
end
% kappa(N/10+1:3*N/10)      = -kappa(N/10+1:3*N/10);
% kappa(5*N/10+1:7*N/10) = -kappa(5*N/10+1:7*N/10);
% kappa(9*N/10+1:end)    = - kappa(9*N/10+1:end);
% 
 ind = ((nfold-1)*N/nfold+1):(N+1);
% 
 kappa(ind) = -kappa(ind);

 
kappa = kappa';

%-- unknot 
% A_unknot = 12.330674378995440; 
% kappa = kappa/max(kappa)*A_unknot;

%kappa(236:702,1) = -kappa(236:702,1);
% %--- evaluating at 2nd point using integration
i = 1;
rx(i+1,1) = rx(i,1) + h*tx(i,1);
ry(i+1,1) = ry(i,1) + h*ty(i,1);
rz(i+1,1) = rz(i,1) + h*tz(i,1);


tx(i+1,1) = tx(i,1) + kappa(i,1)*h*nx(i,1);
ty(i+1,1) = ty(i,1) + kappa(i,1)*h*ny(i,1);
tz(i+1,1) = tz(i,1) + kappa(i,1)*h*nz(i,1);


nx(i+1,1) = nx(i,1) + (-kappa(i,1)*tx(i,1) + tau*bx(i,1))*h ;
ny(i+1,1) = ny(i,1) + (-kappa(i,1)*ty(i,1) + tau*by(i,1))*h ;
nz(i+1,1) = nz(i,1) + (-kappa(i,1)*tz(i,1) + tau*bz(i,1))*h ;

bx(i+1,1) = bx(i,1) - tau*nx(i,1)*h;
by(i+1,1) = by(i,1) - tau*ny(i,1)*h;
bz(i+1,1) = bz(i,1) - tau*nz(i,1)*h;



for i = 2:N
    
    rx(i,1) = rx(i-1,1)+h*tx(i-1,1);
    ry(i,1) = ry(i-1,1)+h*ty(i-1,1);
    rz(i,1) = rz(i-1,1)+h*tz(i-1,1);
    
    
    tx(i+1,1) = tx(i-1,1) + kappa(i,1)*2*h*nx(i,1);
    ty(i+1,1) = ty(i-1,1) + kappa(i,1)*2*h*ny(i,1);
    tz(i+1,1) = tz(i-1,1) + kappa(i,1)*2*h*nz(i,1);
    
    
    nx(i+1,1) = nx(i-1,1) + (-kappa(i,1)*tx(i,1) + tau*bx(i,1))*2*h ;
    ny(i+1,1) = ny(i-1,1) + (-kappa(i,1)*ty(i,1) + tau*by(i,1))*2*h ;
    nz(i+1,1) = nz(i-1,1) + (-kappa(i,1)*tz(i,1) + tau*bz(i,1))*2*h ;
    
    bx(i+1,1) = bx(i-1,1) - tau*nx(i,1)*2*h;
    by(i+1,1) = by(i-1,1) - tau*ny(i,1)*2*h;
    bz(i+1,1) = bz(i-1,1) - tau*nz(i,1)*2*h;
    
 
    
    
end
 
i = N+1;

rx(i,1) = rx(i-1,1)+h*tx(i-1,1);
ry(i,1) = ry(i-1,1)+h*ty(i-1,1);
rz(i,1) = rz(i-1,1)+h*tz(i-1,1);


f_r = [rx(1)-rx(end) ry(1)-ry(end) rz(1)-rz(end)];
f_t = [tx(1)-tx(end) ty(1)-ty(end) tz(1)-tz(end)];
f_n = [nx(1)- or*nx(end) ny(1)- or*ny(end) nz(1)- or*nz(end)];
f_b = [bx(1)- or*bx(end) by(1)- or*by(end) bz(1)- or*bz(end)];

%
% %fval = sqrt(sum(f_r.^2+ f_t.^2 + f_n.^2  ) )
f1 =   sqrt(sum(f_r.^2));
f2 =   sqrt(sum(f_b.^2));
f3 =   sqrt(sum(f_t.^2));
% %f2 =   sqrt(sum(f_b.^2));
%
fval = [f1 f2] ;
%

  end
