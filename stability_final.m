clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold id tau fac

format longE

col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;


N = 105;
N1 = N-1;
h = 1/N;


sig = 0.01;

  
tau = 16;




branch =2;

path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
strtau = ['tau_branch_' num2str(branch) '.txt'];
tau1 = load([path strtau]);
%==========================================================================
 
st = zeros(1,length(tau1));
%  st = zeros(1,p2);

tau0 = 16;
[temp,p1] =  min(abs(tau1-tau0)) ;   % find index closest to tau = tau0;

tau = tau1(p1);
n = tau/2/pi;

fac =  (asinh(n*pi*sig))/(4*pi^3*n^3);




str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;

ht = 0.002;
p1 = 2;
 str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_h_' num2str(round(10^5*ht)) '_step_' num2str(p1-1) '.txt'];
    
 
 
f= load([path str0]);



bx(2:N,1)     = f(1:3:3*N1-2,1)    ;
by(2:N,1)     = f(2:3:3*N1-1,1)    ;
bz(2:N,1)     = f(3:3:3*N1,1)      ;
%--- boundary points --

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;



id = [1 0 0;0 1 0;0 0 1];


%============ generating tangent and midline =========
i = 1:N;


tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);


tx(N+1,1) = tx(1,1);
ty(N+1,1) = ty(1,1);
tz(N+1,1) = tz(1,1);

tx = tx/h;
ty = ty/h;
tz = tz/h;

%--------------------

% %--- position vector using integration of tangent
% % initialization

rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;


for i = 1:N
    
    rx(i+1,1) = rx(i,1) + h*tx(i,1);
    ry(i+1,1) = ry(i,1) + h*ty(i,1);
    rz(i+1,1) = rz(i,1) + h*tz(i,1);
    
end


[fval,jac] = fun_jacobian6(f);

num_eq = 3*N1   ;    % number of equilibrium equations
num_C  = 2*N1+4   ;   % number of constraints
%
A = (jac(num_eq+1:end,1:num_eq))';   % num_C x num_eq constraint gradient matrix

%--- Qr decomposition of the constraint matrix to obtain null space

[Q R] = qr(A)                ;
Z     = Q(:,num_C+1:num_eq)  ;

%==== checking the true dimension of the null space ======
count = 0;
tol = 10^-10;

[U,S,V2] = svd(A);


%-- finding rank of A

s = svd(A);
rank = length(s) - length(find(s<tol));

Z2  = U(:,rank+1:end);

for i = 1:num_eq
    
    temp(i) = sqrt(sum((Q(:,i)'*A).^2));
    
    
    if(temp(i)<tol)
        count = count + 1;
        
        Z1(:,count) = Q(:,i);
    end
end




%   %---- Projected hessian ----
%
%Ac = Z'*jac(1:num_eq,1:num_eq)*Z     ;
Ac = Z'*jac(1:num_eq,1:num_eq)*Z     ;

%
%   %--- eigsen values and eigsen vectors ----
%
[V,D] = eigs(Ac,5,'smallestabs') ;
%
[B, I ] = sort(diag(D),'ascend');

B

%-- perturbed shape with respect to negative eigen value ---

i11 = 1;

x1 = V(:,I(i11));

x  = Z*x1;

%
dbx(1,1) = 0;
dby(1,1) = 0;
dbz(1,1) = 0;




dbx(2:N,1)   = x(1:3:3*N1-2,1)    ;
dby(2:N,1)   = x(2:3:3*N1-1,1)    ;
dbz(2:N,1)   = x(3:3:3*N1,1)     ;

dbx(N+1,1) = 0;
dby(N+1,1) = 0;
dbz(N+1,1) = 0;

ep = .5;

bx1 = bx + ep*dbx;
by1 = by + ep*dby;
bz1 = bz + ep*dbz;

%==== saving the perturbed binormal shapes ===

% path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Constant torsion/Full_discrete_formulation/orientable/' ;
%
% strb = ['b_eig' num2str(i11) '.txt'];
%
% fileID = fopen(strb,'w');
% fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx1'    ;by1'    ; bz1'    ] );
% fclose(fileID);


ds1 = sqrt((bx1(1:end-1,1)-bx1(2:end,1)).^2 + (by1(1:end-1,1)-by1(2:end,1)).^2 + (bz1(1:end-1,1)-bz1(2:end,1)).^2);

l2 = sum(ds1);

figure(1)
plotbrowser on
title('binormal')
plot_sphere()
hold on
plot3(bx,by,bz,'color',col{1},'LineWidth',2)
hold on
plot3(bx1,by1,bz1,'r','LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d
%  plot3(bx(1,1),by(1,1),bz(1,1),'ok')


%--- -perturbed midline

%----------------- Post processing -----
h =1;% mean(len);
%---- Tangent ti = bi \times bi+1

i = 1:N;


tx(i,1) = by1(i,1).*bz1(i+1,1) - bz1(i,1).*by1(i+1,1);
ty(i,1) = bz1(i,1).*bx1(i+1,1) - bx1(i,1).*bz1(i+1,1);
tz(i,1) = bx1(i,1).*by1(i+1,1) - by1(i,1).*bx1(i+1,1);


tx(N+1,1) = tx(1,1);
ty(N+1,1) = ty(1,1);
tz(N+1,1) = tz(1,1);

%--------------------

% %--- position vector using integration of tangent
% % initialization

rx1(1,1) = 0; ry1(1,1) = 0; rz1(1,1) = 0;



for i = 1:N
    
    rx1(i+1,1) = rx1(i,1) + h*tx(i,1);
    ry1(i+1,1) = ry1(i,1) + h*ty(i,1);
    rz1(i+1,1) = rz1(i,1) + h*tz(i,1);
    
end

figure(2)
plotbrowser on
title('midline')
plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
hold on
plot3(rx1,ry1,rz1,'r','LineWidth',1)
hold on
plot3(rx(1,1),ry(1,1),rz(1,1),'ok','LineWidth',2)
% daspect([1,1,1]);
axis off
set(gcf,'color',colw);
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d
% %----------------------------------------------------
% %           N        Lk        L                          alpha
% %
% %          240       4         1.167879653094892e+01     1.167994912472664e+01
% %          240       6         1.618789163145617e+01     1.619096178086449e+01
% %          240       8         2.023333699809010e+01     2.023933375306768e+01
% %          240       8         2.396846060643186e+01     2.397843242686577e+01
% %          480       8         2.395742609593649e+01     2.395991351234742e+01
% %          240      12fold     2.747534830559159e+01     2.749037408464824e+01
% 
% %=== curvature calculation ===
% 
% %=== base curve ====
% 
% rxg = [rx(N,1) rx'  rx(2,1)]' ;
% ryg = [ry(N,1)  ry'  ry(2,1)]' ;
% rzg = [rz(N,1)  rz'  rz(2,1)]' ;
% 
% rxp = rxg(3:end,1) - rxg(2:end-1,1) ;
% ryp = ryg(3:end,1) - ryg(2:end-1,1) ;
% rzp = rzg(3:end,1) - rzg(2:end-1,1) ;
% 
% 
% rxpp = rxg(3:end,1) + rxg(1:end-2,1)-2*rxg(2:end-1,1);
% rypp = ryg(3:end,1) + ryg(1:end-2,1)-2*ryg(2:end-1,1);
% rzpp = rzg(3:end,1) + rzg(1:end-2,1)-2*rzg(2:end-1,1);
% 
% k20  = rxpp.^2 + rypp.^2 + rzpp.^2;
% 
% %==============
% 
% rx1g = [rx1(N,1) rx1'  rx1(2,1)]' ;
% ry1g = [ry1(N,1)  ry1'  ry1(2,1)]' ;
% rz1g = [rz1(N,1)  rz1'  rz1(2,1)]' ;
% 
% rx1p = rx1g(3:end,1) - rx1g(2:end-1,1) ;
% ry1p = ry1g(3:end,1) - ry1g(2:end-1,1) ;
% rz1p = rz1g(3:end,1) - rz1g(2:end-1,1) ;
% 
% 
% rx1pp = rx1g(3:end,1) + rx1g(1:end-2,1)-2*rx1g(2:end-1,1);
% ry1pp = ry1g(3:end,1) + ry1g(1:end-2,1)-2*ry1g(2:end-1,1);
% rz1pp = rz1g(3:end,1) + rz1g(1:end-2,1)-2*rz1g(2:end-1,1);
% 
% k21  = rx1pp.^2 + ry1pp.^2 + rz1pp.^2;
% 
% 
% %=== curvature of the binormal vector ===
% 
% bxg = [-bx(N,1) bx'  -bx(2,1)]' ;
% byg = [-by(N,1)  by'  -by(2,1)]' ;
% bzg = [-bz(N,1)  bz'  -bz(2,1)]' ;
% 
% 
% bxp = bxg(3:end,1) - bxg(2:end-1,1) ;
% byp = byg(3:end,1) - byg(2:end-1,1) ;
% bzp = bzg(3:end,1) - bzg(2:end-1,1) ;
% 
% bxpp = bxg(3:end,1) + bxg(1:end-2,1)-2*bxg(2:end-1,1);
% bypp = byg(3:end,1) + byg(1:end-2,1)-2*byg(2:end-1,1);
% bzpp = bzg(3:end,1) + bzg(1:end-2,1)-2*bzg(2:end-1,1);
% 
% k2b  = bxpp.^2 + bypp.^2 + bzpp.^2;
% 
% l = tau;
% temp = (k2b/h^4)/l^2-l^2;
% 
% 
% %==== checking if constraints are satisfied ===
% 
% %=== unit b ===
% % . b\cdot \delta b = 0
% 
% du = bx.*dbx + by.*dby + bz.*dbz;
% 
% %=== unit speed ===
% % b' \cdot \delta b' = 0
% 
% dbxp = dbx(2:N+1,1)-dbx(1:N);
% dbyp = dby(2:N+1,1)-dby(1:N);
% dbzp = dbz(2:N+1,1)-dbz(1:N);
% 
% 
% dup = bxp(1:N,1).*dbxp + byp(1:N,1).*dbyp + bzp(1:N,1).*dbzp ;
% 
% 
% %=== closure constraint ===
% 
% b_dpx = by(1:N,1).*dbzp - bz(1:N,1).*dbyp;
% b_dpy = bz(1:N,1).*dbxp - bx(1:N,1).*dbzp;
% b_dpz = bx(1:N,1).*dbyp - by(1:N,1).*dbxp;
% 
% db_bxp = dby(1:N,1).*bzp(1:N,1) - dbz(1:N,1).*byp(1:N,1);
% db_byp = dbz(1:N,1).*bxp(1:N,1) - dbx(1:N,1).*bzp(1:N,1);
% db_bzp = dbx(1:N,1).*byp(1:N,1) - dby(1:N,1).*bxp(1:N,1);
% 
% 
% closure = u*(b_dpx + db_bxp) + v*(b_dpy + db_byp) + w*(b_dpz + db_bzp) ;
% 
% 
% clx = (dby(1:N,1).*bz(2:N+1,1) -  dbz(1:N,1).*by(2:N+1,1)) + (by(1:N,1).*dbz(2:N+1,1) -  bz(1:N,1).*dby(2:N+1,1)) ;
% cly = (dbz(1:N,1).*bx(2:N+1,1) -  dbx(1:N,1).*bz(2:N+1,1)) + (bz(1:N,1).*dbx(2:N+1,1) -  bx(1:N,1).*dbz(2:N+1,1)) ;
% clz = (dbx(1:N,1).*by(2:N+1,1) -  dby(1:N,1).*bx(2:N+1,1)) + (bx(1:N,1).*dby(2:N+1,1) -  by(1:N,1).*dbx(2:N+1,1)) ;
% 
% 
% clx2 = (dby(1:N,1).*dbz(2:N+1,1) -  dbz(1:N,1).*dby(2:N+1,1)) + (dby(1:N,1).*dbz(2:N+1,1) -  dbz(1:N,1).*dby(2:N+1,1)) ;
% cly2 = (dbz(1:N,1).*dbx(2:N+1,1) -  dbx(1:N,1).*dbz(2:N+1,1)) + (dbz(1:N,1).*dbx(2:N+1,1) -  dbx(1:N,1).*dbz(2:N+1,1)) ;
% clz2 = (dbx(1:N,1).*dby(2:N+1,1) -  dby(1:N,1).*dbx(2:N+1,1)) + (dbx(1:N,1).*dby(2:N+1,1) -  dby(1:N,1).*dbx(2:N+1,1)) ;

