clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold id tau sig

format longE

  col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
  col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
  col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
  col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
  col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
  col6 = [1, 0.929, 0.278];
  col = {col1,col2,col3,col4,col5,col6} ;
  colw = [41,44,52]/255;

branch = 2;

if(branch==1)
     %      tau_3 = [8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
     %          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
     %      tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
     %          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
     %
     %      tau_7 = [13.224  13.225 13.23 13.3:.1:13.9 14.5:.5:21 21.5:.5:45];
     %
     %
     %
     %      tau1 =      [tau_3 tau_5 tau_7];
     
   
%         
%           tau_3 = [ 8.0957344694550208E+00  8.09574 8.095745 8.09575 8.0958  8.09585 8.096 8.097 8.09975 8.1:.05:8.5 8.6:.1:9.5 9.6:.1:11.9 12:.1:12.5 12.55 12.552 12.554 12.555 12.556 12.56 12.57 12.6:.1:13.2 13.21 13.215 ...
%          13.2165 13.21675 13.216865 13.21686525 13.2168653 13.2168654 13.21686545 13.21686546 13.21686547 13.2168654725    ] ;
%        tau_5 = [13.21686547250001  13.21686547255 13.2168654726 13.21686547275 13.216865473 13.216865474355 13.216865474 13.2168655 13.21686575 13.216866...
%          13.21687 13.2169 13.217 13.2175 13.22 13.2225 13.2235 13.2236 13.223675 13.22375 13.2238];
%      
%        tau_7 = [13.224  13.225 13.23 13.3 13.31 13.32  13.33 13.3 13.34 13.35  13.36 13.37 13.38 13.39 13.395  13.4:.1:13.9 14.5:.5:21 21.5:.5:33.5 33.6:.1:33.9 33.95 33.96 33.965 33.966 33.9665 33.96675 33.96676 33.96677 ];
% %      
%     
%      tau1 =      [tau_3 tau_5 tau_7];
%      
      


     N = 72;
     N1 = N-1;
     h = 1/N;
     nfold=3;
     
    str1 =  '3fold_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_4/';
    
    strtau = 'tau_branch_1.txt';

    tau1 = load([path strtau]);

 elseif(branch==2)
     %================ Branch 2 =====================
     
     %tau1 = [ 11.984235 11.98425  11.9845 11.985 12:.1:25 25.5:.5:34 34.1 34.15 34.155 34.157  34.158 34.159];
     
     
     N = 120;
     N1 = N-1;
     h = 1/N;
     nfold=5;
     
    str1 = '5pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_4/';
    
     strtau = 'tau_branch_2.txt';

    tau1 = load([path strtau]);
    
     % =============== Branch 3 =====================
     
 else
     tau1 = [15.45:.1:25 25.1:.1:40];
     
     N = 84;
     N1 = N-1;
     h = 1/N;
     
     nfold=7;
    str1 = '7pi_knot_N' ;
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
 
 end


p1 = 7;
tau = tau1(p1);

%=== Unknot solutions ===
 str0 = [str1  num2str(N) '_tau_' num2str((10^10*tau)) '.txt'];
%str0 = ['3fold_N' num2str(N) '_tau_' num2str(10^10*tau) '_sig_' num2str(sig) '.txt'];

if(nfold==7)
    str0 = [str1  num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
end

% str0 = '5pi_knot_N120_tau_1197871.txt';
% tau = 11.97871;

strb =  [path 'b_' str0] ;
strlm = [path 'lm_' str0];
strrho = [path 'rho_' str0];
struvw = [path 'uvw_' str0];
strbext = [path 'bext_' str0];

%==== load b ===========
temp = load(strb);

% N = N/2;
bx = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);



b(1:3:3*N-2) = bx(1:N,1)  ;
b(2:3:3*N-1) = by(1:N,1)  ;
b(3:3:3*N)   = bz(1:N,1)  ;
%====================


temp =     load(struvw);
%
u  =    temp(1,1);
v  =    temp(2,1);
w =     temp(3,1);
 % %===============================

%=== load rho ==========
rho =      load(strrho);
%=== load lm  ==========
lm =   load(strlm);

f = [b(1:3*N)  rho' lm' u v w  ]'  ;         % for fun_jacobian6

 

%--- boundary points --

 
% bz(2,1) = 0;
% bz(N,1) = 0;

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


[fval,jac] = fun_jacobian8(f);

num_eq = 3*N     ;    % number of equilibrium equations
num_C  = 2*N+3   ;   % number of constraints
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
Ac = Z1'*jac(1:num_eq,1:num_eq)*Z1     ;

%
%   %--- eigsen values and eigsen vectors ----
%
[V,D] = eigs(Ac,15,'smallestabs') ;
%
[B, I ] = sort(diag(D),'ascend');

B

%-- perturbed shape with respect to negative eigen value ---

i11 = 1;

x1 = V(:,I(i11));

x  = Z*x1;

%
 
 

dbx(1:N,1)   = x(1:3:3*N-2,1)    ;
dby(1:N,1)   = x(2:3:3*N-1,1)    ;
dbz(1:N,1)   = x(3:3:3*N,1)      ;

dbx(N+1,1) = -dbx(1,1);
dby(N+1,1) = -dby(1,1);
dbz(N+1,1) = -dbz(1,1);

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
%----------------------------------------------------
%           N        Lk        L                          alpha
%
%          240       4         1.167879653094892e+01     1.167994912472664e+01
%          240       6         1.618789163145617e+01     1.619096178086449e+01
%          240       8         2.023333699809010e+01     2.023933375306768e+01
%          240       8         2.396846060643186e+01     2.397843242686577e+01
%          480       8         2.395742609593649e+01     2.395991351234742e+01
%          240      12fold     2.747534830559159e+01     2.749037408464824e+01

%=== curvature calculation ===

 %=== base curve ====
 
 rxg = [rx(N,1) rx'  rx(2,1)]' ;
 ryg = [ry(N,1)  ry'  ry(2,1)]' ;
 rzg = [rz(N,1)  rz'  rz(2,1)]' ;
 
 rxp = rxg(3:end,1) - rxg(2:end-1,1) ;
 ryp = ryg(3:end,1) - ryg(2:end-1,1) ;
 rzp = rzg(3:end,1) - rzg(2:end-1,1) ;
 
 
 rxpp = rxg(3:end,1) + rxg(1:end-2,1)-2*rxg(2:end-1,1);
 rypp = ryg(3:end,1) + ryg(1:end-2,1)-2*ryg(2:end-1,1);
 rzpp = rzg(3:end,1) + rzg(1:end-2,1)-2*rzg(2:end-1,1);

 k20  = rxpp.^2 + rypp.^2 + rzpp.^2;
 
 %==============
 
 rx1g = [rx1(N,1) rx1'  rx1(2,1)]' ;
 ry1g = [ry1(N,1)  ry1'  ry1(2,1)]' ;
 rz1g = [rz1(N,1)  rz1'  rz1(2,1)]' ;
 
 rx1p = rx1g(3:end,1) - rx1g(2:end-1,1) ;
 ry1p = ry1g(3:end,1) - ry1g(2:end-1,1) ;
 rz1p = rz1g(3:end,1) - rz1g(2:end-1,1) ;
 
 
 rx1pp = rx1g(3:end,1) + rx1g(1:end-2,1)-2*rx1g(2:end-1,1);
 ry1pp = ry1g(3:end,1) + ry1g(1:end-2,1)-2*ry1g(2:end-1,1);
 rz1pp = rz1g(3:end,1) + rz1g(1:end-2,1)-2*rz1g(2:end-1,1);

 k21  = rx1pp.^2 + ry1pp.^2 + rz1pp.^2;
 
 
 %=== curvature of the binormal vector ===
 
 bxg = [-bx(N,1) bx'  -bx(2,1)]' ;
 byg = [-by(N,1)  by'  -by(2,1)]' ;
 bzg = [-bz(N,1)  bz'  -bz(2,1)]' ;
 
 
 bxp = bxg(3:end,1) - bxg(2:end-1,1) ;
 byp = byg(3:end,1) - byg(2:end-1,1) ;
 bzp = bzg(3:end,1) - bzg(2:end-1,1) ;
 
 bxpp = bxg(3:end,1) + bxg(1:end-2,1)-2*bxg(2:end-1,1);
 bypp = byg(3:end,1) + byg(1:end-2,1)-2*byg(2:end-1,1);
 bzpp = bzg(3:end,1) + bzg(1:end-2,1)-2*bzg(2:end-1,1);

 k2b  = bxpp.^2 + bypp.^2 + bzpp.^2;
 
 l = tau;
 temp = (k2b/h^4)/l^2-l^2;
 
 
 %==== checking if constraints are satisfied ===
 
 %=== unit b ===
 % . b\cdot \delta b = 0
 
 du = bx.*dbx + by.*dby + bz.*dbz;
 
 %=== unit speed === 
 % b' \cdot \delta b' = 0
 
 dbxp = dbx(2:N+1,1)-dbx(1:N);
 dbyp = dby(2:N+1,1)-dby(1:N);
 dbzp = dbz(2:N+1,1)-dbz(1:N);

 
 dup = bxp(1:N,1).*dbxp + byp(1:N,1).*dbyp + bzp(1:N,1).*dbzp ;

  
 %=== closure constraint ===
 
 b_dpx = by(1:N,1).*dbzp - bz(1:N,1).*dbyp;
 b_dpy = bz(1:N,1).*dbxp - bx(1:N,1).*dbzp;
 b_dpz = bx(1:N,1).*dbyp - by(1:N,1).*dbxp;
 
 db_bxp = dby(1:N,1).*bzp(1:N,1) - dbz(1:N,1).*byp(1:N,1);
 db_byp = dbz(1:N,1).*bxp(1:N,1) - dbx(1:N,1).*bzp(1:N,1);
 db_bzp = dbx(1:N,1).*byp(1:N,1) - dby(1:N,1).*bxp(1:N,1);
 
 
 closure = u*(b_dpx + db_bxp) + v*(b_dpy + db_byp) + w*(b_dpz + db_bzp) ;
 
 
 clx = (dby(1:N,1).*bz(2:N+1,1) -  dbz(1:N,1).*by(2:N+1,1)) + (by(1:N,1).*dbz(2:N+1,1) -  bz(1:N,1).*dby(2:N+1,1)) ;
 cly = (dbz(1:N,1).*bx(2:N+1,1) -  dbx(1:N,1).*bz(2:N+1,1)) + (bz(1:N,1).*dbx(2:N+1,1) -  bx(1:N,1).*dbz(2:N+1,1)) ;
 clz = (dbx(1:N,1).*by(2:N+1,1) -  dby(1:N,1).*bx(2:N+1,1)) + (bx(1:N,1).*dby(2:N+1,1) -  by(1:N,1).*dbx(2:N+1,1)) ;


 clx2 = (dby(1:N,1).*dbz(2:N+1,1) -  dbz(1:N,1).*dby(2:N+1,1)) + (dby(1:N,1).*dbz(2:N+1,1) -  dbz(1:N,1).*dby(2:N+1,1)) ;
 cly2 = (dbz(1:N,1).*dbx(2:N+1,1) -  dbx(1:N,1).*dbz(2:N+1,1)) + (dbz(1:N,1).*dbx(2:N+1,1) -  dbx(1:N,1).*dbz(2:N+1,1)) ;
 clz2 = (dbx(1:N,1).*dby(2:N+1,1) -  dby(1:N,1).*dbx(2:N+1,1)) + (dbx(1:N,1).*dby(2:N+1,1) -  dby(1:N,1).*dbx(2:N+1,1)) ;

 