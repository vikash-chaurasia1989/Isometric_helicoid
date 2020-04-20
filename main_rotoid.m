clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 t x y z  a b c p

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

branch = 1;

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
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    
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
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
    
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



bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

%---- derivative calculation for b2 --- bN  ---

ig = 2:N;


h =1/N;% mean(len);


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


%----------------- Post processing -----
%---- Tangent ti = bi \times bi+1

i = 1:N+1   ;

%
tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);

% %
tx = tx/tau;
ty = ty/tau;
tz = tz/tau;

rx(1) = 0; ry(1) = 0; rz(1) = 0;

 

for i = 1:N

     rx(i+1) = rx(i) + h*tx(i+1);
     ry(i+1) = ry(i) + h*ty(i+1);
     rz(i+1) = rz(i) + h*tz(i+1);

end

l = sum(sqrt((rx(1:N)-rx(2:N+1)).^2 + (ry(1:N)-ry(2:N+1)).^2 + (rz(1:N)-rz(2:N+1)).^2));

rx = rx/l;
ry = ry/l;
rz = rz/l;

rx = rx - mean(rx);
ry = ry - mean(ry);
rz = rz - mean(rz);

 ind = [1 N/nfold+1 2*N/nfold+1];
  
 x1 = rx(ind );
 y1 = ry(ind );
 z1 = rz(ind );
 
n1 = [x1(1) y1(1) z1(1)];
n2 = [x1(2) y1(2) z1(2)];
n3 = [x1(3) y1(3) z1(3)];

normal = cross(n1 - n2, n1 - n3);
 
norm = sum(sqrt(normal.*normal));

a = normal(1)/norm;
b = normal(2)/norm;
c = normal(3)/norm;

%-- rotating the mid-plane to align with x-y plane 
%--- axis of rotation u = (u1, u2, u3)
%  u = n x k  = (b,-a,0)
u1 =  b;
u2 = -a;

th = acos(c);

Ru = [(cos(th) + u1^2*(1-cos(th))) u1*u2*(1-cos(th))                u2*sin(th);
       u1*u2*(1-cos(th))           (cos(th) + u2^2*(1-cos(th)))    -u1*sin(th);
       -u2*sin(th)                 u1*sin(th)                       cos(th)   ] ;
   
  
%-- rotating the mid-line such that the midplane is x-y plane 
for i = 1:N+1
    
temp =  Ru*[rx(i);ry(i);rz(i)] ;
 
rx(i) = temp(1);
ry(i) = temp(2);
rz(i) = temp(3);

end

 


a = .9*mean(sqrt(rx.^2+ry.^2+rz.^2));
b = .5*(max(rz)+ abs(min(rz)));
c = .9*b;
p =  nfold ;

var_in = [a,b,c];


t = linspace(0,2*pi,N+1);

%=== initialization for the parameteric curve 
x = 0*t;
y = 0*t;
z = 0*t;



  options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','off', 'TOlFun',10^(-20),'TOlX',10^(-20),'MaxFunEvals',69500  ) ;
  options.MaxIter     = 50000  ;
    %   [x2,fval,exitflag,output,qd1] =  fsolve(@fun_rotoid    ,var_in,options)          ;

   fun_rotoid(var_in)  
   
  
%--- rigid body rotation of rotoid   --

[maxv1,ind1] = max(rz);
[maxv2,ind2] = max(z);


th = acos((rx(ind1)*x(ind2) + ry(ind1)*y(ind2))/(sqrt(rx(ind1)^2+ry(ind1)^2)*sqrt(x(ind2)^2+y(ind2)^2)));

R = [cos(th) -sin(th);
     sin(th) cos(th)];
 
for i = 1:N+1 
    temp = R*[x(i);y(i)];
    
    x(i) = temp(1);
    y(i) = temp(2);
end


   
figure(1)

plot3(rx,ry,rz,'color',col{1},'LineWidth',.5)
hold on 

plot3(x,y,z,'color',col{4},'LineWidth',.5)
hold on  