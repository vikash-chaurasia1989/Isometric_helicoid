 function fval =   fun_frenet2( var_in)
  
 global k0 kappa  k or s  h N   A   f_t f_b f_n   tau   f_r tx ty tz nx ny nz bx by bz  rx ry rz al w p th B tau1 s1

A = var_in(1);
B = var_in(2);
 %tau1 = var_in(3);
 
%s = linspace(0,6*K/r,N+1);
s = linspace(0,1,N+1);

s1= linspace(0,1,N+1);
 
 
 
 opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
 

in_rk = [0 -0.1276 0 0.9280 0 -0.3726 0.3726 0 0.9280 0 1 0] ;  % initial point for rk4 integration 
[t y] = ode45(@fun_rk4, s,in_rk,opts);

rx = y(:,1);
ry = y(:,2);
rz = y(:,3);

tx = y(:,4);
ty = y(:,5);
tz = y(:,6);


nx = y(:,7);
ny = y(:,8);
nz = y(:,9);

bx = y(:,10);
by = y(:,11);
bz = y(:,12);
  
 
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
 fval = [f1 f2 f3] ;
%       
 
 end
 
 