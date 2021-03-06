      function fval =   fun_rk4(s, var_in)
  
      global     tau    A B s1 N

 %-- tangent      
 
 %tau1 = var_in(3);

al1 = A +sqrt(A^2+B^2);
al2 = 0;
al3 = -A +sqrt(A^2+B^2);
 
p = sqrt((al3-al2)/(al3+al1));
q = sqrt(1-al2/al3);
r = 1/2*sqrt(al3+al1);


[K,E] = ellipke(p);
%s = linspace(0,6*K/r,N+1);
 %tau = tau1;

[sn,cn,dn] = ellipj(r*s*6*K/r,p);


kappa =  sqrt(al3*(1-q^2*sn.^2));
% 
% kappa(N/6+1:(N/6+N/3)) = - kappa(N/6+1:(N/6+N/3));
% kappa(N-N/6+1:N+1) = -kappa(N-N/6+1:N+1);

if((s1(N/6+1)<s<s1(N/6+N/3))||(s1(N-N/6+1)<s<s1(N+1)))
    kappa = -kappa;
end
 %--- position vector 
 
 r1x = var_in(1);
 r1y = var_in(2);
 r1z = var_in(3);

  
 %--- tangent  
  
 t1x = var_in(4);
 t1y = var_in(5);
 t1z = var_in(6);
 
 
 %--- normal 
 
 n1x = var_in(7);
 n1y = var_in(8);
 n1z = var_in(9);
 
 
 b1x = var_in(10);
 b1y = var_in(11);
 b1z = var_in(12);
 
 
 
%  b1x = t1y*n1z - t1z*n1y;
%  b1y = t1z*n1x - t1x*n1z;
%  b1z = t1x*n1y - t1y*n1x;
%  
 
 
 %-- ode 
 
 
 f_r1x = t1x  ; 
 f_r1y = t1y  ;
 f_r1z = t1z  ;
 
 f_t1x = kappa*n1x  ; 
 f_t1y = kappa*n1y  ;
 f_t1z = kappa*n1z  ;
 
 f_n1x = -kappa*t1x + tau*b1x  ;
 f_n1y = -kappa*t1y + tau*b1y  ;
 f_n1z = -kappa*t1z + tau*b1z  ;
 
 
 f_b1x = -tau*n1x  ; 
 f_b1y = -tau*n1y  ;
 f_b1z = -tau*n1z  ;
 
 
 fval = [f_r1x f_r1y f_r1z f_t1x f_t1y f_t1z f_n1x f_n1y f_n1z f_b1x f_b1y f_b1z]'  ;
      
      end
 
 