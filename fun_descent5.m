  function [F,jac] = fun_descent5(x)


global  N  N1  u v w    err   h    lm x0 gamma bt
global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p bext   k  v0 ht  

  %  x = initial_descent5();
  

 
%----- reading input and constructing bx by bz ----
% 
% bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
% by(2:N,1)     = x(2:3:3*N1-1,1)    ;
% bz(2:N,1)     = x(3:3:3*N1,1)      ;
% 
% 
% rho           = x(3*N1+1:4*N1,1)   ;
% lm            = x(4*N1+1:5*N1+1,1) ;
% 
% u             = x(5*N1+2 ,1)       ;
% V             = x(5*N1+3,1)        ;
% w             = x(5*N1+4,1)        ;

x1 = x(1:5*N1+4,1);   % array storing bx by bz and the lagrange multipliers 

%==== initializations ===
temp = 0;
V = zeros(5*N1+4,k);        % array storing eigenvectors in k columns 

Ex   = zeros(5*N1+4,1)      ;
Expv = zeros(5*N1+4,1)      ;
Exmv = zeros(5*N1+4,1)      ;
Hx   = zeros(5*N1+4,5*N1+4) ;
Hxpv = zeros(5*N1+4,5*N1+4) ;
Hxmv = zeros(5*N1+4,5*N1+4) ;

jac = zeros((1+k)*(5*N1+4),(1+k)*(5*N1+4));

for i = 1:k
    
    indi = (5*N1+4)*(i-1)+1:(5*N1+4)*i;
    V(:,i) = x(5*N1+4+indi,1);
    temp = temp + V(:,i)'*V(:,i);
end

prefac_x = (1-2*temp);

 
%=== Evaluating gradient and hessian at x  without k-HISOD method =========

[Ex,Hx] = fun_jacobian6(x1);


%===== Gradient and Hessian for k-HISOD method with respect to x

f_b(1:3*N1,1)        = bt*prefac_x*Ex(1:3*N1,:) + (x1(1:3*N1,1)-x0(1:3*N1,1))/ht   ;
f_b(3*N1+1:5*N1+4,1) = Ex(3*N1+1:5*N1+4,1)                                        ;


jac(1:3*N1,1:3*N1)        = bt*prefac_x*Hx(1:3*N1,1:3*N1) + eye(3*N1)/ht         ;
jac(1:3*N1,3*N1+1:5*N1+4) = bt*prefac_x*Hx(1:3*N1,3*N1+1:5*N1+4)                 ;
 
jac(3*N1+1:5*N1+4,1:5*N1+4)    = Hx(3*N1+1:5*N1+4,:);


%===== Gradient and Hessian for k-HISOD method with respect to V

hv = 0.001;

for i = 1:k
    
    indi = (5*N1+4)*(i-1)+1:(5*N1+4)*i;
      %=== d x_dot/dV
    jac(1:3*N1,5*N1+4+indi) = -4*bt*Ex(1:3*N1,1)*V(:,i)';
    
end

%===========================  Steps for eValuating V_dot in k--HIOSD
%METHOD


%==== Projected hessian calculation for V_dot =========

i = 1;

indi = (5*N1+4)*(i-1)+1:(5*N1+4)*i;

 

[Expv,HxpV] = fun_jacobian6(x1+hv*V(:,i));   %   gradE(x+hv)
[Exmv,HxmV] = fun_jacobian6(x1-hv*V(:,i));   %   gradE(x-hv)

f_b(5*N1+4+indi,1) = gamma*(1-V(:,i)'*V(:,i))*(Expv-Exmv)/2/hv  + (V(:,1)-v0(:,1))/ht ;

%== jacobian ===

jac(5*N1+4+indi,1:5*N1+4) = gamma*(1-V(:,i)'*V(:,i))*(HxpV-HxmV)/2/hv; % -- with respect
% -- lagrange multipliers

jac(5*N1+4+indi,5*N1+4+indi) = - gamma*2*(Expv-Exmv)/2/hv*V(:,i)' ...
    + gamma*(1-V(:,i)'*V(:,i))*(HxpV-HxmV)/2  + eye(5*N1+4)/ht;


for i = 2:k
    
    indi = (5*N1+4)*(i-1)+1:(5*N1+4)*i;
    
    temp = 0;
    for j = 1:i-1
     %   temp = temp + V(:,j)'*V(:,j)*(Expv-Exmv +Expv-Exmv')/2/hv;
        temp = temp + V(:,j)'*V(:,j);
    end
    
    
    [Expv,HxpV] = fun_jacobian6(x1+hv*V(:,i));   %   gradE(x+hv)
    [Exmv,HxmV] = fun_jacobian6(x1-hv*V(:,i));   %   gradE(x-hv)
    
    f_b(5*N1+4+indi,1) = gamma*(1-V(:,i)'*V(:,i)-2*temp)*(Expv-Exmv)/2/hv  + (V(:,i)-v0(:,i))/ht ;
    
    %== jacobian
    
    jac(5*N1+4+indi,1:5*N1+4) = gamma*(1-V(:,i)'*V(:,i)-2*temp)*(HxpV-HxmV)/2/hv; % -- with respect
    % -- lagrange multipliers
    
    jac(5*N1+4+indi,5*N1+4+indi) = -gamma*2*(Expv-Exmv)/2/hv*V(:,i)' ...
        + gamma*(1-V(:,i)'*V(:,i)-2*temp)*(HxpV-HxmV)/2 + eye(5*N1+4)/ht; 
    
    for j = 1:i-1
        indj = (5*N1+4)*(j-1)+1:(5*N1+4)*j;
        
        jac(5*N1+4+indi,5*N1+4+indj) = -gamma*4*(Expv-Exmv)/2/hv*V(:,j)'   ;
        
    end
    
end




F = f_b;

err = sqrt(sum(f_b.^2)) ;
  end