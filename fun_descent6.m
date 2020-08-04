function F = fun_descent6(x)


global  N  N1      err    x0 gamma bt l
global   k  v0 ht

%  x = initial_descent5();



x1 = x(1:3*N1,1);   % array storing bx by bz

%==== initializations ===
temp = 0;
V = zeros(3*N1,k);        % array storing eigenvectors in k columns  of dimension 3*N1 x 1

Ex   = zeros(3*N1,1)      ;

Expv = zeros(3*N1,1)      ;
Exmv = zeros(3*N1,1)      ;


for i = 1:k
    
    indi = (3*N1)*(i-1)+1:(3*N1)*i;
    V(:,i) = x(3*N1+indi,1);
    temp = temp + V(:,i)'*V(:,i);
end

prefac_x = (1-2*temp);


%=== Evaluating gradient and hessian at x  without k-HISOD method =========

[Ex,gradG] = fun_projection(x1);


%===== Gradient and Hessian for k-HISOD method with respect to x

f_b(1:3*N1,1)  = bt*prefac_x*Ex(1:3*N1,:) + (x1(1:3*N1,1)-x0(1:3*N1,1))/ht   ;



%===== Gradient and Hessian for k-HISOD method with respect to V

hv = 0.001;


%===========================  Steps for eValuating V_dot in k--HIOSD
%METHOD


%==== Projected hessian calculation for V_dot =========

i = 1;

indi = (3*N1)*(i-1)+1:(3*N1)*i;

[Expv,gradGp] = fun_projection(x1+hv*V(:,i));   %   gradE(x+hv)
[Exmv,gradGm] = fun_projection(x1-hv*V(:,i));   %   gradE(x-hv)


gradGbt =  gradG'*(l*bt/gamma*((gradG*gradG'))\(gradGp-gradGm)/2/hv*(Ex-2*(V(:,i)'*Ex)*V(:,i)));

f_b(3*N1+indi,1) = gamma*(1-V(:,i)'*V(:,i))*(Expv-Exmv)/2/hv/l  + (V(:,1)-v0(:,1))/ht - gradGbt ;



for i = 2:k
    
    indi = (3*N1)*(i-1)+1:(3*N1)*i;
    
    temp = 0;
    for j = 1:i-1
        temp = temp + V(:,j)'*V(:,j);
    end
    
    
    [Expv,gradGp] = fun_projection(x1+hv*V(:,i));   %   gradE(x+hv)
    [Exmv,gradGm] = fun_projection(x1-hv*V(:,i));   %   gradE(x-hv)
    
    gradGbt =  gradG'*(l*bt/gamma*((gradG*gradG'))\(gradGp-gradGm)/2/hv*(Ex-2*(V(:,i)'*Ex)*V(:,i)));

    f_b(3*N1+indi,1) = gamma*(1-V(:,i)'*V(:,i)-2*temp)*(Expv-Exmv)/2/hv/l  + (V(:,i)-v0(:,i))/ht - gradGbt  ;
    
    
end




F = f_b ;

err = sqrt(sum(f_b.^2))  
end