function f = fun_rotoid(var_in)

global rx ry rz   x y z t a b c p N

a = var_in(1);
b = var_in(2);
c = var_in(3);
%p = var_in(4);




x = (a-b*cos(p*t)).*cos(t) - c*sin(t) - (a-b)  ;
y = (a-b*cos(p*t)).*sin(t) + c*cos(t)-c  ;
z =  b*sin(p*t) +rz(1);



x = x-mean(x);
y = y-mean(y);
z = z-mean(z);


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


f = sum(sqrt((rx-x).^2 +  (ry-y).^2 + (rz-z).^2));


end

