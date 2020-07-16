clear all
clc

format longE

branch = 1;
 path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
 path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

%--- midline
temp = load([path 'branch_1_N105_tau_80940900000.txt']);

N = 105;%length(temp)-1;
N1=N-1;
h = 1/N;

bx(2:N,1)     = temp(1:3:3*N1-2,1)    ;
by(2:N,1)     = temp(2:3:3*N1-1,1)    ;
bz(2:N,1)     = temp(3:3:3*N1,1)      ;



bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;


% N = length(temp)-1;
% N1=N-1;
% h = 1/N;
% bx = temp(:,1);
% by = temp(:,2);
% bz = temp(:,3);
% 


%---- derivative calculation for b2 --- bN  ---

bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

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

i = 1:N+1   ;
%
tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
 
%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1) = 0; ry(1) = 0; rz(1) = 0;



for i = 1:N
    
    rx(i+1) = rx(i) + h*tx(i+1);
    ry(i+1) = ry(i) + h*ty(i+1);
    rz(i+1) = rz(i) + h*tz(i+1);
    
end


%================================================

 
%--- symmetry


ind = [1 N/3+1 2*N/3+1 ];


%--- symmetry point ---

x  = rx(ind);
y  = ry(ind);
z  = rz(ind);

xc  = mean(x);
yc  = mean(y);
zc  = mean(z);

%----------- transforming the mid-line to the origin and in alignment with

%-- shifting to the origin --
rx = rx-mean(rx(ind));
ry = ry-mean(ry(ind));
rz = rz-mean(rz(ind));
 
 
p1 = [rx(ind(1)) ry(ind(1)) rz(ind(1))];
p2 = [rx(ind(2)) ry(ind(2)) rz(ind(2))];
p3 = [rx(ind(3)) ry(ind(3)) rz(ind(3))];

normal = cross(p1 - p2, p1 - p3);

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


%------ Rotation matrix that will transform each point
R1u = [(cos(th) + u1^2*(1-cos(th))) u1*u2*(1-cos(th))                u2*sin(th);
    u1*u2*(1-cos(th))           (cos(th) + u2^2*(1-cos(th)))    -u1*sin(th);
    -u2*sin(th)                 u1*sin(th)                       cos(th)   ] ;
R1u = R1u/det(R1u)^(1/3);
 
%-- rotating the mid-line such that the midplane is x-y plane
for i = 1:N+1 
    
    temp =  R1u*[bx(i,1);by(i,1);bz(i,1)] ;
    
    bx(i,1) = temp(1);
    by(i,1) = temp(2);
    bz(i,1) = temp(3);
    
    
end


%---- derivative calculation for b2 --- bN  ---

bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];

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

i = 1:N+1   ;
%
tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
 
%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1) = 0; ry(1) = 0; rz(1) = 0;



for i = 1:N
    
    rx(i+1) = rx(i) + h*tx(i+1);
    ry(i+1) = ry(i) + h*ty(i+1);
    rz(i+1) = rz(i) + h*tz(i+1);
    
end
dl = sqrt((bx(1:N,1)-bx(2:N+1,1)).^2 + (by(1:N,1)-by(2:N+1,1)).^2+(bz(1:N,1)-bz(2:N+1,1)).^2);
tau = sum(dl);
 path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
  path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

      %  str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
        
        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_symmetry.txt'];
        
       %str0 = ['coil_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];

x(1:3:3*N1-2,1) = bx(2:N,1) ;
x(2:3:3*N1-1,1) = by(2:N,1) ;
x(3:3:3*N1,1)   = bz(2:N,1) ;
%                
        fileID = fopen([path str0],'w');
        fprintf(fileID,'%30.16E   \r\n',x );
        fclose(fileID);