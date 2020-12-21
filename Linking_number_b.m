%---
clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err qd   E bext

format longE


%-- This is the final file for calculating Linking number
%-- Here, I have used codes developed by Zin Arai and modified it for my
%purpose. The accuracy is impressive. Way better than the double integrals
%I was using to calculate Writhe and then the linking number using the
%Caligerano theorem  Lk = 2*(\tau/(2\pi) + Wr);

%---- read r and b file and create ruling data ----
col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [0.341, 0.505, 0.819]/2  ;         % -- color for mode n = 6
col7 = [0.831, 0.545, 0.247]/2  ;      % -- color for mode n = 3


col = {col1,col2,col3,col4,col5,col6,col7} ;

%temp = load('05_fold.txt');

N = 120;

N = 600;
branch = 15;
N1 = N-1;
h = 1/N;


tau = 78.2;
%str0 =  'branch_unknot_5pi_N105_tau_200000000000.txt';
str0 = 'branch_unknot_5pi_N120_tau_150000000000_unstable.txt';
str0 = 'branch_unknot_5pi_N120_tau_100000000000.txt';

str0 = 'branch_15_N600_tau_782000000000_unknot1.txt';

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_unknot_5pi/' ;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

rx = zeros(N+1,1);
ry = zeros(N+1,1);
rz = zeros(N+1,1);


bx = zeros(N+1,1);
by = zeros(N+1,1);
bz = zeros(N+1,1);



x = load([path str0]);


bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;

%--- boundary points --

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

bext(1:3,1) = -1*[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -1*[bx(2,1);by(2,1);bz(2,1)];

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
tx = tx/tau;
ty = ty/tau;
tz = tz/tau;

%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;




rx(1,1) = 0;
ry(1,1) = 0;
rz(1,1) = 0;

for i=1:N
    
    rx(i+1,1) = rx(i,1) + h*tx(i,1);
    ry(i+1,1) = ry(i,1) + h*ty(i,1);
    rz(i+1,1) = rz(i,1) + h*tz(i,1);
end



w = 0.002;


nx = by.*tz - bz.*ty;
ny = bx.*tz - bz.*tx;
nz = bx.*ty - by.*tx;


%


%-- edge curve and its tangent

rx0 = rx + w*bx;
ry0 = ry + w*by;
rz0 = rz + w*bz;



tx0 = tx -w*tau*nx;
ty0 = ty -w*tau*ny;
tz0 = tz -w*tau*nz;


%-- creating array for linking number subroutine --
gamma0 = zeros(3,N+1);
gamma1 = zeros(3,N+1);
for i = 1:N+1
    
    gamma0(1,i) = rx0(i,1);
    gamma0(2,i) = ry0(i,1);
    gamma0(3,i) = rz0(i,1);
    
    gamma1(1,i) = rx(i,1);
    gamma1(2,i) = ry(i,1);
    gamma1(3,i) = rz(i,1);
end
n= 0;


% gamma0 = gamma1;

length_0 = size(gamma0, 2);
length_1 = size(gamma1, 2);

% gamma0 = intval(gamma0);
% gamma1 = intval(gamma1);

for i = 1:length_0
    for j = 1:length_1
        a = gamma1(:, j)                    + 0     - gamma0(:, i);
        b = gamma1(:, j)                    - gamma0(:, mod(i, length_0) + 1);
        c = gamma1(:, mod(j, length_1) + 1) - gamma0(:, mod(i, length_0) + 1);
        d = gamma1(:, mod(j, length_1) + 1) + 0     - gamma0(:, i);
        n = n + solid_angle_quadrilateral(a, b, c, d);
        
        
    end
    
    
    
end


%===================================================
%n = n / (4 * midrad(pi, 1e-14))

n = n/(4*pi)*2;



Lk  = round(n);



