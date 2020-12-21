%=== This file plots curvature of the midline for given torsion and branch




clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac kappa bxp byp bzp

global bx by bz    N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2  bxp byp bzp

format longE
tic


%---- read r and b file and create ruling data ----
col1 = [0,     0.976, 0.968]  ;      % -- color for m = 3
col2 = [0.831, 0.545, 0.247]  ;      % -- color for m = 5
col3 = [0.341, 0.505, 0.819]  ;      % -- color for m = 7
col4 = [0.705, 0.701, 0.070]  ;      % -- color for m = 9
col5 = [0.301, 0.811, 0.498]  ;      % -- color for m = 11

clr = linspecer(15);

col6  = [0.98, 0.325, 0];%[0.1705, 0.2525,  0.4095]   ;    % -- color for m = 13
col7  = col2/2;% [2.941176470588235e-01     5.447058823529412e-01     7.494117647058823e-01]  ;    % -- color for m = 15
col8  = [3.717647058823529e-01     7.176470588235294e-01     3.611764705882353e-01]/2 ;    % -- color for m = 17
col9  = [ 1.000000000000000e+00     5.482352941176470e-01     9.999999999999998e-02]  ;    % -- color for m = 19
col10 = col1/1.5;%[8.650000000000000e-01     8.109999999999999e-01     4.330000000000000e-01]   ;    % -- color for m = 21
col11 = [0.533, 0.671, 0.525];%col2/1.5;%[6.858823529411765e-01     4.035294117647059e-01     2.411764705882353e-01]   ;    % -- color for m = 23
col12 = [9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]   ;    % -- color for m = 25
col13 = [0.549, 0.035, 0.486];%[0.196, 0.522, 0.788];%col3/1.5;%[9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]/2 ;    % -- color for m = 27
col14 = [0.996, 0.773, 0.243];%col5/1.5;%[     4.350351799515555e-01     7.800722978266937e-01     6.461779788838612e-01]  ;  % -- color for m = 29                                                                        % -- color for m = 29
col15 = [0.11, 0.18, 1];%[0.988, 0.427, 0.412];%col6/1.5;%[     4.350351799515555e-01     7.800722978266937e-01     6.461779788838612e-01]/3  ;  % -- color for m = 31                                                                      % -- color for m = 29




col = {col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12, col13,col14,col15} ;
%col ={clr(1,:),clr(2,:),clr(3,:),clr(4,:),clr(5,:),clr(6,:),clr(7,:),clr(8,:),clr(9,:),clr(10,:),clr(11,:),clr(12,:),clr(13,:),clr(14,:),clr(15,:)};


%col = {col1,col2,col3,col4,col5,clr(1,:),clr(2,:),clr(3,:),clr(4,:),clr(5,:),clr(6,:)} ;
colw = [41,44,52]/255;






%==========================================================================



sig = .01;


branch =  15;

tau = 78.35855  ;%  8.09365;
tau = 78.33;%
%tau = 83.6;
%tau = 46.8125;
nfold = 25*2;


path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%  path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

% strtau = ['tau_branch_' num2str(branch) '.txt'];
% strLk = ['Linking_branch_' num2str(branch) '.txt'];
% tau1 = load([path strtau]);
% Lk =  load([path strLk]);
%
%
% [temp,ind] = min(abs(tau0-tau1));  % index of the saved data nearest to the given tau value
%
% tau  = tau1(ind);
%
%
%
%


if (branch==1||branch==2||branch==3)
    N=105;
end

if(branch==4||branch==9)
    N = 180; %
end
if(branch==5)
    N = 140;
end

if(branch== 6)
    N = 90;
end

if(branch== 7)
    N = 150;
end


if(branch== 8)
    N = 150;
end

if(branch== 10)
    N = 220;
end

if(branch==11||branch==14||branch==15)
    N=170;
end

if(branch== 12)
    N = 190;
end


if(branch==13)
    N=210;
    
    %N=300;
    
    
end

N = 600;%225;
N1 = N-1;
h = 1/N;




%=== reading data==
%str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
str0 =  'branch_15_N600_tau_782000000000_unknot1.txt';
str0 = 'branch_15_N600_tau_783300000000_cos.txt';
% str0 = 'branch_15_N600_tau_783300000000_cos.txt';
% str0 =  'branch_15_N600_tau_782000000000_cos.txt';
% str0 =   'branch_15_N600_tau_783585500000_cos.txt';    % this is the final data used for curvature of m=25 in the SI final version
%str0 = 'branch_8_N405_tau_844150000000.txt';
temp = load([path str0]);

x = temp(1:3*N1,1);


bx = zeros(N+1,1);
by = zeros(N+1,1);
bz = zeros(N+1,1);


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





E = energy_b(bx,by,bz);  % == this also gives the curvature as global parameter

i = 1:N+1   ;


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
rx(1) = 0; ry(1) = 0; rz(1) = 0;



for i = 1:N
    
    rx(i+1) = rx(i) + h*tx(i+1);
    ry(i+1) = ry(i) + h*ty(i+1);
    rz(i+1) = rz(i) + h*tz(i+1);
    
end



figure(3)
plotbrowser on
title('midline')
hold on
plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d





kappa = sqrt(kappa);

%== for unknotted_ m25_1
% ind = 192;
% kappa(ind) = .5*(kappa(ind-1)+kappa(ind+1));
% ind = 218;
% kappa(ind) = .5*(kappa(ind-1)+kappa(ind+1));
% ind = 384;
% kappa(ind) = .5*(kappa(ind-1)+kappa(ind+1));
% ind = 410;
% kappa(ind) = .5*(kappa(ind-1)+kappa(ind+1));



temp = N/nfold:3*N/nfold;
for i = 1:((nfold/2-1)/2)
    
    ind = 4*(i-1)*N/nfold+1 + temp;%((2*i-1)*N/nfold + 1):(((2*i+1)*N/nfold));
    
    kappa(ind) = -kappa(ind);
    
end
% kappa(N/10+1:3*N/10)      = -kappa(N/10+1:3*N/10);
% kappa(5*N/10+1:7*N/10) = -kappa(5*N/10+1:7*N/10);
% kappa(9*N/10+1:end)    = - kappa(9*N/10+1:end);
%
ind = ((nfold-1)*N/nfold+1):(N+1);
%
kappa(ind) = -kappa(ind);

% %=== for p3 ---
%  kappa(ind) = -kappa(ind);
%
%  %kappa(41:61) = -kappa(41:61) ;
%
%  kappa(22:40) = -kappa(22:40) ;
%
%
%  kappa(62:73) = -kappa(62:73) ;
%
%  kappa(169:178) = -kappa(169:178) ;
%
%  kappa(201:219) = -kappa(201:219) ;
% %--------------------------

kappa = kappa';

s = linspace(0,1,N+1);

figure(5)
plotbrowser on
plot(s,kappa,'-o')