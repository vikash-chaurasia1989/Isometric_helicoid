clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac E

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2 algo

format longE
tic


col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;


%==== loading data ====
branch = 1;

path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];


N = 240;
N1 = N-1;
h = 1/N;


str0 = 'branch_1_N240_tau_80936500000.txt';


x = load([path str0]);



bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;



%--------------------------------------------------------------------------

%------------------------- Plotting  --------------------------------------

%----------------- Post processing -----
h = 1;% mean(len);
%---- Tangent ti = bi \times bi+1

i = 1:N;


tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);


tx(N+1,1) = tx(1,1);
ty(N+1,1) = ty(1,1);
tz(N+1,1) = tz(1,1);

tx = -tx/h;
ty = -ty/h;
tz = -tz/h;

%--------------------

% %--- position vector using integration of tangent
% % initialization

rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;


for i = 1:N
    
    rx(i+1,1) = rx(i,1) + h*tx(i,1);
    ry(i+1,1) = ry(i,1) + h*ty(i,1);
    rz(i+1,1) = rz(i,1) + h*tz(i,1);
    
end




figure(2)
plotbrowser on
title('binormal')
plot_sphere()
hold on
plot3(bx,by,bz,'color',col{1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d



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

%
% figure(4)
% plotbrowser on
% title('tangent')
% hold on
% plot3(tx,ty,tz,'color',col{1},'LineWidth',2)
% hold on
% set(gca, 'DataAspectRatio',[1,1,1]);
% daspect([1,1,1]);
% axis off
% set(gcf,'color',colw);
% axis vis3d





