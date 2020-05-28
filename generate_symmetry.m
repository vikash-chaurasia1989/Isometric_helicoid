%=== This file reads a data and transform it such that bx(1,1)=1 and
%bx(N+1,1) = -1
% Also, it resizes the data. I do it for modifying the number of points
% Used to create initial guess for searching for symmetric data

clear all
clc

col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;



path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Constant torsion/For_hannes/'    ;
strb =  '05_fold_binormal_N240.txt';

temp = load([path strb]);

N = round((length(temp)-1)/2);

bx = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);


th = atan(-by(1,1)/bx(1,1));
R = [cos(th) -sin(th);
    sin(th)  cos(th)];

%
for i = 1:N+1
    temp = R*[bx(i,1) ; by(i,1)] ;
    bx(i,1) = temp(1);
    by(i,1) = temp(2);
end
%




%=== Midline

ig = 1:N;

tx(ig,1) = by(ig,1).*bz(ig+1,1) - bz(ig,1).*by(ig+1,1);
ty(ig,1) = bz(ig,1).*bx(ig+1,1) - bx(ig,1).*bz(ig+1,1);
tz(ig,1) = bx(ig,1).*by(ig+1,1) - by(ig,1).*bx(ig+1,1);

tx(N+1,1) = tx(1,1);
ty(N+1,1) = ty(1,1);
tz(N+1,1) = tz(1,1);

h = 1/N;

rx(1,1) = 0; ry(1,1) = 0; rz(1,1) = 0;



for i = 1:N
    
    rx(i+1,1) = rx(i,1) + h*tx(i+1,1);
    ry(i+1,1) = ry(i,1) + h*ty(i+1,1);
    rz(i+1,1) = rz(i,1) + h*tz(i+1,1);
    
end




figure(3)
plotbrowser on
title('midline')
hold on
plot_sphere()
hold on
plot3(bx,by,bz,'color',col{1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
%axis off
%set(gcf,'color',colw);
axis vis3d



figure(4)
plotbrowser on
title('midline')
hold on
plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
%axis off
%set(gcf,'color',colw);
axis vis3d


ds = sqrt((bx(2:end,1)-bx(1:end-1,1)).^2 + (by(2:end,1)-by(1:end-1,1)).^2 + (bz(2:end,1)-bz(1:end-1,1)).^2);

l = sum(ds);   

branch = 2;


path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

str0 = ['b_branch_' num2str(branch) '_N' num2str(N) '_symmetric.txt'];

fileID = fopen([path str0],'w');
fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
fclose(fileID);

