clear all
clc


%---- read r and b file and create ruling data ----
col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6


% col = {'','','r','b','g','k','m'} ;
col = {col1,col2,col3,col4,col5} ;
%--- midline ---

%===== tau for 3 fold ====
tau1 = [8.1:.1:9.5 10:1:19];

tau1 = [8.1:.1:9.5 10:1:11.9 11.97871 11.98071 11.98571 11.99071 11.99571  12:.1:13.9 14:.5:21];

N = 72;
h = 1/N;
%
% N1 = N-1;
%==================

%===== tau for 5 fold ====
% % % %
%     tau1 = [12:.1:13.9 14:.5:25];
% % %
%     tau1 = [11.97871 11.98071 11.98571 11.99071 11.99571  tau1];
%
%     tau1 = 12.727;
%     tau1 = 1.480000000000000e+01;
% % %
%  %   tau1 =  12.1;%11.97871 ;
% % % %
%       N = 120      ;
% % %==================

h = 1;

% N = length(temp)-1;
%
%  N = N/2;
% h = 1/N;


p1 = 1;

tau = tau1(p1);



%-- 3 fold --

%-- 3 fold --

str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '.txt'];

if(tau>11.9)
    str0 = ['3fold_N' num2str(N) '_tau_' num2str(10*tau) '_1.txt'];
end



%     %-- 5 fold ----
%       str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '_1.txt'];
% %
%       if(tau>13.9)
%               str0 = ['5pi_knot_N120_tau_' num2str(10*tau) '.txt'];
%      end
% %


strb = ['b_' str0];

strb = 'b_3fold_N288_z0_tau81.txt';
%--- load b
temp = load(strb);
N = length(temp)-1;

bx  = temp(1:N+1,1);
by = temp(1:N+1,2);
bz = temp(1:N+1,3);


%--- rotating about x axis by pi/2 --
R = [ cos(pi/2) -sin(pi/2);
    sin(pi/2)   cos(pi/2)];


for i = 1:N+1
    
    %-- Rotation about x axis by pi/2
    temp = R*[by(i,1);
        bz(i,1)];
    
    by(i,1) = temp(1);
    bz(i,1) = temp(2);
    
    %-- Rotation about z axis by pi/2
    temp = R'*[bx(i,1);
        by(i,1)];
    
    bx(i,1) = temp(1);
    by(i,1) = temp(2);
end


dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);

tau = sum(dl);


i = 1:N;

tx = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
ty = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
tz = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);



tx(N+1,1) = tx(1,1);
ty(N+1,1) = ty(1,1);
tz(N+1,1) = tz(1,1);

tx = tx/tau/h;
ty = ty/tau/h;
tz = tz/tau/h;


rx(1,1) = 0;
ry(1,1) = 0;
rz(1,1) = 0;

for i=1:N
    
    rx(i+1,1) = rx(i,1) + h*tx(i,1);
    ry(i+1,1) = ry(i,1) + h*ty(i,1);
    rz(i+1,1) = rz(i,1) + h*tz(i,1);
end

nfold = 3;
ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];



rx = rx -  mean(rx(ind,1));
ry = ry  -  mean(ry(ind,1));
rz = rz -   mean(rz(ind,1));




%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%                 Creating edges
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--- parameter for ruling --
%--- 5pi twist
wd = .005;

%-- edges---

rx0 = rx + wd*bx;
ry0 = ry + wd*by;
rz0 = rz + wd*bz;


%-- edges---

rx1 = rx -wd*bx;
ry1 = ry -wd*by;
rz1 = rz -wd*bz;


str12 = ['3pi_N_' num2str(N) '_tau_' num2str(10*tau) '_edge.txt'];




for i = 1:N
    
    p1 = 1+ 4*(i-1):4*i;
    
    %=====================
    data(p1(1),:) = [rx0(i,1) ry0(i,1) rz0(i,1)];
    data(p1(2),:) = [rx1(i,1) ry1(i,1) rz1(i,1)];
    
    data(p1(3),:) = [rx1(i+1,1) ry1(i+1,1) rz1(i+1,1)];
    data(p1(4),:) = [rx0(i+1,1) ry0(i+1,1) rz0(i+1,1)];
    
    
    
end



plot3(data(:,1),data(:,2),data(:,3))

fileID = fopen(str12,'w');
fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[data(:,1)'    ;data(:,2)'    ; data(:,3)'    ] );
fclose(fileID);
