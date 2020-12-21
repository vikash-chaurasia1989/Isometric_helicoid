clear all
clc

global    N1  tau   h tx ty tz  s nx ny nz k0 al w r

global bx by bz   N     rx ry rz   bext or kappa A B tau1 nfold
set(0,'defaultlinelinewidth',2)


col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;



or =    -1;

nfold = 25*2;

N = 1000;% nfold*50;
N1 = N-1;
h = 1/N;
 
 
if (nfold==3*2)
    var_in = [ 4.043808837982855e+01,  1.678602921926909e+04] ;   % for 3pi solution
    tau = 8.09361 ;
end
    
if (nfold==5*2)
    var_in = [ 1.280530390782753e+01, 1.837823905891308e+05 ];      % for 5 pi solution 
    tau = 11.98423030321598 ;
end   

if (nfold==7*2)
    var_in = [ -8.732663887726436e+01, 7.988414727595540e+05]; % for 7pi solution 
    tau = 15.37971 ;
end 

if (nfold==9*2)
    var_in = [-2.664560362414188e+02, 2.325042020521571e+06 ];   % for 9pi solution 
    tau = 18.47472 ;
end 

if (nfold==11*2)
    var_in = [-5.284567289865532e+02, 5.388830397947116e+06];    % for 11pi solution 
    tau = 21.35104 ;
end 

if (nfold==25*2)
   % var_in = [-8.846833560597079e+03, 6.617454204512132e+08];    % for 11pi solution 
   
 var_in = [ 3.676534655636211e+03     9.711831243309331e+07];  % knot 1 -- N = 1000
 %var_in = [3.963756950488346e+03  9.478849000952253e+07];      % knot 1 -- N 300;
 
  var_in = [-9.752624690199411e+03     6.632285173077818e+08];  % knot 2 
  var_in = [-7.310643750503460e+03  6.448980633271666e+08];
  var_in = [-38.263037790275611,-38.263037790275611,pi/5];
    tau = 78.33 ;
    var_in = 12.330674378995440;%[12.33,tau];
end 
%tau = 78.5;
tau = 78.33;
tau1 = tau;
%tau1 = 78.4;
%var_in = [ 4.037960852507271e+01,  tau^4] ;   
  %var_in = [15  ,tau1^4];


%--------------------------- Solver ---------------------------------------
%==========================================================================
%
rx = zeros(N+1,1);
ry = zeros(N+1,1);
rz = zeros(N+1,1);
tx = zeros(N+1,1);
ty = zeros(N+1,1);
tz = zeros(N+1,1);
nx = zeros(N+1,1);
ny = zeros(N+1,1);
nz = zeros(N+1,1);
bx = zeros(N+1,1);
by = zeros(N+1,1);
bz = zeros(N+1,1);


%== initial point same as the saved data

tx(1,1) = 0;
ty(1,1) = 1;
tz(1,1) = 0;

bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

nx(1,1) = by(1,1)*tz(1,1) - bz(1,1)*ty(1,1);
ny(1,1) = bz(1,1)*tx(1,1) - bx(1,1)*tz(1,1);
nz(1,1) = bx(1,1)*ty(1,1) - by(1,1)*tx(1,1);

algo = 3;

stralgo = {'levenberg-marquardt' ,'trust-region-reflective' ,'trust-region-dogleg'};

options             = optimset('Display','iter', 'Algorithm',stralgo{algo},'Jacobian','off', 'TOlFun',10^(-32),'TOlX',10^(-32),'MaxFunEvals',695000  ) ;
options.MaxIter     = 50000  ;
[x,fval,exitflag,output,qd1] =  fsolve(@fun_coskappa  ,var_in,options)          ;
 
 

%==== constructing band from frenet frame data =====

N = length(kappa)-1;
h = 1/N;

% bz1 = x(:,3);

bext(1:3,1) = -1*[bx(N,1);by(N,1);bz(N,1)];
bext(4:6,1) = -1*[bx(2,1);by(2,1);bz(2,1)];

%---- derivative calculation for b2 --- bN  ---

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

title(' ')

%== saving the data ===


for i = 2:N+1
    bx(i,1) = bx(i-1,1) -tau*nx(i-1,1)*h;
    by(i,1) = by(i-1,1) -tau*ny(i-1,1)*h;
    bz(i,1) = bz(i-1,1) -tau*nz(i-1,1)*h;
    
end

branch = 15;


sv=0;

if(sv==1)
  path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
   % str0 = ['b_' num2str(nfold/2) '_knot_N' num2str(N) '_tau_' num2str(10^10*tau1) '_cos.txt'];
 
   str0 = ['b_' num2str(nfold/2) '_N' num2str(N) '_tau_' num2str(10^10*tau1) '_cos.txt'];

    %
    fileID = fopen([path str0],'w');
    fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
    fclose(fileID);
end
    
    %  A unknot = amp = 12.330674378995440 
 
