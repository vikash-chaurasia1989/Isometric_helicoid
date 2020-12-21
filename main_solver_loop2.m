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




branch = 1;

if(branch==1)   % 8.0957344694550208E+00
    
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    
%    tau2 = load([path strtau]);  % old tau for which we have the data
    
    tau1 = [8.09409  8.1:.1:46.5];% 46.525:.015:47];
    
    
    
elseif(branch==2) % 1.1984230303215979E+01
    %================ Branch 2 =====================
    
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_5pi_3/';
    
    tau2 = load([path strtau]);  % old tau for which we have the data
    
    tau1 =  [12:.1:34.1 34.11 34.12 34.125 34.126 34.13 34.14 34.16 34.17 34.18 34.19 34.3]   ;
    
    tau1 =  [1.1984275004796786e+01 11.9875 11.99 12:.1:34 ]   ;%1.197871304796786e+01 tau for symmetric knot

    % =============== Branch 3 =====================
     
elseif(branch==3) % 1.1984230303215979E+01
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_7pi_3/';
    
    tau2 = load([path strtau]);  % old tau for which we have the data
    
    tau1 = [1.537950616265322e+01    15.4195 15.425 15.45 15.5:.1:46.5]; % 1.5379706162645375E+01 
    
else 
    tau1 =  1:.1:20;
end

   
%==========================================================================
N = 120;%336;%126; %105;
N1 = N-1;
h = 1/N;

branch = 6;

% nfold=3;
% sig = .06   

tau1 =  13:-.1:8.5;% 12:-.1:10;%15:-.1:12;%[15 15.05  15.1:0.1:17 17.1:.1:40 40.1:.1:48];%15:.1:18;
%=== most symmetric 5pi knot ==
%str0 = 'branch_2_N120_tau_119842500480_symmetry.txt';%;  % tau = 11.9842500480
 tau1 = 15.1:.1:48;   % branch 4
 
 tau1 = 21.5:.1:60;
 
 tau1 = 27:.1:60;% [27 26.5];
 
 tau1 = [27.95 28];
 tau1 = 28:.1:65;
 
%  branch = 4;
%  strtau = ['tau_branch_' num2str(branch) '.txt'];
%         path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%  tau1 = load([path strtau]);

%  tau2 = 8.1:.1:46.5;
% %tau1 = [tau2' 48.1:.1:61.3];
% 
% tau1 = [tau2  46.6:.1:85];
% % % % % % %
N = 105;
N1 = N-1;
h = 1/N;

 

tau1 = 21.4:.1:85; 

tau1 = 80.4:.1:85; 

branch = 15  ;
for p1  = 1:1%2:length(tau1)%1:10%446:length(tau1)%388:length(tau1)%2:length(tau1)%53:length(tau1)%53:length(tau1)%308:332%1:1%length(tau1)%length(tau1)%length(tau1)%%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
    sv  = 1;
    tau =  78.35855  ;% 8.09365 ;% tau1(p1);%16;% 2.148057962831568e+01
% tau1(p1);%8.093946633549898e+00
    
   % tau = 11;
   % var_initial = initial_guess_loop2();
    
   % tau = 18.58 ;
   %  tau = 19 ;
   
     sig = 0.01;

tau2 = tau;
 %===
n1 = tau/2/pi;
fac = (asinh(n1*pi*sig)/(4*pi^3*n1^3))  ;

E = .62355;
    
    %count = 2;
    %--------------------------- Solver ---------------------------------------
   stralgo = {'levenberg-marquardt' ,'trust-region-reflective' ,'trust-region-dogleg'}; 
   path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

 % path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
 
    %=== two step algorithm  
    % first we use trust region method for certain steps. Then we use
    % levenberg-marquardt method using the saved solution as initial guess 
        %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
   %=== Trust region method 
%    algo = 3 ;
%     
%      
%    var_initial = initial_guess_evolution2();%initial_orientable();% initial_guess_symmetry();% initial_guess_evolution2();
%   
%    options             = optimset('Display','iter', 'Algorithm',stralgo{algo},'Jacobian','on', 'TOlFun',10^(-15),'TOlX',10^(-15),'MaxFunEvals',6950    ) ;
%    options.MaxIter     = 50         ;
%    [x,fval,exitflag,output,qd1] =  fsolve(@fun_jacobian6     ,var_initial,options)          ;
%    % [x,fval,exitflag,output,qd1] =  fsolve(@fun_orientable     ,var_initial,options)          ;
%      
%    %=== saving the solution ===
%    
%    %branch = 14;
%    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
% 
%   % branch = 13;
%    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
% 
%    fileID = fopen([path str0],'w');
%    fprintf(fileID,'%30.16E   \r\n',x );
%    fclose(fileID);
% % %     
%    
   algo = 	1; 
   var_initial = initial_guess_evolution2();%initial_orientable();% initial_guess_symmetry();% initial_guess_evolution2();
   options             = optimset('Display','iter', 'Algorithm',stralgo{algo},'Jacobian','on', 'TOlFun',10^(-18),'TOlX',10^(-18),'MaxFunEvals',6950    ) ;
   options.MaxIter     = 3000      ;
   [x,fval,exitflag,output,qd1] =  fsolve(@fun_jacobian6   ,var_initial,options)          ;
   % [x,fval,exitflag,output,qd1] =  fsolve(@fun_orientable     ,var_initial,options)          ;
     
   %=== saving the solution ===
     
   str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_cos.txt'];
   
    % x = [x(1:3*N1,1)'  zeros(1,2*N1+4)]';
   
   fileID = fopen([path str0],'w');
   fprintf(fileID,'%30.16E   \r\n',x );
   fclose(fileID);     
        
   err1(p1) = err; 
   
   %--------------------------------------------------------------------------
    
   %------------------------- Plotting  --------------------------------------
      
    toc
        
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
    
     
    
    %
    % figure(2)
    % plotbrowser on
    % title('binormal')
    % plot_sphere()
    % hold on
    % plot3(bx,by,bz,'color',col{1},'LineWidth',2)
    % hold on
    % set(gca, 'DataAspectRatio',[1,1,1]);
    % daspect([1,1,1]);
    % axis off
    % set(gcf,'color',colw);
    % axis vis3d
    %
    %
    
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
    
    
    
    
    
    
%     if(sv==1)
%         
%         path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%        % path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_unknot_5pi/' ;
% 
%     % path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
% 
%           str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
%         
%       %    str0 = ['branch_unknot_5pi_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
% 
%      %   str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_symmetry.txt'];
%         
%        %str0 = ['coil_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
% 
%                 
%         fileID = fopen([path str0],'w');
%         fprintf(fileID,'%30.16E   \r\n',x );
%         fclose(fileID);
%     
%         
%         
%     end
end
%=== saving the torsion values for a given branch
 
if(sv==1)
%     strtau = ['tau_branch_' num2str(branch) '.txt'];
%     fileID = fopen([path strtau],'w');
%     fprintf(fileID,'%30.16E   \r\n',tau1 );
%     fclose(fileID);
end

sv2 = 0;

if(sv2==1)
    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    str0 = ['9pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
    strb = ['b_' str0];
    strr = ['r_' str0];
    
    strlm = ['lm_' str0];
    strrho = ['rho_' str0];
    struvw = ['uvw_' str0];
    
    
    %
    fileID = fopen([path strb],'w');
    fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
    fclose(fileID);
    
    
    
    fileID = fopen([path strlm],'w');
    fprintf(fileID,'%30.16E   \r\n',lm );
    fclose(fileID);
    
    fileID = fopen([path strrho],'w');
    fprintf(fileID,'%30.16E   \r\n',rho );
    fclose(fileID);
    
    
    
    fileID = fopen([path struvw],'w');
    fprintf(fileID,'%30.16E   \r\n',[u v w]' );
    fclose(fileID);
    
end


% kappa2 = sqrt((bx2p.*bx2p + by2p.*by2p + bz2p.*bz2p -tau^4)/tau^2);
% kappa2 = [kappa2(end,1)' kappa2(2:end-1,1)']';
% %kappa2 = [kappa2(end-1,1) kappa2(1:end-1,1)']';
% kappa3 = kappa2-min(kappa2);
% kappa3(N/6+1) = 0;
% kappa3(N/3+1) = kappa3(end);
% 
% kappa3(N/6+1:2*N/6-1) = -kappa3(N/6-1:-1:1);
% kappa3(N/3) = -kappa3(end);
% kappa3(N/3+1:N/3+N/6-1) = kappa3(N/3-1:-1:N/6+1);
% kappa3(N/2+1:5*N/6-1)  = -kappa3(N/2-1:-1:N/6+1);
% kappa3(5*N/6+1:end) =-kappa3(5*N/6-1:-1:2*N/3);
% 
% kappa3 = [-kappa3(end,1) kappa3(1:end,1)']';
% 
% s3 = linspace(0,1,length(kappa3));
% plot(s3,kappa3/1.380983773188140e+01);
% 
% 
% 
% %=== derivative of the curvature === 
% strkp = 'curvature_numerical';% ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
% 
% fileID = fopen([path strkp],'w');
%     fprintf(fileID,'%30.16E   \r\n',kappa3 );
%     fclose(fileID);



% golden ratio = 1.618033988749895e+00

sv2 = 0;


if(sv2==1)
    
    
    b(1:3:3*N1-2) = bx(2:N,1)  ;
    b(2:3:3*N1-1) = by(2:N,1)  ; 
    b(3:3:3*N1)   = bz(2:N,1)  ;
    % %====================
     x =  [b(1:3*N1)  rho' lm'  u  v     w  ]'  ;         % for fun_jacobian6
    
%    x =  [b(1:3*N1)   ]'  ;         % for fun_jacobian6

    path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    % path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_unknot_5pi/' ;
    
 %    path = ['/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
    
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_knot2.txt'];
    
    fileID = fopen([path str0],'w');
    fprintf(fileID,'%30.16E   \r\n',x );
    fclose(fileID);
    
    
    
end

    


    