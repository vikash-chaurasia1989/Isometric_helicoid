clear all
clc

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2

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




branch = 4;

if(branch==1)   % 8.0957344694550208E+00
    
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_3pi_3/';
    
    tau2 = load([path strtau]);  % old tau for which we have the data
    
    tau1 = [8.09409  8.1:.1:46.5];
    
    
    
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
    tau1 = 18.7;%19:.1:46.5;
end

   
%==========================================================================
N = 336;%126; %105;
N1 = N-1;
h = 1/N;

nfold=3;
sig = .06;


%=== most symmetric 5pi knot ==
%str0 = 'branch_2_N120_tau_119842500480_symmetry.txt';%;  % tau = 11.9842500480

for p1  =1:1%308:332%1:1%length(tau1)%length(tau1)%length(tau1)%%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
    sv = 1;
    tau =   tau1(p1);
    
    tau = 33;
    
    var_initial = initial_guess_coil();% initial_guess_symmetry();% initial_guess_evolution2();
    
   % tau = 18.58 ;
   %  tau = 19 ;
    
    
    count = 1;
    %--------------------------- Solver ---------------------------------------
    
    %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
    options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-30),'TOlX',10^(-30),'MaxFunEvals',69500  ) ;
    options.MaxIter     = 50000  ;
    [x,fval,exitflag,output,qd1] =  fsolve(@fun_jacobian6     ,var_initial,options)          ;
   
    err1(p1) = err;
    %--------------------------------------------------------------------------
    
    %------------------------- Plotting  --------------------------------------
  
    
    toc
    
        
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
    %
    % figure(3)
    % plotbrowser on
    % title('midline')
    % hold on
    % plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
    % hold on
    % set(gca, 'DataAspectRatio',[1,1,1]);
    % daspect([1,1,1]);
    % axis off
    % set(gcf,'color',colw);
    % axis vis3d
    
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
    
    
    
    
    
    
    if(sv==1)
        
        path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
       % str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
        
        s%tr0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_symmetry.txt'];
        
       str0 = ['coil_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];

                
        fileID = fopen([path str0],'w');
        fprintf(fileID,'%30.16E   \r\n',x );
        fclose(fileID);
        
        
        
    end
end
%=== saving the torsion values for a given branch
sv = 0;
if(sv==1)
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    
    fileID = fopen([path strtau],'w');
    fprintf(fileID,'%30.16E   \r\n',tau1' );
    fclose(fileID);
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