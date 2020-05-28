clear all
clc


%===== steepest descent method == 

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac steps t1 d21 d13 d23 d N2

global bx by bz bxp byp bzp bx2p by2p bz2p bx4p by4p bz4p  N   len rx ry rz Lk err   nfold branch sig w1 w2 w3 tau2 c ht

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


N = 105;
% ht = 0.005;
% c =  0.2;
sig = 0.01;

 

t1 = 0:ht:3;
tau = 16;





n = tau/2/pi;

fac = (asinh(n*pi*sig))/(4*pi^3*n^3);


N2 = 1000;

pre_descent2(N2);

ht = d21/N2;




%pre_descent();

%==========================================================================

t1 = 0:ht:d21;




for p1  = 2:5%:length(t1)%100%3:length(steps)%length(tau1)%length(tau1)%length(tau1)%%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
    sv = 1;
      
    d = t1(p1);    % incremental distance between the two curves 
    
    var_initial = initial_descent2();
    
   
    
    count = 1;
    %--------------------------- Solver ---------------------------------------
    
    %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
    options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-30),'TOlX',10^(-30),'MaxFunEvals',69500  ) ;
    options.MaxIter     = 50000  ;
    [x,fval,exitflag,output,qd1] =  fsolve(@fun_descent      ,var_initial,options)          ;
   
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
    
    
    
    
    
    
    if(sv==1)
        
%         path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent/'    ;
%         str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) 'step' num2str(p1) '_c_' num2str(round(1000*c)) '_ht_' num2str(round(1000*ht)) '.txt'];
%         
        path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent2/'    ;
     %   str0 = ['branch_' num2str(213) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_h_' num2str(round(10^5*ht)) '_step_' num2str(p1) '.txt'];
        str0 = ['branch_' num2str(21) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(N2) '_' num2str(p1)   '.txt'];

        
        fileID = fopen([path str0],'w');
        fprintf(fileID,'%30.16E   \r\n',x );
        fclose(fileID);
        
        
        
    end
end
%=== saving the torsion values for a given branch

if(sv==1)
    strtau = ['tau_branch_' num2str(branch) '.txt'];
    
    fileID = fopen([path strtau],'w');
    fprintf(fileID,'%30.16E   \r\n',tau1' );
    fclose(fileID);
end

% sv2 = 1;
% 
% if(sv2==1)
%     path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
%     str0 = ['9pi_knot_N' num2str(N) '_tau_' num2str(10^10*tau) '.txt'];
%     strb = ['b_' str0];
%     strr = ['r_' str0];
%     
%     strlm = ['lm_' str0];
%     strrho = ['rho_' str0];
%     struvw = ['uvw_' str0];
%     
%     
%     %
%     fileID = fopen([path strb],'w');
%     fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
%     fclose(fileID);
%     
%     
%     
%     fileID = fopen([path strlm],'w');
%     fprintf(fileID,'%30.16E   \r\n',lm );
%     fclose(fileID);
%     
%     fileID = fopen([path strrho],'w');
%     fprintf(fileID,'%30.16E   \r\n',rho );
%     fclose(fileID);
%     
%     
%     
%     fileID = fopen([path struvw],'w');
%     fprintf(fileID,'%30.16E   \r\n',[u v w]' );
%     fclose(fileID);
%     
% end
