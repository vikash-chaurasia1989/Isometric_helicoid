clear all
clc


%===== steepest descent method ==

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac steps t1  N2 bt

global bx by bz bxp byp bzp    N   len rx ry rz Lk err    branch sig   ht beta gamma x2

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


parameters5();

 

pre_descent5();

 





for p1  = 1:length(t1)%100%3:length(steps)%length(tau1)%length(tau1)%length(tau1)%%262:length(tau1)%1:length(tau1)%3:length(tau1)%48:length(tau1)
    
    sv = 1;
    
    d = t1(p1);    % incremental distance between the two curves
    
    var_initial = initial_descent5();
    
    
    
    count = 1;
    %--------------------------- Solver ---------------------------------------
    
    %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
    options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-9),'TOlX',10^(-9),'MaxFunEvals',695  ) ;
    options.MaxIter     = 500  ;
    [x,fval,exitflag,output,qd1] =  fsolve(@fun_descent5      ,var_initial,options)          ;
    
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
        
        path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent7/'    ;
        
     %   path = '/Users/rtodres/Documents/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_descent5/';
        str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(p1)   '.txt'];
       
        
        fileID = fopen([path str0],'w');
        fprintf(fileID,'%30.16E   \r\n',x );
        fclose(fileID);
        
        
        
    end
end
 
 
