clear all
clc


%===== steepest descent method ==

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1 fac steps t1  N2 bt k

global bx by bz bxp byp bzp    N   len rx ry rz Lk err    branch sig   ht beta gamma x2 nstep path dsum

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

 

pre_descent6();

 


t_start=0;

 t_end=  nstep*ht;

  opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','on');
 %tspan =  t_start:0.001:1*0.001 ;

tspan = t_start:ht:t_end;

dsum = 0;
for p1  =1:500%:500%1:500%1:2000%1600%:1700%1:nstep%300:500%1:length(tspan)%100%3: 
    sv = 1;
    
    p1
      
    var_initial = initial_descent6();
    
     x = my_rk4(var_initial);
    
%     [tode45,x]=ode45(@fun_descent7,tspan, var_initial,opts);
%     x=x';
    count = 1;
    %--------------------------- Solver ---------------------------------------
    
%     %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
%     options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','off', 'TOlFun',10^(-9),'TOlX',10^(-9),'MaxFunEvals',300  ) ;
%     options.MaxIter     = 200  ;
%     [x,fval,exitflag,output,qd1] =  fsolve(@fun_descent6      ,var_initial,options)          ;
%     
%     err1(p1) = err;
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
    
    
    sz = size(x);
    str = '';
    for i =1:sz(2)
        str = ['%30.16E ' str ];
    end
    str = [str ' ' '\r\n'];
    
    
    if(sv==1)
        
          str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' num2str(p1)   '.txt'];
        
      % str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' '_ht_' num2str(ht) '_nstep_' num2str(nstep)   '.txt'];

        
        fileID = fopen([path str0],'w');
                  fprintf(fileID,'%30.16E   \r\n',x );
                 fclose(fileID);
        %
%         
%         for i = 1:sz(1)
%             fprintf(fileID,str,x(i,:) );
%         end
%         fclose(fileID);
%         
    end
end
 
%  sz = size(x);
%  
%  for p =1:sz(2)
%      
%     bx(2:N,1)     = x(1:3:3*N1-2,p)    ;
%     by(2:N,1)     = x(2:3:3*N1-1,p)    ;
%     bz(2:N,1)     = x(3:3:3*N1,p)      ;
%     
%       
%     toc
%     
%     
%       
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
%     
%     %
%     % figure(2)
%     % plotbrowser on
%     % title('binormal')
%     % plot_sphere()
%     % hold on
%     % plot3(bx,by,bz,'color',col{1},'LineWidth',2)
%     % hold on
%     % set(gca, 'DataAspectRatio',[1,1,1]);
%     % daspect([1,1,1]);
%     % axis off
%     % set(gcf,'color',colw);
%     % axis vis3d
%     %
%     %
%     %
%     figure(3)
%     plotbrowser on
%     title('midline')
%     hold on
%     plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
%     hold on
%     set(gca, 'DataAspectRatio',[1,1,1]);
%     daspect([1,1,1]);
%     axis off
%     set(gcf,'color',colw);
%     axis vis3d
%      
%  end