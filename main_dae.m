clear all
clc


%===== steepest descent method ==

global    lm  N1 p1 count tau  rho f_b  h tx ty tz  u v w tau1 p1  t1  N2 id

global bx by bz bxp byp bzp    N   rx ry rz  err    branch   ht   c tspan nnt

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


parameters();

branch =2;

path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];
strtau = ['tau_branch_' num2str(branch) '.txt'];
tau1 = load([path strtau]);
[val,ind] = min(abs(tau1-tau));
tau2 = tau1(ind);               % nearest value to input tau for which we have the data
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau2)) '.txt'];
x0 = load([path str0]);  % initial point for pseudodynamics

%===== a perturbed initial point from stationary point with nonzero
%gradient
path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution2/' ;
strb = 'b_branch_213_N_105_t_5_tau_160000000000.txt';
temp = load([path strb]);

%== inserting binormals from the perturbed data into x0

x0(1:3:3*N1-2,1) = temp(2:N,1);
x0(2:3:3*N1-1,1) = temp(2:N,2);
x0(3:3:3*N1,1)   = temp(2:N,3);

%x0(3*N1+1:5*N1+4,1) = 0;
% 



id = eye(3);
bx(1,1) = 1;
by(1,1) = 0;
bz(1,1) = 0;

bx(N+1,1) = -1;
by(N+1,1) =  0;
bz(N+1,1) =  0;

x_out = dae2('fun_dae',tspan,x0,nnt);

%--------------------------- Solver ---------------------------------------
%
% %options.Algorithm   = 'levenberg-marquardt'  ;    %'trust-region-reflective' ;
% options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','off', 'TOlFun',10^(-9),'TOlX',10^(-9),'MaxFunEvals',69500  ) ;
% options.MaxIter     = 50000  ;
% [x,fval,exitflag,output,qd1] =  fsolve(@fun_descent4      ,var_initial,options)          ;
%
%
err1(p1) = err;
%--------------------------------------------------------------------------

sv = 1;

sz = size(x_out);

str = '';
for i =1:sz(2)
    str = ['%30.16E ' str ];
end
str = [str ' ' '\r\n'];

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_dae/'    ;
str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '_step_' '_ht_' num2str(ht) '_tspan_' num2str(length(tspan))   '.txt'];

%
fileID = fopen([path str0],'w');



for i = 1:sz(1)
    fprintf(fileID,str,x_out(i,:)' );
end
fclose(fileID);

%=== post processing ===

% for p1 = 1:sz(2)
%     x = x_out(:,p1);
%     
%     %----- reading input and constructing bx by bz ----
%     
%     bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
%     by(2:N,1)     = x(2:3:3*N1-1,1)    ;
%     bz(2:N,1)     = x(3:3:3*N1,1)      ;
%     
%     
%     %------------------------- Plotting  --------------------------------------
%     
%     
%     toc
%     
%     
%     i = 1:N+1   ;
%     
%     %
%     tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
%     ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
%     tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
%     tx = tx/tau;
%     ty = ty/tau;
%     tz = tz/tau;
%     
%     %--------------------
%     
%     % %--- position vector using integration of tangent
%     % % initialization
%     %h=1;
%     rx(1) = 0; ry(1) = 0; rz(1) = 0;
%     
%     
%     
%     for i = 1:N
%         
%         rx(i+1) = rx(i) + h*tx(i+1);
%         ry(i+1) = ry(i) + h*ty(i+1);
%         rz(i+1) = rz(i) + h*tz(i+1);
%         
%     end
%     
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
%     hold on
%     
%     %
%     % figure(4)
%     % plotbrowser on
%     % title('tangent')
%     % hold on
%     % plot3(tx,ty,tz,'color',col{1},'LineWidth',2)
%     % hold on
%     % set(gca, 'DataAspectRatio',[1,1,1]);
%     % daspect([1,1,1]);
%     % axis off
%     % set(gcf,'color',colw);
%     % axis vis3d
%     
% end

sv = 0;

if (sv==1)



%=== loading stability results ===
fid = figure;

str = ['animation_' num2str(21)];

%----------Set up the movie--------.
writerObj = VideoWriter([path str]); % Name it.
writerObj.FrameRate = 5; % How many frames per second.

open(writerObj);

for p1 = 1:sz(2)
    
    
rx = zeros(N+1,1);
ry = zeros(N+1,1);
rz = zeros(N+1,1);


bx = zeros(N+1,1);
by = zeros(N+1,1);
bz = zeros(N+1,1);

x = x_out(:,p1);

bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
by(2:N,1)     = x(2:3:3*N1-1,1)    ;
bz(2:N,1)     = x(3:3:3*N1,1)      ;
%     

    %--- boundary points --
    
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
 
    ig = 2:N;
    
    
    
    bext(1:3,1) = -[bx(N,1);by(N,1);bz(N,1)];
    bext(4:6,1) = -[bx(2,1);by(2,1);bz(2,1)];
    
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
 
    tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
    ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
    tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
    
  
    tx = tx/tau;
    ty = ty/tau;
    tz = tz/tau;
    
    
    
    rx(1,1) = 0;
    ry(1,1) = 0;
    rz(1,1) = 0;
    
    for i=1:N
        
        rx(i+1,1) = rx(i,1) + h*tx(i,1);
        ry(i+1,1) = ry(i,1) + h*ty(i,1);
        rz(i+1,1) = rz(i,1) + h*tz(i,1);
    end
    
    rx = rx-mean(rx);
    ry = ry-mean(ry);
    rz = rz-mean(rz);
    
    
    
    %== rescaling the midline
    
    l = sum(sqrt((rx(1:end-1,1)-rx(2:end,1)).^2 + (ry(1:end-1,1)-ry(2:end,1)).^2+(rz(1:end-1,1)-rz(2:end,1)).^2));
    
    rx = rx/l;
    ry = ry/l;
    rz = rz/l;
     
    
    %================================================
    
    
    figure(fid)
    
   %  subplot(1,2,1)
    %subplot_tight(1,2,1,[0,0])
    
    % set(fid, 'Position', [100, 100, 1049, 895]);
    
   % set(fid, 'Position', [100, 100, 2560, 1274]);
    %set(fid, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    %plotbrowser on
    title('midline')
    %hold on
    plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
  
    
    set(gca, 'DataAspectRatio',[1,1,1]);
    daspect([1,1,1]);
    axis off
    set(gcf,'color',colw);
    
     xlim auto
     ylim auto
     zlim auto
   
    axis vis3d
    %view(normal)
    view([1,1,1])
     
    hold off
    
    
    
    pause(.1);
    
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
    
    
    
    
end
close(writerObj);

end

