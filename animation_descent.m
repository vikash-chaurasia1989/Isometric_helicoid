%---
clear all
clc
global    lm  N1 p1   tau     h tx ty tz   fac  

global bx by bz bxp byp bzp   N  rx ry rz Lk err    branch sig   ht beta gamma x2 nstep path  
format longE


 
%---- read r and b file and create ruling data ----
col1 = [0,     0.976, 0.968]  ;      % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]  ;      % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]  ;      % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]  ;      % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]  ;      % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1,col2,col3,col4,col5,col6} ;
colw = [41,44,52]/255;

parameters5() 

%=== loading stability results ===
fid = figure;

str = 'animation_descent_213c';% ['animation_' num2str(213)];
str = 'animation_descent_1-unstable';% ['animation_' num2str(213)];


% Set up the movie.
writerObj = VideoWriter([path str]); % Name it.
writerObj.FrameRate = 5; % How many frames per second.

open(writerObj);




rx = zeros(N+1,1);
ry = zeros(N+1,1);
rz = zeros(N+1,1);


bx = zeros(N+1,1);
by = zeros(N+1,1);
bz = zeros(N+1,1);


nfold = 5;
 

 branch = 2;
 
 num1 = 3255;
 num3 = 2244;
         
for p1= 1:397%(num1+num3) %1:253%length(t1)
    
    %==== saving data =====
    
    % str0 = ['branch_213_N_' num2str(N) '_step_' num2str(p1) '_tau_' num2str(10^10*tau) '.txt'];
     
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(p1-1)   '.txt'];

%     if(p1>num1)
%         branch =3;
%         str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*16)) '_step_' num2str(num1+num3-p1)   '.txt'];
%     end
%     

    x = load([path str0]);

    bx(2:N,1)     = x(1:3:3*N1-2,1)    ;
    by(2:N,1)     = x(2:3:3*N1-1,1)    ;
    bz(2:N,1)     = x(3:3:3*N1,1)      ;

    
    %--- boundary points --
    
    bx(1,1) = 1;
    by(1,1) = 0;
    bz(1,1) = 0;
    
    bx(N+1,1) = -1;
    by(N+1,1) =  0;
    bz(N+1,1) =  0;
    
    E(p1) = energy_b(bx,by,bz);
    
    %dl = sqrt((bx(1:end-1,1)-bx(2:end,1)).^2 + (by(1:end-1,1)-by(2:end,1)).^2 +(bz(1:end-1,1)-bz(2:end,1)).^2);
    
    %tau = sum(dl);
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
    
    % %
    %    tx(i,1) = by(i,1).*bz(i+1,1) - bz(i,1).*by(i+1,1);
    %    ty(i,1) = bz(i,1).*bx(i+1,1) - bx(i,1).*bz(i+1,1);
    %    tz(i,1) = bx(i,1).*by(i+1,1) - by(i,1).*bx(i+1,1);
    %
    tx(i,1) = by(i,1).*bzp(i,1) - bz(i).*byp(i,1);
    ty(i,1) = bz(i,1).*bxp(i,1) - bx(i).*bzp(i,1);
    tz(i,1) = bx(i,1).*byp(i,1) - by(i).*bxp(i,1);
    
    %
    %    tx(N+1,1) = tx(1,1);
    %    ty(N+1,1) = ty(1,1);
    %    tz(N+1,1) = tz(1,1);
    % % %
    % %
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
    %==== reorienting the midline for animation
    
    
    ind = [1 1 + (N/nfold:N/nfold:((nfold-1)*N/nfold))];
    
    x1 = rx(ind );
    y1 = ry(ind );
    z1 = rz(ind );
    
    xc = mean(x1);
    yc = mean(y1);
    zc = mean(z1);
    
    
    rx = rx-xc;
    ry = ry-yc;
    rz = rz-zc;
    
    x1 = x1-xc;
    y1 = y1-yc;
    z1 = z1-zc;
    
    
    n1 = [x1(1) y1(1) z1(1)];
    n2 = [x1(2) y1(2) z1(2)];
    n3 = [x1(3) y1(3) z1(3)];
    
    normal = cross(n1 - n2, n1 - n3);
    
    norm = sum(sqrt(normal.*normal));
    
    a = normal(1)/norm;
    b = normal(2)/norm;
    c = normal(3)/norm;
    
    %-- rotating the mid-plane to align with x-y plane
    %--- axis of rotation u = (u1, u2, u3)
    %  u = n x k  = (b,-a,0)
    u1 =  b;
    u2 = -a;
    
    th = acos(c);
    
    Ru = [(cos(th) + u1^2*(1-cos(th))) u1*u2*(1-cos(th))                u2*sin(th);
        u1*u2*(1-cos(th))           (cos(th) + u2^2*(1-cos(th)))    -u1*sin(th);
        -u2*sin(th)                 u1*sin(th)                       cos(th)   ] ;
    
    
    %-- rotating the mid-line such that the midplane is x-y plane
%     for i = 1:N+1
%         
%         temp =  Ru*[rx(i);ry(i);rz(i)] ;
%         
%         rx(i) = temp(1);
%         ry(i) = temp(2);
%         rz(i) = temp(3);
%         
%     end
%     
    
    %================================================
    
    
    figure(fid)
    
    subplot(1,2,1)
   % subplot_tight(1,2,1,[0,0])
    
    % set(fid, 'Position', [100, 100, 1049, 895]);
    
    set(fid, 'Position', [100, 100, 2560, 1274]);
    %set(fid, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    %plotbrowser on
    title('midline')
    %hold on
    plot3(rx,ry,rz,'color',col{1},'LineWidth',2)
    
%     if(t(p1)==0||t(p1)==0.5||t(p1)==1)
%         plot3(rx,ry,rz,'r','LineWidth',2)
%         
%     end
%     
     set(gca, 'DataAspectRatio',[1,1,1]);
     daspect([1,1,1]);
    axis off
    set(gcf,'color',colw);
    
    % xlim auto
    % ylim auto
    % zlim auto
    lm = 0.25;
    % xlim([min(rx) max(rx)])
    %ylim([min(ry) max(ry)])
    %zlim([min(rz) max(rz)])
    
%     xlim([-lm lm])
%     ylim([-lm lm])
%     zlim([-lm lm])
   % axis vis3d
     %view(normal)
   view([1,1,1])
    
    
       hold off
%     
   subplot(1,2,2)
  %   subplot_tight(1,2,2,[0,0])
     plot(p1,E(p1),'--o','color','w','LineWidth',.5)
%     if(t(p1)==0||t(p1)==0.5||t(p1)==1)
%         plot(t(p1),Et(p1),'--or', 'LineWidth',.5)
%         
%     end
    xlim([1 500])
    ylim([.5 22])
    
    set(gca,'FontSize',25,'LineWidth',.5)
    set(gcf,'color', colw );
    set(gca,'XColor','w');
    %
    set(gca,'color',colw);
    set(gca,'XColor','w');
    set(gca,'YColor','w');%hold on
    %=== labels ===
    xlabel('$t$','FontSize',25,'Interpreter','latex');
    ylabel('$F$','FontSize',25,'Interpreter','latex');
    title('Bending energy for $\sigma=0.01$','FontSize',25,'Interpreter','latex','color','w');

    
    grid on
    hold on
    
    
    
    pause(.05);
    
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
    
    
    
    
end
close(writerObj);

