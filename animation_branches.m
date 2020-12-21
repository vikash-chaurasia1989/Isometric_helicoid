%---
clear all
clc
global    lm  N1 p1   tau     h tx ty tz   fac

global bx by bz bxp byp bzp   N  rx ry rz Lk err    branch sig   ht beta gamma x2 nstep path sig
format longE



%---- read r and b file and create ruling data ----

%---- read r and b file and create ruling data ----
col1 = [0,     0.976, 0.968]  ;      % -- color for m = 3
col2 = [0.831, 0.545, 0.247]  ;      % -- color for m = 5
col3 = [0.341, 0.505, 0.819]  ;      % -- color for m = 7
col4 = [0.705, 0.701, 0.070]  ;      % -- color for m = 9
col5 = [0.301, 0.811, 0.498]  ;      % -- color for m = 11

clr = linspecer(15);

col6  = [0.98, 0.325, 0];%[0.1705, 0.2525,  0.4095]   ;    % -- color for m = 13
col7  = col2/2;% [2.941176470588235e-01     5.447058823529412e-01     7.494117647058823e-01]  ;    % -- color for m = 15 
col8  = [3.717647058823529e-01     7.176470588235294e-01     3.611764705882353e-01]/2 ;    % -- color for m = 17
col9  = [ 1.000000000000000e+00     5.482352941176470e-01     9.999999999999998e-02]  ;    % -- color for m = 19
col10 = col1/1.5;%[8.650000000000000e-01     8.109999999999999e-01     4.330000000000000e-01]   ;    % -- color for m = 21
col11 = [0.533, 0.671, 0.525];%col2/1.5;%[6.858823529411765e-01     4.035294117647059e-01     2.411764705882353e-01]   ;    % -- color for m = 23
col12 = [9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]   ;    % -- color for m = 25
col13 = [0.549, 0.035, 0.486];%[0.196, 0.522, 0.788];%col3/1.5;%[9.717647058823530e-01     5.552941176470587e-01     7.741176470588236e-01]/2 ;    % -- color for m = 27
col14 = [0.996, 0.773, 0.243];%col5/1.5;%[     4.350351799515555e-01     7.800722978266937e-01     6.461779788838612e-01]  ;  % -- color for m = 29                                                                        % -- color for m = 29
col15 = [0.11, 0.18, 1];%[0.988, 0.427, 0.412];%col6/1.5;%[     4.350351799515555e-01     7.800722978266937e-01     6.461779788838612e-01]/3  ;  % -- color for m = 31                                                                      % -- color for m = 29




col = {col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12, col13,col14,col15} ;


% col1 = [1,1,1];
% col = {col1,col1,col1,col1,col1,col1,col1,col1,col1,col1,col1,col1, col1,col1,col1} ;

sig = 0.01;

colw = [41,44,52]/255;

%===== Loading saved energy =====
str = 'Energy_branches.txt';
path = '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch1/';

E1 = load([path str]);
E1 = E1';

%=== loading stability results ===
fid = figure;
figure(fid)

subplot(1,2,1)
set(fid, 'Position', [100, 100, 2560, 1274]);
hold on 

h1 =subplot(1,2,2)
strfig = 'energy_for_animation.fig';

addpath('/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch1/')

g = openfig(strfig,'invisible');
 
copyobj(allchild(get(g,'CurrentAxes')),h1);

    

branch =12;
nfold = 10;
path = ['/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_branch' num2str(branch) '/'];

strtau = ['tau_branch_' num2str(branch) '.txt'];
strLk = ['Linking_branch_' num2str(branch) '.txt'];
tau1 = load([path strtau]);
Lk = load([path strLk]);

if (branch==1||branch==2||branch==3)
    N=105;
end

if(branch==4||branch==9)
    N = 180; %
end
if(branch==5)
    N = 140;
end

if(branch== 6)
    N = 90;
end

if(branch== 7)
    N = 150;
end


if(branch== 8)
    N = 150;
end

if(branch== 10)
    N = 220;
end

if(branch==11||branch==14||branch==15)
    N=170;
end

if(branch== 12)
    N = 190;
end


if(branch==13)
    N=210;
end
        

N1 = N-1;
h = 1/N;

     set(gca,'FontSize',25,'LineWidth',.5)
    set(gcf,'color', colw );
    set(gca,'XColor','w');
    
    set(gca,'color',colw);
    set(gca,'XColor','w');
    set(gca,'YColor','w');%hold on
   % === labels ===
    xlabel('$n$','FontSize',25,'Interpreter','latex');
    ylabel('$F_b$','FontSize',25,'Interpreter','latex');
    title(['Bending energy for $\sigma=0.01$ and branch ' num2str(branch)],'FontSize',25,'Interpreter','latex','color','w');
    box on
    
    grid on
    hold on    
    
    
    
    
    
    


str =   ['animation_branch_' num2str(branch)];


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

 

for p1=  1:length(tau1)
    
    tau = tau1(p1);
    str0 = ['branch_' num2str(branch) '_N' num2str(N) '_tau_' num2str(round(10^10*tau)) '.txt'];
    
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
        for i = 1:N+1
    
            temp =  Ru*[rx(i);ry(i);rz(i)] ;
    
            rx(i) = temp(1);
            ry(i) = temp(2);
            rz(i) = temp(3);
    
        end
    
    
    %================================================
    Lk1 = Lk(p1);
    icol =  round((Lk1-1)/2);
    
    figure(fid)
    subplot(1,2,1)

    %plotbrowser on
    title('midline')
    %hold on
    
    
    plot3(rx,ry,rz,'color',col{icol},'LineWidth',2)
    
  
    set(gca, 'DataAspectRatio',[1,1,1]);
    daspect([1,1,1]);
    axis off
    set(gcf,'color',colw);
    
    % xlim auto
    % ylim auto
    % zlim auto
    
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
    %plot(tau1(p1)/2/pi,E(branch,p1),'--o','color',col{icol},'LineWidth',.5)
        plot(tau1(p1)/2/pi,E(p1),'--o','color',col{icol},'LineWidth',.5)

    
    %     if(t(p1)==0||t(p1)==0.5||t(p1)==1)
    %         plot(t(p1),Et(p1),'--or', 'LineWidth',.5)
    %
    %     end
    xlim([8.4 85]/2/pi)
    ylim([.3 14])
    
%      
%     txt = {['m =' num2str(Lk1)]};
%     text(2,13,txt,'color', 'w' ,'FontSize',25)
%     
    
    
%     set(gca,'FontSize',25,'LineWidth',.5)
%     set(gcf,'color', colw );
%     set(gca,'XColor','w');
%     %
%     set(gca,'color',colw);
%     set(gca,'XColor','w');
%     set(gca,'YColor','w');%hold on
%     %=== labels ===
%     xlabel('$t$','FontSize',25,'Interpreter','latex');
%     ylabel('$F$','FontSize',25,'Interpreter','latex');
%     title('Bending energy for $\sigma=0.01$','FontSize',25,'Interpreter','latex','color','w');
%     
%     
%     grid on
    hold on
    
    
    
   
    
    pause(.05);
    
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
    
    
    
    
end
close(writerObj);

