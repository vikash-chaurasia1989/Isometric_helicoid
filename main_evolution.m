%===  In this file, we evolve from one binormal curve to another

global    p1   tau   f_b  p1   d p q

global bx by bz bxp byp bzp bx2p by2p bz2p   N     rx ry rz Lk err   del  E fac 


path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';



% steps
steps = 10 ;


N = 105;
tau = 16;
p = 1;
q = 2;


d1 = distance_curves(p,q,tau);

d = d1/steps; 
 
sv=0;
for p1 = 1:1%steps
    
    
    %var_in = initial_evolution();
    var_in = initial_evolution3();

    
    options             = optimset('Display','iter', 'Algorithm','levenberg-marquardt','Jacobian','on', 'TOlFun',10^(-10),'TOlX',10^(-10),'MaxFunEvals',69500  ) ;
    options.MaxIter     = 50000  ;
    %     [x,fval,exitflag,output,qd1] =  fsolve(@fun_curve2 ,var_initial,options)          ;
    [x,fval,exitflag,output,qd1] =  fsolve(@fun_evolution        ,var_in,options)          ;
    
    
    
    
    
   %----------------- Post processing -----
h =1/N;% mean(len);
%---- Tangent ti = bi \times bi+1

i = 1:N+1   ;

 
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

%--------------------

% %--- position vector using integration of tangent
% % initialization
%h=1;
rx(1) = 0; ry(1) = 0; rz(1) = 0;

% i = 1 ;
%
%     rx(i+1) = rx(i) + h*tx(i);
%     ry(i+1) = ry(i) + h*ty(i);
%     rz(i+1) = rz(i) + h*tz(i);
%
% for i = 2:N-1
%     rx(i+1) = rx(i-1) + 2*h*tx(i);
%     ry(i+1) = ry(i-1) + 2*h*ty(i);
%     rz(i+1) = rz(i-1) + 2*h*tz(i);
% end
%
 

for i = 1:N

     rx(i+1) = rx(i) + h*tx(i+1);
     ry(i+1) = ry(i) + h*ty(i+1);
     rz(i+1) = rz(i) + h*tz(i+1);

end
rx = rx-mean(rx);
ry = ry-mean(ry);
rz = rz-mean(rz);

 

%== rescaling the midline 

l = sum(sqrt((rx(1:end-1)-rx(2:end)).^2 + (ry(1:end-1)-ry(2:end)).^2+(rz(1:end-1)-rz(2:end)).^2));

rx = rx/l;
ry = ry/l;
rz = rz/l;
 
% 
figure(2)
plotbrowser on
title('binormal')
plot_sphere()
hold on
plot3(bx,by,bz,'color',col{1},'LineWidth',2)
hold on
set(gca, 'DataAspectRatio',[1,1,1]);
daspect([1,1,1]);
axis off
set(gcf,'color',colw);
axis vis3d



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
 
    
    
    
    %==== saving data =====
    
   str0 = ['branch_' num2str(p) num2str(q) '_N_' num2str(N) '_step_' num2str(p1) '_tau_' num2str(10^10*tau) '.txt'];
    
   %== for evolution ==
strb = ['b_' str0];
strr = ['r_' str0];

strlm = ['lm_' str0];
strrho = ['rho_' str0];
struvw = ['uvw_' str0];


if(sv==1)

 % 
fileID = fopen([path strb],'w');
fprintf(fileID,'%30.16E %30.16E %30.16E \r\n',[bx'    ;by'    ; bz'    ] );
fclose(fileID);

 

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


 

    
    
    
end
    
    
     
    
    
    
    