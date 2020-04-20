function a = plot_sphere()
 
SPHERE_RES = 100; % resolution for your spheres
SPHERE_RAD = 1.0; % radius of spheres
 % Make base sphere data
 [xb, yb, zb] = sphere(SPHERE_RES);
% Make a figure
 
% Plot spheres
 
     surf(SPHERE_RAD*xb , SPHERE_RAD*yb , SPHERE_RAD*zb , 'facecolor', 'Y', 'edgealpha', 0,'FaceAlpha',.3);
 
 
% Make sure they're smooth and shaded
light;
 lighting gouraud;
 az = 45;
 al = 20;
 view(az,al)
  

set(0,'defaultlinelinewidth',4)
set(gca,'Fontsize',30)
set(gca,'XMinorTick','on','YMinorTick','on')
xlhand = get(gca,'xlabel');
set(xlhand,'string','$x$ ','Interpreter','LaTex','fontsize',30)
ylhand = get(gca,'ylabel');
set(ylhand,'string','$y$','Interpreter','LaTex','fontsize',30)
zlhand = get(gca,'zlabel');
set(zlhand,'string','$z$','Interpreter','LaTex','fontsize',30)

title('')
 
%ylim([0 max(y)])
%grid on
a = 1;
 %axis square
 %axis off
end