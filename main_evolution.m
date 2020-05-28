%===  In this file, we evolve from one binormal curve to another
clear all
clc
global p1 tau f_b p1 d p q N1 h

global bx by bz bxp byp bzp bx2p by2p bz2p N rx ry rz Lk err del E fac sig tau fac

col1 = [0, 0.976, 0.968]; % -- color for mode n = 2
col2 = [0.831, 0.545, 0.247]; % -- color for mode n = 3
col3 = [0.341, 0.505, 0.819]; % -- color for mode n = 4
col4 = [0.705, 0.701, 0.070]; % -- color for mode n = 5
col5 = [0.301, 0.811, 0.498]; % -- color for mode n = 6
col6 = [1, 0.929, 0.278];
col = {col1, col2, col3, col4, col5, col6};
colw = [41, 44, 52] / 255;

path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution/';

%N = 72;
tau = 16;

N = 105;
h = 1 / N
N1 = N - 1;
sig = 0.01;

tau2 = tau;
 %===
n1 = tau/2/pi;
fac = (asinh(n1*pi*sig)/(4*pi^3*n1^3))  ;


str1 = ['5pi_knot_N' num2str(N) '_tau_' num2str(10^10 * tau2) '.txt'];
str2 = ['3fold_N' num2str(N) '_tau_' num2str(10^10 * tau2) '.txt'];
str3 = ['7pi_knot_N' num2str(N) '_tau_' num2str(10^10 * tau2) '.txt'];
strb1 = [path 'b_' str1];
strb2 = [path 'b_' str2];
strb3 = [path 'b_' str3];
%==== load b ===========
temp = load(strb1);

bx1 = temp(:, 1);
by1 = temp(:, 2);
bz1 = temp(:, 3);

%==== load b ===========
temp = load(strb2);

bx2 = temp(:, 1);
by2 = temp(:, 2);
bz2 = temp(:, 3);

%==== load b ===========
temp = load(strb3);

bx3 = temp(:, 1);
by3 = temp(:, 2);
bz3 = temp(:, 3);

%=== Isotopy between 1 and 2 ====
%==== Isotopy =====
t1 = linspace(-0.05, 0, 20);
t2 = linspace(0, .5, 50);
t3 = linspace(.5, 1, 50);
t4 = linspace(1, 1.05, 20);

t = [t1 t2 t3 t4];

t = -.1:.005:1.1;

E1 = energy_b(bx1, by1, bz1);
E2 = energy_b(bx2, by2, bz2);
E3 = energy_b(bx3, by3, bz3);

%Et = (1-t).*(0.5-t)/0.5*energy_b(bx1,by1,bz1)  + (1-t).*t/0.25*energy_b(bx2,by2,bz2) - t.*(0.5-t)/0.5*energy_b(bx3,by3,bz3);

%=== Creating energy profile
x = 0.5;

f1 = x^5/5 - 3/8 * x^4 + x^3/6;

f2 = x^4/4 - x^3/2 + x^2/4;

x = 1;
f3 = x^5/5 - 3/8 * x^4 + x^3/6;

f4 = x^4/4 - x^3/2 + x^2/4;

B = (E3 - E1) / (E2 - E1);

d = (f3 - B * f1) / (f4 - B * f2);

A = (E2 - E1) / (f1 - d * f2);

%Et = A*(t.^5/5 - 3*t.^4/8 + t.^3/6 -t.^2/4  -d*(t.^4/4-t.^3/2 +t.^2/4-t/2)) + E1;
Et = A/2*(2/5 * t.^5 -(3 + 2 * d) / 4 * t.^4 + (3 * d + 1) / 3 * t.^3 -d / 2 * t.^2) + E1;

 
sv = 1;

t1 = -.1:.02:1.1;

for p1 =22:22%113:length(t)%:110%-1:70%-1:75%111%:61%steps

    bx = (1 - t(p1)) * (0.5 - t(p1)) / 0.5 * bx1 + (1 - t(p1)) * t(p1) / 0.25 * bx2 - t(p1) * (0.5 - t(p1)) / 0.5 * bx3;
    by = (1 - t(p1)) * (0.5 - t(p1)) / 0.5 * by1 + (1 - t(p1)) * t(p1) / 0.25 * by2 - t(p1) * (0.5 - t(p1)) / 0.5 * by3;
    bz = (1 - t(p1)) * (0.5 - t(p1)) / 0.5 * bz1 + (1 - t(p1)) * t(p1) / 0.25 * bz2 - t(p1) * (0.5 - t(p1)) / 0.5 * bz3;

      %  bx = bx3;
       % by = by3;
      %  bz = bz3;
    
        
        
        %== find nearest value in the old t 
        
        [val,ind] = min(abs(t1-t(p1)));
        
        t2 = t1(ind(1));
       
   str0 = ['branch_213_N_' num2str(N) '_t_' num2str(round(1000*t2)) '_tau_' num2str(round(10^10*tau)) '.txt'];
   
    str0 = ['branch_213_N_' num2str(N) '_t_' num2str(round(1000*t(p1+1))) '_tau_' num2str(round(10^10*tau)) '.txt'];

      
    path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution2/';
  %  str0 = ['branch_213_N_' num2str(N) '_step_' num2str(1) '_tau_' num2str(round(10^10 * tau)) '.txt'];

 
       %== for evolution ==
    strb = ['b_' str0];
    
    temp = load([path strb]);
    
    bx = temp(:,1);
    by = temp(:,2);
    bz = temp(:,3);
%     

    E = Et(p1);
    %var_in = initial_evolution();
    % var_in = initial_evolution3();
    var_in(1:3:3 * N1 - 2, 1) = bx(2:N, 1);
    var_in(2:3:3 * N1 - 1, 1) = by(2:N, 1);
    var_in(3:3:3 * N1, 1)     = bz(2:N, 1);

    options = optimset('Display', 'iter',  'Algorithm', 'levenberg-marquardt', 'Jacobian', 'on',  'TOlFun', 10^(-10), 'TOlX', 10^(-10), 'MaxFunEvals', 69500);
    options.MaxIter = 50000;
    %     [x,fval,exitflag,output,qd1] =  fsolve(@fun_curve2 ,var_initial,options)          ;
    [x, fval, exitflag, output, qd1] = fsolve(@fun_evolution, var_in, options);

    err1(p1) = err;

    %----------------- Post processing -----
    h = 1 / N; % mean(len);
    %---- Tangent ti = bi \times bi+1

    i = 1:N + 1;

    %
    tx(i, 1) = by(i, 1) .* bzp(i, 1) - bz(i) .* byp(i, 1);
    ty(i, 1) = bz(i, 1) .* bxp(i, 1) - bx(i) .* bzp(i, 1);
    tz(i, 1) = bx(i, 1) .* byp(i, 1) - by(i) .* bxp(i, 1);

    %
    %    tx(N+1,1) = tx(1,1);
    %    ty(N+1,1) = ty(1,1);
    %    tz(N+1,1) = tz(1,1);
    % % %
    % %
    tx = tx / tau;
    ty = ty / tau;
    tz = tz / tau;

    %--------------------

    % %--- position vector using integration of tangent
    % % initialization
    %h=1;
    rx(1) = 0; ry(1) = 0; rz(1) = 0;

 

    for i = 1:N

        rx(i + 1) = rx(i) + h * tx(i + 1);
        ry(i + 1) = ry(i) + h * ty(i + 1);
        rz(i + 1) = rz(i) + h * tz(i + 1);

    end

    rx = rx - mean(rx);
    ry = ry - mean(ry);
    rz = rz - mean(rz);

    %== rescaling the midline

    l = sum(sqrt((rx(1:end - 1) - rx(2:end)).^2 + (ry(1:end - 1) - ry(2:end)).^2 + (rz(1:end - 1) - rz(2:end)).^2));

    rx = rx / l;
    ry = ry / l;
    rz = rz / l;

    % %
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

    %==== saving data =====

    str0 = ['branch_213_N_' num2str(N) '_step_' num2str(p1) '_tau_' num2str(10^10 * tau) '.txt'];
    
   str0 = ['branch_213_N_' num2str(N) '_t_' num2str(round(1000*t(p1))) '_tau_' num2str(round(10^10*tau)) '.txt'];
 path =  '/Users/vikashchaurasia/OneDrive/Vikash_Documents/Isometric_deformation/Matlab_files/fixed_rotation_final/data_evolution2/';

    %== for evolution ==
    strb = ['b_' str0];
    

    if (sv == 1)

        %
        fileID = fopen([path strb], 'w');
        fprintf(fileID, '%30.16E %30.16E %30.16E \r\n', [bx'; by'; bz']);
        fclose(fileID);

    end

end


if(sv==1)
strt  = ['t_branch_' num2str(213) '_tau_' num2str(round(10^10*tau)) '.txt'];

fileID = fopen([path strt],'w');
fprintf(fileID,'%30.16E   \r\n',t' ); 
fclose(fileID);
end

