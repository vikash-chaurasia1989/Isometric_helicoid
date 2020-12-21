clear all
clc


%=== In this file, we save the valley and peak values of the lower
%envelope.  ====


%valleys ==

valley = [ 1.289155039044352e+00     4.245388624242500e-01;
    2.387324146378430e+00     4.060921822279959e-01;
    3.421831276475750e+00     3.951753406387553e-01;
    4.440422912263880e+00     3.954620667861414e-01;
    5.427183559433631e+00     3.884348090486561e-01;
    6.445775195221762e+00     3.937379554008180e-01;
    7.432535842391513e+00     3.878716811985634e-01;
    8.435211983870452e+00     3.832641079440240e-01;
    9.406057136731015e+00     3.826964716265375e-01;
    1.040873327820996e+01     3.778989645809290e-01;
    1.136366293676133e+01     3.722611910979878e-01;
    1.236633907824027e+01     3.711484112899749e-01;
    1.330535324248245e+01     3.548818205885949e-01];


figure(1)
plotbrowser on
n = valley(:,1);
E = valley(:,2);

plot(n,E,'--ok')
set(gca,'FontSize',25,'LineWidth',.5)

xlim([min(valley(:,1)), max(valley(:,1))]);
ylim([min(valley(:,2)), max(valley(:,2))]);

%=== labels ===
xlabel('$n_m$','FontSize',25,'Interpreter','latex');
ylabel('$F_b$','FontSize',25,'Interpreter','latex');
grid on
title('Bending energy vs $n_m$ for $\sigma=0.01$','FontSize',25,'Interpreter','latex');





%===== peak points ====

peak =    [2.100845248813018e+00     1.219087130531306e+00;
    3.262676333383855e+00     1.156702066668591e+00;
    4.329014452099553e+00     1.128521891222190e+00;
    5.283944110650926e+00     1.117270206919960e+00;
    6.350282229366624e+00     1.079768368869583e+00;
    7.352958370845565e+00     1.100480613717072e+00;
    8.339719018015316e+00     1.079342437071986e+00;
    9.358310653803446e+00     1.077173737178754e+00;
    1.032915580666401e+01     1.053744141370867e+00;
    1.134774744245214e+01     1.053811969846629e+00;
    1.228676160669432e+01     1.024513831517825e+00;
    1.330535324248245e+01     1.027488610374204e+00; ];

figure(2)
plotbrowser on
n = peak(:,1);
E = peak(:,2);

plot(n,E,'--ok')
set(gca,'FontSize',25,'LineWidth',.5)

xlim([min(peak(:,1)), max(peak(:,1))]);
ylim([min(peak(:,2)), max(peak(:,2))]);
%=== labels ===
xlabel('$p_m$','FontSize',25,'Interpreter','latex');
ylabel('$F_b$','FontSize',25,'Interpreter','latex');
grid on
title('Bending energy vs $p_m$ for $\sigma=0.01$','FontSize',25,'Interpreter','latex');

 

%==== intersection values ====
inter  = [1.687042396774091e+00     1.067712756496647e+00;
    2.785211504108168e+00     9.493622856626490e-01;
    3.851549622823867e+00     9.393732938703463e-01;
    4.870141258611998e+00     9.238836344178447e-01;
    5.888732894400127e+00     9.217635588082717e-01;
    6.875493541569879e+00     8.936960462887049e-01;
    7.894085177358010e+00     8.994036862015496e-01;
    8.880845824527759e+00     8.838838503601212e-01;
    9.883521966006702e+00     8.894946644450717e-01;
    1.085436711886726e+01     8.618977740531107e-01;
    1.185704326034620e+01     8.670164808859114e-01;
    1.281197291889758e+01     8.396116448533957e-01; ];


figure(3)
plotbrowser on
n = inter(:,1);
E = inter(:,2);

plot(n,E,'--ok')
set(gca,'FontSize',25,'LineWidth',.5)

xlim([min(inter(:,1)), max(inter(:,1))]);
ylim([min(inter(:,2)), max(inter(:,2))]);     

%=== labels ===
xlabel('$i_m$','FontSize',25,'Interpreter','latex');
ylabel('$F_b$','FontSize',25,'Interpreter','latex');
grid on
title('Bending energy vs $i_m$ for $\sigma=0.01$','FontSize',25,'Interpreter','latex');

 


%==== Table for supplementary information ==========
%
%   3  &   1.289155039044352e+00  &   1.687042396774091e+00   &  2.100845248813018e+00 \\
%   5  &   2.387324146378430e+00  &   2.785211504108168e+00   &  3.262676333383855e+00 \\
%   7  &   3.421831276475750e+00  &   3.851549622823867e+00   &  4.329014452099553e+00 \\
%   9  &   4.440422912263880e+00  &   4.870141258611998e+00   &  5.283944110650926e+00 \\
%   11 &   5.427183559433631e+00  &   5.888732894400127e+00   &  6.350282229366624e+00 \\
%   13 &   6.445775195221762e+00  &   6.875493541569879e+00   &  7.352958370845565e+00 \\
%   15 &   7.432535842391513e+00  &   7.894085177358010e+00   &  8.339719018015316e+00 \\
%   17 &   8.435211983870452e+00  &   8.880845824527759e+00   &  9.358310653803446e+00 \\
%   19 &   9.406057136731015e+00  &   9.883521966006702e+00   &  1.032915580666401e+01 \\
%   21 &   1.040873327820996e+01  &   1.085436711886726e+01   &  1.134774744245214e+01 \\
%   23 &   1.136366293676133e+01  &   1.185704326034620e+01   &  1.228676160669432e+01 \\
%   25 &   1.236633907824027e+01  &   1.281197291889758e+01   &  1.330535324248245e+01 \\
%   27 &   1.330535324248245e+01  &           --              &         --
%

