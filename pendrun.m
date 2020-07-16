% run the differential algebraic equation solver dae4, dae2 or dae4o
% on the 5 equation model for the dynamics of a pendulum:
% dx/dt=u
% dy/dt=v
% du/dt=-x*s
% dv/dt=-y*s+1
% 0=x^2+y^2-1
% where (x,y) is the location of the bob, (u,v) is its velocity
% and s is proportional to tension in the wire of the pendulum.
% Uses penddae.m for the coded DAEs, peng.m for some graphics,
% dae4o.m, dae4.m or dae2.m to do the integration.
%
% Tony Roberts, 1 August 1998, aroberts@usq.edu.au 
format compact;
w0=[1;0;0;0;0];
ts=linspace(0,4,21);
w=dae2('penddae',ts,w0,2,'pendg');
% w=dae4('penddae',ts,w0,1,'pendg');
% w=dae4o('penddae',ts,w0,1,'pendg');
plot(ts,w);
legend(' x',' y',' u',' v',' s')
xlabel(' t')
