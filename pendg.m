function pendg(t,w,wd)
% Used by pendrun.m
% Tony Roberts, 24 July 1998, aroberts@usq.edu.au
x=w(1); y=w(2); u=w(3); v=w(4); s=w(5);
plot([0 -x],[0 -y],-x,-y,'o',-x-[0 u]/8,-y-[0 v]/8)
axis([-1 1 -1 1])
axis('equal')
axis('off')
drawnow
