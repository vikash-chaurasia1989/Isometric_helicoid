function [f,k,m]=pendae(t,w,wd)
% Used by pendrun.m
% Tony Roberts, 21 July 1998, aroberts@usq.edu.au
x=w(1); y=w(2); u=w(3); v=w(4); s=w(5);
xd=wd(1); yd=wd(2); ud=wd(3); vd=wd(4); 
size(w)
size(wd)
f=[-xd+u
   -yd+v
   -ud-x*s
   -vd-y*s+1
       x^2+y^2-1];
k=[0 0 1 0 0
   0 0 0 1 0
   -s 0 0 0 -x
   0 -s 0 0 -y
   2*x 2*y 0 0 0];
m=-diag([1,1,1,1,0]);
