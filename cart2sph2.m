function [th,ph] = cart2sph2(x,y,z)

%=== cartesian to spherical 


% if (z<0)
%     th = pi-acos(-z);
% else
%     th = acos(z);
% end

 th = acos(z);
if(x>0&&y>0)
   ph =  atan(y/x);
elseif(x<0&&y>0)
   ph = pi - atan(-y/x);
elseif(x<0&&y<0)
   ph = pi +atan(y/x);
else
    ph = pi + pi/2 + atan(-y/x);
    %  ph = -atan(-y/x);

  
end

    
if(x==0)
    if(y>0)
        ph = pi/2;
    else
        ph = pi+pi/2;
    end
end

if(y==0)
    if(x>0)
        ph = 0;
    else
        ph = pi;
    end
end
    
end
 