 function [indi,indj] = intersection_2curve(rx,ry,rx0,ry0)
%function [indi,indj] = intersection_2curve(rx,ry,rz,rx0,ry0,rz0)

sz = length(rx);


for i = 1:sz
  
        ind = (i-1)*sz + (1:sz);
        
      %  dr(ind) = sqrt((rx(i,1)-rx0(1:sz,1)).^2 + (ry(i,1)-ry0(1:sz,1)).^2 + (rz(i,1)-rz0(1:sz,1)).^2) ;
        dr(ind) = sqrt((rx(i,1)-rx0(1:sz,1)).^2 + (ry(i,1)-ry0(1:sz,1)).^2  ) ;

        
end

% %======   nearest point between the two curves =====
 
  [dr2,ind] = sort(dr);
% 

%--- first num_s points that are nearest to each other --
  num_s = 30;
  
  ind1 = ind(1:num_s);
  
  %--- i is index of rx and j is index of rx0---
  % from ind1, we extract values of i and j 
  
  % j index
  indj = rem(ind1,sz);
  
  indi = (ind1 -rem(ind1,sz))/sz + 1;
   
  for m1 = 1:length(indj)
    
    if(indj(m1)==0)
        indj(m1) = 1;
        indi(m1) =  indi(m1)-1;
    end
  end
  %==== Now, we eliminate the neighborhood points of the actual point of
  %intersections 
  
  %-- index of rx0 and rx in increasing order 
[indj2,indx]  = sort(indj);

indi2 = indi(indx);

tempi = indi2;
tempj = indj2;

ne_tol = 4;

for i = 1:length(indj2)-1
    
    %--- distance between rx at indi2(i) and rx0 at indj2(i) ------
   % dsi = sqrt((rx(indi2(i),1)-rx0(indj2(i),1))^2 + (ry(indi2(i),1)-ry0(indj2(i),1))^2 +  (rz(indi2(i),1)-rz0(indj2(i),1))^2);
    dsi = sqrt((rx(indi2(i),1)-rx0(indj2(i),1))^2 + (ry(indi2(i),1)-ry0(indj2(i),1))^2  );

    for j = (i+1):length(indj2)
        if(indj2(j)-indj2(i)<ne_tol )    % if the next point is in the neighborhood
            
            %-- comparing distance between the pair(indi2(i),indj2(i)) and (indi2(i),indj2(i))
        %    dsj = sqrt((rx(indi2(j),1)-rx0(indj2(j),1))^2 + (ry(indi2(j),1)-ry0(indj2(j),1))^2 +  (rz(indi2(j),1)-rz0(indj2(j),1))^2);
           dsj = sqrt((rx(indi2(j),1)-rx0(indj2(j),1))^2 + (ry(indi2(j),1)-ry0(indj2(j),1))^2 );
        
            if(dsi<=dsj)
%                 tempi(j) = tempi(i);
%                 tempj(j) = tempj(i);
                
                 tempi(i:j) = ones(j-i+1,1)*tempi(i);  % replacing index i:j with i:i
                tempj(i:j) = ones(j-i+1,1)*tempj(i);
            else
                tempi(i:j) = ones(j-i+1,1)*tempi(j);   % replacing index i:j with j:j
                tempj(i:j) = ones(j-i+1,1)*tempj(j);
            end
            
        end
    end
end
 
indi = unique(tempi);
indj = unique(tempj);

