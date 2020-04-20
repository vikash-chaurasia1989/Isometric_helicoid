function [indi,indj] = intersection_1curve(rx,ry)

sz = length(rx) ;


for i = 1:sz
  
        
        ind = (i-1)*sz + (1:sz);
        
        dr(ind) = sqrt((rx(i,1)-rx(1:sz,1)).^2 + (ry(i,1)-ry(1:sz,1)).^2 ) ;
        
        tempind = find(dr==0);
        
        dr(tempind) = 100;
        
end

% %======   nearest point between the two curves =====
 
  [dr2,ind] = sort(dr);
% 

%--- first num_s points that are nearest to each other --
  num_s = 20;
  
  ind1 = ind(1:num_s);
  
  %--- i is index of rx and j is index of rx---
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
  
  %-- index of rx and rx in increasing order 
[indj2,indx]  = sort(indj);

indi2 = indi(indx);

tempi = indi2;
tempj = indj2;

ne_tol = 8;

for i = 1:length(indj2)-1
    
    %--- distance between rx at indi2(i) and rx at indj2(i) ------
    dsi = sqrt((rx(indi2(i),1)-rx(indj2(i),1))^2 + (ry(indi2(i),1)-ry(indj2(i),1))^2);

    for j = (i+1):length(indj2)
        if(indj2(j)-indj2(i)<ne_tol )    % if the next point is in the neighborhood
            
            %-- comparing distance between the pair(indi2(i),indj2(i)) and (indi2(i),indj2(i))
            dsj = sqrt((rx(indi2(j),1)-rx(indj2(j),1))^2 + (ry(indi2(j),1)-ry(indj2(j),1))^2);
            
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

