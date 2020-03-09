function cen = cen3(ne,ng,c)
                                                                                
%=================================
% Evaluation of the element-to-node
% connectivity matrix cen
%=================================

for i=1:ng
                                                                                
   cen(i,1) = 0;
   Nodeplace = 1;
                                                                                
   for j=1:ne
     for k=1:3
      if(c(j,k) == i)        % the ith node and the test node
        cen(i,1)=cen(i,1)+1;     % are identical
        Nodeplace = Nodeplace+1;
        cen(i,Nodeplace) = j;
      end;
     end
   end
                                                                                
end

%-----
% done
%-----
                                                                                
return
