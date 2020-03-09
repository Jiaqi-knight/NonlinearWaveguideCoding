function cen = cen6(ne,ng,c)

%===================================
% Evaluation of the element-to-node
% connectivity matrix: cen
%
% cen(i,1) is the number of elements sharing
% the ith global node, where i=1, 2, ..., N_G.
%
% cen(i,j), for j=2, 3, ..., ce(i,1)+1
% are the corresponding element labels.
%===================================

  for i=1:ng

   cen(i,1) = 0;
   Icount = 1;

   for j=1:ne
     for k=1:6
      if(c(j,k)==i)   %  the ith node and the test node are identical
        cen(i,1)=cen(i,1)+1; 
        Icount = Icount+1;
        cen(i,Icount) = j;
      end
     end
   end

  end

%-----
% done
%-----

return
