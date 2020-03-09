function ces = ces3(ne,ng,c)
                                                                                
%===================================
% Evaluation of the element-to-sides
% connectivity matrix ces
%===================================
                                                                                
%--------------------
% wrap the first node
%--------------------

  for i=1:ne
    c(i,4) = c(i,1)
  end
                                                                                
%----------------------
% run over the elements
%----------------------
                                                                                
  for i=1:ne
   for j=1:3
                                                                                
    ces(i,j) = 0;

    for k=1:ne
      if(k ~= i)  % skip the self-element
        for l=1:3
          if(c(k,l)==c(i,j) &  c(k,l+1)==c(i,j+1) )
           ces(i,j) = k;
          end
          if(c(k,l)==c(i,j+1) &  c(k,l+1)==c(i,j) )
           ces(i,j) = k;
          end
        end
      end
    end
                                                                                
   end
  end
return

%-----
% done
%-----
