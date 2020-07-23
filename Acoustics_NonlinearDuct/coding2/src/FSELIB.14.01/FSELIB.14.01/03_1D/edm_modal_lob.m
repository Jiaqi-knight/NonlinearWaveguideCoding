function elm_dm =  edm_modal_lob(m,h)

%---------------------------------
% element difussion matrix for the
% Lobatto modal expansion
%---------------------------------

   for i=1:m+1
     for j=1:m+1
       elm_dm(i,j)= 0;
     end
   end

   for i=2:m
    elm_dm(i,i) = 4*(i-1)^2 * i^2/((2*i-1)*h);
   end

   elm_dm(1,1)  = 1/h; elm_dm(1,m+1)  =-1/h;
   elm_dm(m+1,1)=-1/h; elm_dm(m+1,m+1)= 1/h;

%---
% done
%---

return
