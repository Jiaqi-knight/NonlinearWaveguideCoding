function elm_mm = emm_modal_lob(m,h)

%==================================
% element mass matrix (emm) for the
% Lobatto modal expansion
%==================================

   for i=1:m+1
     for j=1:m+1
       elm_mm(i,j) = 0;
     end
   end

   for i=2:m
    elm_mm(i,i) = h*2*(i-1)^2*i^2/((2*i-3)*(2*i-1)*(2*i+1));
   end
   for i=2:m-1
    elm_mm(i,i+2) = - h*(i-1)*i*(i+1)*(i+2)/((2*i-1)*(2*i+1)*(2*i+3));
   end
   for i=3:m
    elm_mm(i,i-2) = - h* (i-3)*(i-2)*(i-1)*i/((2*i-5)*(2*i-3)*(2*i-1));
   end

   elm_mm(1,1) = h/3;
   elm_mm(1,2) = h/3;
   elm_mm(2,1) = h/3;
   elm_mm(1,3) = - h/5;
   elm_mm(3,1) = - h/5;

   elm_mm(1,m+1) = h/6;
   elm_mm(2,m+1) = h/3;
   elm_mm(3,m+1) = h/5;
   elm_mm(m+1,1) = h/6;
   elm_mm(m+1,2) = h/3;
   elm_mm(m+1,3) = h/5;
   elm_mm(m+1,m+1) = h/3;

%---
% done
%---

return
