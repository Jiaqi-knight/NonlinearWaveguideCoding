function vdm = vdm_modal_lob(m,xi)

%===============================
% compute the vandermonde matrix
% for a given nodal set xi
% for the Lobatto modal expansion
%===============================

  vdm = zeros(m+1:m+1);

  vdm(1,1)=1;
  vdm(m+1,m+1)=1;

  if(m==1)
   return
  end

%---
   for j=2:m
%---

   x = xi(j);

   vdm(1,j)   = (1-x)/2;
   vdm(2,j)   =  1-x^2;
   vdm(m+1,j) = (1+x)/2;
   if(m>2)
     vdm(3,j) = 3*(1-x^2)*x;
   end
   if(m>3)
     vdm(4,j)=3*(1-x^2)*(5*x^2-1)/2;
   end
   if(m>4)
     vdm(5,j)=5*(1-x^2)*(7*x^2-3)*x/2;
   end
   if(m>5)
     vdm(6,j)=15*(1-x^2)*(21*x^4-14*x^2+1)/8;
   end
   if(m>6)
     vdm(7,j)=(1-x^2)*(693*x^4-630*x^2+105)*x/8;
   end
   if(m>7)
     vdm(8,j)=(1-x^2)*(3003*x^6-3465*x^4+945*x^2-35)/16;
   if(m>8)
     disp(' -->');
     error(' vdm_modal_lob: Sorry this high order not yet implemented');
   end
  end

%---
   end % over j
%---

%---
% done
%---

return
