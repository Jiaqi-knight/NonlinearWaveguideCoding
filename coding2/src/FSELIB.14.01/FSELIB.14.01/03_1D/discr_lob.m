function [xe,xen,xien,xg,c,ng] = discr_lob (x1,x2,ne,ratio,np)

%========================================
% Discretize the interval (x1, x2)
% into ne elements,
% and generate the np-order
% Lobatto element interpolation nodes
%
% xe:  element end-nodes
% xen: element interpolation nodes
% xg:  unique global nodes
% c:   connectivity matrix
% ng:  number of unique global nodes
%========================================

%-----------------------------
% define the element end-nodes
%-----------------------------

xe = elm_line1(x1,x2,ne,ratio);

%---------------------------------------------
% define the element interpolation nodes (xen)
%---------------------------------------------

for l=1:ne

  m = np(l);

  xi(1) = -1.0;
  if(m>1)
   [tL, wL] = lobatto(m-1);
   for j=2:m
    xi(j) = tL(j-1);
   end
  end
  xi(m+1) = 1.0;

  mx = 0.5*(xe(l+1)+xe(l));
  dx = 0.5*(xe(l+1)-xe(l));

  for j=1:m+1
    xen(l,j) = mx + xi(j)*dx;
  end

  for j=1:m+1
    xien(l,j) =  xi(j);
  end

end

%-------------------------------
% define the connectivity matrix
% and the global nodes
%-------------------------------

Ic = 2;      % global node counter

for l=1:ne
  Ic = Ic-1;
  for j=1:np(l)+1
    c(l,j) = Ic;
    xg(Ic) = xen(l,j);
    Ic = Ic+1;
  end
end

%----------------------
% total number of nodes
%----------------------

ng = Ic-1;

%-----
% done
%-----

return;
