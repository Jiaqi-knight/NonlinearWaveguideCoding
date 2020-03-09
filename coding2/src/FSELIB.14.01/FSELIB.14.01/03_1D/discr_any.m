function [xe,xen,xien,xg,c,ng] = discr_any (x1,x2,ne,ratio,np,idistr)

%========================================
% Discretize the interval (x1, x2) into "ne" elements
% and generate the np-order element interpolation nodes
%
% xe:   element end-nodes
% xen:  element x interpolation nodes
% xien: element xi interpolation nodes
% xg:   unique global nodes
% c:    connectivity matrix
% ng:   number of unique global nodes
%========================================

%-----------------------------
% define the element end-nodes
%-----------------------------

xe = elm_line1 (x1,x2,ne,ratio);

%---------------------------------------------
% define the element interpolation nodes (xen)
%---------------------------------------------

for l=1:ne

  m = np(l);
  
  if(idistr==1)  % uniform

    for j=1:m+1
       xi(j) = -1.0+2.0*(j-1)/m;
    end

  elseif(idistr==2)  % lobatto

    xi(1) = -1.0;

    if(m >1)
     [tL, wL] = lobatto(m-1);
     for j=2:m
      xi(j) = tL(j-1);
     end
    end

    xi(m+1) = 1.0;

  else

    error('unknown option in discr_any')

  end

  for j=1:m+1
    xix(j) = 0.5*(xe(l+1)+xe(l)) + 0.5*(xe(l+1)-xe(l)) * xi(j);
  end

  for j=1:m+1
    xien(l,j) =  xi(j);
     xen(l,j) = xix(j);
  end

end

%-------------------------------
% define the connectivity matrix
% and the global nodes
%-------------------------------

Ic = 2;      % global node counter

for l=1:ne

  Ic = Ic-1;
  for je=1:np(l)+1
    c(l,je) = Ic;
    xg(Ic) = xen(l,je);
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
