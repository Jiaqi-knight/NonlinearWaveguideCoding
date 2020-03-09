function [ebsm, rhs, arel] = HCT_ebsm (x1,y1, x2,y2, x3,y3, NQ, w0)

%=================================================
% Evaluation of the element bending stiffness matrix
% and right-hand side of the global system
% for the HCT triangle, using a Gauss
% integration quadrature
%
% ebsm: element bending stiffness matrix
% rhs:  right-hand side
% arel: element area 
%===============================================

%---------------------------------
% compute the area of the triangle
%---------------------------------

d12x = x1-x2; d12y = y1-y2;
d31x = x3-x1; d31y = y3-y1;

arel = 0.5*(d31x*d12y - d31y*d12x);

%--------------
% interior node
%--------------

x7 = (x1+x2+x3)/3.0;
y7 = (y1+y2+y3)/3.0;

%-----------------------------
% read the triangle quadrature
%-----------------------------

[xi, eta, w] = gauss_trgl(NQ);

%---------
% initialize the element bending stiffness matrix
% and right-hand side
%---------

for k=1:12
 for l=1:12
  ebsm(k,l) = 0.0;
 end
 rhs(k) = 0.0;
end

%------------------------------
% compute the 30 cubic coefficients
% and put them in the vector: am
%------------------------------

for mode=1:12

  dof=zeros(1,12); dof(mode)=1.0;

  [a] = HCT_sys (x1,y1, x2,y2, x3,y3, dof);

  for i=1:30
    am(i,mode) = a(i);
  end

end

%-----------------------
% perform the quadrature
%-----------------------

for pass=1:3  % run over the 3 sub-triangles

 for i=1:NQ    % quadrature

%----
% modes and laplacian:
%---

%--
 for mode=1:12
 
   for j=1:30
     c(j) = am(j,mode);
   end

   if(pass==1)       % triangle 1-2-7

      x = x1 + (x2-x1)*xi(i) + (x7-x1)*eta(i);
      y = y1 + (y2-y1)*xi(i) + (y7-y1)*eta(i);

      psi(mode) = c(1) + c(2)*x   +c(3)*y ...
                       + c(4)*x^2 +c(5)*x*y   + c(6)*y^2 ...
                       + c(7)*x^3 +c(8)*x^2*y + c(9)*x*y^2+c(10)*y^3;
      lpsi(mode) = 2.0*c(4)+2.0*c(6) ...
               +6.0*c(7)*x + 2.0*c(8)*y + 2.0*c(9)*x + 6.0*c(10)*y;

     elseif(pass==2)  % triangle 2-3-7

      x = x2 + (x3-x2)*xi(i) + (x7-x2)*eta(i);
      y = y2 + (y3-y2)*xi(i) + (y7-y2)*eta(i);

      psi(mode) = c(11)+c(12)*x+c(13)*y+c(14)*x^2+c(15)*x*y+c(16)*y^2 ...
                 +c(17)*x^3+c(18)*x^2*y+c(19)*x*y^2+c(20)*y^3;

      lpsi(mode) = 2.0*c(14)+2.0*c(16) ...
               +6.0*c(17)*x+2.0*c(18)*y+2.0*c(19)*x+6.0*c(20)*y;

     else            % triangle 3-1-7

      x = x3 + (x1-x3)*xi(i) + (x7-x3)*eta(i);
      y = y3 + (y1-y3)*xi(i) + (y7-y3)*eta(i);

      psi(mode) = c(21)+c(22)*x+c(23)*y+c(24)*x^2+c(25)*x*y+c(26)*y^2 ...
                 +c(27)*x^3+c(28)*x^2*y+c(29)*x*y^2+c(30)*y^3;

      lpsi(mode) = 2.0*c(24)+2.0*c(26) ...
               +6.0*c(27)*x+2.0*c(28)*y+2.0*c(29)*x+6.0*c(30)*y;
     end

%    if(mode == 1) disp(lpsi(mode)); end
%    if(mode == 4) disp(lpsi(mode)); end
%    if(mode == 1) disp(psi(mode)); 
%       end
%    if(mode == 4) disp(psi(mode)); 
%     end

    end   % end of loop over modes
%--

   cf = arel*w(i)/3.0;

   for k=1:12
    for l=k:12
     ebsm(k,l) = ebsm(k,l) + lpsi(k)*lpsi(l)*cf;
    end
    wload = 0.10*(1.0-x^2-y^2) ;   % parabolic load
    wload = x-1.0;   % linear
    wload = w0;                    % uniform load
    rhs(k) = rhs(k) + wload*psi(k)*cf;
   end

 end   % end of quadrature

end % end of run over sub-triangles

%-------------------------------
% fill in lower diagonal of ebsm
%-------------------------------

 for k=1:12
  for l=1:k-1
   ebsm(k,l) = ebsm(l,k);
  end
 end

%-----
% done
%-----

return;
