function [edm, V] = edm_t4 (x1,y1,z1,x2,y2,z2 ...
                           ,x3,y3,z3,x4,y4,z4)

%===========================================
% Evaluation of the element diffusion matrix
% for a 4-node tetrahedron
%===========================================

%--------
% prepare
%--------
                                                                                
x21 = x2-x1; y21 = y2-y1; z21 = z2-z1;
x31 = x3-x1; y31 = y3-y1; z31 = z3-z1;
x41 = x4-x1; y41 = y4-y1; z41 = z4-z1;

x32 = x3-x2; y32 = y3-y2; z32 = z3-z2;
x42 = x4-x2; y42 = y4-y2; z42 = z4-z2;

%----------------
% Jacobian matrix
%----------------

Jac(1,1) = x21; Jac(1,2) = x31; Jac(1,3) = x41;
Jac(2,1) = y21; Jac(2,2) = y31; Jac(2,3) = y41;
Jac(3,1) = z21; Jac(3,2) = z31; Jac(3,3) = z41;

%-------
% volume
%-------
V = det(Jac)/6.0;

V6 = 6.0*V;

% first gradient

gx(1) = (-y32*z42 + z32*y42)/V6;
gy(1) = ( x32*z42 - z32*x42)/V6;
gz(1) = (-x32*y42 + y32*x42)/V6;

% second gradient
                                                                                
gx(2) = ( y31*z41 - z31*y41)/V6;
gy(2) = (-x31*z41 + z31*x41)/V6;
gz(2) = ( x31*y41 - y31*x41)/V6;

% third gradient
                                                                                
gx(3) = (-y21*z41 + z21*y41)/V6;
gy(3) = ( x21*z41 - z21*x41)/V6;
gz(3) = (-x21*y41 + y21*x41)/V6;

% fourth gradient

gx(4) = ( y21*z31 - z21*y31)/V6;
gy(4) = (-x21*z31 + z21*x31)/V6;
gz(4) = ( x21*y31 - y21*x31)/V6;

%=========
% alternative
%=========

ido=1;
ido=0;

if(ido==1)

 rhs=[-1,-1,-1];sln=rhs/Jac;
 gx(1)=sln(1); gy(1)=sln(2);gz(1)=sln(3);
 rhs=[1,0,0];sln=rhs/Jac;
 gx(2)=sln(1); gy(2)=sln(2);gz(2)=sln(3);
 rhs=[0,1,0];sln=rhs/Jac;
 gx(3)=sln(1); gy(3)=sln(2);gz(3)=sln(3);
 rhs=[0,0,1];sln=rhs/Jac;
 gx(4)=sln(1); gy(4)=sln(2);gz(4)=sln(3);

end

%---
% element diffusion matrix
%---

for i=1:4
 for j=1:4
  edm(i,j) = V*( gx(i)*gx(j) + gy(i)*gy(j) + gz(i)*gz(j) );
 end
end

%sum1 = edm(1,1)+edm(1,2)+edm(1,3)+edm(1,4);
%sum2 = edm(2,1)+edm(2,2)+edm(2,3)+edm(2,4);
%sum3 = edm(3,1)+edm(3,2)+edm(3,3)+edm(3,4);
%sum4 = edm(4,1)+edm(4,2)+edm(4,3)+edm(4,4);

%edm
%[sum1 sum2 sum3 sum4]
%det(edm)
%V
%pause

%-----
% done
%-----
                                                                                
return;
