function emm = emm_t4  (x1,y1,z1,x2,y2,z2 ...
                       ,x3,y3,z3,x4,y4,z4);

%======================================
% Evaluation of the element mass matrix
% for a 4-node tetrahedron
%======================================

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
                                                                                
Jac(1,1)=x21; Jac(1,2)=x31; Jac(1,3)=x41;
Jac(2,1)=y21; Jac(2,2)=y31; Jac(2,3)=y41;
Jac(3,1)=z21; Jac(3,2)=z31; Jac(3,3)=z41;

%-------
% volume
%-------
                                                                                
V = det(Jac)/6.0;

%--------------------
% element mass matrix
%--------------------

fc = V/20;
fo = 2.0*fc;

emm(1,1) = fo; emm(1,2) = fc; emm(1,3) = fc; emm(1,4) = fc;
emm(2,1) = fc; emm(2,2) = fo; emm(2,3) = fc; emm(2,4) = fc;
emm(3,1) = fc; emm(3,2) = fc; emm(3,3) = fo; emm(3,4) = fc;
emm(4,1) = fc; emm(4,2) = fc; emm(4,3) = fc; emm(4,4) = fo;
                                                                                
%-----
% done
%-----
                                                                                
return;
