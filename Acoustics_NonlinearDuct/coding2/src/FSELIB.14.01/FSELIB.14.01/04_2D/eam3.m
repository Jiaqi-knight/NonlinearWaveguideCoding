function [eam] = eam3 (ux,uy,x1,y1,x2,y2,x3,y3)

%===========================================
% Evaluation of the element advection matrix
% for a 3-node triangle.
%===========================================

 d32x = x3-x2; d32y = y3-y2;
 d13x = x1-x3; d13y = y1-y3;
 d21x = x2-x1; d21y = y2-y1;

 eam(1,1) = (-ux*d32y + uy*d32x)/6.0;
 eam(1,2) = (-ux*d13y + uy*d13x)/6.0;
 eam(1,3) = (-ux*d21y + uy*d21x)/6.0;

 eam(2,1) = eam(1,1);
 eam(2,2) = eam(1,2);
 eam(2,3) = eam(1,3);

 eam(3,1) = eam(1,1);
 eam(3,2) = eam(1,2);
 eam(3,3) = eam(1,3);

%-----
% done
%-----

return;
