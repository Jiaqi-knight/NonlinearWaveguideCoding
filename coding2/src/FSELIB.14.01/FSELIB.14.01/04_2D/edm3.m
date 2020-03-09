function [edm] = edm3 (x1,y1, x2,y2, x3,y3)

%===========================================
% Evaluation of the element diffusion matrix
% for a 3-node triangle
%===========================================

%--------
% prepare
%--------

d23x = x2-x3; d23y = y2-y3;
d31x = x3-x1; d31y = y3-y1;
d12x = x1-x2; d12y = y1-y2;

A = 0.5*(d31x*d12y-d31y*d12x);  % element area

A4 = 4.0*A;

%---------
% evaluate
%---------

edm(1,1) = (d23x*d23x+d23y*d23y)/A4;
edm(1,2) = (d23x*d31x+d23y*d31y)/A4;
edm(1,3) = (d23x*d12x+d23y*d12y)/A4;

edm(2,1) = (d31x*d23x+d31y*d23y)/A4;
edm(2,2) = (d31x*d31x+d31y*d31y)/A4;
edm(2,3) = (d31x*d12x+d31y*d12y)/A4;

edm(3,1) = (d12x*d23x+d12y*d23y)/A4;
edm(3,2) = (d12x*d31x+d12y*d31y)/A4;
edm(3,3) = (d12x*d12x+d12y*d12y)/A4;

%-----
% done
%-----

return;
