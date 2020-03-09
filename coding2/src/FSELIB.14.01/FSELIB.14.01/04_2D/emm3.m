function [emm, area] = emm3 (x1,y1,x2,y2,x3,y3)

%======================================
% Evaluation of the element mass matrix
% from the coordinates of the three
% vertices of a 3-node triangle.
%
% area: triangle area
%======================================

 d23x = x2-x3; d23y = y2-y3;
 d31x = x3-x1; d31y = y3-y1;
 d12x = x1-x2; d12y = y1-y2;

 area = 0.5*(d31x*d12y - d31y*d12x);
 fc = area/12;

 emm(1,1) = fc*2; emm(1,2) = fc;   emm(1,3) = fc;
 emm(2,1) = fc;   emm(2,2) = fc*2; emm(2,3) = fc;
 emm(3,1) = fc;   emm(3,2) = fc;   emm(3,3) = fc*2;

%---
% done
%---

return;
