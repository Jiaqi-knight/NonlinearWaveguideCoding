function[Fgeom]=computegeometricforcing3D(N,Q,secondChristoffelsymbol,scheme)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 30th, 2014
%    Last update: May 30th, 2014
%
%    Description: 
%          Input: number of lattice points N
%                 number of lattice velocities Q
%         Output: 

%%

Fgeom = zeros(N,3*Q);

xi1 = scheme(:,1);
xi2 = scheme(:,2);
xi3 = scheme(:,3);

Gamma111 = secondChristoffelsymbol(:,1);
Gamma112 = secondChristoffelsymbol(:,2);
Gamma113 = secondChristoffelsymbol(:,3);
Gamma121 = secondChristoffelsymbol(:,4);
Gamma122 = secondChristoffelsymbol(:,5);
Gamma123 = secondChristoffelsymbol(:,6);
Gamma131 = secondChristoffelsymbol(:,7);
Gamma132 = secondChristoffelsymbol(:,8);
Gamma133 = secondChristoffelsymbol(:,9);
Gamma211 = secondChristoffelsymbol(:,10);
Gamma212 = secondChristoffelsymbol(:,11);
Gamma213 = secondChristoffelsymbol(:,12);
Gamma221 = secondChristoffelsymbol(:,13);
Gamma222 = secondChristoffelsymbol(:,14);
Gamma223 = secondChristoffelsymbol(:,15);
Gamma231 = secondChristoffelsymbol(:,16);
Gamma232 = secondChristoffelsymbol(:,17);
Gamma233 = secondChristoffelsymbol(:,18);
Gamma311 = secondChristoffelsymbol(:,19);
Gamma312 = secondChristoffelsymbol(:,20);
Gamma313 = secondChristoffelsymbol(:,21);
Gamma321 = secondChristoffelsymbol(:,22);
Gamma322 = secondChristoffelsymbol(:,23);
Gamma323 = secondChristoffelsymbol(:,24);
Gamma331 = secondChristoffelsymbol(:,25);
Gamma332 = secondChristoffelsymbol(:,26);
Gamma333 = secondChristoffelsymbol(:,27);

for i=1:Q
    Fgeom(:,3*(i-1)+1) = -(Gamma111*xi1(i,1)*xi1(i,1) + Gamma112*xi1(i,1)*xi2(i,1) + Gamma113*xi1(i,1)*xi3(i,1) + Gamma121*xi2(i,1)*xi1(i,1) + Gamma122*xi2(i,1)*xi2(i,1) + Gamma123*xi2(i,1)*xi3(i,1) + Gamma131*xi3(i,1)*xi1(i,1) + Gamma132*xi3(i,1)*xi2(i,1) + Gamma133*xi3(i,1)*xi3(i,1));
    Fgeom(:,3*(i-1)+2) = -(Gamma211*xi1(i,1)*xi1(i,1) + Gamma212*xi1(i,1)*xi2(i,1) + Gamma213*xi1(i,1)*xi3(i,1) + Gamma221*xi2(i,1)*xi2(i,1) + Gamma222*xi2(i,1)*xi2(i,1) + Gamma223*xi2(i,1)*xi3(i,1) + Gamma231*xi3(i,1)*xi1(i,1) + Gamma232*xi3(i,1)*xi2(i,1) + Gamma233*xi3(i,1)*xi3(i,1));
    Fgeom(:,3*(i-1)+3) = -(Gamma311*xi1(i,1)*xi1(i,1) + Gamma312*xi1(i,1)*xi2(i,1) + Gamma313*xi1(i,1)*xi3(i,1) + Gamma321*xi2(i,1)*xi1(i,1) + Gamma322*xi2(i,1)*xi2(i,1) + Gamma323*xi2(i,1)*xi3(i,1) + Gamma331*xi3(i,1)*xi1(i,1) + Gamma332*xi3(i,1)*xi2(i,1) + Gamma333*xi3(i,1)*xi3(i,1));
end

return