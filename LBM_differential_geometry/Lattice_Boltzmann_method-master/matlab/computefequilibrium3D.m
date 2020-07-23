function[feq]=computefequilibrium3D(N,Q,rho,u,reciprocalmetriccoefficients,scheme,invcssq)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 30th, 2014
%    Last update: June 3rd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

feq = zeros(N,Q);

G11 = reciprocalmetriccoefficients(:,1);
G22 = reciprocalmetriccoefficients(:,2);
G33 = reciprocalmetriccoefficients(:,3);
G12 = reciprocalmetriccoefficients(:,4);
G13 = reciprocalmetriccoefficients(:,5);
G23 = reciprocalmetriccoefficients(:,6);

for i=1:Q
    feq(:,i) = scheme(i,4).*rho(:,1).*(0.5*5 +...
    2*(scheme(i,1).*u(:,1)+scheme(i,2).*u(:,2)+scheme(i,3).*u(:,3)).*invcssq +...
    0.5*(scheme(i,1).*G11.*scheme(i,1) + scheme(i,1).*G12.*scheme(i,2) + scheme(i,1).*G13.*scheme(i,3) + scheme(i,2).*G12.*scheme(i,1) + scheme(i,2).*G22.*scheme(i,2) + scheme(i,2).*G23.*scheme(i,3) + scheme(i,3).*G13.*scheme(i,1) + scheme(i,3).*G23.*scheme(i,2) + scheme(i,3).*G33.*scheme(i,3)).*invcssq +...
    -0.5*(scheme(i,1).*scheme(i,1)+scheme(i,2).*scheme(i,2)+scheme(i,3).*scheme(i,3)).*invcssq +...
    0.5*((scheme(i,1).*u(:,1)+scheme(i,2).*u(:,2)+scheme(i,3).*u(:,3)).^2).*invcssq.*invcssq +...
    -0.5*(G11+G22+G33)+...
    -0.5*(u(:,1).*u(:,1)+u(:,2).*u(:,2)+u(:,3).*u(:,3)).*invcssq+...
    ((scheme(i,1).*u(:,1)+scheme(i,2).*u(:,2)+scheme(i,3).*u(:,3)).^3).*invcssq.*invcssq.*invcssq./6 +...
    0.5*(scheme(i,1).*u(:,1)+scheme(i,2).*u(:,2)+scheme(i,3).*u(:,3)).*(scheme(i,1).*G11*scheme(i,1) + scheme(i,1).*G12*scheme(i,2) + scheme(i,1).*G13*scheme(i,3) + scheme(i,2).*G12*scheme(i,1) + scheme(i,2).*G22*scheme(i,2) + scheme(i,2).*G23*scheme(i,3) + scheme(i,3).*G13*scheme(i,1) + scheme(i,3).*G23*scheme(i,2) + scheme(i,3).*G33*scheme(i,3) -scheme(i,1).*scheme(i,1)-scheme(i,2).*scheme(i,2)-scheme(i,3).*scheme(i,3)).*invcssq.*invcssq +...
    -0.5*(scheme(i,1).*u(:,1)+scheme(i,2).*u(:,2)+scheme(i,3).*u(:,3)).*(u(:,1).*u(:,1)+u(:,2).*u(:,2)+u(:,3).*u(:,3)).*invcssq.*invcssq +...
    -0.5*(scheme(i,1).*u(:,1)+scheme(i,2).*u(:,2)+scheme(i,3).*u(:,3)).*(G11+G22+G33 - 3).*invcssq +...
    -(u(:,1).*G11.*scheme(i,1) + u(:,1).*G12.*scheme(i,2) + u(:,1).*G13.*scheme(i,3) + u(:,2).*G12.*scheme(i,1) + u(:,2).*G22.*scheme(i,2) + u(:,2).*G23.*scheme(i,3) + u(:,3).*G13.*scheme(i,1) + u(:,3).*G23.*scheme(i,2) + u(:,3).*G33.*scheme(i,3)).*invcssq);
end

return