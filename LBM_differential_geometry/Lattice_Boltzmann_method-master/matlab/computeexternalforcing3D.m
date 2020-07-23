function[Fext]=computeexternalforcing3D(N,Q,contravariantbase,F)

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

G1 = contravariantbase(:,1:3);
G2 = contravariantbase(:,4:6);
G3 = contravariantbase(:,7:9);

Nforces = round(size(F,2)/3);

Fext = zeros(N,3*Q);

for i=1:Q
    Fext(:,3*(i-1)+1) = 0;
    Fext(:,3*(i-1)+2) = 0;
    Fext(:,3*(i-1)+3) = 0;
    for j=1:Nforces
        Fext(:,3*(i-1)+1) = Fext(:,3*(i-1)+1) + sum(F(:,3*(j-1)+1:3*(j-1)+3).*G1,2);
        Fext(:,3*(i-1)+2) = Fext(:,3*(i-1)+2) + sum(F(:,3*(j-1)+1:3*(j-1)+3).*G2,2);
        Fext(:,3*(i-1)+3) = Fext(:,3*(i-1)+3) + sum(F(:,3*(j-1)+1:3*(j-1)+3).*G3,2);
    end
end

return