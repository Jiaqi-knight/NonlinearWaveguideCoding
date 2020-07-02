function[lattice]=ellipticgridgen2D(Nx,N,lattice,deltaq,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 28th, 2014
%    Last update: July 17th, 2014
%
%          Input: meshed computational domain compdomain
%         Output: mesh in the physical domain

%%

boundaryindices = [indicesE1;indicesE2;indicesE3;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4];

it = 0;
err = 1;

while it<=itmax && err>=tol
    covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);
    [metriccoefficients,J] = computemetriccoefficients2D(covariantbase);
    %J = covariantbase(:,1).*covariantbase(:,4) - covariantbase(:,3).*covariantbase(:,2);
    bx = zeros(N,1);
    by = zeros(N,1);
    bx(boundaryindices,:) = lattice(boundaryindices,3);
    by(boundaryindices,:) = lattice(boundaryindices,4);
    Pvec = P(lattice);
    Qvec = Q(lattice);
    bx(indicesbulk,:) = -(J(indicesbulk,:).^2).*(Pvec(indicesbulk,:).*covariantbase(indicesbulk,1)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,3));
    by(indicesbulk,:) = -(J(indicesbulk,:).^2).*(Pvec(indicesbulk,:).*covariantbase(indicesbulk,2)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,4));
    clear covariantbase Pvec Qvec J
    A = eye(N);
    A(indicesbulk,indicesbulk) = diag(-2*(metriccoefficients(indicesbulk,2)./(deltaq(1).^2)+metriccoefficients(indicesbulk,1)./(deltaq(2).^2)));
    A(indicesbulk,indicesbulk-1) = A(indicesbulk,indicesbulk-1) + diag(metriccoefficients(indicesbulk,2)./(deltaq(1).^2));
    A(indicesbulk,indicesbulk+1) = A(indicesbulk,indicesbulk+1) + diag(metriccoefficients(indicesbulk,2)./(deltaq(1).^2));
    A(indicesbulk,indicesbulk-Nx) = A(indicesbulk,indicesbulk-Nx) + diag(metriccoefficients(indicesbulk,1)./(deltaq(2).^2));
    A(indicesbulk,indicesbulk+Nx) = A(indicesbulk,indicesbulk+Nx) + diag(metriccoefficients(indicesbulk,1)./(deltaq(2).^2));
    A(indicesbulk,indicesbulk-1-Nx) = A(indicesbulk,indicesbulk-1-Nx) + diag(-0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));
    A(indicesbulk,indicesbulk+1-Nx) = A(indicesbulk,indicesbulk+1-Nx) + diag(0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));
    A(indicesbulk,indicesbulk-1+Nx) = A(indicesbulk,indicesbulk-1+Nx) + diag(0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));
    A(indicesbulk,indicesbulk+1+Nx) = A(indicesbulk,indicesbulk+1+Nx) + diag(-0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));   
%     lattice(:,5) = gepp(A,bx);
%     lattice(:,6) = gepp(A,by);
    lattice(:,5) = A\bx;
    lattice(:,6) = A\by;
    err = sqrt(sum([bx-A*lattice(:,5);by-A*lattice(:,6)].^2));
    it = it + 1;
    clear A bx by metriccoefficients
end

return