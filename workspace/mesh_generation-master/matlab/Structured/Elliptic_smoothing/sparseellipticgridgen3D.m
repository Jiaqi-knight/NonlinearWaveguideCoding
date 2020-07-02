function[lattice]=sparseellipticgridgen3D(Nx,Ny,N,lattice,deltaq,indicesbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,firstdevneighbours,itmax,tol,spyflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 23rd, 2014
%    Last update: July 24th, 2014
%
%          Input: meshed computational domain compdomain
%         Output: mesh in the physical domain

%%

boundaryindices = [indicesF1;indicesF2;indicesF3;indicesF4;indicesF5;indicesF6;indicesE1;indicesE2;indicesE3;indicesE4;indicesE5;indicesE6;indicesE7;indicesE8;indicesE9;indicesE10;indicesE11;indicesE12;indicesC1;indicesC2;indicesC3;indicesC4;indicesC5;indicesC6;indicesC7;indicesC8];

it = 0;
err = 1;

if spyflag
    covariantbase = computecovariantbase3D(N,deltaq,lattice,firstdevneighbours);
    [metriccoefficients,J,sqrtg] = computemetriccoefficients3D(covariantbase);
    clear metriccoefficients J
    reciprocalmetriccoefficients = computereciprocalmetriccoefficients3D(computecontravariantbase3D(covariantbase,sqrtg));   
    %J = covariantbase(:,1).*(covariantbase(:,5).*covariantbase(:,9) - covariantbase(:,8).*covariantbase(:,6)) + covariantbase(:,4).*(covariantbase(:,8).*covariantbase(:,3) - covariantbase(:,2).*covariantbase(:,9)) + covariantbase(:,7).*(covariantbase(:,2).*covariantbase(:,6) - covariantbase(:,5).*covariantbase(:,3));
    clear covariantbase
    G11 = reciprocalmetriccoefficients(:,1);
    G22 = reciprocalmetriccoefficients(:,2);
    G33 = reciprocalmetriccoefficients(:,3);
    G12 = reciprocalmetriccoefficients(:,4);
    G13 = reciprocalmetriccoefficients(:,5);
    G23 = reciprocalmetriccoefficients(:,6);
    clear reciprocalmetriccoefficients
    A = sparse([boundaryindices;                 indicesbulk;                                                                                             indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk],...
               [boundaryindices;                 indicesbulk;                                                                                             indicesbulk-1;                    indicesbulk+1;                    indicesbulk-Nx;                   indicesbulk+Nx;                   indicesbulk-Nx*Ny;                indicesbulk+Nx*Ny;                indicesbulk-1-Nx;                              indicesbulk+1-Nx;                              indicesbulk-1+Nx;                              indicesbulk+1+Nx;                              indicesbulk-1-Nx*Ny;                           indicesbulk+1-Nx*Ny;                           indicesbulk-1+Nx*Ny;                           indicesbulk+1+Nx*Ny;                           indicesbulk-Nx-Nx*Ny;                          indicesbulk+Nx-Nx*Ny;                          indicesbulk-Nx+Nx*Ny;                          indicesbulk+Nx+Nx*Ny],...
               [ones(size(boundaryindices,1),1); -2*(G11(indicesbulk)./(deltaq(1).^2)+G22(indicesbulk)./(deltaq(2).^2)+G33(indicesbulk)./(deltaq(3).^2)); G11(indicesbulk)./(deltaq(1).^2); G11(indicesbulk)./(deltaq(1).^2); G22(indicesbulk)./(deltaq(2).^2); G22(indicesbulk)./(deltaq(2).^2); G33(indicesbulk)./(deltaq(3).^2); G33(indicesbulk)./(deltaq(3).^2); +0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); -0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); -0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); +0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); +0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); -0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); -0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); +0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); +0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3)); -0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3)); -0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3)); +0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3))],N,N); 
    clear G11 G22 G33 G12 G13 G23
    figure;
    spy(A)
    hold on
    grid on
    title('Structure of solving matrix for 3D elliptic grid generation')
    clear A
end

while it<=itmax && err>=tol
    covariantbase = computecovariantbase3D(N,deltaq,lattice,firstdevneighbours);
    [metriccoefficients,J,sqrtg] = computemetriccoefficients3D(covariantbase);
    clear metriccoefficients J
    reciprocalmetriccoefficients = computereciprocalmetriccoefficients3D(computecontravariantbase3D(covariantbase,sqrtg));   
    %J = covariantbase(:,1).*(covariantbase(:,5).*covariantbase(:,9) - covariantbase(:,8).*covariantbase(:,6)) + covariantbase(:,4).*(covariantbase(:,8).*covariantbase(:,3) - covariantbase(:,2).*covariantbase(:,9)) + covariantbase(:,7).*(covariantbase(:,2).*covariantbase(:,6) - covariantbase(:,5).*covariantbase(:,3));
    Pvec = sparse(P(lattice));
    Qvec = sparse(Q(lattice));
    Rvec = sparse(R(lattice));
    bx = sparse([boundaryindices;indicesbulk],[ones(length(boundaryindices),1);ones(length(indicesbulk),1)],[lattice(boundaryindices,4);-(Pvec(indicesbulk,:).*covariantbase(indicesbulk,1)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,4)+Rvec(indicesbulk,:).*covariantbase(indicesbulk,7))],N,1);
    by = sparse([boundaryindices;indicesbulk],[ones(length(boundaryindices),1);ones(length(indicesbulk),1)],[lattice(boundaryindices,5);-(Pvec(indicesbulk,:).*covariantbase(indicesbulk,2)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,5)+Rvec(indicesbulk,:).*covariantbase(indicesbulk,8))],N,1);
    bz = sparse([boundaryindices;indicesbulk],[ones(length(boundaryindices),1);ones(length(indicesbulk),1)],[lattice(boundaryindices,6);-(Pvec(indicesbulk,:).*covariantbase(indicesbulk,3)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,6)+Rvec(indicesbulk,:).*covariantbase(indicesbulk,9))],N,1);
    clear covariantbase Pvec Qvec
    G11 = reciprocalmetriccoefficients(:,1);
    G22 = reciprocalmetriccoefficients(:,2);
    G33 = reciprocalmetriccoefficients(:,3);
    G12 = reciprocalmetriccoefficients(:,4);
    G13 = reciprocalmetriccoefficients(:,5);
    G23 = reciprocalmetriccoefficients(:,6);
    clear reciprocalmetriccoefficients
    A = sparse([boundaryindices;                 indicesbulk;                                                                                             indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                      indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk;                                   indicesbulk],...
               [boundaryindices;                 indicesbulk;                                                                                             indicesbulk-1;                    indicesbulk+1;                    indicesbulk-Nx;                   indicesbulk+Nx;                   indicesbulk-Nx*Ny;                indicesbulk+Nx*Ny;                indicesbulk-1-Nx;                              indicesbulk+1-Nx;                              indicesbulk-1+Nx;                              indicesbulk+1+Nx;                              indicesbulk-1-Nx*Ny;                           indicesbulk+1-Nx*Ny;                           indicesbulk-1+Nx*Ny;                           indicesbulk+1+Nx*Ny;                           indicesbulk-Nx-Nx*Ny;                          indicesbulk+Nx-Nx*Ny;                          indicesbulk-Nx+Nx*Ny;                          indicesbulk+Nx+Nx*Ny],...
               [ones(size(boundaryindices,1),1); -2*(G11(indicesbulk)./(deltaq(1).^2)+G22(indicesbulk)./(deltaq(2).^2)+G33(indicesbulk)./(deltaq(3).^2)); G11(indicesbulk)./(deltaq(1).^2); G11(indicesbulk)./(deltaq(1).^2); G22(indicesbulk)./(deltaq(2).^2); G22(indicesbulk)./(deltaq(2).^2); G33(indicesbulk)./(deltaq(3).^2); G33(indicesbulk)./(deltaq(3).^2); +0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); -0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); -0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); +0.5*G12(indicesbulk)./(deltaq(1).*deltaq(2)); +0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); -0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); -0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); +0.5*G13(indicesbulk)./(deltaq(1).*deltaq(3)); +0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3)); -0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3)); -0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3)); +0.5*G23(indicesbulk)./(deltaq(2).*deltaq(3))],N,N); 
    clear G11 G22 G33 G12 G13 G23
    lattice(:,7) = full(A\bx);
    lattice(:,8) = full(A\by);
    lattice(:,9) = full(A\bz);
    err = sqrt(sum([bx-A*lattice(:,7);by-A*lattice(:,8);bz-A*lattice(:,9)].^2));
    it = it + 1;
    clear A bx by
end

return