function[lattice]=transfiniteinterpolation3D(Nx,Ny,Nz,compdomain,xi_min,xi_max,eta_min,eta_max,zeta_min,zeta_max,interpolanttype,f1,f2,f3,f4,f5,f6,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,c1,c2,c3,c4,c5,c6,c7,c8)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 23rd, 2014
%    Last update: July 28th, 2014
%
%          Input: meshed computational domain compdomain
%                 numerical flag to choose the interpolant type
%                 f1 = e1(xi,eta_min,zeta) = f1(xi) = (x(xi),y(xi))
%                 f2 = e2(xi_min,eta,zeta) = f2(eta) = (x(eta),y(eta))
%                 f3 = e3(xi,eta_max,zeta) = f3(xi) = (x(xi),y(xi))
%                 f4 = e4(xi_max,eta,zeta) = f4(eta) = (x(eta),y(eta))
%                 f5 = e1(xi,eta_min,zeta) = f5(xi) = (x(xi),y(xi))
%                 f6 = e2(xi_min,eta,zeta) = f6(eta) = (x(eta),y(eta))
%                 e1 = e1(xi,eta_min,zeta) = e1(xi) = (x(xi),y(xi))
%                 e2 = e2(xi_min,eta,zeta) = e2(eta) = (x(eta),y(eta))
%                 e3 = e3(xi,eta_max,zeta) = e3(xi) = (x(xi),y(xi))
%                 e4 = e4(xi_max,eta,zeta) = e4(eta) = (x(eta),y(eta))
%                 e5 = e1(xi,eta_min,zeta) = e5(xi) = (x(xi),y(xi))
%                 e6 = e2(xi_min,eta,zeta) = e6(eta) = (x(eta),y(eta))
%                 e7 = e3(xi,eta_max,zeta) = e7(xi) = (x(xi),y(xi))
%                 e8 = e4(xi_max,eta,zeta) = e8(eta) = (x(eta),y(eta))
%                 e9 = e1(xi,eta_min,zeta) = e9(xi) = (x(xi),y(xi))
%                 e10 = e2(xi_min,eta,zeta) = e10(eta) = (x(eta),y(eta))
%                 e11 = e3(xi,eta_max,zeta) = e11(xi) = (x(xi),y(xi))
%                 e12 = e4(xi_max,eta,zeta) = e12(eta) = (x(eta),y(eta))
%                 c1 = c1(xi_min,eta_min) = (x1,y1)
%                 c2 = c2(xi_max,eta_min) = (x2,y2)
%                 c3 = c3(xi_max,eta_max) = (x3,y3)
%                 c4 = c4(xi_min,eta_max) = (x4,y4)
%                 c5 = c1(xi_min,eta_min) = (x1,y1)
%                 c6 = c2(xi_max,eta_min) = (x2,y2)
%                 c7 = c3(xi_max,eta_max) = (x3,y3)
%                 c8 = c4(xi_min,eta_max) = (x4,y4)
%                 f1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 f2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 f3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 f4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 f5 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 f6 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e5 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e6 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e7 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e8 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e9 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e10 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e11 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e12 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c5 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c6 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c7 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c8 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%
%    Interpolant: 1 --> Lagrange
%                 2 --> Hermite
%
%         Output: mesh in the physical domain

%%

N = Nx*Ny*Nz;

mesh = zeros(N,3);

switch interpolanttype
    case 1 %--------------------------- Lagrange
        for s=1:Nz
            for r=1:Ny
                for p=1:Nx
                    xi   = compdomain(p + (r-1)*Nx + (s-1)*Nx*Ny,1);
                    eta  = compdomain(p + (r-1)*Nx + (s-1)*Nx*Ny,2);
                    zeta = compdomain(p + (r-1)*Nx + (s-1)*Nx*Ny,3);
                    mesh(p + (r-1)*Nx + (s-1)*Nx*Ny,1) = Lagrangeinterps1D([xi_min;xi_max],[f5(r+(s-1)*Ny,1);f3(r+(s-1)*Ny,1)],xi) + Lagrangeinterps1D([eta_min;eta_max],[f2(p+(s-1)*Nx,1);f4(p+(s-1)*Nx,1)],eta) + Lagrangeinterps1D([zeta_min;zeta_max],[f1(p+(r-1)*Nx,1);f6(p+(r-1)*Nx,1)],zeta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[e5(s,1) e8(s,1);e6(s,1) e7(s,1)],xi,eta) - Lagrangeinterps2D([xi_min;xi_max],[zeta_min;zeta_max],[e4(r,1) e12(r,1);e2(r,1) e10(r,1)],xi,zeta) - Lagrangeinterps2D([eta_min;eta_max],[zeta_min;zeta_max],[e1(p,1) e9(p,1);e3(p,1) e11(p,1)],eta,zeta) + Lagrangeinterps3D([xi_min;xi_max],[eta_min;eta_max],[zeta_min;zeta_max],cat(3,[c1(1,1) c4(1,1);c5(1,1) c8(1,1)],[c2(1,1) c3(1,1);c6(1,1) c7(1,1)]),xi,eta,zeta);
                    mesh(p + (r-1)*Nx + (s-1)*Nx*Ny,2) = Lagrangeinterps1D([xi_min;xi_max],[f5(r+(s-1)*Ny,2);f3(r+(s-1)*Ny,2)],xi) + Lagrangeinterps1D([eta_min;eta_max],[f2(p+(s-1)*Nx,2);f4(p+(s-1)*Nx,2)],eta) + Lagrangeinterps1D([zeta_min;zeta_max],[f1(p+(r-1)*Nx,2);f6(p+(r-1)*Nx,2)],zeta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[e5(s,2) e8(s,2);e6(s,2) e7(s,2)],xi,eta) - Lagrangeinterps2D([xi_min;xi_max],[zeta_min;zeta_max],[e4(r,2) e12(r,2);e2(r,2) e10(r,2)],xi,zeta) - Lagrangeinterps2D([eta_min;eta_max],[zeta_min;zeta_max],[e1(p,2) e9(p,2);e3(p,2) e11(p,2)],eta,zeta) + Lagrangeinterps3D([xi_min;xi_max],[eta_min;eta_max],[zeta_min;zeta_max],cat(3,[c1(1,2) c4(1,2);c5(1,2) c8(1,2)],[c2(1,2) c3(1,2);c6(1,2) c7(1,2)]),xi,eta,zeta);
                    mesh(p + (r-1)*Nx + (s-1)*Nx*Ny,3) = Lagrangeinterps1D([xi_min;xi_max],[f5(r+(s-1)*Ny,3);f3(r+(s-1)*Ny,3)],xi) + Lagrangeinterps1D([eta_min;eta_max],[f2(p+(s-1)*Nx,3);f4(p+(s-1)*Nx,3)],eta) + Lagrangeinterps1D([zeta_min;zeta_max],[f1(p+(r-1)*Nx,3);f6(p+(r-1)*Nx,3)],zeta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[e5(s,3) e8(s,3);e6(s,3) e7(s,3)],xi,eta) - Lagrangeinterps2D([xi_min;xi_max],[zeta_min;zeta_max],[e4(r,3) e12(r,3);e2(r,3) e10(r,3)],xi,zeta) - Lagrangeinterps2D([eta_min;eta_max],[zeta_min;zeta_max],[e1(p,3) e9(p,3);e3(p,3) e11(p,3)],eta,zeta) + Lagrangeinterps3D([xi_min;xi_max],[eta_min;eta_max],[zeta_min;zeta_max],cat(3,[c1(1,3) c4(1,3);c5(1,3) c8(1,3)],[c2(1,3) c3(1,3);c6(1,3) c7(1,3)]),xi,eta,zeta);
                end
            end
        end
    case 2 %--------------------------- Hermite
        for s=1:Nz
            for r=1:Ny
                for p=1:Nx
                    xi   = compdomain(p + (r-1)*Nx + (s-1)*Nx*Ny,1);
                    eta  = compdomain(p + (r-1)*Nx + (s-1)*Nx*Ny,2);
                    zeta = compdomain(p + (r-1)*Nx + (s-1)*Nx*Ny,3);
                    mesh(p + (r-1)*Nx + (s-1)*Nx*Ny,1) = Hermiteinterps1D([xi_min;xi_max],[e2(r,1);e4(r,1)],[e2(r,3);e4(r,3)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,1);e3(p,1)],[e2(r,4);e4(r,4)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,1) c4(1,1);c2(1,1) c3(1,1)],[c1(1,3) c4(1,3);c2(1,3) c3(1,3)],[c1(1,4) c4(1,4);c2(1,4) c3(1,4)],[c1(1,7) c4(1,7);c2(1,7) c3(1,7)],xi,eta);
                    mesh(p + (r-1)*Nx + (s-1)*Nx*Ny,2) = Hermiteinterps1D([xi_min;xi_max],[e2(r,2);e4(r,2)],[e2(r,5);e4(r,5)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,2);e3(p,2)],[e2(r,6);e4(r,6)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,2) c4(1,2);c2(1,2) c3(1,2)],[c1(1,5) c4(1,5);c2(1,5) c3(1,5)],[c1(1,6) c4(1,6);c2(1,6) c3(1,6)],[c1(1,8) c4(1,8);c2(1,8) c3(1,8)],xi,eta);
                    mesh(p + (r-1)*Nx + (s-1)*Nx*Ny,3) = Hermiteinterps1D([xi_min;xi_max],[e2(r,2);e4(r,2)],[e2(r,5);e4(r,5)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,2);e3(p,2)],[e2(r,6);e4(r,6)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,2) c4(1,2);c2(1,2) c3(1,2)],[c1(1,5) c4(1,5);c2(1,5) c3(1,5)],[c1(1,6) c4(1,6);c2(1,6) c3(1,6)],[c1(1,8) c4(1,8);c2(1,8) c3(1,8)],xi,eta);
                end
            end
        end
end

[indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8] = getindices3D(Nx,Ny,Nz);

f1 = f1(Nx+1:end-Nx,:);
f1temp = zeros((Nx-2)*(Ny-2),3);
for i=1:Ny-2
    f1temp((i-1)*(Nx-2)+1:i*(Nx-2),1:3) = f1((i-1)*Nx+2:i*Nx-1,1:3);
end
f1 = f1temp;
clear f1temp
mesh(indicesF1,:) = f1(:,:);

f6 = f6(Nx+1:end-Nx,:);
f6temp = zeros((Nx-2)*(Ny-2),3);
for i=1:Ny-2
    f6temp((i-1)*(Nx-2)+1:i*(Nx-2),1:3) = f6((i-1)*Nx+2:i*Nx-1,1:3);
end
f6 = f6temp;
clear f6temp
mesh(indicesF6,:) = f6(:,:);

f2 = f2(Nx+1:end-Nx,:);
f2temp = zeros((Nx-2)*(Nz-2),3);
for i=1:Nz-2
    f2temp((i-1)*(Nx-2)+1:i*(Nx-2),1:3) = f2((i-1)*Nx+2:i*Nx-1,1:3);
end
f2 = f2temp;
clear f2temp
mesh(indicesF2,:) = f2(:,:);

f4 = f4(Nx+1:end-Nx,:);
f4temp = zeros((Nx-2)*(Nz-2),3);
for i=1:Nz-2
    f4temp((i-1)*(Nx-2)+1:i*(Nx-2),1:3) = f4((i-1)*Nx+2:i*Nx-1,1:3);
end
f4 = f4temp;
clear f4temp
mesh(indicesF4,:) = f4(:,:);

f3 = f3(Ny+1:end-Ny,:);
f3temp = zeros((Ny-2)*(Nz-2),3);
for i=1:Nz-2
    f3temp((i-1)*(Ny-2)+1:i*(Ny-2),1:3) = f3((i-1)*Ny+2:i*Ny-1,1:3);
end
f3 = f3temp;
clear f3temp
mesh(indicesF3,:) = f3(:,:);

f5 = f5(Ny+1:end-Ny,:);
f5temp = zeros((Ny-2)*(Nz-2),3);
for i=1:Nz-2
    f5temp((i-1)*(Ny-2)+1:i*(Ny-2),1:3) = f5((i-1)*Ny+2:i*Ny-1,1:3);
end
f5 = f5temp;
clear f5temp
mesh(indicesF5,:) = f5(:,:);

mesh(indicesC1,1) = c1(1,1);
mesh(indicesC1,2) = c1(1,2);
mesh(indicesC1,3) = c1(1,3);

mesh(indicesC2,1) = c2(1,1);
mesh(indicesC2,2) = c2(1,2);
mesh(indicesC2,3) = c2(1,3);

mesh(indicesC3,1) = c3(1,1);
mesh(indicesC3,2) = c3(1,2);
mesh(indicesC3,3) = c3(1,3);

mesh(indicesC4,1) = c4(1,1);
mesh(indicesC4,2) = c4(1,2);
mesh(indicesC4,3) = c4(1,3);

mesh(indicesC5,1) = c5(1,1);
mesh(indicesC5,2) = c5(1,2);
mesh(indicesC5,3) = c5(1,3);

mesh(indicesC6,1) = c6(1,1);
mesh(indicesC6,2) = c6(1,2);
mesh(indicesC6,3) = c6(1,3);

mesh(indicesC7,1) = c7(1,1);
mesh(indicesC7,2) = c7(1,2);
mesh(indicesC7,3) = c7(1,3);

mesh(indicesC8,1) = c8(1,1);
mesh(indicesC8,2) = c8(1,2);
mesh(indicesC8,3) = c8(1,3);

mesh(indicesE1,1) = e1(2:end-1,1);
mesh(indicesE1,2) = e1(2:end-1,2);
mesh(indicesE1,3) = e1(2:end-1,3);

mesh(indicesE2,1) = e2(2:end-1,1);
mesh(indicesE2,2) = e2(2:end-1,2);
mesh(indicesE2,3) = e2(2:end-1,3);

mesh(indicesE3,1) = e3(2:end-1,1);
mesh(indicesE3,2) = e3(2:end-1,2);
mesh(indicesE3,3) = e3(2:end-1,3);

mesh(indicesE4,1) = e4(2:end-1,1);
mesh(indicesE4,2) = e4(2:end-1,2);
mesh(indicesE4,3) = e4(2:end-1,3);

mesh(indicesE5,1) = e5(2:end-1,1);
mesh(indicesE5,2) = e5(2:end-1,2);
mesh(indicesE5,3) = e5(2:end-1,3);

mesh(indicesE6,1) = e6(2:end-1,1);
mesh(indicesE6,2) = e6(2:end-1,2);
mesh(indicesE6,3) = e6(2:end-1,3);

mesh(indicesE7,1) = e7(2:end-1,1);
mesh(indicesE7,2) = e7(2:end-1,2);
mesh(indicesE7,3) = e7(2:end-1,3);

mesh(indicesE8,1) = e8(2:end-1,1);
mesh(indicesE8,2) = e8(2:end-1,2);
mesh(indicesE8,3) = e8(2:end-1,3);

mesh(indicesE9,1) = e9(2:end-1,1);
mesh(indicesE9,2) = e9(2:end-1,2);
mesh(indicesE9,3) = e9(2:end-1,3);

mesh(indicesE10,1) = e10(2:end-1,1);
mesh(indicesE10,2) = e10(2:end-1,2);
mesh(indicesE10,3) = e10(2:end-1,3);

mesh(indicesE11,1) = e11(2:end-1,1);
mesh(indicesE11,2) = e11(2:end-1,2);
mesh(indicesE11,3) = e11(2:end-1,3);

mesh(indicesE12,1) = e12(2:end-1,1);
mesh(indicesE12,2) = e12(2:end-1,2);
mesh(indicesE12,3) = e12(2:end-1,3);

lattice = [compdomain mesh mesh];

return