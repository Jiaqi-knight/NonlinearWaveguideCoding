function[indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8]=getindices3D(Nx,Ny,Nz)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Z眉rich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: June 27th, 2014
%    Last update: June 27th, 2014
%
%    Description: 
%          Input: 
%         Output: indicesbulk:去掉外边界的其他内点，再扣掉bound
%                 indicesinternalbulk:再扣掉一层的内点
%                 indicesF1:面点
%                 indicesinternalF2:面点+里面一层，但去掉边上的点
%                 indicesE1:edge上的点-12条边
%                 indicesinternalE1:外边+里面一层，但去掉边上的点
%                 indicesC1:角点

%

%%

% ---> compute indices

% bulk fluid

indicesbulk = zeros((Nx-2)*(Ny-2)*(Nz-2),1);

for k=1:Nz-2
    for j=1:Ny-2
        for i=1:Nx-2
            indicesbulk(i+(j-1)*(Nx-2)+(k-1)*(Nx-2)*(Ny-2),1) = (i+1)+j*Nx+k*Nx*Ny;
        end
    end
end

indicesinternalbulk = zeros((Nx-4)*(Ny-4)*(Nz-4),1);

for k=1:Nz-4
    for j=1:Ny-4
        for i=1:Nx-4
            indicesinternalbulk(i+(j-1)*(Nx-4)+(k-1)*(Nx-4)*(Ny-4),1) = (i+2)+(j+1)*Nx+(k+1)*Nx*Ny;
        end
    end
end

% faces

indicesF1 = zeros((Nx-2)*(Ny-2),1); % xy

for j=1:Ny-2
    for i=1:Nx-2
        k = 0;
        indicesF1(i+(j-1)*(Nx-2),1) = (i+1)+j*Nx+k*Nx*Ny;
    end
end

indicesF2 = zeros((Nx-2)*(Nz-2),1); % xz

for k=1:Nz-2
    for i=1:Nx-2
        j = 0;
        indicesF2(i+(k-1)*(Nx-2),1) = (i+1)+j*Nx+k*Nx*Ny;
    end
end

indicesF3 = zeros((Ny-2)*(Nz-2),1); % yz

for k=1:Nz-2
    for j=1:Ny-2
        i = Nx;
        indicesF3(j+(k-1)*(Ny-2),1) = i+j*Nx+k*Nx*Ny;
    end
end

indicesF4 = zeros((Nx-2)*(Nz-2),1); % xz

for k=1:Nz-2
    for i=1:Nx-2
        j = Ny-1;
        indicesF4(i+(k-1)*(Nx-2),1) = (i+1)+j*Nx+k*Nx*Ny;
    end
end

indicesF5 = zeros((Ny-2)*(Nz-2),1); % yz

for k=1:Nz-2
    for j=1:Ny-2
        i = 1;
        indicesF5(j+(k-1)*(Ny-2),1) = i+j*Nx+k*Nx*Ny;
    end
end

indicesF6 = zeros((Nx-2)*(Ny-2),1); % xy

for j=1:Ny-2
    for i=1:Nx-2
        k = Nz-1;
        indicesF6(i+(j-1)*(Nx-2),1) = (i+1)+j*Nx+k*Nx*Ny;
    end
end

indicesinternalF1 = zeros((Nx-4)*(Ny-4),1); % xy

for j=1:Ny-4
    for i=1:Nx-4
        k = 2;
        indicesinternalF1(i+(j-1)*(Nx-4),1) = (i+2)+(j+1)*Nx+(k-1)*Nx*Ny;
    end
end

indicesinternalF1 = [indicesinternalF1;indicesinternalF1-Nx*Ny];

indicesinternalF2 = zeros((Nx-4)*(Nz-4),1); % xz

for k=1:Nz-4
    for i=1:Nx-4
        j = 2;
        indicesinternalF2(i+(k-1)*(Nx-4),1) = (i+2)+(j-1)*Nx+(k+1)*Nx*Ny;
    end
end

indicesinternalF2 = [indicesinternalF2;indicesinternalF2-Nx];

indicesinternalF3 = zeros((Ny-4)*(Nz-4),1); % yz

for k=1:Nz-4
    for j=1:Ny-4
        i = Nx-1;
        indicesinternalF3(j+(k-1)*(Ny-4),1) = i+(j+1)*Nx+(k+1)*Nx*Ny;
    end
end

indicesinternalF3 = [indicesinternalF3;indicesinternalF3+1];

indicesinternalF4 = zeros((Nx-4)*(Nz-4),1); % xz

for k=1:Nz-4
    for i=1:Nx-4
        j = Ny-1;
        indicesinternalF4(i+(k-1)*(Nx-4),1) = (i+2)+(j-1)*Nx+(k+1)*Nx*Ny;
    end
end

indicesinternalF4 = [indicesinternalF4;indicesinternalF4+Nx];

indicesinternalF5 = zeros((Ny-4)*(Nz-4),1); % yz

for k=1:Nz-4
    for j=1:Ny-4
        i = 2;
        indicesinternalF5(j+(k-1)*(Ny-4),1) = i+(j+1)*Nx+(k+1)*Nx*Ny;
    end
end

indicesinternalF5 = [indicesinternalF5;indicesinternalF5-1];

indicesinternalF6 = zeros((Nx-4)*(Ny-4),1); % xy

for j=1:Ny-4
    for i=1:Nx-4
        k = Nz-1;
        indicesinternalF6(i+(j-1)*(Nx-4),1) = (i+2)+(j+1)*Nx+(k-1)*Nx*Ny;
    end
end

indicesinternalF6 = [indicesinternalF6;indicesinternalF6+Nx*Ny];

% edges

indicesE1 = zeros(Nx-2,1); % x

for i=1:Nx-2
    j = 1;
    k = 1;
    indicesE1(i,1) = (i+1)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesE2 = zeros(Ny-2,1); % y

for j=1:Ny-2
    i = Nx;
    k = 1;
    indicesE2(j,1) = i+j*Nx+(k-1)*Nx*Ny;
end

indicesE3 = zeros(Nx-2,1); % x

for i=1:Nx-2
    j = Ny;
    k = 1;
    indicesE3(i,1) = (i+1)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesE4 = zeros(Ny-2,1); % y

for j=1:Ny-2
    i = 1;
    k = 1;
    indicesE4(j,1) = i+j*Nx+(k-1)*Nx*Ny;
end

indicesE5 = zeros(Nz-2,1); % z

for k=1:Nz-2
    i = 1;
    j = 1;
    indicesE5(k,1) = i+(j-1)*Nx+k*Nx*Ny;
end

indicesE6 = zeros(Nz-2,1); % z

for k=1:Nz-2
    i = Nx;
    j = 1;
    indicesE6(k,1) = i+(j-1)*Nx+k*Nx*Ny;
end

indicesE7 = zeros(Nz-2,1); % z

for k=1:Nz-2
    i = Nx;
    j = Ny;
    indicesE7(k,1) = i+(j-1)*Nx+k*Nx*Ny;
end

indicesE8 = zeros(Nz-2,1); % z

for k=1:Nz-2
    i = 1;
    j = Ny;
    indicesE8(k,1) = i+(j-1)*Nx+k*Nx*Ny;
end

indicesE9 = zeros(Nx-2,1); % x

for i=1:Nx-2
    j = 1;
    k = Nz;
    indicesE9(i,1) = (i+1)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesE10 = zeros(Ny-2,1); % y

for j=1:Ny-2
    i = Nx;
    k = Nz;
    indicesE10(j,1) = i+j*Nx+(k-1)*Nx*Ny;
end

indicesE11 = zeros(Nx-2,1); % x

for i=1:Nx-2
    j = Ny;
    k = Nz;
    indicesE11(i,1) = (i+1)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesE12 = zeros(Ny-2,1); % y

for j=1:Ny-2
    i = 1;
    k = Nz;
    indicesE12(j,1) = i+j*Nx+(k-1)*Nx*Ny;
end

indicesinternalE1 = zeros(Nx-4,1); % x

for i=1:Nx-4
    j = 2;
    k = 2;
    indicesinternalE1(i,1) = (i+2)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE1 = [indicesinternalE1;indicesinternalE1-Nx*Ny;indicesinternalE1-Nx-Nx*Ny;indicesinternalE1-Nx];

indicesinternalE2 = zeros(Ny-4,1); % y

for j=1:Ny-4
    i = Nx-1;
    k = 2;
    indicesinternalE2(j,1) = i+(j+1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE2 = [indicesinternalE2;indicesinternalE2-Nx*Ny;indicesinternalE2+1-Nx*Ny;indicesinternalE2+1];

indicesinternalE3 = zeros(Nx-4,1); % x

for i=1:Nx-4
    j = Ny-1;
    k = 2;
    indicesinternalE3(i,1) = (i+2)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE3 = [indicesinternalE3;indicesinternalE3-Nx*Ny;indicesinternalE3-Nx*Ny+Nx;indicesinternalE3+Nx];

indicesinternalE4 = zeros(Ny-4,1); % y

for j=1:Ny-4
    i = 2;
    k = 2;
    indicesinternalE4(j,1) = i+(j+1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE4 = [indicesinternalE4;indicesinternalE4-Nx*Ny;indicesinternalE4-Nx*Ny-1;indicesinternalE4-1];

indicesinternalE5 = zeros(Nz-4,1); % z

for k=1:Nz-4
    i = 2;
    j = 2;
    indicesinternalE5(k,1) = i+(j-1)*Nx+(k+1)*Nx*Ny;
end

indicesinternalE5 = [indicesinternalE5;indicesinternalE5-1;indicesinternalE5-Nx;indicesinternalE5-1-Nx];

indicesinternalE6 = zeros(Nz-4,1); % z

for k=1:Nz-4
    i = Nx-1;
    j = 2;
    indicesinternalE6(k,1) = i+(j-1)*Nx+(k+1)*Nx*Ny;
end

indicesinternalE6 = [indicesinternalE6;indicesinternalE6+1;indicesinternalE6-Nx;indicesinternalE6+1-Nx];

indicesinternalE7 = zeros(Nz-4,1); % z

for k=1:Nz-4
    i = Nx-1;
    j = Ny-1;
    indicesinternalE7(k,1) = i+(j-1)*Nx+(k+1)*Nx*Ny;
end

indicesinternalE7 = [indicesinternalE7;indicesinternalE7+1;indicesinternalE7+Nx;indicesinternalE7+1+Nx];

indicesinternalE8 = zeros(Nz-4,1); % z

for k=1:Nz-4
    i = 2;
    j = Ny-1;
    indicesinternalE8(k,1) = i+(j-1)*Nx+(k+1)*Nx*Ny;
end

indicesinternalE8 = [indicesinternalE8;indicesinternalE8-1;indicesinternalE8+Nx;indicesinternalE8-1+Nx];

indicesinternalE9 = zeros(Nx-4,1); % x

for i=1:Nx-4
    j = 2;
    k = Nz-1;
    indicesinternalE9(i,1) = (i+2)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE9 = [indicesinternalE9;indicesinternalE9+Nx*Ny;indicesinternalE9-Nx;indicesinternalE9-Nx+Nx*Ny];

indicesinternalE10 = zeros(Ny-4,1); % y

for j=1:Ny-4
    i = Nx-1;
    k = Nz-1;
    indicesinternalE10(j,1) = i+(j+1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE10 = [indicesinternalE10;indicesinternalE10+Nx*Ny;indicesinternalE10+1;indicesinternalE10+1+Nx*Ny];

indicesinternalE11 = zeros(Nx-4,1); % x

for i=1:Nx-4
    j = Ny-1;
    k = Nz-1;
    indicesinternalE11(i,1) = (i+2)+(j-1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE11 = [indicesinternalE11;indicesinternalE11+Nx*Ny;indicesinternalE11+Nx;indicesinternalE11+Nx+Nx*Ny];

indicesinternalE12 = zeros(Ny-4,1); % y

for j=1:Ny-4
    i = 2;
    k = Nz-1;
    indicesinternalE12(j,1) = i+(j+1)*Nx+(k-1)*Nx*Ny;
end

indicesinternalE12 = [indicesinternalE12;indicesinternalE12+Nx*Ny;indicesinternalE12-1;indicesinternalE12-1+Nx*Ny];

% corners

indicesC1 = 1;

indicesC2 = Nx;

indicesC3 = Nx + Nx*(Ny-1);

indicesC4 = 1 + Nx*(Ny-1);

indicesC5 = 1 + Nx*Ny*(Nz-1);

indicesC6 = Nx + Nx*Ny*(Nz-1);

indicesC7 = Nx + Nx*(Ny-1) + Nx*Ny*(Nz-1);

indicesC8 = 1 + Nx*(Ny-1) + Nx*Ny*(Nz-1);

%A = -1-Nx-Nx*Ny;
B =   -Nx-Nx*Ny;
%C = +1-Nx-Nx*Ny;
D = -1   -Nx*Ny;
E =      -Nx*Ny;
F = +1   -Nx*Ny;
%G = -1+Nx-Nx*Ny;
H =   +Nx-Nx*Ny;
%I = +1+Nx-Nx*Ny;
L = -1-Nx      ;
M =   -Nx      ;
N = +1-Nx      ;
O = -1         ;
P = +1         ;
Q = -1+Nx      ;
R =   +Nx      ;
S = +1+Nx      ;
%T = -1-Nx+Nx*Ny;
U =   -Nx+Nx*Ny;
%V = +1-Nx+Nx*Ny;
Z = -1   +Nx*Ny;
K =      +Nx*Ny;
J = +1   +Nx*Ny;
%Y = -1+Nx+Nx*Ny;
X =   +Nx+Nx*Ny;
%W = +1+Nx+Nx*Ny;

l = 2 + Nx + Nx*Ny;
indicesinternalC1 = [l;l+B;l+D;l+E;l+L;l+M;l+O];

l = (Nx-1) + Nx + Nx*Ny;
indicesinternalC2 = [l;l+B;l+E;l+F;l+M;l+N;l+P];

l = Nx*Ny - Nx - 1 + Nx*Ny;
indicesinternalC3 = [l;l+E;l+F;l+H;l+P;l+R;l+S];

l = 2 + Nx*(Ny-2) + Nx*Ny;
indicesinternalC4 = [l;l+D;l+E;l+H;l+O;l+Q;l+R];

l = 2 + Nx + Nx*Ny + Nx*Ny*(Nz-2);
indicesinternalC5 = [l;l+L;l+M;l+O;l+U;l+K;l+Z];

l = (Nx-1) + Nx + Nx*Ny*(Nz-2);
indicesinternalC6 = [l;l+M;l+N;l+P;l+U;l+K;l+J];

l = Nx*Ny - Nx - 1 + Nx*Ny*(Nz-2);
indicesinternalC7 = [l;l+P;l+R;l+S;l+K;l+J;l+X];

l = 2 + Nx*(Ny-2) + Nx*Ny*(Nz-2);
indicesinternalC8 = [l;l+O;l+Q;l+R;l+Z;l+K;l+X];

return