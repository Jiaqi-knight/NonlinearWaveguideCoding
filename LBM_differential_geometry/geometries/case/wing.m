function[wing,halfwing,upperhalfwing,lowerhalfwing]=wing(Nx,Ny,L,x1func,x2func,x3_4func,cfunc,ctip,croot,meanlinefunc)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 14th, 2014
%    Last update: July 14th, 2014
%
%          Input: 
%         Output: 

%%

y = (0:L/(Ny-1):L)';

upperhalfwing = [];
lowerhalfwing = [];
wing = [];

for i=1:Ny
    x1 = x1func(y(i));
    x2 = x2func(y(i));
    x3_4 = x3_4func(y(i));
    c = cfunc(y(i),ctip,croot,L);
    coords = meanlinefunc(y(i));
    [airfoil,upperairfoil,lowerairfoil] = NACAxxxx(Nx,x1,x2,x3_4,c);
    upperdim = size(upperairfoil,1);
    lowerdim = size(lowerairfoil,1);
    dim = size(airfoil,1);
    upperhalfwing = [upperhalfwing; coords(1)*ones(upperdim,1)+upperairfoil(:,1) coords(2)*ones(upperdim,1) coords(3)*ones(upperdim,1)+upperairfoil(:,2)];
    lowerhalfwing = [lowerhalfwing; coords(1)*ones(lowerdim,1)+lowerairfoil(:,1) coords(2)*ones(lowerdim,1) coords(3)*ones(lowerdim,1)+lowerairfoil(:,2)];
    wing = [wing; coords(1)*ones(dim,1)+airfoil(:,1) coords(2)*ones(dim,1) coords(3)*ones(dim,1)+airfoil(:,2)];
end

halfwing = wing;
wing = [wing; wing(:,1) -wing(:,2) wing(:,3)];

return 