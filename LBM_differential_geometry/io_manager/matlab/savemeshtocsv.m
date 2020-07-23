function[]=savemeshtocsv(date,outfilename,latticefolder,meshfolder,N,Nx,Ny,Nz,lattice)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 29th, 2014
%    Last update: July 29th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

baseoutfilename = strcat(outfilename,'_compdomain');

flags = zeros(37,1);
flags(1) = 1;

savetovtk(date,0,0,baseoutfilename,latticefolder,meshfolder,flags,N,Nx,Ny,Nz,lattice,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

baseoutfilename = strcat(outfilename,'_physdomain');

flags(1) = 0;
flags(2) = 1;

savetovtk(date,0,0,baseoutfilename,latticefolder,meshfolder,flags,N,Nx,Ny,Nz,lattice,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

return