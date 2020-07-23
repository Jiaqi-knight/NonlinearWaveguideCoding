function[indicesbulk,indicesinternalbulk,...
         indicesE1,indicesE2,indicesE3,indicesE4,...
         indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,...
         indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,...
         indicesC1,indicesC2,indicesC3,indicesC4,...
         indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4]=getindices2D(logfullfile,Nx,Ny)

%%
%==============================================================================
% Copyright (c) 2016 - 2017 Universit? de Lorraine & Lule? tekniska universitet
% Author: Luca Di Stasio <luca.distasio@gmail.com>
%                        <luca.distasio@ingpec.eu>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 
% Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% Neither the name of the Universit? de Lorraine or Lule? tekniska universitet
% nor the names of its contributors may be used to endorse or promote products
% derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================
%
%  DESCRIPTION
%  
%  A function to
%
%%

writeToLogFile(logfullfile,'In function: getindices2D\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

indicesbulk = zeros((Nx-2)*(Ny-2),1);


for j=1:Ny-2
    for i=1:Nx-2
        indicesbulk(i+(j-1)*(Nx-2),1) = (i+1)+j*Nx;
    end
end


indicesinternalbulk = zeros((Nx-4)*(Ny-4),1);


for j=1:Ny-4
    for i=1:Nx-4
        indicesinternalbulk(i+(j-1)*(Nx-4),1) = (i+2)+(j+1)*Nx;
    end
end

% edges

indicesE1 = zeros(Nx-2,1); % x…œ±ﬂ

for i=1:Nx-2
    j = 1;
    indicesE1(i,1) = (i+1)+(j-1)*Nx;
end

indicesE2 = zeros(Ny-2,1); % y”“±ﬂ

for j=1:Ny-2
    i = Nx;
    indicesE2(j,1) = i+j*Nx;
end

indicesE3 = zeros(Nx-2,1); % xœ¬±ﬂ

for i=1:Nx-2
    j = Ny;
    indicesE3(i,1) = (i+1)+(j-1)*Nx;
end

indicesE4 = zeros(Ny-2,1); % y◊Û±ﬂ

for j=1:Ny-2
    i = 1;
    indicesE4(j,1) = i+j*Nx;
end

indicesinternalE1 = zeros(Nx-4,1); % x;…œ-2≈≈

for i=1:Nx-4
    j = 2;
    indicesinternalE1(i,1) = (i+2)+(j-1)*Nx;
end

indicesinternalE2 = zeros(Ny-4,1); % y;”“-2≈≈

for j=1:Ny-4
    i = Nx-1;
    indicesinternalE2(j,1) = i+(j+1)*Nx;
end

indicesinternalE3 = zeros(Nx-4,1); % x;œ¬-2≈≈

for i=1:Nx-4
    j = Ny-1;
    indicesinternalE3(i,1) = (i+2)+(j-1)*Nx;
end

indicesinternalE4 = zeros(Ny-4,1); % y;◊Û-2≈≈

for j=1:Ny-4
    i = 2;
    indicesinternalE4(j,1) = i+(j+1)*Nx;
end

indicesexternalE1 = indicesE1(indicesE1(indicesE1~=2)~=Nx-1);

indicesexternalE2 = indicesE2(indicesE2(indicesE2~=2*Nx)~=Nx+(Ny-2)*Nx);

indicesexternalE3 = indicesE3(indicesE3(indicesE3~=Nx*Ny-1)~=2+(Ny-1)*Nx);

indicesexternalE4 = indicesE4(indicesE4(indicesE4~=1+Nx)~=1+(Ny-2)*Nx);

% corners

indicesC1 = 1;

indicesC2 = Nx;

indicesC3 = Nx + Nx*(Ny-1);

indicesC4 = 1 + Nx*(Ny-1);

indicesinternalC1 = [2 + Nx + Nx*Ny;2;1+Nx];

indicesinternalC2 = [(Nx-1) + Nx + Nx*Ny;Nx-1;2*Nx];

indicesinternalC3 = [Nx*Ny - Nx - 1 + Nx*Ny;Nx+(Ny-2)*Nx;Nx*Ny-1];

indicesinternalC4 = [2 + Nx*(Ny-2) + Nx*Ny;1+(Ny-2)*Nx;2+(Ny-1)*Nx];

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: getindices2D\n')

return