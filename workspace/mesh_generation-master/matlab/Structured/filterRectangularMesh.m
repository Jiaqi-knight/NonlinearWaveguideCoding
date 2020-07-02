function[nodes,elements,edges,...
         nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
         edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
         elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside]=filterRectangularMesh(logfullfile,...
                                                                                                                                                                          elType,elOrder,Nx,Ny,NxEquiv,NyEquiv,...
                                                                                                                                                                          baseNodes,isCircular,circularity)
%%
%==============================================================================
% Copyright (c) 2016 - 2017 Université de Lorraine & Luleå tekniska universitet
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
% Neither the name of the Université de Lorraine or Luleå tekniska universitet
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
%  A function to mesh a simply connected 2D rectangular geometry with elements of
%  shape and order of choice
%
%  Available:
%  - 1st and 2nd order quads
%  - 1st and 2nd order tris
%
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         lx - [M x 1] vector - Side length in each one of the M mesh regions in x-direction
%         ly - [N x 1] vector - Side length in each one of the N mesh regions in y-direction
%         Nx - [M x 1] vector - Number of ELEMENTS in each one of the M mesh regions in x-direction
%         Ny - [N x 1] vector - Number of ELEMENTS in each one of the N mesh regions in y-direction
%
%         For
%            --> 1st order quadrilaterals
%                    NxEquiv = Nx
%                    NyEquiv = Ny
%            --> 2nd order quadrilaterals
%                    NxEquiv = 2*Nx
%                    NyEquiv = 2*Ny
%            --> 1st order triangles
%                    NxEquiv = Nx
%                    NyEquiv = Ny
%            --> 2nd order triangles
%                    NxEquiv = 2*Nx
%                    NyEquiv = 2*Ny
%%


writeToLogFile(logfullfile,'In function: filterRectangularMesh\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;
% filter nodes if elements are not linear quadrilaterals
writeToLogFile(logfullfile,'Filtering nodes if elements are not linear quadrilaterals and building connectivity of edges and elements ...\n')
try
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      if isCircular
        if strcomp(circularity,'EW') || strcomp(circularity,'E-W') || strcomp(circularity,'EastWest') || strcomp(circularity,'East-West')
          filteredNodes = zeros(sum(Nx)*(sum(Ny)+1),size(baseNodes,2));
          edges = zeros((2*sum(Nx))*sum(Ny)+sum(Nx),2);
          elements = zeros(sum(Nx)*sum(Ny),4);
          for j=1:sum(Ny)
            filteredNodes((j-1)*sum(Nx)+1:j*sum(Nx),:) = baseNodes((j-1)*(sum(Nx)+1)+1:(j-1)*(sum(Nx)+1)+sum(Nx),:);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,2) = ((j-1)*sum(Nx)+2:j*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,3) = (j*sum(Nx)+2:(j+1)*sum(Nx))';
            elements(j*sum(Nx),2) = (j-1)*sum(Nx)+1;
            elements(j*sum(Nx),3) = j*sum(Nx)+1;
            elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          end
          j = sum(Ny)+1;
          filteredNodes((j-1)*sum(Nx)+1:j*sum(Nx),:) = baseNodes((j-1)*(sum(Nx)+1)+1:(j-1)*(sum(Nx)+1)+sum(Nx),:);
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
          nodes = filteredNodes;
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = (1:sum(Nx))'; %except corners
          nodesEASTside = []; %except corners
          nodesNORTHside = sum(Ny)*sum(Nx)+(1:sum(Nx))'; %except corners
          nodesWESTside = []; %except corners
          edgesSOUTHside = (1:sum(Nx))';
          edgesEASTside = [];
          edgesNORTHside = ((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx))';
          edgesWESTside = [];
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = (1:sum(Nx))'; %except corners
          elementsEASTside = []; %except corners
          elementsNORTHside = sum(Ny)*sum(Nx)+(1:sum(Nx))'; %except corners
          elementsWESTside = []; %except corners
        else % North-South circularity
          filteredNodes = zeros((sum(Nx)+1)*sum(Ny),size(baseNodes,2));
          edges = zeros((2*sum(Nx)+1)*sum(Ny),2);
          elements = zeros(sum(Nx)*sum(Ny),4);
          for j=1:sum(Ny)-1
            filteredNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),:) = baseNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),:);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
            edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
            edges(2*j*sum(Nx)+1,1) = j*sum(Nx)+1;
            edges(2*j*sum(Nx)+1,2) = (j+1)*sum(Nx)+1;
          end
          j = sum(Ny);
          filteredNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),:) = baseNodes((j-1)*(sum(Nx)+1)+1:j*(sum(Nx)+1),:);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (2:sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (1:sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges(2*j*sum(Nx)+1,1) = j*sum(Nx)+1;
          edges(2*j*sum(Nx)+1,2) = (j+1)*sum(Nx)+1;
          nodes = filteredNodes;
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = []; %except corners
          nodesEASTside = (sum(Nx)+1:sum(Nx)+1:(sum(Ny)-1)*(sum(Nx)+1)+sum(Nx)+1)'; %except corners
          nodesNORTHside = []; %except corners
          nodesWESTside = (1:sum(Nx)+1:(sum(Ny)-1)*(sum(Nx)+1)+1)'; %except corners
          edgesSOUTHside = [];
          edgesEASTside = (2*sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+2*sum(Nx)+1)';
          edgesNORTHside = [];
          edgesWESTside = (sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+sum(Nx)+1)';
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = []; %except corners
          elementsEASTside = (sum(Nx):sum(Nx):(sum(Ny)-1)*sum(Nx)+sum(Nx))'; %except corners
          elementsNORTHside = []; %except corners
          elementsWESTside = (1:sum(Nx):(sum(Ny)-1)*sum(Nx)+1)'; %except corners
        end
      else
        edges = zeros((2*sum(Nx)+1)*sum(Ny)+sum(Nx),2);
        elements = zeros(sum(Nx)*sum(Ny),4);
        nodes = baseNodes;
        for j=1:sum(Ny)
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*sum(Nx)+2:j*sum(Nx)+1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*sum(Nx)+1:j*sum(Nx))';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
          edges(2*j*sum(Nx)+1,1) = j*sum(Nx)+1;
          edges(2*j*sum(Nx)+1,2) = (j+1)*sum(Nx)+1;
        end
        j = sum(Ny)+1;
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = (j*sum(Nx)+1:(j+1)*sum(Nx))';
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = (j*sum(Nx)+2:(j+1)*sum(Nx)+1)';
        nodesSWcorner = 1;
        nodesSEcorner = sum(Nx)+1;
        nodesNEcorner = (sum(Nx)+1)*(sum(Ny)+1);
        nodesNWcorner = (sum(Nx)+1)*sum(Ny)+1;
        nodesSOUTHside = (2:sum(Nx))'; %except corners
        nodesEASTside = (2*(sum(Nx)+1):sum(Nx)+1:sum(Ny)*(sum(Nx)+1))'; %except corners
        nodesNORTHside = ((sum(Nx)+1)*sum(Ny)+2:(sum(Nx)+1)*(sum(Ny)+1)-1)'; %except corners
        nodesWESTside = (sum(Nx)+1+1:sum(Nx)+1:(sum(Ny)-1)*(sum(Nx)+1)+1)'; %except corners
        edgesSOUTHside = (1:sum(Nx))';
        edgesEASTside = (2*sum(Nx)+1:2*sum(Nx)+1:sum(Ny)*(2*sum(Nx)+1))';
        edgesNORTHside = ((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx))';
        edgesWESTside = (sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+sum(Nx)+1)';
        elementsSWcorner = 1;
        elementsSEcorner = sum(Nx);
        elementsNEcorner = sum(Nx)*sum(Ny);
        elementsNWcorner = sum(Nx)*(sum(Ny)-1)+1;
        elementsSOUTHside = (2:sum(Nx)-1)'; %except corners
        elementsEASTside = (2*sum(Nx):sum(Nx):sum(Nx)*(sum(Ny)-1))'; %except corners
        elementsNORTHside = (sum(Nx)*(sum(Ny)-1)+2:sum(Nx)*sum(Ny)-1)'; %except corners
        elementsWESTside = (sum(Nx)+1:sum(Nx):(sum(Ny)-2)*sum(Nx)+1)'; %except corners
      end
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      if isCircular
        if strcomp(circularity,'EW') || strcomp(circularity,'E-W') || strcomp(circularity,'EastWest') || strcomp(circularity,'East-West')
          filteredNodes = zeros((3*sum(Nx))*sum(Ny)+2*sum(Nx),size(baseNodes,2));
          edges = zeros((2*sum(Nx))*sum(Ny)+sum(Nx),3);
          elements = zeros(sum(Nx)*sum(Ny),8);
          for j=1:sum(Ny)
            filteredNodes((j-1)*(3*sum(Nx))+1:(j-1)*(3*sum(Nx))+2*sum(Nx),:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv),:);
            filteredNodes((j-1)*(3*sum(Nx))+2*sum(Nx)+1:j*(3*sum(Nx)),:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1+1:2:2*j*(sum(NxEquiv)+1)-2,:);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*(3*sum(Nx))+1:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,2) = ((j-1)*(3*sum(Nx))+3:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,3) = (j*(3*sum(Nx))+3:2:j*(3*sum(Nx))+2*sum(Nx)-1)';
            elements(j*sum(Nx),2) = (j-1)*(3*sum(Nx))+1;
            elements(j*sum(Nx),3) = j*(3*sum(Nx))+1;
            elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*(3*sum(Nx))+1:2:j*(3*sum(Nx))+2*sum(Nx)-1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),5) = ((j-1)*(3*sum(Nx))+2:2:(j-1)*(3*sum(Nx))+2*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx)-1,6) = (j-1)*(3*sum(Nx))+2*sum(Nx)+2:2:j*(3*sum(Nx))';
            elements(j*sum(Nx),6) = (j-1)*(3*sum(Nx))+2*sum(Nx)+1;
            elements((j-1)*sum(Nx)+1:j*sum(Nx),7) = (j*(3*sum(Nx))+2:2:j*(3*sum(Nx))+2*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),8) = (j-1)*(3*sum(Nx))+2*sum(Nx)+1:2:j*(3*sum(Nx))';
            edges((j-1)*2*sum(Nx)+1:(j-1)*2*sum(Nx)+sum(Nx),1) = ((j-1)*(3*sum(Nx))+1:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
            edges((j-1)*2*sum(Nx)+1:(j-1)*2*sum(Nx)+sum(Nx),2) = ((j-1)*(3*sum(Nx))+2:2:(j-1)*(3*sum(Nx))+2*sum(Nx))';
            edges((j-1)*2*sum(Nx)+1:(j-1)*2*sum(Nx)+sum(Nx),3) = ((j-1)*(3*sum(Nx))+3:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
            edges((j-1)*2*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*(3*sum(Nx))+1:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
            edges((j-1)*2*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j-1)*(3*sum(Nx))+2*sum(Nx)+1:2:j*(3*sum(Nx))';
            edges((j-1)*2*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),3) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          end
          j = sum(Ny)+1;
          filteredNodes((j-1)*(3*sum(Nx))+1:(j-1)*(3*sum(Nx))+2*sum(Nx),:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv),:);
          edges((j-1)*2*sum(Nx)+1:(j-1)*2*sum(Nx)+sum(Nx),1) = ((j-1)*(3*sum(Nx))+1:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
          edges((j-1)*2*sum(Nx)+1:(j-1)*2*sum(Nx)+sum(Nx),2) = ((j-1)*(3*sum(Nx))+2:2:(j-1)*(3*sum(Nx))+2*sum(Nx))';
          edges((j-1)*2*sum(Nx)+1:(j-1)*2*sum(Nx)+sum(Nx),3) = ((j-1)*(3*sum(Nx))+3:2:(j-1)*(3*sum(Nx))+2*sum(Nx)-1)';
          nodes = filteredNodes;
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = (1:2*sum(Nx))'; %except corners
          nodesEASTside = []; %except corners
          nodesNORTHside = (sum(Ny)*(3*sum(Nx))+1:sum(Ny)*(3*sum(Nx))+2*sum(Nx))'; %except corners
          nodesWESTside = []; %except corners
          edgesSOUTHside = (1:sum(Nx))';
          edgesEASTside = [];
          edgesNORTHside = (sum(Ny)*2*sum(Nx)+1:sum(Ny)*2*sum(Nx)+sum(Nx))';
          edgesWESTside = [];
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = (1:sum(Nx))'; %except corners
          elementsEASTside = []; %except corners
          elementsNORTHside = ((sum(Ny)-1)*sum(Nx)+1:(sum(Ny)*sum(Nx))'; %except corners
          elementsWESTside = []; %except corners
% ===============> circularity and boundaries completed on 18/08/2017 at 14:38
        else % North-South circularity
          filteredNodes = zeros((3*sum(Nx)+2)*sum(Ny),size(baseNodes,2));
          edges = zeros((2*sum(Nx)+1)*sum(Ny),3);
          elements = zeros(sum(Nx)*sum(Ny),8);
          for j=1:sum(Ny)-1
            filteredNodes((j-1)*(3*sum(Nx)+2)+1:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1,:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
            filteredNodes((j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:j*(3*sum(Nx)+2),:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1+1:2:2*j*(sum(NxEquiv)+1),:);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*(3*sum(Nx)+2)+3:2:j*(3*sum(Nx)+2)+2*sum(Nx)+1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),5) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),6) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+2:2:j*(3*sum(Nx)+2);
            elements((j-1)*sum(Nx)+1:j*sum(Nx),7) = (j*(3*sum(Nx)+2)+2:2:j*(3*sum(Nx)+2)+2*sum(Nx))';
            elements((j-1)*sum(Nx)+1:j*sum(Nx),8) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
            edges((j-1)*(2*sum(Nx)+1)+1:(j-1)*(2*sum(Nx)+1)+sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
            edges((j-1)*(2*sum(Nx)+1)+1:(j-1)*(2*sum(Nx)+1)+sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
            edges((j-1)*(2*sum(Nx)+1)+1:(j-1)*(2*sum(Nx)+1)+sum(Nx),3) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
            edges((j-1)*(2*sum(Nx)+1)+sum(Nx)+1:(j-1)*(2*sum(Nx)+1)+2*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
            edges((j-1)*(2*sum(Nx)+1)+sum(Nx)+1:(j-1)*(2*sum(Nx)+1)+2*sum(Nx),2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
            edges((j-1)*(2*sum(Nx)+1)+sum(Nx)+1:(j-1)*(2*sum(Nx)+1)+2*sum(Nx),3) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
            edges((j-1)*(2*sum(Nx)+1)+2*sum(Nx)+1,1) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
            edges((j-1)*(2*sum(Nx)+1)+2*sum(Nx)+1,2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx);
            edges((j-1)*(2*sum(Nx)+1)+2*sum(Nx)+1,3) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
          end
          j = sum(Ny);
          filteredNodes((j-1)*(3*sum(Nx)+2)+1:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1,:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:j*(3*sum(Nx)+2),:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1+1:2:2*j*(sum(NxEquiv)+1),:);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (3:2:2*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (1:2:2*sum(Nx)-1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),5) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),6) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+2:2:j*(3*sum(Nx)+2);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),7) = (2:2:2*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),8) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
          edges((j-1)*(2*sum(Nx)+1)+1:(j-1)*(2*sum(Nx)+1)+sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*(2*sum(Nx)+1)+1:(j-1)*(2*sum(Nx)+1)+sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
          edges((j-1)*(2*sum(Nx)+1)+1:(j-1)*(2*sum(Nx)+1)+sum(Nx),3) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          edges((j-1)*(2*sum(Nx)+1)+sum(Nx)+1:(j-1)*(2*sum(Nx)+1)+2*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*(2*sum(Nx)+1)+sum(Nx)+1:(j-1)*(2*sum(Nx)+1)+2*sum(Nx),2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
          edges((j-1)*(2*sum(Nx)+1)+sum(Nx)+1:(j-1)*(2*sum(Nx)+1)+2*sum(Nx),3) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*(2*sum(Nx)+1)+2*sum(Nx)+1,1) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
          edges((j-1)*(2*sum(Nx)+1)+2*sum(Nx)+1,2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx);
          edges((j-1)*(2*sum(Nx)+1)+2*sum(Nx)+1,3) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
          j = sum(Ny)+1;
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = (1:2:2*sum(Nx)-1)';
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = (2:2:2*sum(Nx))';
          edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),3) = (1:2:2*sum(Nx)-1)';
          nodes = filteredNodes;
        end
        nodesSWcorner = -1;
        nodesSEcorner = -1;
        nodesNEcorner = -1;
        nodesNWcorner = -1;
        nodesSOUTHside = []; %except corners
        nodesEASTside = (2*sum(Nx)+1:(3*sum(Nx)+2):(sum(Ny)-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)'; %except corners
        nodesNORTHside = []; %except corners
        nodesWESTside = (1:(3*sum(Nx)+2):(sum(Ny)-1)*(3*sum(Nx)+2)+1)'; %except corners
        edgesSOUTHside = [];
        edgesEASTside = (2*sum(Nx)+1:2*sum(Nx)+1:sum(Ny)*(2*sum(Nx)+1))';
        edgesNORTHside = [];
        edgesWESTside = (sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+sum(Nx)+1)';
        elementsSWcorner = -1;
        elementsSEcorner = -1;
        elementsNEcorner = -1;
        elementsNWcorner = -1;
        elementsSOUTHside = []; %except corners
        elementsEASTside = (sum(Nx):sum(Nx):(sum(Ny)*sum(Nx))'; %except corners
        elementsNORTHside = []; %except corners
        elementsWESTside = (1:sum(Nx):(sum(Ny)-1)*sum(Nx)+1)'; %except corners
% ===============> circularity and boundaries completed on 18/08/2017 at 14:46
      else
        filteredNodes = zeros((3*sum(Nx)+2)*sum(Ny)+2*sum(Nx)+1,size(baseNodes,2));
        edges = zeros((2*sum(Nx)+1)*sum(Ny)+sum(Nx),3);
        elements = zeros(sum(Nx)*sum(Ny),8);
        for j=1:sum(Ny)
          filteredNodes((j-1)*(3*sum(Nx)+2)+1:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1,:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:j*(3*sum(Nx)+2),:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1+1:2:2*j*(sum(NxEquiv)+1),:);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),3) = (j*(3*sum(Nx)+2)+3:2:j*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),4) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),5) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),6) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+2:2:j*(3*sum(Nx)+2);
          elements((j-1)*sum(Nx)+1:j*sum(Nx),7) = (j*(3*sum(Nx)+2)+2:2:j*(3*sum(Nx)+2)+2*sum(Nx))';
          elements((j-1)*sum(Nx)+1:j*sum(Nx),8) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
          edges((j-1)*sum(Nx)+1:(j-1)*sum(Nx)+sum(Nx),3) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1+1:2:j*(3*sum(Nx)+1);
          edges((j-1)*sum(Nx)+sum(Nx)+1:2*j*sum(Nx),3) = (j*(3*sum(Nx)+2)+1:2:j*(3*sum(Nx)+2)+2*sum(Nx)-1)';
          edges(2*j*sum(Nx)+1,1) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
          edges(2*j*sum(Nx)+1,2) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx);
          edges(2*j*sum(Nx)+1,3) = (j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1;
        end
        j = sum(Ny)+1;
        filteredNodes((j-1)*(3*sum(Nx)+2)+1:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1,:) = baseNodes(2*(j-1)*(sum(NxEquiv)+1)+1:2*(j-1)*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),1) = ((j-1)*(3*sum(Nx)+2)+1:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)-1)';
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),2) = ((j-1)*(3*sum(Nx)+2)+2:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx))';
        edges((2*sum(Nx)+1)*sum(Ny)+1:(2*sum(Nx)+1)*sum(Ny)+sum(Nx),3) = ((j-1)*(3*sum(Nx)+2)+3:2:(j-1)*(3*sum(Nx)+2)+2*sum(Nx)+1)';
        nodes = filteredNodes;
        nodesSWcorner = 1;
        nodesSEcorner = 2*sum(Nx)+1;
        nodesNEcorner = sum(Ny)*(3*sum(Nx)+2)+2*sum(Nx);
        nodesNWcorner = sum(Ny)*(3*sum(Nx)+2)+1;
        nodesSOUTHside = (2:2*sum(Nx))'; %except corners
        nodesEASTside = (3*sum(Nx)+2:(3*sum(Nx)+2):(sum(Ny)-1)*(3*sum(Nx)+2))'; %except corners
        nodesNORTHside = ((sum(Ny)-1)*(3*sum(Nx)+2)+2:(sum(Ny)-1)*(3*sum(Nx)+2)+2*sum(Nx))'; %except corners
        nodesWESTside = (2*sum(Nx)+2:(3*sum(Nx)+2):(sum(Ny)-2)*(3*sum(Nx)+2)+2*sum(Nx)+2)'; %except corners
        edgesSOUTHside = (1:sum(Nx))';
        edgesEASTside = (2*sum(Nx)+1:2*sum(Nx)+1:sum(Ny)*(2*sum(Nx)+1))';
        edgesNORTHside = (sum(Ny)*(2*sum(Nx)+1)+1:sum(Ny)*(2*sum(Nx)+1)+sum(Nx))';
        edgesWESTside = (sum(Nx)+1:2*sum(Nx)+1:(sum(Ny)-1)*(2*sum(Nx)+1)+sum(Nx)+1)';
        elementsSWcorner = 1;
        elementsSEcorner = sum(Nx);
        elementsNEcorner = sum(Ny)*sum(Nx);
        elementsNWcorner = (sum(Ny)-1)*sum(Nx);
        elementsSOUTHside = (2:sum(Nx)-1)'; %except corners
        elementsEASTside = (2*sum(Nx):sum(Nx):(sum(Ny)-1)*sum(Nx))'; %except corners
        elementsNORTHside = ((sum(Ny)-1)*sum(Nx)+2:sum(Ny)*sum(Nx)-1)'; %except corners
        elementsWESTside = (sum(Nx)+1:sum(Nx):(sum(Ny)-2)*sum(Nx)+1)'; %except corners
      end
    end
% ===============> circularity and boundaries completed on 18/08/2017 at 15:31
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      if isCircular
        if strcomp(circularity,'EW') || strcomp(circularity,'E-W') || strcomp(circularity,'EastWest') || strcomp(circularity,'East-West')
          NxTri = Nx/2;
          NyTri = Ny/2;
          filteredNodes = zeros((2*sum(NxTri))*sum(NyTri)+sum(NxTri),size(baseNodes,2));
          edges = zeros((6*sum(NxTri))*sum(NyTri)+sum(NxTri),2);
          elements = zeros(sum(Nx)*sum(Ny),3);
% ===============> OK
          for j=1:sum(NyTri)
            filteredNodes((j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri),:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv),:);
            filteredNodes((j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri)),:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1+2:2:j*2*(sum(NxEquiv)+1)-1,:);
% ===============> OK
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri)-1,2) = (j-1)*(2*sum(NxTri))+2:(j-1)*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(2*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri)-1,1) = (j-1)*(2*sum(NxTri))+2:(j-1)*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri)-1,2) = j*(2*sum(NxTri))+2:j*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(2*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri)-1,2) = j*(2*sum(NxTri))+2:j*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(2*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(2*sum(NxTri))+1:j*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(2*sum(NxTri))+1:j*(2*sum(NxTri))+sum(NxTri);
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri),1)                  = (j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri)-1,2)                = (j-1)*(2*sum(NxTri))+2:(j-1)*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+sum(NxTri),2) = (j-1)*(2*sum(NxTri))+1;
            edges((j-1)*(6*sum(NxTri))+sum(NxTri)+1:(j-1)*(6*sum(NxTri))+2*sum(NxTri),1)     = (j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+sum(NxTri)+1:(j-1)*(6*sum(NxTri))+2*sum(NxTri),2)     = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            edges((j-1)*(6*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+3*sum(NxTri),1)   = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            edges((j-1)*(6*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+3*sum(NxTri)-1,2) = (j-1)*(2*sum(NxTri))+2:(j-1)*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+3*sum(NxTri),2) = (j-1)*(2*sum(NxTri))+1;
            edges((j-1)*(6*sum(NxTri))+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+4*sum(NxTri),1)   = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            edges((j-1)*(6*sum(NxTri))+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+4*sum(NxTri)-1,2) = j*(2*sum(NxTri))+2:j*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+4*sum(NxTri),2) = j*(2*sum(NxTri))+1;
            edges((j-1)*(6*sum(NxTri))+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+5*sum(NxTri),1)   = j*(2*sum(NxTri))+1:j*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+5*sum(NxTri),2)   = (j-1)*(2*sum(NxTri))+sum(NxTri)+1:j*(2*sum(NxTri));
            edges((j-1)*(6*sum(NxTri))+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+6*sum(NxTri),1)   = (j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+6*sum(NxTri),2)   = j*(2*sum(NxTri))+1:j*(2*sum(NxTri))+sum(NxTri);
% ===============> OK
          end
          j = sum(NyTri)+1;
          filteredNodes((j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri),:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv),:);
          edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri),1) = (j-1)*(2*sum(NxTri))+1:(j-1)*(2*sum(NxTri))+sum(NxTri);
          edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri)-1,2) = (j-1)*(2*sum(NxTri))+2:(j-1)*(2*sum(NxTri))+sum(NxTri);
          edges((j-1)*(6*sum(NxTri))+sum(NxTri),2) = (j-1)*(2*sum(NxTri))+1;
          nodes = filteredNodes;
% ===============> OK
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = (1:sum(NxTri))'; %except corners
          nodesEASTside = []; %except corners
          nodesNORTHside = (sum(NyTri)*2*sum(NxTri)+1:sum(NyTri)*2*sum(NxTri)+sum(NxTri))'; %except corners
          nodesWESTside = []; %except corners
          edgesSOUTHside = (1:sum(NxTri))';
          edgesEASTside = [];
          edgesNORTHside = (sum(NyTri)*(6*sum(NxTri))+1:sum(NyTri)*(6*sum(NxTri))+sum(NxTri))';
          edgesWESTside = [];
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = (1:sum(NxTri))'; %except corners
          elementsEASTside = []; %except corners
          elementsNORTHside = ((sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri)+1:(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri))'; %except corners
          elementsWESTside = []; %except corners
% ===============> circularity and boundaries completed on 18/08/2017 at 18:16
        else % South-North circularity
          NxTri = Nx/2;
          NyTri = Ny/2;
          filteredNodes = zeros((2*sum(NxTri)+1)*sum(NyTri),size(baseNodes,2));
          edges = zeros((6*sum(NxTri)+1)*sum(NyTri),2);
          elements = zeros(sum(Nx)*sum(Ny),3);
% ===============> OK
          for j=1:sum(NyTri)-1
            filteredNodes((j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1,:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
            filteredNodes((j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1),:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1+2:2:j*2*(sum(NxEquiv)+1)-1,:);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
            edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = j*(2*sum(NxTri)+1)+sum(NxTri)+1;
          end
% ===============> OK
          j = sum(NyTri);
          filteredNodes((j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1,:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1),:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1+2:2:j*2*(sum(NxEquiv)+1)-1,:);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = 2:sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = 2:sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = 1:sum(NxTri);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = 1:sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = 2:sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = 1:sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = 1:sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = sum(NxTri)+1;
          nodes = filteredNodes;
% ===============> OK
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = []; %except corners
          nodesEASTside = (sum(NxTri)+1:2*sum(NxTri)+1:(sum(NyTri)-1)*(2*sum(NxTri)+1)+sum(NxTri)+1)'; %except corners
          nodesNORTHside = []; %except corners
          nodesWESTside = (1:2*sum(NxTri)+1:(sum(NyTri)-1)*(2*sum(NxTri)+1)+1)'; %except corners
          edgesSOUTHside = [];
          edgesEASTside = (6*sum(NxTri)+1:6*sum(NxTri)+1:sum(NyTri)*(6*sum(NxTri)+1))';
          edgesNORTHside = [];
          edgesWESTside = (5*sum(NxTri)+1:6*sum(NxTri)+1:(sum(NyTri)-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1)';
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = []; %except corners
          elementsEASTside = (2*sum(NxTri):4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri))'; %except corners
          elementsNORTHside = []; %except corners
          elementsWESTside = (3*sum(NxTri)+1:4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri)+1)'; %except corners
        end
% ===============> circularity and boundaries completed on 18/08/2017 at 19:12
      else
        NxTri = Nx/2;
        NyTri = Ny/2;
        filteredNodes = zeros((2*sum(NxTri)+1)*sum(NyTri)+sum(NxTri)+1,size(baseNodes,2));
        edges = zeros((6*sum(NxTri)+1)*sum(NyTri)+sum(NxTri),2);
        elements = zeros(sum(Nx)*sum(Ny),3);
        for j=1:sum(NyTri)
          filteredNodes((j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1,:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1),:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1+2:2:j*2*(sum(NxEquiv)+1)-1,:);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = j*(2*sum(NxTri)+1)+2:j*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1+1:j*(2*sum(NxTri)+1);
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = j*(2*sum(NxTri)+1)+1:j*(2*sum(NxTri)+1)+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = j*(2*sum(NxTri)+1)+sum(NxTri)+1;
        end
        j = sum(NyTri)+1;
        filteredNodes((j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1,:) = baseNodes((j-1)*2*(sum(NxEquiv)+1)+1:2:(j-1)*2*(sum(NxEquiv)+1)+sum(NxEquiv)+1,:);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(2*sum(NxTri)+1)+1:(j-1)*(2*sum(NxTri)+1)+sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(2*sum(NxTri)+1)+2:(j-1)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        nodes = filteredNodes;
        nodesSWcorner = 1;
        nodesSEcorner = sum(NxTri)+1;
        nodesNEcorner = sum(NyTri)*(2*sum(NxTri)+1)+sum(NxTri)+1;
        nodesNWcorner = sum(NyTri)*(2*sum(NxTri)+1)+1;
        nodesSOUTHside = (2:sum(NxTri))'; %except corners
        nodesEASTside = (3*sum(NxTri)+2:2*sum(NxTri)+1:(sum(NyTri)-1)*(2*sum(NxTri)+1)+sum(NxTri)+1)'; %except corners
        nodesNORTHside = (sum(NyTri)*(2*sum(NxTri)+1)+2:sum(NyTri)*(2*sum(NxTri)+1)+sum(NxTri))'; %except corners
        nodesWESTside = ((2*sum(NxTri)+1)+1:(2*sum(NxTri)+1):(sum(NyTri)-1)*(2*sum(NxTri)+1)+1)'; %except corners
        edgesSOUTHside = (1:sum(NxTri))';
        edgesEASTside = (6*sum(NxTri)+1:6*sum(NxTri)+1:)';
        edgesNORTHside = (sum(NyTri)*(6*sum(NxTri)+1)+1:sum(NyTri)*(6*sum(NxTri)+1)+sum(NxTri))';
        edgesWESTside = (5*sum(NxTri)+1:6*sum(NxTri)+1:(sum(NyTri)-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1)';
        elementsSWcorner = -1;
        elementsSEcorner = -1;
        elementsNEcorner = -1;
        elementsNWcorner = -1;
        elementsSOUTHside = (1:sum(NxTri))'; %except corners
        elementsEASTside = (3*sum(NxTri)+1:4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri)+1)'; %except corners
        elementsNORTHside = ((sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri)+1:(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri))'; %except corners
        elementsWESTside = (2*sum(NxTri):4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri))'; %except corners
      end
% ===============> circularity and boundaries completed on 18/08/2017 at 20:01
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      if isCircular
        if strcomp(circularity,'EW') || strcomp(circularity,'E-W') || strcomp(circularity,'EastWest') || strcomp(circularity,'East-West')
          NxTri = Nx/2;
          NyTri = Ny/2;
          filteredNodes = zeros((8*sum(NxTri))*sum(NyTri)+2*sum(NxTri),size(baseNodes,2));
          edges = zeros((6*sum(NxTri))*sum(NyTri)+sum(NxTri),3);
          elements = zeros(sum(Nx)*sum(Ny),6);
% ===============> OK
          for j=1:sum(NyTri)
            filteredNodes((j-1)*(8*sum(NxTri))+1:(j-1)*(8*sum(NxTri))+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv),:);
            filteredNodes((j-1)*(8*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(8*sum(NxTri))+4*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv),:);
            filteredNodes((j-1)*(8*sum(NxTri))+4*sum(NxTri)+1:(j-1)*(8*sum(NxTri))+6*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv),:);
            filteredNodes((j-1)*(8*sum(NxTri))+6*sum(NxTri)+1:(j-1)*(8*sum(NxTri))+8*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv),:);
% ===============> OK
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri)-1,2) = (j-1)*(8*sum(NxTri))+3:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(8*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),4) = (j-1)*(8*sum(NxTri))+2:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),5) = (j-1)*(8*sum(NxTri))+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),6) = (j-1)*(8*sum(NxTri))+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)+2*sum(NxTri)-1;
% ===============> OK
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri)-1,1) = (j-1)*(8*sum(NxTri))+3:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri)-1,2) = j*(8*sum(NxTri))+3:2:j*(8*sum(NxTri))+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(8*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri)-1,4) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+3:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            elements(:(j-1)*4*sum(NxTri)+2*sum(NxTri),4) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),5) = (j-1)*(8*sum(NxTri))+6*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+7*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),6) = (j-1)*(8*sum(NxTri))+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+4*sum(NxTri);
% ===============> OK
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri)-1,2) = j*(8*sum(NxTri))+3:2:j*(8*sum(NxTri))+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(8*sum(NxTri))+1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(8*sum(NxTri))+1:2:j*(8*sum(NxTri))+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),4) = (j-1)*(8*sum(NxTri))+6*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+7*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),5) = j*(8*sum(NxTri))+2:2:j*(8*sum(NxTri))+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),6) = (j-1)*(8*sum(NxTri))+6*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+7*sum(NxTri)-1;
% ===============> OK
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(8*sum(NxTri))+1:2:j*(8*sum(NxTri))+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),4) = (j-1)*(8*sum(NxTri))+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+4*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),5) = (j-1)*(8*sum(NxTri))+6*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+7*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),6) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri)-1;
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri),2) = (j-1)*(8*sum(NxTri))+2:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri)-1,3) = (j-1)*(8*sum(NxTri))+3:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+sum(NxTri),3) = (j-1)*(8*sum(NxTri))+1;
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+sum(NxTri)+1:(j-1)*(6*sum(NxTri))+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+sum(NxTri)+1:(j-1)*(6*sum(NxTri))+2*sum(NxTri),2) = (j-1)*(8*sum(NxTri))+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+4*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+sum(NxTri)+1:(j-1)*(6*sum(NxTri))+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+3*sum(NxTri),2) = (j-1)*(8*sum(NxTri))+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+4*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+3*sum(NxTri),3) = (j-1)*(8*sum(NxTri))+3:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+3*sum(NxTri),3) = (j-1)*(8*sum(NxTri))+1;
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri))+6*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+7*sum(NxTri);
            edges((j-1)*(6*sum(NxTri))+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+4*sum(NxTri)-1,3) = j*(8*sum(NxTri))+3:2:j*(8*sum(NxTri))+2*sum(NxTri);
            edges((j-1)*(6*sum((j-1)*(6*sum(NxTri))+4*sum(NxTri),3) = j*(8*sum(NxTri))+1;
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+5*sum(NxTri),1) = j*(8*sum(NxTri))+1:2:j*(8*sum(NxTri))+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+5*sum(NxTri),2) = (j-1)*(8*sum(NxTri))+6*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+7*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+5*sum(NxTri),3) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri);
% ===============> OK
            edges((j-1)*(6*sum(NxTri))+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+6*sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+6*sum(NxTri),2) = (j-1)*(8*sum(NxTri))+4*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri))+6*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri))+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri))+6*sum(NxTri),3) = j*(8*sum(NxTri))+1:2:j*(8*sum(NxTri))+2*sum(NxTri)-1;
% ===============> OK
          end
          j = sum(NyTri)+1;
          filteredNodes((j-1)*(8*sum(NxTri))+1:(j-1)*(8*sum(NxTri))+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv),:);
          edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri),1) = (j-1)*(8*sum(NxTri))+1:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri),2) = (j-1)*(8*sum(NxTri))+2:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri)-1,3) = (j-1)*(8*sum(NxTri))+3:2:(j-1)*(8*sum(NxTri))+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri))+sum(NxTri),3) = (j-1)*(8*sum(NxTri))+1;
          nodes = filteredNodes;
  % ===============> OK
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = (1:2*sum(NxTri))'; %except corners
          nodesEASTside = []; %except corners
          nodesNORTHside = (sum(NyTri)*(8*sum(NxTri))+1:sum(NyTri)*(8*sum(NxTri))+2*sum(NxTri))'; %except corners
          nodesWESTside = []; %except corners
          edgesSOUTHside = (1:sum(NxTri))';
          edgesEASTside = [];
          edgesNORTHside = ((j-1)*(6*sum(NxTri))+1:(j-1)*(6*sum(NxTri))+sum(NxTri))';
          edgesWESTside = [];
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = (1:sum(NxTri))'; %except corners
          elementsEASTside = []; %except corners
          elementsNORTHside = ((sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri)+1:(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri))'; %except corners
          elementsWESTside = []; %except corners
% ===============> circularity and boundaries completed on 21/08/2017 at 18:34
        else % South-North circularity
          NxTri = Nx/2;
          NyTri = Ny/2;
          filteredNodes = zeros((8*sum(NxTri)+2)*sum(NyTri),size(baseNodes,2));
          edges = zeros((6*sum(NxTri)+1)*sum(NyTri),3);
          elements = zeros(sum(Nx)*sum(Ny),6);
% ===============> OK
          for j=1:sum(NyTri)-1
            filteredNodes((j-1)*(8*sum(NxTri)+2)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1,:);
            filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv),:);
            filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1,:);
            filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv),:);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),5) = j*(8*sum(NxTri)+2)+2:2:j*(8*sum(NxTri)+2)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
            elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),3) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
            edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
            edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
            edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,3) = j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          end
% ===============> OK
          j = sum(NyTri);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv),:);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv),:);
% ===============> OK
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
% ===============> OK
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = 3:2:2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
% ===============> OK
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = 3:2:2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = 1:2:2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),5) = 2:2:2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
% ===============> OK
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = 1:2:2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),3) = 3:2:2*sum(NxTri)+1;
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = 1:2:2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),3) = 1:2:2*sum(NxTri)-1;
% ===============> OK
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,3) = 2*sum(NxTri)+1;
% ===============> OK
          j = sum(NyTri)+1;
          nodes = filteredNodes;
% ===============> OK
          nodesSWcorner = -1;
          nodesSEcorner = -1;
          nodesNEcorner = -1;
          nodesNWcorner = -1;
          nodesSOUTHside = []; %except corners
          nodesEASTside = (2*sum(NxTri)+1:8*sum(NxTri)+2:(sum(NyTri)-1)*(8*sum(NxTri)+2)+6*sum(NxTri)+2)'; %except corners
          nodesNORTHside = []; %except corners
          nodesWESTside = (1:8*sum(NxTri)+2:(sum(NyTri)-1)*(8*sum(NxTri)+2)+4*sum(NxTri)+1+1)'; %except corners
          edgesSOUTHside = [];
          edgesEASTside = (6*sum(NxTri)+1:6*sum(NxTri)+1:sum(NyTri)*(6*sum(NxTri)+1))';
          edgesNORTHside = [];
          edgesWESTside = (5*sum(NxTri)+1:6*sum(NxTri)+1:(sum(NyTri)-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1)';
          elementsSWcorner = -1;
          elementsSEcorner = -1;
          elementsNEcorner = -1;
          elementsNWcorner = -1;
          elementsSOUTHside = []; %except corners
          elementsEASTside = (2*sum(NxTri):4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri))'; %except corners
          elementsNORTHside = []; %except corners
          elementsWESTside = (3*sum(NxTri)+1:4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri)+1)'; %except corners
        end
% ===============> circularity and boundaries completed on 22/08/2017 at 10:41
      else
        NxTri = Nx/2;
        NyTri = Ny/2;
        filteredNodes = zeros((8*sum(NxTri)+2)*sum(NyTri)+2*sum(NxTri)+1,size(baseNodes,2));
        edges = zeros((6*sum(NxTri)+1)*sum(NyTri)+sum(NxTri),3);
        elements = zeros(sum(Nx)*sum(Ny),6);
        for j=1:sum(NyTri)
          filteredNodes((j-1)*(8*sum(NxTri)+2)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv),:);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1,:);
          filteredNodes((j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2*sum(NxTri),:) = baseNodes((j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+2:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv)+1+sum(NxEquiv),:);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+1:(j-1)*4*sum(NxTri)+sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),2) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+sum(NxTri)+1:(j-1)*4*sum(NxTri)+2*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),2) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),5) = j*(8*sum(NxTri)+2)+2:2:j*(8*sum(NxTri)+2)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+2*sum(NxTri)+1:(j-1)*4*sum(NxTri)+3*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),4) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),5) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
          elements((j-1)*4*sum(NxTri)+3*sum(NxTri)+1:(j-1)*4*sum(NxTri)+4*sum(NxTri),6) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+2*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+2*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+3*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+3*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+4*sum(NxTri),3) = j*(8*sum(NxTri)+2)+3:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),1) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1+sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+4*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+5*sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri);
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1:(j-1)*(6*sum(NxTri)+1)+6*sum(NxTri),3) = j*(8*sum(NxTri)+2)+1:2:j*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,1) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,2) = (j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1+2*sum(NxTri)+2*sum(NxTri)+1;
          edges((j-1)*(6*sum(NxTri)+1)+6*sum(NxTri)+1,3) = j*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        end
% ===============> OK
        j = sum(NyTri)+1;
        filteredNodes((j-1)*(8*sum(NxTri)+2)+1:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1,:) = baseNodes((j-1)*4*sum(NxEquiv)+1:2:(j-1)*4*sum(NxEquiv)+sum(NxEquiv)+1,:);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),1) = (j-1)*(8*sum(NxTri)+2)+1:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)-1;
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),2) = (j-1)*(8*sum(NxTri)+2)+2:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri);
        edges((j-1)*(6*sum(NxTri)+1)+1:(j-1)*(6*sum(NxTri)+1)+sum(NxTri),3) = (j-1)*(8*sum(NxTri)+2)+3:2:(j-1)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        nodes = filteredNodes;
% ===============> OK
        nodesSWcorner = 1;
        nodesSEcorner = 2*sum(NxTri)+1;
        nodesNEcorner = sum(NyTri)*(8*sum(NxTri)+2)+2*sum(NxTri)+1;
        nodesNWcorner = sum(NyTri)*(8*sum(NxTri)+2)+1;
% ===============> OK
        nodesSOUTHside = (2:2*sum(NxTri))'; %except corners
        nodesEASTside = (6*sum(NxTri)+2:8*sum(NxTri)+2:(sum(NyTri)-1)*(8*sum(NxTri)+2)+6*sum(NxTri)+2)'; %except corners
        nodesNORTHside = (sum(NyTri)*(8*sum(NxTri)+2)+2:sum(NyTri)*(8*sum(NxTri)+2)+2*sum(NxTri))'; %except corners
        nodesWESTside = (4*sum(NxTri)+1+1:8*sum(NxTri)+2:(sum(NyTri)-1)*(8*sum(NxTri)+2)+(4*sum(NxTri)+1+1)'; %except corners
% ===============> OK
        edgesSOUTHside = (1:sum(NxTri))';
        edgesEASTside = (6*sum(NxTri)+1:(6*sum(NxTri)+1):sum(NyTri)*(6*sum(NxTri)+1))';
        edgesNORTHside = (sum(NyTri)*(6*sum(NxTri)+1)+1:sum(NyTri)*(6*sum(NxTri)+1)+sum(NxTri))';
        edgesWESTside = (5*sum(NxTri)+1:(6*sum(NxTri)+1):(sum(NyTri)-1)*(6*sum(NxTri)+1)+5*sum(NxTri)+1)';
% ===============> OK
        elementsSWcorner = 1;
        elementsSEcorner = -1;
        elementsNEcorner = -1;
        elementsNWcorner = -1;
% ===============> OK
        elementsSOUTHside = (1:sum(NxTri))'; %except corners
        elementsEASTside = (2*sum(NxTri):4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri))'; %except corners
        elementsNORTHside = ((sum(NyTri)-1)*4*sum(NxTri)+2*sum(NxTri)+1:(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri))'; %except corners
        elementsWESTside = (3*sum(NxTri)+1:4*sum(NxTri):(sum(NyTri)-1)*4*sum(NxTri)+3*sum(NxTri)+1)'; %except corners
      end
% ===============> circularity and boundaries completed on 22/08/2017 at 11:38
    end
  end
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier,'\n'])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end
writeToLogFile(logfullfile,['... done.','\n'])

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: filterRectangularMesh\n')

return
