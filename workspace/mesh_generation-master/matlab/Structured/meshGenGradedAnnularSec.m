function[nodes,elements,edges,...
    nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
    edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
    elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside]=meshGenGradedAnnularSec(logfullfile,elType,elOrder,isClosed,x0,y0,Rin,theta0,ltheta,lR,deltas,NR)
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
% Neither the name of the Université de Lorraine & Luleå tekniska universitet
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
%  A function to mesh a non-simply connected, i.e. with rectangular holes,
%  rectangular geometry with rectangular elements
%
%  Available:
%  - 1st and 2nd order quads
%  - 1st and 2nd order tris
%
%  Input: x0   - scalar - x-coordinate of center
%         y0   - scalar - y-coordinate of center
%         Rin  - scalar - Radius
%         isClosed - scalar - Equal to 1 or True if the annular section is closed (ring)
%         ltheta - [M x 1] vector - Angular aperture of each one of the M mesh regions along the tangential direction
%         lR     - [N x 1] vector - Length of the N mesh regions along the radial direction
%         deltas - [M x 1] vector - Angular aperture of ELEMENTS in each one of the M mesh regions along the tangential direction
%         NR     - [M x 1] vector - Number of ELEMENTS in each one of the N mesh regions along the radial direction
%
%%

writeToLogFile(logfullfile,'In function: meshGenGradedAnnularSec\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

writeToLogFile(logfullfile,['Creating mesh in computational space ...','\n'])
try
  Nx = ltheta./deltas;
  Ny = NR;
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order quadrilaterals','\n'])
      NxEquiv = Nx;
      NyEquiv = Ny;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order quadrilaterals','\n'])
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
    end
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order triangles','\n'])
      NxEquiv = Nx;
      NyEquiv = Ny;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order triangles','\n'])
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
    end
  end
  writeToLogFile(logfullfile,['    Calling function ', 'gradedRectangle',' ...\n']);
  baseMesh = gradedRectangle(logfullfile,theta0+0.5*sum(ltheta),Rin+0.5*sum(lR),ltheta,lR,NxEquiv,NyEquiv);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

writeToLogFile(logfullfile,['Creating mesh in physical space ...','\n'])
try
  mesh = zeros(length(baseMesh),2);
  for i=1:length(mesh)
    mesh(i,1) = x0 + baseMesh(i,2)*cos(baseMesh(i,1));
    mesh(i,2) = y0 + baseMesh(i,2)*sin(baseMesh(i,1));
  end
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

% filter nodes if elements are not linear quadrilaterals
writeToLogFile(logfullfile,'Filtering nodes if elements are not linear quadrilaterals ...\n')
try
  writeToLogFile(logfullfile,['Calling function ', 'filterRectangularMesh',' ...\n']);
  if isCircular
    [nodes,elements,edges,...
        nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
        edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
        elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside] = filterRectangularMesh(logfullfile,...
                                                                                                                                                                        elType,elOrder,Nx,Nx,NxEquiv,NxEquiv,...
                                                                                                                                                                        mesh,1,'EW')
  else
    [nodes,elements,edges,...
      nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
      edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
      elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside] = filterRectangularMesh(logfullfile,...
                                                                                                                                                                      elType,elOrder,Nx,Nx,NxEquiv,NxEquiv,...
                                                                                                                                                                      mesh,0,'none')
  end
  writeToLogFile(logfullfile,['... done.','\n'])                                                                                                                                                                                                                                                                          baseNodes,isCircular,circularity)
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier,'\n'])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: meshGenGradedAnnularSec\n')

return
