function[nodes,elements,edges,...
         nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
         edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
         elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside]=higherOrderGradedRectangle(logfullfile,latexFolder,elType,elOrder,isCircular,circularity,x0,y0,lx,ly,Nx,Ny)
%%
%==============================================================================
% Copyright (c) 2016-2017 Universit� de Lorraine & Lule� tekniska universitet
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
% Neither the name of the Universit� de Lorraine or Lule� tekniska universitet
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
%         holes - [H x 4] matrix - H is the number of holes; for each hole the following data must pe provided:
%                                  xC - scalar - x-coordinate of hole's center
%                                  yC - scalar - y-coordinate of hole's center
%                                  xL - scalar - half-length of side parallel to x-axis
%                                  yL - scalar - half-length of side parallel to y-axis
%%

writeToLogFile(logfullfile,'In function: higherOrderGradedRectangle\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

% calculate equivalent number of linear quadrilateral elements depending on type and order of elements chosen
writeToLogFile(logfullfile,'Calculating equivalent number of linear quadrilateral elements based on type and order of elements chosen ...\n')
try
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      NxEquiv = Nx;
      NyEquiv = Ny;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
    end
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      NxEquiv = Nx;
      NyEquiv = Ny;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
    end
  end
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier,'\n'])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end
writeToLogFile(logfullfile,['... done.','\n'])

% create equivalent mesh with linear quadrilateral elements
writeToLogFile(logfullfile,['Creating equivalent mesh with linear quadrilateral elements ...','\n'])
try
  baseNodes = gradedRectangle(logfullfile,x0,y0,lx,ly,NxEquiv,NyEquiv)
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end
writeToLogFile(logfullfile,['... done.','\n'])

% filter nodes if elements are not linear quadrilaterals
writeToLogFile(logfullfile,'Filtering nodes if elements are not linear quadrilaterals ...\n')
writeToLogFile(logfullfile,['Calling function ', 'filterRectangularMesh',' ...\n']);
try
  nodes,elements,edges,...
  nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
  edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
  elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside=filterRectangularMesh(logfullfile,latexFolder,...
                                                                                                                                                                   elType,elOrder,Nx,Ny,NxEquiv,NyEquiv,...
                                                                                                                                                                   baseNodes,isCircular,circularity)
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier,'\n'])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end
writeToLogFile(logfullfile,['... done.','\n'])


elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: higherOrderGradedRectangle\n')



return
