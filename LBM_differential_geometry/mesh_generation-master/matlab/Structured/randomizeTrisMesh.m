function[baseNodes,baseEdges]=randomizeTrisMesh(logfullfile,baseNodes,iterations,baseEdges,fixedNodes)
%%
%==============================================================================
% Copyright (c) 2016 Universit� de Lorraine & Lule� tekniska universitet
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
%  A function to mesh a non-simply connected, i.e. with rectangular holes,
%  rectangular geometry with rectangular elements
%
%  Input: baseNodes - N x 3 matrix - for each node, an integer label, x and y coordinates are given
%
%%

writeToLogFile(logfullfile,'In function: randomizeTrisMesh\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

if length(baseEdges(0,:))==2
  for n=1:iterations
    for i=1:length(baseNodes)
      if isempty(find(fixedNodes==baseNodes(i,1)))
        [edgeIndeces,vertexIndexes]=find(baseEdges==baseNodes(i,1));
        vertexIndexes = mod(vertexIndexes,2) + 1;
        R = -1;
        for j=1:length(edgeIndeces)
          r = sqrt((baseNodes(baseNodes(i,1),2)-baseNodes(baseEdges(edgeIndeces(j),vertexIndexes(j)),2))^2+(baseNodes(baseNodes(i,1),3)-baseNodes(baseEdges(edgeIndeces(j),vertexIndexes(j)),3))^2);
          if R<0
            R = r;
          elseif r<R
            R = r;
          end
        end
        rFac = rand();
        thetaFac = rand();
        baseNodes(i,2) = baseNodes(i,2) + rFac*R*cos(thetaFac*2*pi);
        baseNodes(i,3) = baseNodes(i,3) + rFac*R*sin(thetaFac*2*pi);
      end
    end
  end
else
  for n=1:iterations
    for i=1:length(baseNodes)
      if isempty(find(fixedNodes==baseNodes(i,1)))
        [edgeIndeces,vertexIndexes]=find(baseEdges==baseNodes(i,1));
        vertexIndexes = mod(vertexIndexes,2) + 1;
        R = -1;
        for j=1:length(edgeIndeces)
          for k=1:length(baseEdges(edgeIndeces(j),:))
            if k~=vertexIndexes(j)
              r = sqrt((baseNodes(baseNodes(i,1),2)-baseNodes(baseEdges(edgeIndeces(j),k),2))^2+(baseNodes(baseNodes(i,1),3)-baseNodes(baseEdges(edgeIndeces(j),k),3))^2);
              if R<0
                R = r;
              elseif r<R
                R = r;
              end
            end
          end
        end
        rFac = rand();
        thetaFac = rand();
        baseNodes(i,2) = baseNodes(i,2) + rFac*R*cos(thetaFac*2*pi);
        baseNodes(i,3) = baseNodes(i,3) + rFac*R*sin(thetaFac*2*pi);
      end
    end
  end
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: randomizeTrisMesh\n')

return
