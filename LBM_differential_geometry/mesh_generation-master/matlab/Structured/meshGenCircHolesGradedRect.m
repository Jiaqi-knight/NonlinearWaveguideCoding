function[nodes,elements,edges,...
         nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,...
         nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
         edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
         elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,...
         elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside,...
         nodesInternalBoundaries,elementsInternalBoundaries]=meshGenCircHolesGradedRect(logfullfile,elType,elOrder,x0,y0,lx,ly,Nx,Ny,holes)
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
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         lx - [M x 1] vector - Side length in each one of the M mesh regions in x-direction
%         ly - [N x 1] vector - Side length in each one of the N mesh regions in y-direction
%         Nx - [M x 1] vector - Number of ELEMENTS in each one of the M mesh regions in x-direction
%         Ny - [N x 1] vector - Number of ELEMENTS in each one of the N mesh regions in y-direction
%         holes - [H x 3] matrix - H is the number of holes; for each hole the following data must be provided:
%                                  xC - scalar - x-coordinate of hole's center
%                                  yC - scalar - y-coordinate of hole's center
%                                  R - scalar - Radius
%%

writeToLogFile(logfullfile,'In function: meshGenCircHolesGradedRect\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

writeToLogFile(logfullfile,['Creating base rectangular mesh ...','\n'])
try
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order quadrilaterals','\n'])
      NxEquiv = Nx;
      NyEquiv = Ny;
      startX = 2;
      jumpX = 1;
      endX = sum(NxEquiv);
      startY = 2;
      jumpY = 1;
      endY = sum(NyEquiv);
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order quadrilaterals','\n'])
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
      startX = 3;
      jumpX = 2;
      endX = sum(NxEquiv)-1;
      startY = 3;
      jumpY = 2;
      endY = sum(NyEquiv)-1;
    end
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order triangles','\n'])
      NxEquiv = Nx;
      NyEquiv = Ny;
      startX = 2;
      jumpX = 1;
      endX = sum(NxEquiv);
      startY = 2;
      jumpY = 1;
      endY = sum(NyEquiv);
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order triangles','\n'])
      NxEquiv = 2*Nx;
      NyEquiv = 2*Ny;
      startX = 3;
      jumpX = 2;
      endX = sum(NxEquiv)-1;
      startY = 3;
      jumpY = 2;
      endY = sum(NyEquiv)-1;
    end
  end
  writeToLogFile(logfullfile,['    Calling function ', 'gradedRectangle',' ...\n']);
  baseMesh = gradedRectangle(logfullfile,x0,y0,lx,ly,NxEquiv,NyEquiv);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

writeToLogFile(logfullfile,['Create holes and deform internal boundaries ...','\n'])
try
  baseMesh = [baseMesh (1:length(baseMesh))' zeros(length(baseMesh),1) ones(length(baseMesh),1) 2*ones(length(baseMesh),1) zeros(length(baseMesh),1) zeros(length(baseMesh),1)]; % 1:x 2:y 3:<node index in full mesh> 4:<node index in reduced mesh> 5:isInBulk(1, else 0) 6:isInHole(0)/isOnBoundary(1)/isInBulk(2) 7:<internal boundary topology> 8:<hole number for interface nodes>
  defMesh = baseMesh;
  DeltaX = sum(lx)/sum(Nx);
  DeltaY = sum(ly)/sum(Ny);
  x1 = x0 - 0.5*sum(lx) + 10*DeltaX;
  x2 = x0 + 0.5*sum(lx) - 10*DeltaX;
  y1 = y0 - 0.5*sum(ly) + 10*DeltaY;
  y2 = y0 + 0.5*sum(ly) - 10*DeltaY;
  for k=1:length(holes)
    xC = holes(k,1);
    yC = holes(k,2);
    R = holes(k,3);
    if x1>xC+R || x2<xC+R || x1>xC-R || x2<xC-R
      writeToLogFile(logfullfile,['Hole number ',num2str(k),' is at least partially outside the domain. ','\n'])
      writeToLogFile(logfullfile,['Skipping it. ','\n'])
      continue
    end
    A.center = [xC yC];
    A.radius = R + 5*max([DeltaX DeltaY]);
    for h=1:k-1
      B.center = [holes(h,1) holes(h,2)];
      B.radius = holes(h,3) + 5*max([DeltaX DeltaY]);
      if checkCircleIntersections(A,B)>0
        writeToLogFile(logfullfile,['Hole number ',num2str(k),' intersects a preceding one. ','\n'])
        writeToLogFile(logfullfile,['Skipping it. ','\n'])
        continue
      end
    end
    L = R/sqrt(2);
    xH1 = x0 - L;
    xH2 = x0 + L;
    yH1 = y0 - L;
    yH2 = y0 + L;
    for j=startY:jumpY:endY
      for i=startX:jumpX:endX
        xP = baseMesh((j-1)*(sum(NxEquiv)+1)+i,1);    % current point
        yP = baseMesh((j-1)*(sum(NxEquiv)+1)+i,2);
        xRi = baseMesh((j-1)*(sum(NxEquiv)+1)+i+jumpX,1); % neighbour at the right
        yRi = baseMesh((j-1)*(sum(NxEquiv)+1)+i+jumpX,2);
        xLe = baseMesh((j-1)*(sum(NxEquiv)+1)+i-jumpX,1); % neighbour at the left
        yLe = baseMesh((j-1)*(sum(NxEquiv)+1)+i-jumpX,2);
        xUp = baseMesh((j-1+jumpY)*(sum(NxEquiv)+1)+i,1);       % upper neighbour
        yUp = baseMesh((j-1+jumpY)*(sum(NxEquiv)+1)+i,2);
        xLo = baseMesh((j-1-jumpY)*(sum(NxEquiv)+1)+i,1);   % lower neighbour
        yLo = baseMesh((j-1-jumpY)*(sum(NxEquiv)+1)+i,2);
        xLoLe = baseMesh((j-1-jumpY)*(sum(NxEquiv)+1)+i-jumpX,1); % lower left neighbour
        yLoLe = baseMesh((j-1-jumpY)*(sum(NxEquiv)+1)+i-jumpX,2);
        xLoRi = baseMesh((j-1-jumpY)*(sum(NxEquiv)+1)+i+jumpX,1); % lower right neighbour
        yLoRi = baseMesh((j-1-jumpY)*(sum(NxEquiv)+1)+i+jumpX,2);
        xUpRi = baseMesh((j-1+jumpY)*(sum(NxEquiv)+1)+i+jumpX,1); % upper right neighbour
        yUpRi = baseMesh((j-1+jumpY)*(sum(NxEquiv)+1)+i+jumpX,2);
        xUpLe = baseMesh((j-1+jumpY)*(sum(NxEquiv)+1)+i-jumpX,1); % upper left neighbour
        yUpLe = baseMesh((j-1+jumpY)*(sum(NxEquiv)+1)+i-jumpX,2);
        if (xH1==xP || xH2==xP) && (yH1==yP || yH2==yP)     % on one of the four boundary corners
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,8) = k;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,8) = k;
          if xH1==xP && yH1==yP                             % find which corner
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 1;
          elseif xH2==xP && yH1==yP
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 2;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 2;
          elseif xH2==xP && yH2==yP
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 3;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 3;
          else
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 4;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 4;
          end
        elseif xH1<xP && xH2>xP && (yH1==yP || yH2==yP)     % either on the south or north hole's boundary line
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
          defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,8) = k;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,8) = k;
          if yH1==yP                                        % on south line
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 7;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 7;
          else
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 5;       % on north line
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 5;
          %  if strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
          %    baseMesh((j-1-1)*(sum(NxEquiv)+1)+i,5) = 0;
         %     defMesh((j-1-1)*(sum(NxEquiv)+1)+i,5) = 0;
          %    baseMesh((j-1-1)*(sum(NxEquiv)+1)+i,6) = 0;
          %    defMesh((j-1-1)*(sum(NxEquiv)+1)+i,6) = 0;
          %  end
          end
        elseif (xH1==xP || xH2==xP) && yH1<yP && yH2>yP % either on the left or right hole's boundary line
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
          defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,8) = k;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,8) = k;
          if xH1==xP                                   % on left line
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 6;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 6;
          %  if strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
          %    baseMesh((j-1)*(sum(NxEquiv)+1)+i+1,5) = 0;
          %    defMesh((j-1)*(sum(NxEquiv)+1)+i+1,5) = 0;
          %    baseMesh((j-1)*(sum(NxEquiv)+1)+i+1,6) = 0;
          %    defMesh((j-1)*(sum(NxEquiv)+1)+i+1,6) = 0;
          %  end
          else
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 8;  % on right line
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 8;
          end
        elseif xH1<xP && xH2>xP && yH1<yP && yH2>yP     % inside hole
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 0;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 0;
          baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 0;
          defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 0;
          if strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
            baseMesh((j-1)*(sum(NxEquiv)+1)+i+1,5) = 0;  % east neighbour
            defMesh((j-1)*(sum(NxEquiv)+1)+i+1,5) = 0;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i+1,6) = 0;
            defMesh((j-1)*(sum(NxEquiv)+1)+i+1,6) = 0;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i-1,5) = 0;  % west neighbour
            defMesh((j-1)*(sum(NxEquiv)+1)+i-1,5) = 0;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i-1,6) = 0;
            defMesh((j-1)*(sum(NxEquiv)+1)+i-1,6) = 0;
            baseMesh((j-1-1)*(sum(NxEquiv)+1)+i,5) = 0;  % south neighbour
            defMesh((j-1-1)*(sum(NxEquiv)+1)+i,5) = 0;
            baseMesh((j-1-1)*(sum(NxEquiv)+1)+i,6) = 0;
            defMesh((j-1-1)*(sum(NxEquiv)+1)+i,6) = 0;
            baseMesh((j-1+1)*(sum(NxEquiv)+1)+i,5) = 0;  % north neighbour
            defMesh((j-1+1)*(sum(NxEquiv)+1)+i,5) = 0;
            baseMesh((j-1+1)*(sum(NxEquiv)+1)+i,6) = 0;
            defMesh((j-1+1)*(sum(NxEquiv)+1)+i,6) = 0;
            baseMesh((j-1+1)*(sum(NxEquiv)+1)+i+1,5) = 0;  % north-east neighbour
            defMesh((j-1+1)*(sum(NxEquiv)+1)+i+1,5) = 0;
            baseMesh((j-1+1)*(sum(NxEquiv)+1)+i+1,6) = 0;
            defMesh((j-1+1)*(sum(NxEquiv)+1)+i+1,6) = 0;
            baseMesh((j-1+1)*(sum(NxEquiv)+1)+i-1,5) = 0;  % north-west neighbour
            defMesh((j-1+1)*(sum(NxEquiv)+1)+i-1,5) = 0;
            baseMesh((j-1+1)*(sum(NxEquiv)+1)+i-1,6) = 0;
            defMesh((j-1+1)*(sum(NxEquiv)+1)+i-1,6) = 0;
            baseMesh((j-1-1)*(sum(NxEquiv)+1)+i+1,5) = 0;  % south-east neighbour
            defMesh((j-1-1)*(sum(NxEquiv)+1)+i+1,5) = 0;
            baseMesh((j-1-1)*(sum(NxEquiv)+1)+i+1,6) = 0;
            defMesh((j-1)*(sum(NxEquiv)+1)+i+1,6) = 0;
            baseMesh((j-1-1)*(sum(NxEquiv)+1)+i-1,5) = 0;  % south-west neighbour
            defMesh((j-1-1)*(sum(NxEquiv)+1)+i-1,5) = 0;
            baseMesh((j-1-1)*(sum(NxEquiv)+1)+i-1,6) = 0;
            defMesh((j-1-1)*(sum(NxEquiv)+1)+i-1,6) = 0;
          end
        else                                            % outside hole
          if xH1<xLoLe && xH2>xLoLe && yH1<yLoLe && yH2>yLoLe       % lower left neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 1;
          elseif xH1<xLoRi && xH2>xLoRi && yH1<yLoRi && yH2>yLoRi   % lower right neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 2;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 2;
          elseif xH1<xUpRi && xH2>xUpRi && yH1<yUpRi && yH2>yUpRi   % upper right neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 3;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 3;
          elseif xH1<xUpLe && xH2>xUpLe && yH1<yUpLe && yH2>yUpLe   % upper left neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 4;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 4;
          elseif xH1<xRi && xH2>xRi && yH1<yRi && yH2>yRi   % right neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 6;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 6;
          elseif xH1<xLe && xH2>xLe && yH1<yLe && yH2>yLe   % left neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 8;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 8;
          elseif xH1<xUp && xH2>xUp && yH1<yUp && yH2>yUp   % upper neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 7;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 7;
          elseif xH1<xLo && xH2>xLo && yH1<yLo && yH2>yLo   % lower neighbour inside hole
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 1;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,1) = xC + R*cos(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            defMesh((j-1)*(sum(NxEquiv)+1)+i,2) = yC + R*sin(atan2((baseMesh((j-1)*(sum(NxEquiv)+1)+i,2)-yC),(baseMesh((j-1)*(sum(NxEquiv)+1)+i,1)-xC)));
            baseMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 5;
            defMesh((j-1)*(sum(NxEquiv)+1)+i,7) = 5;
          %else                                              % bulk
          %  baseMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
          %  defMesh((j-1)*(sum(NxEquiv)+1)+i,5) = 1;
          %  baseMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 2;
          %  defMesh((j-1)*(sum(NxEquiv)+1)+i,6) = 2;
          end
        end
      end
    end
  end
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

writeToLogFile(logfullfile,['Deforming the mesh ...','\n'])
try
  writeToLogFile(logfullfile,['    Creating reduced mesh ...','\n'])
  reducedMesh = defMesh(defMesh(:,5)>0,:); % select only nodes in bulk or on boundary
  reducedMesh(:,4) = (1:length(reducedMesh))'; % create node indeces for the reduced mesh
  defMesh(reducedMesh(:,3),4) = reducedMesh(:,4); % update node indeces of reduced mesh in full mesh
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Extracting node indeces and type of internal boundaries ...','\n'])
  internalboundaries = reducedMesh(reducedMesh(:,6)==1,4);
  boundarytype = reducedMesh(reducedMesh(:,6)==1,7);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'getindices2D',' ...\n']);
  [indicesbulk,indicesinternalbulk,...
   indicesE1,indicesE2,indicesE3,indicesE4,...
   indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,...
   indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,...
   indicesC1,indicesC2,indicesC3,indicesC4,...
   indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,sum(NxEquiv)+1,sum(NyEquiv)+1);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'build_neighbourhoods2D',' ...\n']);
  [temp1,temp2,temp3,firstdevneighbours] = build_neighbourhoods2D((sum(NxEquiv)+1)*(sum(NyEquiv)+1),sum(NxEquiv)+1,...
                                                                  0,0,...
                                                                  1,internalboundaries,boundarytype,...
                                                                  indicesbulk,indicesinternalbulk,...
                                                                  indicesE1,indicesE2,indicesE3,indicesE4,...
                                                                  indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,...
                                                                  indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,...
                                                                  indicesC1,indicesC2,indicesC3,indicesC4,...
                                                                  indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'sparseellipticgridgen2Dinternalboundaries',' ...\n']);
  itmax = 10;
  tol = 10^-8;
  smoothedMesh = sparseellipticgridgen2Dinternalboundaries(logfullfile,sum(NxEquiv)+1,(sum(NxEquiv)+1)*(sum(NyEquiv)+1),...
                                            reducedMesh(:,1:2),[sum(lx)/sum(NxEquiv) sum(ly)/sum(NyEquiv)],internalboundaries,0,0,...
                                            indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,...
                                            indicesC1,indicesC2,indicesC3,indicesC4,...
                                            firstdevneighbours,itmax,tol,0)
  reducedMesh(:,1:2) = smoothedMesh;
  defMesh(reducedMesh(:,3),1:2) = smoothedMesh;
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

writeToLogFile(logfullfile,['Filtering the mesh ...','\n'])
try
  % UR = UnReduced
  [URnodes,URelements,URedges,...
    URnodesSWcorner,URnodesSEcorner,URnodesNEcorner,URnodesNWcorner,...
    URnodesSOUTHside,URnodesEASTside,URnodesNORTHside,URnodesWESTside,...
    URedgesSOUTHside,URedgesEASTside,URedgesNORTHside,URedgesWESTside,...
    URelementsSWcorner,URelementsSEcorner,URelementsNEcorner,URelementsNWcorner,...
    URelementsSOUTHside,URelementsEASTside,URelementsNORTHside,URelementsWESTside]=filterRectangularMesh(logfullfile,...
                                                                                                 elType,elOrder,Nx,Ny,NxEquiv,NyEquiv,...
                                                                                                 defMesh(:,1:2),0,'none')
  % URnodes => 1:x 2:y 3:<node index in full mesh> 4:<node index in reduced mesh> 5:isInBulk(1, else 0) 6:isInHole(0)/isOnBoundary(1)/isInBulk(2) 7:<internal boundary topology> 8:<hole number for interface nodes>
  URnodes = [URnodes(:,1:4) (1:length(URnodes))' zeros(length(URnodes),1) URnodes(:,5:end)];
  % URnodes => 1:x 2:y 3:<node index in full mesh> 4:<node index in reduced mesh> 5:<node index in filtered mesh> 6:<node index in final mesh> 7:isInBulk(1, else 0) 8:isInHole(0)/isOnBoundary(1)/isInBulk(2) 9:<internal boundary topology> 10:<hole number for interface nodes>
  nodes = [];
  for i=1:length(URnodes)
    if URnodes(i,7)>0  % if node is in bulk
      index = length(nodes) + 1;
      URnodes(i,6) = index;
      nodes = [nodes;URnodes(i,:)];
    end
  end
  nodesSWcorner = [];
  for i=1:length(URnodesSWcorner)
    if URnodes(URnodesSWcorner(i),7)>0  % if node is in bulk
      index = URnodes(URnodesSWcorner(i),6);
      nodesSWcorner = [nodesSWcorner;index];
    end
  end
  nodesSEcorner = [];
  for i=1:length(URnodesSEcorner)
    if URnodes(URnodesSEcorner(i),7)>0  % if node is in bulk
      index = URnodes(URnodesSEcorner(i),6);
      nodesSEcorner = [nodesSEcorner;index];
    end
  end
  nodesNEcorner = [];
  for i=1:length(URnodesNEcorner)
    if URnodes(URnodesNEcorner(i),7)>0  % if node is in bulk
      index = URnodes(URnodesNEcorner(i),6);
      nodesNEcorner = [nodesNEcorner;index];
    end
  end
  nodesNWcorner = [];
  for i=1:length(URnodesNWcorner)
    if URnodes(URnodesNWcorner(i),7)>0  % if node is in bulk
      index = URnodes(URnodesNWcorner(i),6);
      nodesNWcorner = [nodesNWcorner;index];
    end
  end
  nodesSOUTHside = [];
  for i=1:length(URnodesSOUTHside)
    if URnodes(URnodesSOUTHside(i),7)>0  % if node is in bulk
      index = URnodes(URnodesSOUTHside(i),6);
      nodesSOUTHside = [nodesSOUTHside;index];
    end
  end
  nodesEASTside = [];
  for i=1:length(URnodesEASTside)
    if URnodes(URnodesEASTside(i),7)>0  % if node is in bulk
      index = URnodes(URnodesEASTside(i),6);
      nodesEASTside = [nodesEASTside;index];
    end
  end
  nodesNORTHside = [];
  for i=1:length(URnodesNORTHside)
    if URnodes(URnodesNORTHside(i),7)>0  % if node is in bulk
      index = URnodes(URnodesNORTHside(i),6);
      nodesNORTHside = [nodesNORTHside;index];
    end
  end
  nodesWESTside = [];
  for i=1:length(URnodesWESTside)
    if URnodes(URnodesWESTside(i),7)>0  % if node is in bulk
      index = URnodes(URnodesWESTside(i),6);
      nodesWESTside = [nodesWESTside;index];
    end
  end
  numElNodes = size(URelements,2);
  URelements = [URelements (1:length(URelements))' zeros(length(URelements),1) 2*ones(length(URelements),1) zeros(length(URelements),1)];
  % -4:<element index in filtered mesh> -3:<element index in final mesh> -2:<isOnInternalBoundary (2:inBulk, 1:onBoundary, else 0)> -1:<on which hole>
  elements = [];
  for i=1:length(URelements)
    isElInBulk = 1;
    for j=1:numElNodes
      if URnodes(URelements(i,j),7)==0  % if node is not in bulk
        isElInBulk = 0;
        URelements(i,-2) = 0;
        break
      end
    end
    if isElInBulk
      index = length(elements) + 1;
      URelements(i,-3) = index;
      for j=1:numElNodes
        if URnodes(URelements(i,j),8)==1  % if node is on internal boundary
          URelements(i,-1) = URnodes(URelements(i,j),10);
          URelements(i,-2) = 1;
          break
        end
      end
      for j=1:numElNodes
        URelements(i,j) = URnodes(URelements(i,j),6);
      end
      elements = [elements; URelements(i,:)];
    end
  end
  elementsSWcorner = [];
  for i=1:length(URelementsSWcorner)
    if URelements(URelementsSWcorner(i),-2)>0  % if element is in bulk
      index = URelements(URelementsSWcorner(i),-3);
      elementsSWcorner = [elementsSWcorner;index];
    end
  end
  elementsSEcorner = [];
  for i=1:length(URelementsSWcorner)
    if URelements(URelementsSWcorner(i),-2)>0  % if element is in bulk
      index = URelements(URelementsSWcorner(i),-3);
      elementsSWcorner = [elementsSWcorner;index];
    end
  end
  elementsNEcorner = [];
  for i=1:length(URelementsNEcorner)
    if URelements(URelementsNEcorner(i),-2)>0  % if element is in bulk
      index = URelements(URelementsNEcorner(i),-3);
      elementsNEcorner = [elementsNEcorner;index];
    end
  end
  elementsNWcorner = [];
  for i=1:length(URelementsNWcorner)
    if URelements(URelementsNWcorner(i),-2)>0  % if element is in bulk
      index = URelements(URelementsNWcorner(i),-3);
      elementsNWcorner = [elementsNWcorner;index];
    end
  end
  elementsSOUTHside = [];
  for i=1:length(URelementsSOUTHside)
    if URelements(URelementsSOUTHside(i),-2)>0  % if element is in bulk
      index = URelements(URelementsSOUTHside(i),-3);
      elementsSOUTHside = [elementsSOUTHside;index];
    end
  end
  elementsEASTside = [];
  for i=1:length(URelementsEASTside)
    if URelements(URelementsEASTside(i),-2)>0  % if element is in bulk
      index = URelements(URelementsEASTside(i),-3);
      elementsEASTside = [elementsEASTside;index];
    end
  end
  elementsNORTHside = [];
  for i=1:length(URelementsNORTHside)
    if URelements(URelementsNORTHside(i),-2)>0  % if element is in bulk
      index = URelements(URelementsNORTHside(i),-3);
      elementsNORTHside = [elementsNORTHside;index];
    end
  end
  elementsWESTside = [];
  for i=1:length(URelementsWESTside)
    if URelements(URelementsWESTside(i),-2)>0  % if element is in bulk
      index = URelements(URelementsWESTside(i),-3);
      elementsWESTside = [elementsWESTside;index];
    end
  end
  numEdgeNodes = size(URedges,2);
  URedges = [URedges (1:length(URelements))' zeros(length(URelements),1) ones(length(URelements),1)];
  % -3:<edge index in filtered mesh> -2:<edge index in final mesh> -1:<inInBulk (1, else 0)>
  edges = [];
  for i=1:length(URedges)
    isEdgeInBulk = 1;
    for j=1:numEdgeNodes
      if URnodes(URedges(i,j),7)==0  % if node is not in bulk
        isEdgeInBulk = 0;
        break
      end
    end
    if isEdgeInBulk
      index = length(edges) + 1;
      URedges(i,-2) = index;
      URedges(i,-1) = 0;
      for j=1:numEdgeNodes
        URedges(i,j) = URnodes(URedges(i,j),6);
      end
      edges = [edges; URedges(i,:)];
    end
  end
  edgesSOUTHside = [];
  for i=1:length(URedgesSOUTHside)
    if URedges(URedgesSOUTHside(i),-1)>0  % if edge is in bulk
      index = URedges(URedgesSOUTHside(i),-2);
      edgesSOUTHside = [edgesSOUTHside;index];
    end
  end
  edgesEASTside = [];
  for i=1:length(URedgesEASTside)
    if URedges(URedgesEASTside(i),-1)>0  % if edge is in bulk
      index = URedges(URedgesEASTside(i),-2);
      edgesEASTside = [edgesEASTside;index];
    end
  end
  edgesNORTHside = [];
  for i=1:length(URedgesNORTHside)
    if URedges(URedgesNORTHside(i),-1)>0  % if edge is in bulk
      index = URedges(URedgesNORTHside(i),-2);
      edgesNORTHside = [edgesNORTHside;index];
    end
  end
  edgesWESTside = [];
  for i=1:length(URedgesWESTside)
    if URedges(URedgesWESTside(i),-1)>0  % if edge is in bulk
      index = URedges(URedgesWESTside(i),-2);
      edgesWESTside = [edgesWESTside;index];
    end
  end
  nodesInternalBoundaries = {};
  elementsInternalBoundaries = {};
  for k=1:length(holes)
    nodesInternalBoundaries{k} = nodes(nodes(:,10)==k,6);
    elementsInternalBoundaries{k} = elements(elements(:,-1)==k,-3);
  end
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: meshGenCircHolesGradedRect\n')

return
