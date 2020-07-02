function[nodes,elements,edges,...
    nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
    edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
    elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside]=meshGenGradedCircSec(logfullfile,elType,elOrder,x0,y0,R,thetas,deltas)
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
%  A function to mesh a simply connected 2D rectangular geometry with elements of
%  shape and order of choice
%
%  Available:
%  - 1st and 2nd order quads
%  - 1st and 2nd order tris
%
%  Input: x0 - scalar - x-coordinate of center
%         y0 - scalar - y-coordinate of center
%         R  - scalar - Radius
%         thetas - [M x 1] vector - Angular aperture of each one of the M mesh regions
%                                   in each 90° section with wich the circle is mapped to a square,
%                                   the last one included; i.e. if the 90° region is not divided
%                                   the length of this vector is 1
%         deltas - [M x 1] vector - Angular aperture of ELEMENTS in each one of the M mesh regions in x-direction
%%

writeToLogFile(logfullfile,'In function: meshGenGradedCircSec\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

try
  writeToLogFile(logfullfile,['Compute boundaries of mesh in real and computational space ...','\n'])
  A = [-R/sqrt(2);-R/sqrt(2)]; % cos(45°)=sin(45°)=1/sqrt(2)
  B = [R/sqrt(2);-R/sqrt(2)];
  C = [R/sqrt(2);R/sqrt(2)];
  D = [-R/sqrt(2);R/sqrt(2)];
  writeToLogFile(logfullfile,['    A = (',num2str(A(1)),',',num2str(A(2)),')','\n'])
  writeToLogFile(logfullfile,['    B = (',num2str(B(1)),',',num2str(B(2)),')','\n'])
  writeToLogFile(logfullfile,['    C = (',num2str(C(1)),',',num2str(C(2)),')','\n'])
  writeToLogFile(logfullfile,['    D = (',num2str(D(1)),',',num2str(D(2)),')','\n'])
  A1 = A; % point on rectangular mesh starting configuration
  B1 = B;
  C1 = C;
  D1 = D;
  writeToLogFile(logfullfile,['    A1 = (',num2str(A1(1)),',',num2str(A1(2)),')','\n'])
  writeToLogFile(logfullfile,['    B1 = (',num2str(B1(1)),',',num2str(B1(2)),')','\n'])
  writeToLogFile(logfullfile,['    C1 = (',num2str(C1(1)),',',num2str(C1(2)),')','\n'])
  writeToLogFile(logfullfile,['    D1 = (',num2str(D1(1)),',',num2str(D1(2)),')','\n'])
  As = [A];
  Bs = [B];
  Cs = [D];
  Ds = [A];
  A1s = [A1];
  B1s = [B1];
  C1s = [D1];
  D1s = [A1];
  AB = [];
  DC = [];
  BC = [];
  AD = [];
  AB1 = [];
  DC1 = [];
  BC1 = [];
  AD1 = [];
  writeToLogFile(logfullfile,['    Number of mesh regions in 90° aperture = ',num2str(length(thetas)),'\n'])
  if length(thetas)>1
    writeToLogFile(logfullfile,['    Entering if in branch for number of regions greater than 1','\n'])
    for i=1:length(theta)
        Ntheta = thetas(i)/deltas(i);
        writeToLogFile(logfullfile,['    ****','\n'])
        writeToLogFile(logfullfile,['    Region number ',num2str(i),'\n'])
        writeToLogFile(logfullfile,['    theta = ',num2str(thetas(i)),'[rad]','\n'])
        writeToLogFile(logfullfile,['            ',num2str(thetas(i)*180/pi),'[°]','\n'])
        writeToLogFile(logfullfile,['    delta = ',num2str(deltas(i)),'[rad]','\n'])
        writeToLogFile(logfullfile,['            ',num2str(deltas(i)*180/pi),'[°]','\n'])
        writeToLogFile(logfullfile,['    Number of elements = ',num2str(Ntheta),'\n'])
        if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
          if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
            writeToLogFile(logfullfile,['    Type ad order of elements : First order quadrilaterals','\n'])
            NthetaEquiv = Ntheta;
          elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
            writeToLogFile(logfullfile,['    Type ad order of elements : Second order quadrilaterals','\n'])
            NthetaEquiv = 2*Ntheta;
          end
        elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
          if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
            writeToLogFile(logfullfile,['    Type ad order of elements : First order triangles','\n'])
            NthetaEquiv = Ntheta;
          elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
            writeToLogFile(logfullfile,['    Type ad order of elements : Second order triangles','\n'])
            NthetaEquiv = 2*Ntheta;
          end
        end
        writeToLogFile(logfullfile,['    Equivalent number of elements = ',num2str(NthetaEquiv),'\n'])
        writeToLogFile(logfullfile,['    **','\n'])
        for j=1:i
          writeToLogFile(logfullfile,['    thetas(',num2str(j),') = ',num2str(thetas(j)),' [rad]','\n'])
        end
        writeToLogFile(logfullfile,['    -----------------------------------------------------','\n'])
        writeToLogFile(logfullfile,['    sum(thetas(1:',num2str(i),')) = ',num2str(sum(thetas(1:i,1))),' [rad]','\n'])
        writeToLogFile(logfullfile,['    **','\n'])
        for j=1:i
          writeToLogFile(logfullfile,['    thetas(',num2str(j),') = ',num2str(thetas(j)*180/pi),' [°]','\n'])
        end
        writeToLogFile(logfullfile,['    -----------------------------------------------------','\n'])
        writeToLogFile(logfullfile,['    sum(thetas(1:',num2str(i),')) = ',num2str(sum(thetas(1:i,1).*180./pi)),' [°]','\n'])
        writeToLogFile(logfullfile,['    **','\n'])
        alpha = -0.75*pi + sum(theta(1:i,1));
        beta = -0.25*pi + sum(theta(1:i,1));
        if i>1
          alpha0 = -0.75*pi + sum(theta(1:i-1,1));
          beta0 = -0.25*pi + sum(theta(1:i-1,1));
        else
          alpha0 = -0.75*pi;
          beta0 = -0.25*pi;
        end
        writeToLogFile(logfullfile,['    alpha = ',num2str(alpha),' [rad]','\n'])
        writeToLogFile(logfullfile,['    alpha = ',num2str(alpha*180/pi),' [°]','\n'])
        writeToLogFile(logfullfile,['    beta = ',num2str(beta),' [rad]','\n'])
        writeToLogFile(logfullfile,['    beta = ',num2str(beta*180/pi),' [°]','\n'])
        writeToLogFile(logfullfile,['    alpha0 = ',num2str(alpha0),' [rad]','\n'])
        writeToLogFile(logfullfile,['    alpha0 = ',num2str(alpha0*180/pi),' [°]','\n'])
        writeToLogFile(logfullfile,['    beta0 = ',num2str(beta0),' [rad]','\n'])
        writeToLogFile(logfullfile,['    beta0 = ',num2str(beta0*180/pi),' [°]','\n'])
        E = [R*cos(alpha);R*sin(alpha)];
        F = [R*cos(beta);R*sin(beta)];
        G = [R*cos(alpha);-R*sin(alpha)];
        H = [-R*cos(beta);R*sin(beta)];
        writeToLogFile(logfullfile,['    E = (',num2str(E(1)),',',num2str(E(2)),')','\n'])
        writeToLogFile(logfullfile,['    F = (',num2str(F(1)),',',num2str(F(2)),')','\n'])
        writeToLogFile(logfullfile,['    G = (',num2str(G(1)),',',num2str(G(2)),')','\n'])
        writeToLogFile(logfullfile,['    H = (',num2str(H(1)),',',num2str(H(2)),')','\n'])
        E1 = [sign(E(1))*abs(A(2))/tan(0.25*pi+sum(theta(1:i,1))),A(2)];
        F1 = [B(1),sign(F(2))*B(1)*abs(tan(-0.25*pi+sum(theta(1:i,1))))];
        G1 = [E1(1),-E1(2)];
        H1 = [-F1(1),F1(2)];
        writeToLogFile(logfullfile,['    E1 = (',num2str(E1(1)),',',num2str(E1(2)),')','\n'])
        writeToLogFile(logfullfile,['    F1 = (',num2str(F1(1)),',',num2str(F1(2)),')','\n'])
        writeToLogFile(logfullfile,['    G1 = (',num2str(G1(1)),',',num2str(G1(2)),')','\n'])
        writeToLogFile(logfullfile,['    H1 = (',num2str(H1(1)),',',num2str(H1(2)),')','\n'])
        As = [As;E];
        Bs = [Bs;F];
        Cs = [Cs;G];
        Ds = [Ds;H];
        A1s = [A1s;E1];
        B1s = [B1s;F1];
        C1s = [C1s;G1];
        D1s = [D1s;H1];
        writeToLogFile(logfullfile,['    Adding points to arc AB','\n'])
        angles = (alpha0:(alpha-alpha0)/NthetaEquiv:alpha)';
        angles = angles(2:end-1);
        ABnew = [R.*cos(angles) R.*sin(angles)];
        AB = [AB;ABnew];
        if i<length(thetas)
          AB = [AB;E];
        end
        writeToLogFile(logfullfile,['    Adding points to cord AB1','\n'])
        AB1new = [sign(ABnew(:,1)).*abs(A(2))./tan(pi-abs(angles)),A(2)*ones(length(ABnew),1)];
        AB1 = [AB1;AB1new];
        if i<length(thetas)
          AB1 = [AB1;E1];
        end
        writeToLogFile(logfullfile,['    Adding points to arc DC (by symmetry of AB with respect to x axis)','\n'])
        DCnew = [ABnew(:,1) -ABnew(:,2)];
        DC = [DC;DCnew];
        if i<length(thetas)
          DC = [DC;G];
        end
        writeToLogFile(logfullfile,['    Adding points to cord DC1 (by symmetry of AB1 with respect to x axis)','\n'])
        DC1new = [AB1new(:,1) -AB1new(:,2)];
        DC1 = [DC1;DC1new];
        if i<length(thetas)
          DC1 = [DC1;G1];
        end
        writeToLogFile(logfullfile,['    Adding points to arc BC','\n'])
        angles = (beta0:(beta-beta0)/NthetaEquiv:beta)';
        angles = angles(2:end-1);
        BCnew = [R.*cos(angles) R.*sin(angles)];
        BC = [BC;BCnew];
        if i<length(thetas)
          BC = [BC;F];
        end
        writeToLogFile(logfullfile,['    Adding points to cord BC1','\n'])
        BC1new = [B(1),sign(BCnew(:,2)).*BCnew(:,1).*abs(tan(angles))];
        BC1 = [BC1;BC1new];
        if i<length(thetas)
          BC1 = [BC1;F1];
        end
        writeToLogFile(logfullfile,['    Adding points to arc AD (by symmetry of BC with respect to y axis)','\n'])
        ADnew = [-BCnew(:,1) BCnew(:,2)];
        AD = [AD;ADnew];
        if i<length(thetas)
          AD = [AD;H];
        end
        writeToLogFile(logfullfile,['    Adding points to cord AD1 (by symmetry of BC1 with respect to y axis)','\n'])
        AD1new = [-BC1new(:,1) BC1new(:,2)];
        AD1 = [AD1;AD1new];
        if i<length(thetas)
          AD1 = [AD1;H1];
        end
    end
  else
    writeToLogFile(logfullfile,['    Entering if in branch for number of regions equal to 1','\n'])
    Ntheta = 0.5*pi/deltas(1);
    writeToLogFile(logfullfile,['    ****','\n'])
    writeToLogFile(logfullfile,['    Region number ',num2str(1),'\n'])
    writeToLogFile(logfullfile,['    theta = ',num2str(0.5*pi),' [rad]','\n'])
    writeToLogFile(logfullfile,['            ',num2str(90),' [°]','\n'])
    writeToLogFile(logfullfile,['    delta = ',num2str(deltas(1)),' [rad]','\n'])
    writeToLogFile(logfullfile,['            ',num2str(deltas(1)*180/pi),' [°]','\n'])
    writeToLogFile(logfullfile,['    Number of elements = ',num2str(Ntheta),'\n'])
    if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
      if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
        writeToLogFile(logfullfile,['    Type ad order of elements : First order quadrilaterals','\n'])
        NthetaEquiv = Ntheta;
      elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
        writeToLogFile(logfullfile,['    Type ad order of elements : Second order quadrilaterals','\n'])
        NthetaEquiv = 2*Ntheta;
      end
    elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
      if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
        writeToLogFile(logfullfile,['    Type ad order of elements : First order triangles','\n'])
        NthetaEquiv = Ntheta;
      elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
        writeToLogFile(logfullfile,['    Type ad order of elements : Second order triangles','\n'])
        NthetaEquiv = 2*Ntheta;
      end
    end
    writeToLogFile(logfullfile,['    Equivalent number of elements = ',num2str(NthetaEquiv),'\n'])
    writeToLogFile(logfullfile,['    **','\n'])
    for j=1:i
      writeToLogFile(logfullfile,['    thetas(',num2str(j),') = ',num2str(thetas(j)),' [rad]','\n'])
    end
    writeToLogFile(logfullfile,['    -----------------------------------------------------','\n'])
    writeToLogFile(logfullfile,['    sum(thetas(1:',num2str(i),')) = ',num2str(sum(thetas(1:i,1))),' [rad]','\n'])
    writeToLogFile(logfullfile,['    **','\n'])
    for j=1:i
      writeToLogFile(logfullfile,['    thetas(',num2str(j),') = ',num2str(thetas(j)*180/pi),' [°]','\n'])
    end
    writeToLogFile(logfullfile,['    -----------------------------------------------------','\n'])
    writeToLogFile(logfullfile,['    sum(thetas(1:',num2str(i),')) = ',num2str(sum(thetas(1:i,1).*180./pi)),' [°]','\n'])
    writeToLogFile(logfullfile,['    **','\n'])
    writeToLogFile(logfullfile,['    Adding points to arc AB','\n'])
    alpha = -0.25*pi;
    alpha0 = -0.75*pi;
    angles = (alpha0:(alpha-alpha0)/NthetaEquiv:alpha)';
    angles = angles(2:end-1);
    ABnew = [R.*cos(angles) R.*sin(angles)];
    AB = [AB;ABnew];
    writeToLogFile(logfullfile,['    Adding points to cord AB1','\n'])
    AB1new = [sign(ABnew(:,1)).*abs(A(2))./tan(pi-abs(angles)),A(2)*ones(length(ABnew),1)];
    AB1 = [AB1;AB1new];
    writeToLogFile(logfullfile,['    Adding points to arc DC (by symmetry of AB with respect to x axis)','\n'])
    DCnew = [ABnew(:,1) -ABnew(:,2)];
    DC = [DC;DCnew];
    writeToLogFile(logfullfile,['    Adding points to cord DC1 (by symmetry of AB1 with respect to x axis)','\n'])
    DC1new = [AB1new(:,1) -AB1new(:,2)];
    DC1 = [DC1;DC1new];
    writeToLogFile(logfullfile,['    Adding points to arc BC','\n'])
    beta = 0.25*pi;
    beta0 = -0.25*pi;
    angles = (beta0:(beta-beta0)/NthetaEquiv:beta)';
    angles = angles(2:end-1);
    BCnew = [R.*cos(angles) R.*sin(angles)];
    BC = [BC;BCnew];
    writeToLogFile(logfullfile,['    Adding points to cord BC1','\n'])
    BC1new = [B(1),sign(BCnew(:,2)).*BCnew(:,1).*abs(tan(angles))];
    BC1 = [BC1;BC1new];
    writeToLogFile(logfullfile,['    Adding points to arc AD (by symmetry of BC with respect to y axis)','\n'])
    ADnew = [-BCnew(:,1) BCnew(:,2)];
    AD = [AD;ADnew];
    writeToLogFile(logfullfile,['    Adding points to cord AD1 (by symmetry of BC1 with respect to y axis)','\n'])
    AD1new = [-BC1new(:,1) BC1new(:,2)];
    AD1 = [AD1;AD1new];
  end
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end



% create equivalent mesh with linear quadrilateral elements
writeToLogFile(logfullfile,['Creating equivalent base mesh with linear quadrilateral elements for mesh generation ...','\n'])
try
  Hline = [A1;AB1;B1];
  Vline = [B1;BC1;C1];
  Xs = Hline(:,1);
  Ys = Vline(:,2);
  baseMesh = zeros(length(Hline)*length(Vline),2);
  for j=1:length(Vline)
    baseMesh((j-1)*length(Hline)+1:j*length(Hline),1) = Xs;
    baseMesh((j-1)*length(Hline)+1:j*length(Hline),2) = Ys(j)*ones(length(Hline),1);
  end
  writeToLogFile(logfullfile,['... done.','\n'])
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

% create structured mesh by mapping
writeToLogFile(logfullfile,['Creating structured mesh by mapping ...','\n'])
try
  writeToLogFile(logfullfile,['    Calling function ', 'transfiniteinterpolation2D',' ...\n']);
  mappedMesh = transfiniteinterpolation2D(logfullfile,length(Hline)*length(Vline),baseMesh,A(1),B(1),A(2),D(2),length(Hline),length(Vline),1,[A1;AB1;B1],[A1;AD1;D1],[D1;DC1;C1],[B1;BC1;C1],A,B,C,D);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'getindices2D',' ...\n']);
  [indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,length(Hline),length(Vline));
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'build_neighbourhoods2D',' ...\n']);
  [temp1,temp2,temp3,firstdevneighbours] = build_neighbourhoods2D(logfullfile,length(Hline)*length(Vline),length(Hline),0,0,0,0,0,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);
  writeToLogFile(logfullfile,['    ... done.','\n'])
  writeToLogFile(logfullfile,['    Calling function ', 'sparseellipticgridgen2D',' ...\n']);
  itmax = 10;
  tol = 10^-8;
  mappedMesh = sparseellipticgridgen2D(logfullfile,length(Hline),length(Hline)*length(Vline),mappedMesh,[(B(1)-A(1))/(length(Hline)-1) (D(2)-A(2))/(length(Vline)-1)],0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0);
  clear temp1 temp2 temp3
  writeToLogFile(logfullfile,['    ... done.','\n'])
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
  Nx = thetas./deltas;
  if strcomp(elType,'quads') || strcomp(elType,'quad') || strcomp(elType,'quadrilaterals') || strcomp(elType,'quadrilateral')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order quadrilaterals','\n'])
      NxEquiv = Nx;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order quadrilaterals','\n'])
      NxEquiv = 2*Nx;
    end
  elseif strcomp(elType,'tris') || strcomp(elType,'tri') || strcomp(elType,'triangles') || strcomp(elType,'triangle')
    if strcomp(elOrder,'first') || strcomp(elOrder,'First') || strcomp(elOrder,'1st') || strcomp(elOrder,'1')
      writeToLogFile(logfullfile,['    Type ad order of elements : First order triangles','\n'])
      NxEquiv = Nx;
    elseif strcomp(elOrder,'second') || strcomp(elOrder,'Second') || strcomp(elOrder,'2nd') || strcomp(elOrder,'2')
      writeToLogFile(logfullfile,['    Type ad order of elements : Second order triangles','\n'])
      NxEquiv = 2*Nx;
    end
  end
  [nodes,elements,edges,...
      nodesSWcorner,nodesSEcorner,nodesNEcorner,nodesNWcorner,nodesSOUTHside,nodesEASTside,nodesNORTHside,nodesWESTside,...
      edgesSOUTHside,edgesEASTside,edgesNORTHside,edgesWESTside,...
      elementsSWcorner,elementsSEcorner,elementsNEcorner,elementsNWcorner,elementsSOUTHside,elementsEASTside,elementsNORTHside,elementsWESTside] = filterRectangularMesh(logfullfile,...
                                                                                                                                                                      elType,elOrder,Nx,Nx,NxEquiv,NxEquiv,...
  nodes(:,1) = x0 + nodes(:,1);
  nodes(:,2) = y0 + nodes(:,2);
  writeToLogFile(logfullfile,['... done.','\n'])                                                                                                                      mappedMesh,0,'none')                                                                                                                                                     baseNodes,isCircular,circularity)
catch ME
  writeToLogFile(logfullfile,['An error occurred: ', ME.identifier,'\n'])
  writeToLogFile(logfullfile,['Terminating program.','\n'])
  exit(2)
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: meshGenGradedCircSec\n')

return
