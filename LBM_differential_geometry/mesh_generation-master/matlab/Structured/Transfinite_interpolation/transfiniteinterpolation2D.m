function[lattice]=transfiniteinterpolation2D(logfullfile,N,compdomain,xi_min,xi_max,eta_min,eta_max,Ndim1,Ndim2,interpolanttype,e1,e2,e3,e4,c1,c2,c3,c4)

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
%          Input: meshed computational domain compdomain
%                 numerical flag to choose the interpolant type
%                 e1 = e1(xi,eta_min) = e1(xi) = (x(xi),y(xi))
%                 e2 = e2(xi_min,eta) = e2(eta) = (x(eta),y(eta))
%                 e3 = e3(xi,eta_max) = e3(xi) = (x(xi),y(xi))
%                 e4 = e4(xi_max,eta) = e4(eta) = (x(eta),y(eta))
%                 c1 = c1(xi_min,eta_min) = (x1,y1)
%                 c2 = c2(xi_max,eta_min) = (x2,y2)
%                 c3 = c3(xi_max,eta_max) = (x3,y3)
%                 c4 = c4(xi_min,eta_max) = (x4,y4)
%                 e1 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e2 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e3 is a Ndim1 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 e4 is a Ndim2 x 6 vector (x,y,dxdxi,dxdeta,dydxi,dydeta)
%                 c1 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c2 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c3 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%                 c4 is a 1 x 8 vector (x,y,dxdxi,dxdeta,dydxi,dydeta,d2xdxideta,d2ydxideta)
%
%    Interpolant: 1 --> Lagrange
%                 2 --> Hermite
%
%         Output: mesh in the physical domain

%%

writeToLogFile(logfullfile,'In function: transfiniteinterpolation2D\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

mesh = zeros(N,2);

switch interpolanttype
    case 1 %--------------------------- Lagrange
        for r=1:Ndim2
            for p=1:Ndim1
                xi  = compdomain(p + (r-1)*Ndim1,1);
                eta = compdomain(p + (r-1)*Ndim1,2);
                mesh(p + (r-1)*Ndim1,1) = Lagrangeinterps1D([xi_min;xi_max],[e4(r,1);e2(r,1)],xi) + Lagrangeinterps1D([eta_min;eta_max],[e1(p,1);e3(p,1)],eta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,1) c4(1,1);c2(1,1) c3(1,1)],xi,eta);
                mesh(p + (r-1)*Ndim1,2) = Lagrangeinterps1D([xi_min;xi_max],[e4(r,2);e2(r,2)],xi) + Lagrangeinterps1D([eta_min;eta_max],[e1(p,2);e3(p,2)],eta) - Lagrangeinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,2) c4(1,2);c2(1,2) c3(1,2)],xi,eta); 
            end
        end
    case 2 %--------------------------- Hermite
        for r=1:Ndim2
            for p=1:Ndim1
                xi  = compdomain(p + (r-1)*Ndim1,1);
                eta = compdomain(p + (r-1)*Ndim1,2);
                mesh(p + (r-1)*Ndim1,1) = Hermiteinterps1D([xi_min;xi_max],[e4(r,1);e2(r,1)],[e4(r,3);e2(r,3)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,1);e3(p,1)],[e4(r,4);e2(r,4)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,1) c4(1,1);c2(1,1) c3(1,1)],[c1(1,3) c4(1,3);c2(1,3) c3(1,3)],[c1(1,4) c4(1,4);c2(1,4) c3(1,4)],[c1(1,7) c4(1,7);c2(1,7) c3(1,7)],xi,eta);
                mesh(p + (r-1)*Ndim1,2) = Hermiteinterps1D([xi_min;xi_max],[e4(r,2);e2(r,2)],[e4(r,5);e2(r,5)],xi) + Hermiteinterps1D([eta_min;eta_max],[e1(p,2);e3(p,2)],[e4(r,6);e2(r,6)],eta) - Hermiteinterps2D([xi_min;xi_max],[eta_min;eta_max],[c1(1,2) c4(1,2);c2(1,2) c3(1,2)],[c1(1,5) c4(1,5);c2(1,5) c3(1,5)],[c1(1,6) c4(1,6);c2(1,6) c3(1,6)],[c1(1,8) c4(1,8);c2(1,8) c3(1,8)],xi,eta); 
            end
        end
end

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,Ndim1,Ndim2);

mesh(indicesC1,1) = c1(1,1);
mesh(indicesC1,2) = c1(1,2);

mesh(indicesC2,1) = c2(1,1);
mesh(indicesC2,2) = c2(1,2);

mesh(indicesC3,1) = c3(1,1);
mesh(indicesC3,2) = c3(1,2);

mesh(indicesC4,1) = c4(1,1);
mesh(indicesC4,2) = c4(1,2);

mesh(indicesE1,1) = e1(2:end-1,1);
mesh(indicesE1,2) = e1(2:end-1,2);

mesh(indicesE2,1) = e2(2:end-1,1);
mesh(indicesE2,2) = e2(2:end-1,2);

mesh(indicesE3,1) = e3(2:end-1,1);
mesh(indicesE3,2) = e3(2:end-1,2);

mesh(indicesE4,1) = e4(2:end-1,1);
mesh(indicesE4,2) = e4(2:end-1,2);

lattice = [compdomain mesh mesh];

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: transfiniteinterpolation2D\n')

return