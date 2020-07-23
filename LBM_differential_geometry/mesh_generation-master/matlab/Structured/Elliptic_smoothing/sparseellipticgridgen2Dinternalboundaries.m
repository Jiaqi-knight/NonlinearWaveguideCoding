function[lattice]=sparseellipticgridgen2Dinternalboundaries(logfullfile,Nx,N,...
                                                            lattice,deltaq,boundary,flagperiodicity,periodicity,...
                                                            indicesbulk,...
                                                            indicesE1,indicesE2,indicesE3,indicesE4,...
                                                            indicesC1,indicesC2,indicesC3,indicesC4,...
                                                            firstdevneighbours,itmax,tol,spy)

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
%%

writeToLogFile(logfullfile,'In function: sparseellipticgridgen2Dinternalboundaries\n')
writeToLogFile(logfullfile,'\nStarting timer\n')
start = tic;

if ~flagperiodicity
    boundaryindices = [indicesE1;indicesE2;indicesE3;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4;boundary];
elseif  any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    boundaryindices = [indicesE2;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4;boundary];
elseif ~any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    boundaryindices = [indicesE1;indicesE3;indicesC1;indicesC2;indicesC3;indicesC4;boundary];    
elseif ~any(periodicity==1) && ~any(periodicity==2) &&  any(periodicity==3) && ~any(periodicity==4)
            
elseif ~any(periodicity==1) && ~any(periodicity==2) && ~any(periodicity==3) &&  any(periodicity==4)
    
elseif  any(periodicity==1) &&  any(periodicity==2) && ~any(periodicity==3) && ~any(periodicity==4)
    boundaryindices = [indicesC1;indicesC2;indicesC3;indicesC4;boundary];
elseif  any(periodicity==1) &&  any(periodicity==2) &&  any(periodicity==3) &&  any(periodicity==4)
    boundaryindices = [boundary];
end

effectivebulk = [];

for i=1:N
    if ~any(boundaryindices==i)
        effectivebulk = [effectivebulk;i];
    end
end

it = 0;
err = 1;

if spy
    covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);
    metriccoefficients = computemetriccoefficients2D(covariantbase); 
    clear covariantbase
    A = sparse([boundaryindices;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk],[boundaryindices;effectivebulk;effectivebulk-1;effectivebulk+1;effectivebulk-Nx;effectivebulk+Nx;effectivebulk-1-Nx;effectivebulk+1-Nx;effectivebulk-1+Nx;effectivebulk+1+Nx],[ones(size(boundaryindices,1),1);-2*(metriccoefficients(effectivebulk,2)./(deltaq(1).^2)+metriccoefficients(effectivebulk,1)./(deltaq(2).^2));metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2))],N,N);   
    figure();
    spy(A);
    hold on
    grid on
    title('Structure of solving matrix for 2D elliptic grid generation')
    clear A metriccoefficients
end

while it<=itmax && err>=tol
    covariantbase = computecovariantbase2D(N,deltaq,lattice,firstdevneighbours);
    [metriccoefficients,J] = computemetriccoefficients2D(covariantbase);
    %J = covariantbase(:,1).*covariantbase(:,4) - covariantbase(:,3).*covariantbase(:,2);
    Pvec = sparse(P(lattice));
    Qvec = sparse(Q(lattice));
%     bx = zeros(N,1);
%     by = zeros(N,1);
%     bx(boundaryindices,:) = lattice(boundaryindices,3);
%     by(boundaryindices,:) = lattice(boundaryindices,4);
%     bx(indicesbulk,:) = -(J(indicesbulk,:).^2).*(Pvec(indicesbulk,:).*covariantbase(indicesbulk,1)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,3));
%     by(indicesbulk,:) = -(J(indicesbulk,:).^2).*(Pvec(indicesbulk,:).*covariantbase(indicesbulk,2)+Qvec(indicesbulk,:).*covariantbase(indicesbulk,4));
    bx = sparse([boundaryindices;effectivebulk],[ones(length(boundaryindices),1);ones(length(effectivebulk),1)],[lattice(boundaryindices,3);-(J(effectivebulk,:).^2).*(Pvec(effectivebulk,:).*covariantbase(effectivebulk,1)+Qvec(effectivebulk,:).*covariantbase(effectivebulk,3))],N,1);
    by = sparse([boundaryindices;effectivebulk],[ones(length(boundaryindices),1);ones(length(effectivebulk),1)],[lattice(boundaryindices,4);-(J(effectivebulk,:).^2).*(Pvec(effectivebulk,:).*covariantbase(effectivebulk,2)+Qvec(effectivebulk,:).*covariantbase(effectivebulk,4))],N,1);
    clear covariantbase Pvec Qvec J
    A = sparse([boundaryindices;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk;effectivebulk],[boundaryindices;effectivebulk;effectivebulk-1;effectivebulk+1;effectivebulk-Nx;effectivebulk+Nx;effectivebulk-1-Nx;effectivebulk+1-Nx;effectivebulk-1+Nx;effectivebulk+1+Nx],[ones(size(boundaryindices,1),1);-2*(metriccoefficients(effectivebulk,2)./(deltaq(1).^2)+metriccoefficients(effectivebulk,1)./(deltaq(2).^2));metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,2)./(deltaq(1).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);metriccoefficients(effectivebulk,1)./(deltaq(2).^2);-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2));-0.5*metriccoefficients(effectivebulk,3)./(deltaq(1).*deltaq(2))],N,N); 
%     A = speye(N);
%     A(indicesbulk,indicesbulk) = diag(-2*(metriccoefficients(indicesbulk,2)./(deltaq(1).^2)+metriccoefficients(indicesbulk,1)./(deltaq(2).^2)));
%     A(indicesbulk,indicesbulk-1) = A(indicesbulk,indicesbulk-1) + diag(metriccoefficients(indicesbulk,2)./(deltaq(1).^2));
%     A(indicesbulk,indicesbulk+1) = A(indicesbulk,indicesbulk+1) + diag(metriccoefficients(indicesbulk,2)./(deltaq(1).^2));
%     A(indicesbulk,indicesbulk-Nx) = A(indicesbulk,indicesbulk-Nx) + diag(metriccoefficients(indicesbulk,1)./(deltaq(2).^2));
%     A(indicesbulk,indicesbulk+Nx) = A(indicesbulk,indicesbulk+Nx) + diag(metriccoefficients(indicesbulk,1)./(deltaq(2).^2));
%     A(indicesbulk,indicesbulk-1-Nx) = A(indicesbulk,indicesbulk-1-Nx) + diag(-0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));
%     A(indicesbulk,indicesbulk+1-Nx) = A(indicesbulk,indicesbulk+1-Nx) + diag(0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));
%     A(indicesbulk,indicesbulk-1+Nx) = A(indicesbulk,indicesbulk-1+Nx) + diag(0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));
%     A(indicesbulk,indicesbulk+1+Nx) = A(indicesbulk,indicesbulk+1+Nx) + diag(-0.5*metriccoefficients(indicesbulk,3)./(deltaq(1).*deltaq(2)));   
    lattice(:,5) = full(A\bx);
    lattice(:,6) = full(A\by);
    err = sqrt(sum([bx-A*lattice(:,5);by-A*lattice(:,6)].^2));
    it = it + 1;
    clear A bx by metriccoefficients
end

elapsed = toc(start);
writeToLogFile(logfullfile,'Timer stopped.\n')
writeToLogFile(logfullfile,['\nELAPSED WALLCLOCK TIME: ', num2str(elapsed),' [s]\n\n'])
writeToLogFile(logfullfile,'Exiting function: sparseellipticgridgen2Dinternalboundaries\n')

return