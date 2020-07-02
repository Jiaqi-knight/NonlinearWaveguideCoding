function[optNodes] = quad4optimize(nodes,elements,boundary,obj,tol,maxIt)
%%
%==============================================================================
% Copyright (c) 2016 Université de Lorraine & Luleå tekniska universitet
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
%  A function to optimize meshes with 4-nodes linear quadrilateral elements
%
% REFERENCES
% [1] Gang Mei, John C. Tipper and Nengxiong Xu, The Modified Direct Method: an Approach for Smoothing Planar and Surface Meshes. Procedia Computer Science 18(2013) 2436-2439
%     Link: http://arxiv.org/pdf/1212.3133.pdf
% [2] Gang Mei, John C. Tipper and Nengxiong Xu, A Hybrid Approach for Optimizing Planar Triangular Meshes. 2012 2nd International Conference on Computer Science and Network Technology (ICCSNT), pp.968,972, 29-31 Dec. 2012
%     Link: https://arxiv.org/ftp/arxiv/papers/1212/1212.6045.pdf
% [3] K. Ho-Le, Finite element mesh generation methods: a review and classification. Computer-Aided Design, 20(1):27–38, 1988
%     Link: http://www.ann.jussieu.fr/frey/papers/meshing/Le%20K.H.,%20Finite%20element%20mesh%20generation%20methods,%20a%20reviewand%20classification.pdf
%%
dof = 2;

N = size(nodes,1);

d = 0.25*ones(dof*N,1);

for i=1:length(boundary)
    d(dof*(boundary(i)-1)+1) = 1;
    d(dof*(boundary(i)-1)+2) = 1;
end

D = spdiags(d,0,dof*N,dof*N);

clear d

k = [0 0 0.5 0.5 0 0 0.5 -0.5;...
     0 0 -0.5 0.5 0 0 0.5 0.5;...
     0.5 -0.5 0 0 0.5 0.5 0 0;...
     0.5 0.5 0 0 -0.5 0.5 0 0;...
     0 0 0.5 -0.5 0 0 0.5 0.5;...
     0 0 0.5 0.5 0 0 -0.5 0.5;...
     0.5 0.5 0 0 0.5 -0.5 0 0;...
     -0.5 0.5 0 0 0.5 0.5 0 0];
 
K = sparse(dof*N,dof*N);
     
for i=1:size(elements,1)
    K = K + nodalmaps(N,elements(i,:),dof,1,1)*k;
end

for i=1:length(boundary)
    for j=1:N
        if j==boundary(i)
            K(dof*(boundary(i)-1)+1,dof*(j-1)+1) = 1;
            K(dof*(boundary(i)-1)+1,dof*(j-1)+2) = 0;
            K(dof*(boundary(i)-1)+2,dof*(j-1)+1) = 0;
            K(dof*(boundary(i)-1)+2,dof*(j-1)+2) = 1;
        else
            K(dof*(boundary(i)-1)+1,dof*(j-1)+1) = 0;
            K(dof*(boundary(i)-1)+1,dof*(j-1)+2) = 0;
            K(dof*(boundary(i)-1)+2,dof*(j-1)+1) = 0;
            K(dof*(boundary(i)-1)+2,dof*(j-1)+2) = 0;
        end
    end
end

fshape = quad4quality(nodes,elements);
it = 1;

solveNodes = zeros(dof*N,1);
for i=1:N
    solveNodes(dof*(i-1)+1) = nodes(i,1);
    solveNodes(dof*(i-1)+2) = nodes(i,2);
end

optNodes = nodes;

while abs(fshape-obj)>=tol && it<=maxIt
    solveNodes = D*(K*solveNodes);
    for i=1:N
        optNodes(i,1) = solveNodes(dof*(i-1)+1);
        optNodes(i,2) = solveNodes(dof*(i-1)+2);
    end
    fshape = quad4quality(optNodes,elements);
    it = it + 1;
end

return

