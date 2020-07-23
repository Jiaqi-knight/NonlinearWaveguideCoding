function[nodes,edges,elements] = tri3optimize(nodes,edges,elements,boundary,obj,tol,maxIt)
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

fshape = tri3quality(nodes,elements);
it = 1;

solveNodes = zeros(dof*N,1);
for i=1:N
    solveNodes(dof*(i-1)+1) = nodes(i,1);
    solveNodes(dof*(i-1)+2) = nodes(i,2);
end

while abs(min(fshape)-obj)>=tol && it<=maxIt
    
    for i=1:size(elements)
        edge1 = [elements(i,1) elements(i,2)];
        edge2 = [elements(i,2) elements(i,3)];
        edge3 = [elements(i,1) elements(i,3)];

        searchEl = elements;
        searchEl(i,:) = [0 0 0];
        
        neigh1 = find(sum((searchEl==edge1(1))+(searchEl==edge1(2)),2)==2);
        neigh2 = find(sum((searchEl==edge2(1))+(searchEl==edge2(2)),2)==2);
        neigh3 = find(sum((searchEl==edge3(1))+(searchEl==edge3(2)),2)==2);
        
        if isempty(neigh1)
            neigh1 = i;
        end
        if isempty(neigh2)
            neigh2 = i;
        end
        if isempty(neigh3)
            neigh3 = i;
        end
        
        localEdge1 = sparse(N,4);
        localEdge2 = sparse(N,4);
        localEdge3 = sparse(N,4);
        
        localEdge1(edge1(1),1) = 1;
        localEdge1(edge1(1),2) = 2;
        find(elements(neigh1,:)==edge1(1))
        localEdge1(edge1(1),3) = find(elements(neigh1,:)==edge1(1));
        find(elements(neigh1,:)==edge1(2))
        localEdge1(edge1(1),4) = find(elements(neigh1,:)==edge1(2));
        
        localEdge2(edge2(1),1) = 2;
        localEdge2(edge2(1),2) = 3;
        find(elements(neigh2,:)==edge2(1))
        localEdge2(edge2(1),3) = find(elements(neigh2,:)==edge2(1));
        find(elements(neigh2,:)==edge2(2))
        localEdge2(edge2(1),4) = find(elements(neigh2,:)==edge2(2));
        
        localEdge3(edge3(2),1) = 1;
        localEdge3(edge3(2),2) = 3;
        find(elements(neigh3,:)==edge3(1))
        localEdge3(edge3(2),3) = find(elements(neigh3,:)==edge3(1));
        find(elements(neigh3,:)==edge3(2))
        localEdge3(edge3(2),4) = find(elements(neigh3,:)==edge3(2));
        
        x1 = nodes(elements([i;neigh1;neigh2;neigh3],1),1);
        x2 = nodes(elements([i;neigh1;neigh2;neigh3],2),1);
        x3 = nodes(elements([i;neigh1;neigh2;neigh3],3),1);
        y1 = nodes(elements([i;neigh1;neigh2;neigh3],1),2);
        y2 = nodes(elements([i;neigh1;neigh2;neigh3],2),2);
        y3 = nodes(elements([i;neigh1;neigh2;neigh3],3),2);

                                                                           % Edges (as vectors in the plane)
        edgeVec1 = [x2-x1 y2-y1];
        edgeVec2 = [x3-x2 y3-y2];
        edgeVec3 = [x3-x1 y3-y1];

        lengths = [sqrt(sum(edgeVec1.^2,2)) sqrt(sum(edgeVec2.^2,2)) sqrt(sum(edgeVec3.^2,2))];
        
                                                                          % Internal angles
        alphas = [acos(sum(edgeVec3.*edgeVec1,2)./(lengths(:,3).*lengths(:,1)))...
                  acos(sum(-edgeVec1.*edgeVec2,2)./(lengths(:,1).*lengths(:,2)))...
                  acos(sum(-edgeVec2.*(-edgeVec3),2)./(lengths(:,2).*lengths(:,3)))].*(180/pi);

        minAlpha = min(alphas,[],2);
        
        isConvex = sum([alphas(1,localEdge1(edge1(1),1))+alphas(2,localEdge1(edge1(1),3)) alphas(1,localEdge1(edge1(1),2))+alphas(2,localEdge1(edge1(1),4));...
                    alphas(1,localEdge2(edge2(1),1))+alphas(3,localEdge2(edge2(1),3)) alphas(1,localEdge2(edge2(1),2))+alphas(3,localEdge2(edge2(1),4));...
                    alphas(1,localEdge3(edge3(2),1))+alphas(4,localEdge3(edge3(2),3)) alphas(1,localEdge3(edge3(2),2))+alphas(4,localEdge3(edge3(2),4))]<180,2)==2;
        
        minBeforeSwap = min([minAlpha(1) minAlpha(2);minAlpha(1) minAlpha(3);minAlpha(1) minAlpha(4)],[],2);
        
        n1 = elements(i,1);
        n2 = elements(i,2);
        n3 = elements(i,3);
        mx = elements(neigh1,find(((elements(neigh1,:)~=edge1(1)) + (elements(neigh1,:)~=edge1(2)))==2));
        px = elements(neigh2,find(((elements(neigh2,:)~=edge2(1)) + (elements(neigh2,:)~=edge2(2)))==2));
        qx = elements(neigh3,find(((elements(neigh3,:)~=edge3(1)) + (elements(neigh3,:)~=edge3(2)))==2));
        
        newelements = [n1 mx n3;...
                       mx n2 n3;...
                       n1 n2 px;...
                       n1 px n3;...
                       n1 n2 qx;...
                       qx n2 n3];
        
        edgesChange = [n1 n2;...
                       mx n3;...
                       n2 n3;...
                       n1 px;...
                       n1 n3;...
                       qx n2];
                   
        x1 = nodes(newelements(:,1),1);
        x2 = nodes(newelements(:,2),1);
        x3 = nodes(newelements(:,3),1);
        y1 = nodes(newelements(:,1),2);
        y2 = nodes(newelements(:,2),2);
        y3 = nodes(newelements(:,3),2);

                                                                           % Edges (as vectors in the plane)
        edgeVec1 = [x2-x1 y2-y1];
        edgeVec2 = [x3-x2 y3-y2];
        edgeVec3 = [x3-x1 y3-y1];

        lengths = [sqrt(sum(edgeVec1.^2,2)) sqrt(sum(edgeVec2.^2,2)) sqrt(sum(edgeVec3.^2,2))];
        
                                                                           % Internal angles
        alphas = [acos(sum(edgeVec3.*edgeVec1,2)./(lengths(:,3).*lengths(:,1)))...
                  acos(sum(-edgeVec1.*edgeVec2,2)./(lengths(:,1).*lengths(:,2)))...
                  acos(sum(-edgeVec2.*(-edgeVec3),2)./(lengths(:,2).*lengths(:,3)))].*(180/pi);

        minAlpha = min(alphas,[],2);
        
        minAfterSwap = min([minAlpha(1) minAlpha(2);minAlpha(3) minAlpha(4);minAlpha(5) minAlpha(6)],[],2);
                
        neighs = [neigh1;...
                  neigh2;...
                  neigh3];
        
        if ~isempty(find(isConvex))
            switch sum(isConvex)
                case 1
                    swapIndex = find(isConvex,1);
                    if minBeforeSwap(swapIndex) < minAfterSwap(swapIndex)
                        elements(i,:) = newelements(2*(swapIndex-1)+1,:);
                        elements(neighs(swapIndex),:) = newelements(2*(swapIndex-1)+2,:);
                        edgeIndex = find(sum(((edges==edgesChange(2*(swapIndex-1)+1,1))+(edges==edgesChange(2*(swapIndex-1)+1,2))),2)==2);
                        edges(edgeIndex,:) = edgesChange(2*(swapIndex-1)+2,:);
                    end
                case 2
                    swapIndeces = find(isConvex);
                    deltaSwap = minAfterSwap(swapIndeces)-minBeforeSwap(swapIndeces);
                    [Y,index] = max(deltaSwap);
                    swapIndex = swapIndeces(index);
                    if minBeforeSwap(swapIndex) < minAfterSwap(swapIndex)
                        elements(i,:) = newelements(2*(swapIndex-1)+1,:);
                        elements(neighs(swapIndex),:) = newelements(2*(swapIndex-1)+2,:);
                        edgeIndex = find(sum(((edges==edgesChange(2*(swapIndex-1)+1,1))+(edges==edgesChange(2*(swapIndex-1)+1,2))),2)==2);
                        edges(edgeIndex,:) = edgesChange(2*(swapIndex-1)+2,:);
                    end
                case 3
                    deltaSwap = minAfterSwap-minBeforeSwap;
                    [Y,swapIndex] = max(deltaSwap);
                    if deltaSwap(swapIndex)>0
                        elements(i,:) = newelements(2*(swapIndex-1)+1,:);
                        elements(neighs(swapIndex),:) = newelements(2*(swapIndex-1)+2,:);
                        edgeIndex = find(sum(((edges==edgesChange(2*(swapIndex-1)+1,1))+(edges==edgesChange(2*(swapIndex-1)+1,2))),2)==2);
                        edges(edgeIndex,:) = edgesChange(2*(swapIndex-1)+2,:);
                    end
            end
        end

    end
    
    d = ones(dof*N,1);

    for i=1:N
        if isempty(find(boundary==i))
            connect = length(find(elements(:,1)==i)) + ...
                      length(find(elements(:,2)==i)) + ...
                      length(find(elements(:,3)==i));
            d(dof*(i-1)+1) = 1/connect;
            d(dof*(i-1)+2) = 1/connect;
        end
    end

    D = spdiags(d,0,dof*N,dof*N);

    clear d

    k = [0           0           0.5          sqrt(3)/2   0.5         -sqrt(3)/2;...
         0           0           -sqrt(3)/2   0.5         sqrt(3)/2   0.5;...
         0.5         -sqrt(3)/2  0            0           0.5         sqrt(3)/2;...
         sqrt(3)/2   0.5         0            0           -sqrt(3)/2  0.5;...
         0.5         sqrt(3)/2   0.5          -sqrt(3)/2  0           0;...
         -sqrt(3)/2  0.5         sqrt(3)/2    0.5         0           0];

    K = sparse(dof*N,dof*N);

    for i=1:size(elements,1)
        K = K + nodalmaps(N,elements(i,:),k,dof,1);
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

    solveNodes = D*(K*solveNodes);
    for i=1:N
        nodes(i,1) = solveNodes(dof*(i-1)+1);
        nodes(i,2) = solveNodes(dof*(i-1)+2);
    end
    fshape = tri3quality(nodes,elements);
    it = it + 1;
end

return

