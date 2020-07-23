function[Aordered,indeces]=mintournamenttree(A,column)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 2.0
%  Creation date: May 14th, 2014
%    Last update: May 20th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

if nargin==1
    column = 1;
end

N = size(A,1);

if N==1
    Aordered = A;
    indeces = 1;
    return
end

if N==2
    if A(1,column)<A(2,column)
        Aordered = A;
        indeces = [1;2];
        return
    else
        Aordered = [A(2,:);A(1,:)];
        indeces = [2;1];
        return
    end
end

indeces = zeros(N,1);
Aordered = zeros(N,size(A,2));

%-------------------------- build binary tree -----------------------------

m = floor(log2(N));
p = N - 2^m;

Nnodes = 2^(m+1)-1+2*p;

% if check
%     if p==0
%         depthmax = m+1;
%     else
%         depthmax = m+2;
%     end
% end

globalindeces = (1:Nnodes)';
levels = floor(log2(globalindeces));
localindeces = globalindeces - 2.^levels + 1;

tree = [globalindeces [0; globalindeces(2:end,:)-localindeces(2:end,:)-(2.^(levels(2:end,:)-1)-1)+ceil(0.5.*localindeces(2:end,:))-1] [2.^(levels(1:Nnodes-N,:)+1)+2.*(localindeces(1:Nnodes-N,:)-1); zeros(N,1)] [2.^(levels(1:Nnodes-N,:)+1)+2.*(localindeces(1:Nnodes-N,:)-1)+1; zeros(N,1)] [zeros(Nnodes-N,1); (1:N)'] [zeros(Nnodes-N,1); A(:,column)] [zeros(Nnodes-N,1); (Nnodes-N+1:Nnodes)']];

% if check
%     ftree = representtree(tree,depthmax,1);
%     pause
% end

%---------------------- play first sorting match --------------------------

sorted = 0;

player1 = Nnodes;
player2 = Nnodes - 1;
play = 1;
while play
    parent = tree(player1,2);
    if tree(player1,6)<tree(player2,6)
        tree(parent,5) = tree(player1,5);
        tree(parent,6) = tree(player1,6);
        tree(parent,7) = tree(player1,7);
    else
        tree(parent,5) = tree(player2,5);
        tree(parent,6) = tree(player2,6);
        tree(parent,7) = tree(player2,7);
    end
    if parent==1
        play = 0;
    else
        player1 = player1 - 2;
        player2 = player2 - 2;
    end
end

indeces(sorted + 1,1) = tree(1,5);
Aordered(sorted + 1,:) = A(tree(1,5),:);
sorted = sorted + 1;
%------------------------------ replay ------------------------------------

while sorted<N-1
    play = 0;
    winnerposition = 1;
    leaf = tree(winnerposition,7);
    tree(leaf,5) = 0;
    tree(leaf,6) = Inf;
    tree(leaf,7) = 0;
%   if check
%       ftree = representtree(tree,depthmax,1);
%       pause
%   end
    rootparentindex = tree(leaf,2);
    if rootparentindex~=1
        if leaf==tree(rootparentindex,3)
            brother = tree(rootparentindex,4);
        else
            brother = tree(rootparentindex,3);
        end
        prune = 1;
        while prune
            if tree(brother,3)==0 && tree(brother,4)==0
                tree(rootparentindex,5) = tree(brother,5);
                tree(rootparentindex,6) = tree(brother,6);
                tree(rootparentindex,7) = rootparentindex;
                tree(brother,5) = 0;
                tree(brother,6) = Inf;
                tree(brother,7) = 0;
                prune = 0;
            elseif tree(brother,6)==Inf
                grandparent = tree(rootparentindex,2);
                if rootparentindex==tree(grandparent,3)
                    uncle = tree(grandparent,4);
                else
                    uncle = tree(grandparent,3);
                end
                if tree(uncle,6)==Inf
                    tree(rootparentindex,5) = 0;
                    tree(rootparentindex,6) = Inf;
                    tree(rootparentindex,7) = 0;
                    rootparentindex = grandparent;
                    brother = uncle;
                else
                    tree(grandparent,5) = tree(uncle,5);
                    tree(grandparent,6) = tree(uncle,6);
                    tree(grandparent,7) = tree(uncle,7);
                    tree(rootparentindex,5) = 0;
                    tree(rootparentindex,6) = Inf;
                    tree(rootparentindex,7) = 0;
                    rootparentindex = grandparent;
                    prune = 0;
                end                      
             else
                tree(rootparentindex,5) = tree(brother,5);
                tree(rootparentindex,6) = tree(brother,6);
                tree(rootparentindex,7) = tree(brother,7);
                prune = 0;
             end
        end
    else
        if leaf==tree(rootparentindex,3)
             rootparentindex = tree(rootparentindex,4);
        else
             rootparentindex = tree(rootparentindex,3);
        end
    end
   if rootparentindex~=1 && rem(rootparentindex,2)==1
        player1 = rootparentindex;
        player2 = rootparentindex - 1;
        play = 1;
   elseif rootparentindex~=1 && rem(rootparentindex,2)==0
        player1 = rootparentindex + 1;
        player2 = rootparentindex;
        play = 1;
   else
        player1 = 3;
        player2 = 2;
   end
   while play
        parent = tree(player1,2);
        if player1~=0 && player2~=0
            if tree(player1,6)<tree(player2,6)
                tree(parent,5) = tree(player1,5);
                tree(parent,6) = tree(player1,6);
                tree(parent,7) = tree(player1,7);
            else
                tree(parent,5) = tree(player2,5);
                tree(parent,6) = tree(player2,6);
                tree(parent,7) = tree(player2,7);
            end
        end
        if parent<=1
            play = 0;
        else
            if rem(parent,2)==1
                player1 = parent;
                player2 = parent - 1;
            else
                player1 = parent + 1;
                player2 = parent;
            end
        end
    end
    indeces(sorted + 1,1) = tree(tree(winnerposition,7),5);
    Aordered(sorted + 1,:) = A(tree(tree(winnerposition,7),5),:);
    sorted = sorted + 1;
end


winnerposition = 1;
leaf = tree(winnerposition,7);
tree(leaf,5) = 0;
tree(leaf,6) = Inf;
tree(leaf,7) = 0;
parent = tree(leaf,2);
if leaf==tree(parent,3)
    child = tree(parent,4);
else
    child = tree(parent,3);
end
searchlast = 1;
while searchlast
    if tree(child,6)~=Inf && tree(child,6)~=A(indeces(sorted,1),column) && tree(child,6)~=A(indeces(sorted-1,1),column)
        tree(winnerposition,5) = tree(child,5);
        tree(winnerposition,6) = tree(child,6);
        tree(winnerposition,7) = child;
        searchlast = 0;
    else
        leaf = parent;
        parent = tree(leaf,2);
        if leaf==tree(parent,3)
            child = tree(parent,4);
        else
            child = tree(parent,3);
        end
    end
end
indeces(sorted + 1,1) = tree(tree(winnerposition,7),5);
Aordered(sorted + 1,:) = A(tree(tree(winnerposition,7),5),:);

return