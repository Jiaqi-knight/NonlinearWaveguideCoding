function[sets]=get_nodesets2D(N,lattice,operator,bounds,functions)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 10th, 2014
%    Last update: July 10th, 2014
%
%    Description: 
%          Input: 
%         Output:
%
%          integer mask for set type:
%          1 --> input
%          2 --> measurements
%          3 --> control
%
%          integer mask for operator:
%          1 --> g() <
%          2 --> g() <=
%          3 --> g() ==
%          4 --> g() >=
%          5 --> g() >
%          6 --> < g() <
%          7 --> <= g() <
%          8 --> < g() <=
%          9 --> <= g() <=

%%

nbounds = length(operator);

sets = zeros(N,nbounds);

for j=1:nbounds
    for i=1:N
        switch operator(j)
            case 1
                if functions{j}(lattice(i,3),lattice(i,4))<bounds(j,1)
                    sets(i,j) = 1;
                end
            case 2
                if functions{j}(lattice(i,3),lattice(i,4))<=bounds(j,1)
                    sets(i,j) = 1;
                end
            case 3
                if functions{j}(lattice(i,3),lattice(i,4))==bounds(j,1)
                    sets(i,j) = 1;
                end
            case 4
                if functions{j}(lattice(i,3),lattice(i,4))>=bounds(j,1)
                    sets(i,j) = 1;
                end
            case 5
                if functions{j}(lattice(i,3),lattice(i,4))>bounds(j,1)
                    sets(i,j) = 1;
                end
            case 6
                if functions{j}(lattice(i,3),lattice(i,4))>bounds(j,1)  && functions{j}(lattice(i,3),lattice(i,4))<bounds(j,2)
                    sets(i,j) = 1;
                end
            case 7
                if functions{j}(lattice(i,3),lattice(i,4))>=bounds(j,1) && functions{j}(lattice(i,3),lattice(i,4))<bounds(j,2)
                    sets(i,j) = 1;
                end
            case 8
                if functions{j}(lattice(i,3),lattice(i,4))>bounds(j,1)  && functions{j}(lattice(i,3),lattice(i,4))<=bounds(j,2)
                    sets(i,j) = 1;
                end
            case 9
                if functions{j}(lattice(i,3),lattice(i,4))>=bounds(j,1) && functions{j}(lattice(i,3),lattice(i,4))<=bounds(j,2)
                    sets(i,j) = 1;
                end
        end
    end
end

return