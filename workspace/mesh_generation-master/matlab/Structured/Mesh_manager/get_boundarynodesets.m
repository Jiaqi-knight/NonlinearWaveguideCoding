function[outmesh]=get_boundarynodesets(inpmesh,type,funcs,operator,bounds,args)

% integer mask for boundary type:
% 1 --> imposed displacement (for solid domain)
% 2 --> imposed velocity (for solid/fluid domain)
% 3 --> imposed flux (for fluid domain)
% 4 --> imposed density (for fluid domain)
% 5 --> imposed pressure (for fluid domain)
% 6 --> imposed distributed load (for solid domain)
% 7 --> imposed concentrated load (for solid domain)
% 8 --> fluid-solid interface (for solid/fluid domain)
% 9 --> fluid-fluid interface (for fluid domain)
% 10 --> free surface (for fluid domain)

% integer mask for operator:
% 1 --> g() <
% 2 --> g() <=
% 3 --> g() ==
% 4 --> g() >=
% 5 --> g() >
% 6 --> < g() >
% 7 --> <= g() >
% 8 --> < g() >=
% 9 --> <= g() >=

outmesh = inpmesh;

nbounds = length(funcs);

outmesh.nodeboundaries = zeros(nbounds,3);
outmesh.boundarynodesets = [];

for i=1:nbounds
    outmesh.nodeboundaries(i,1) = i;
    outmesh.nodeboundaries(i,2) = type(i,1);
    outmesh.nodeboundaries(i,3) = length(outmesh.boundarynodesets) + 1;
    checknodes = [mesh.nodes.id' mesh.nodes.meshcoordinates'];
    foundnodes = [];
    if i==nbounds
        conditionsnum = size(type,1) - type(i,2);
    else
        conditionsnum = type(i+1,2) - type(i,2);
    end
    for j=1:conditionsnum
        g = inline(char(funcs(type(i,2)+(j-1))),char(args(1)),char(args(2)),char(args(3)));
        switch operator(type(i,2)+(j-1))
            case 1
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) < bounds(type(i,2)+(j-1),1)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end
            case 2
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) <= bounds(type(i,2)+(j-1),1)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end  
            case 3
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) == bounds(type(i,2)+(j-1),1)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end 
                end  
            case 4
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) >= bounds(type(i,2)+(j-1),1)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end  
            case 5
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) > bounds(type(i,2)+(j-1),1)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end
            case 6
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) < bounds(type(i,2)+(j-1),1) && g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) > bounds(type(i,2)+(j-1),2)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end
            case 7
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) <= bounds(type(i,2)+(j-1),1) && g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) >= bounds(type(i,2)+(j-1),2)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end
            case 8
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) < bounds(type(i,2)+(j-1),1) && g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) > bounds(type(i,2)+(j-1),2)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end
            case 9
                for k=1:size(checknodes,1)
                    if g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) <= bounds(type(i,2)+(j-1),1) && g(checknodes(k,2),checknodes(k,3),checknodes(k,4)) >= bounds(type(i,2)+(j-1),2)
                        foundnodes = [foundnodes; checknodes(k,:)];
                    end
                end
            otherwise
                disp('The flag provided does not correspond to any operator.');
        end
        checknodes = foundnodes;
        foundnodes = [];
    end
    outmesh.boundarynodesets = [outmesh.boundarynodesets; checknodes(:,1)];  
end

return