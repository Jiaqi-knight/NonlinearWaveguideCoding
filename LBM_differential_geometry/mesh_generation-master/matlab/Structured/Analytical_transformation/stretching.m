function[mesh]=stretching(oldmesh,funcs,args)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: March 27th, 2014
%    Last update: April 1st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

mesh = oldmesh;

switch oldmesh.D
    case 1
        g = inline(char(funcs(1)),char(args(1)));
        if oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g(oldmesh.nodes.meshcoordinates(1,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g(oldmesh.coordlines.dim1free(l,1));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.edges.nodes(2,i)))/2;
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g(oldmesh.nodes.meshcoordinates(1,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g(oldmesh.coordlines.dim1free(l,1));
            end
        end
    case 2
        g1 = inline(char(funcs(1)),char(args(1)),char(args(2)));
        g2 = inline(char(funcs(2)),char(args(1)),char(args(2)));
        if oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.edges.nodes(2,i)))/2;
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(4,i)))/4;
                mesh.faces.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(4,i)))/4;
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.edges.nodes(2,i)))/2;
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2));
            end
        end
    case 3
        g1 = inline(char(funcs(1)),char(args(1)),char(args(2)),char(args(3)));
        g2 = inline(char(funcs(2)),char(args(1)),char(args(2)),char(args(3)));
        g3 = inline(char(funcs(3)),char(args(1)),char(args(2)),char(args(3)));
        if oldmesh.edgeflag && oldmesh.faceflag && oldmesh.cellflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,3) = g3(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,3) = g3(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
            end
            for l=1:mesh.Nlinesdim3
                mesh.coordlines.dim3free(l,1) = g1(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,2) = g2(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,3) = g3(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(3,i) = (mesh.nodes.meshcoordinates(3,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(3,mesh.edges.nodes(2,i)))/2;
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(4,i)))/4;
                mesh.faces.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(4,i)))/4;
                mesh.faces.centroid(3,i) = (mesh.nodes.meshcoordinates(3,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(3,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(3,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(3,mesh.faces.nodes(4,i)))/4;
            end
            for i=1:oldmesh.totC
                mesh.cells.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.cells.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(2,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(3,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(4,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(5,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(6,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(7,i)) + mesh.nodes.meshcoordinates(1,mesh.cells.nodes(8,i)))/8;
                mesh.cells.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.cells.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(2,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(3,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(4,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(5,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(6,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(7,i)) + mesh.nodes.meshcoordinates(2,mesh.cells.nodes(8,i)))/8;
                mesh.cells.centroid(3,i) = (mesh.nodes.meshcoordinates(3,mesh.cells.nodes(1,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(2,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(3,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(4,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(5,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(6,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(7,i)) + mesh.nodes.meshcoordinates(3,mesh.cells.nodes(8,i)))/8;
            end
        elseif oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,3) = g3(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,3) = g3(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
            end
            for l=1:mesh.Nlinesdim3
                mesh.coordlines.dim3free(l,1) = g1(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,2) = g2(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,3) = g3(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(3,i) = (mesh.nodes.meshcoordinates(3,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(3,mesh.edges.nodes(2,i)))/2;
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(1,mesh.faces.nodes(4,i)))/4;
                mesh.faces.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(2,mesh.faces.nodes(4,i)))/4;
                mesh.faces.centroid(3,i) = (mesh.nodes.meshcoordinates(3,mesh.faces.nodes(1,i)) + mesh.nodes.meshcoordinates(3,mesh.faces.nodes(2,i)) + mesh.nodes.meshcoordinates(3,mesh.faces.nodes(3,i)) + mesh.nodes.meshcoordinates(3,mesh.faces.nodes(4,i)))/4;
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,3) = g3(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,3) = g3(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
            end
            for l=1:mesh.Nlinesdim3
                mesh.coordlines.dim3free(l,1) = g1(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,2) = g2(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,3) = g3(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = (mesh.nodes.meshcoordinates(1,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(1,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(2,i) = (mesh.nodes.meshcoordinates(2,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(2,mesh.edges.nodes(2,i)))/2;
                mesh.edges.centroid(3,i) = (mesh.nodes.meshcoordinates(3,mesh.edges.nodes(1,i)) + mesh.nodes.meshcoordinates(3,mesh.edges.nodes(2,i)))/2;
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,1) = g1(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,2) = g2(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
                mesh.coordlines.dim1free(l,3) = g3(oldmesh.coordlines.dim1free(l,1),oldmesh.coordlines.dim1free(l,2),oldmesh.coordlines.dim1free(l,3));
            end
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,1) = g1(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,2) = g2(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
                mesh.coordlines.dim2free(l,3) = g3(oldmesh.coordlines.dim2free(l,1),oldmesh.coordlines.dim2free(l,2),oldmesh.coordlines.dim2free(l,3));
            end
            for l=1:mesh.Nlinesdim3
                mesh.coordlines.dim3free(l,1) = g1(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,2) = g2(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
                mesh.coordlines.dim3free(l,3) = g3(oldmesh.coordlines.dim3free(l,1),oldmesh.coordlines.dim3free(l,2),oldmesh.coordlines.dim3free(l,3));
            end
        end
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return