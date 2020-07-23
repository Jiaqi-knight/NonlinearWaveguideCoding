function[mesh]=analytictransformation(D,oldmesh,funcs,args)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: March 27th, 2014
%    Last update: June 12th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

mesh = oldmesh;

switch D
    case 2
        g1 = inline(char(funcs(1)),char(args(1)),char(args(2)));
        g2 = inline(char(funcs(2)),char(args(1)),char(args(2)));
        if oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = g1(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i));
                mesh.edges.centroid(2,i) = g2(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i));
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = g1(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i));
                mesh.faces.centroid(2,i) = g2(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i));
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = g1(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i));
                mesh.edges.centroid(2,i) = g2(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i));
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i));
            end
        end
         for j=1:mesh.Ndim2
            coord3 = 0;
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,3*(j-1) + 1) = g1(oldmesh.coordlines.dim1free(l,3*(j-1) + 1),oldmesh.coordlines.dim1free(l,3*(j-1) + 2));
                mesh.coordlines.dim1free(l,3*(j-1) + 2) = g2(oldmesh.coordlines.dim1free(l,3*(j-1) + 1),oldmesh.coordlines.dim1free(l,3*(j-1) + 2)); 
                mesh.coordlines.dim1free(l,3*(j-1) + 3) = coord3;
            end
        end
        for k=1:mesh.Ndim1
            coord3 = 0;
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,3*(k-1) + 1) = g1(oldmesh.coordlines.dim1free(l,3*(k-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1) + 2));
                mesh.coordlines.dim2free(l,3*(k-1) + 2) = g2(oldmesh.coordlines.dim1free(l,3*(k-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1) + 2)); 
                mesh.coordlines.dim2free(l,3*(k-1) + 3) = coord3;
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
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = g1(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
                mesh.edges.centroid(2,i) = g2(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
                mesh.edges.centroid(3,i) = g3(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = g1(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i),oldmesh.faces.centroid(3,i));
                mesh.faces.centroid(2,i) = g2(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i),oldmesh.faces.centroid(3,i));
                mesh.faces.centroid(3,i) = g3(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i),oldmesh.faces.centroid(3,i));
            end
            for i=1:oldmesh.totC
                mesh.cells.centroid(1,i) = g1(oldmesh.cells.centroid(1,i),oldmesh.cells.centroid(2,i),oldmesh.cells.centroid(3,i));
                mesh.cells.centroid(2,i) = g2(oldmesh.cells.centroid(1,i),oldmesh.cells.centroid(2,i),oldmesh.cells.centroid(3,i));
                mesh.cells.centroid(3,i) = g3(oldmesh.cells.centroid(1,i),oldmesh.cells.centroid(2,i),oldmesh.cells.centroid(3,i));
            end
        elseif oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = g1(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
                mesh.edges.centroid(2,i) = g2(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
                mesh.edges.centroid(3,i) = g3(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = g1(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i),oldmesh.faces.centroid(3,i));
                mesh.faces.centroid(2,i) = g2(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i),oldmesh.faces.centroid(3,i));
                mesh.faces.centroid(3,i) = g3(oldmesh.faces.centroid(1,i),oldmesh.faces.centroid(2,i),oldmesh.faces.centroid(3,i));
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = g1(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
                mesh.edges.centroid(2,i) = g2(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
                mesh.edges.centroid(3,i) = g3(oldmesh.edges.centroid(1,i),oldmesh.edges.centroid(2,i),oldmesh.edges.centroid(3,i));
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = g1(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(2,i) = g2(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
                mesh.nodes.meshcoordinates(3,i) = g3(oldmesh.nodes.meshcoordinates(1,i),oldmesh.nodes.meshcoordinates(2,i),oldmesh.nodes.meshcoordinates(3,i));
            end
        end
        for i=1:mesh.Ndim3
            for j=1:mesh.Ndim2
                for l=1:mesh.Nlinesdim1
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1) = g1(oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1),oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 2),mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 3));
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 2) = g2(oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1),oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 2),mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 3)); 
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 3) = g3(oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1),oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 2),mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 3));
                end
            end
        end
        for i=1:mesh.Ndim3
            for k=1:mesh.Ndim1
                for l=1:mesh.Nlinesdim2
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 1) = g1(oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2),mesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 3));
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2) = g2(oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2),mesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 3)); 
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 3) = g3(oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2),mesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 3));
                end
            end
        end
        for j=1:mesh.Ndim2
            for k=1:mesh.Ndim1
                for l=1:mesh.Nlinesdim3
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 1) = g1(oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 2),mesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3));
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 2) = g2(oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 2),mesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3)); 
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3) = g3(oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 1),oldmesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 2),mesh.coordlines.dim1free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3));
                end
            end
        end
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return