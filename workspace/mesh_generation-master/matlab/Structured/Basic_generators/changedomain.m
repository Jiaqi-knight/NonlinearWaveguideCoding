function[mesh]=changedomain(D,oldmesh,dim1min,dim1max,dim2min,dim2max,dim3min,dim3max)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: June 12th, 2014
%    Last update: June 12th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

mesh = oldmesh;

switch D
    case 2
        if oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = dim1min + (oldmesh.edges.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.edges.centroid(2,i) = dim2min + (oldmesh.edges.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = dim1min + (oldmesh.faces.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.faces.centroid(2,i) = dim2min + (oldmesh.faces.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = dim1min + (oldmesh.edges.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.edges.centroid(2,i) = dim2min + (oldmesh.edges.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            end
        end
        for j=1:mesh.Ndim2
            coord2 = dim2min + ((oldmesh.dim2min + oldmesh.deltadim2*(j-1))-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
            coord3 = 0;
            for l=1:mesh.Nlinesdim1
                mesh.coordlines.dim1free(l,3*(j-1) + 1) = dim1min + (oldmesh.coordlines.dim1free(l,3*(j-1) + 1)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.coordlines.dim1free(l,3*(j-1) + 2) = coord2; 
                mesh.coordlines.dim1free(l,3*(j-1) + 3) = coord3;
            end
        end
        for k=1:mesh.Ndim1
            coord1 = dim1min + ((oldmesh.dim1min + oldmesh.deltadim1*(k-1))-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
            coord3 = 0;
            for l=1:mesh.Nlinesdim2
                mesh.coordlines.dim2free(l,3*(k-1) + 1) = coord1;
                mesh.coordlines.dim2free(l,3*(k-1) + 2) = dim2min + (oldmesh.coordlines.dim2free(l,3*(k-1) + 2)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min); 
                mesh.coordlines.dim2free(l,3*(k-1) + 3) = coord3;
            end
        end
    case 3
        if oldmesh.edgeflag && oldmesh.faceflag && oldmesh.cellflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.nodes.meshcoordinates(3,i) = dim3min + (oldmesh.nodes.meshcoordinates(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = dim1min + (oldmesh.edges.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.edges.centroid(2,i) = dim2min + (oldmesh.edges.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.edges.centroid(3,i) = dim3min + (oldmesh.edges.centroid(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = dim1min + (oldmesh.faces.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.faces.centroid(2,i) = dim2min + (oldmesh.faces.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.faces.centroid(3,i) = dim3min + (oldmesh.faces.centroid(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
            for i=1:oldmesh.totC
                mesh.cells.centroid(1,i) = dim1min + (oldmesh.cells.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.cells.centroid(2,i) = dim2min + (oldmesh.cells.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.cells.centroid(3,i) = dim3min + (oldmesh.cells.centroid(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
        elseif oldmesh.edgeflag && oldmesh.faceflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.nodes.meshcoordinates(3,i) = dim3min + (oldmesh.nodes.meshcoordinates(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = dim1min + (oldmesh.edges.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.edges.centroid(2,i) = dim2min + (oldmesh.edges.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.edges.centroid(3,i) = dim3min + (oldmesh.edges.centroid(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
            for i=1:oldmesh.totF
                mesh.faces.centroid(1,i) = dim1min + (oldmesh.faces.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.faces.centroid(2,i) = dim2min + (oldmesh.faces.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.faces.centroid(3,i) = dim3min + (oldmesh.faces.centroid(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
        elseif oldmesh.edgeflag
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.nodes.meshcoordinates(3,i) = dim3min + (oldmesh.nodes.meshcoordinates(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
            for i=1:oldmesh.totE
                mesh.edges.centroid(1,i) = dim1min + (oldmesh.edges.centroid(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.edges.centroid(2,i) = dim2min + (oldmesh.edges.centroid(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.edges.centroid(3,i) = dim3min + (oldmesh.edges.centroid(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
        else
            for i=1:oldmesh.totN
                mesh.nodes.meshcoordinates(1,i) = dim1min + (oldmesh.nodes.meshcoordinates(1,i)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                mesh.nodes.meshcoordinates(2,i) = dim2min + (oldmesh.nodes.meshcoordinates(2,i)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                mesh.nodes.meshcoordinates(3,i) = dim3min + (oldmesh.nodes.meshcoordinates(3,i)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
            end
        end
        for i=1:mesh.Ndim3
            for j=1:mesh.Ndim2
                coord2 = dim2min + ((oldmesh.dim2min + oldmesh.deltadim2*(j-1))-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                coord3 = dim3min + ((mesh.dim3min + mesh.deltadim3*(i-1))-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
                for l=1:mesh.Nlinesdim1
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1) = dim1min + (oldmesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 1)-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 2) = coord2; 
                    mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1) + 3) = coord3;
                end
            end
        end
        for i=1:mesh.Ndim3
            for k=1:mesh.Ndim1
                coord1 = dim1min + ((oldmesh.dim1min + oldmesh.deltadim1*(k-1))-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                coord3 = dim3min + ((mesh.dim3min + mesh.deltadim3*(i-1))-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
                for l=1:mesh.Nlinesdim2
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 1) = coord1;
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2) = dim2min + (oldmesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 2)-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min); 
                    mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1) + 3) = coord3;
                end
            end
        end
        for j=1:mesh.Ndim2
            for k=1:mesh.Ndim1
                coord1 = dim1min + ((oldmesh.dim1min + oldmesh.deltadim1*(k-1))-oldmesh.dim1min)*(dim1max-dim1min)/(oldmesh.dim1max-oldmesh.dim1min);
                coord2 = dim2min + ((oldmesh.dim2min + oldmesh.deltadim2*(j-1))-oldmesh.dim2min)*(dim2max-dim2min)/(oldmesh.dim2max-oldmesh.dim2min);
                for l=1:mesh.Nlinesdim3
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 1) = coord1;
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 2) = coord2; 
                    mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3) = dim3min + (oldmesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1) + 3)-oldmesh.dim3min)*(dim3max-dim3min)/(oldmesh.dim3max-oldmesh.dim3min);
                end
            end
        end
    otherwise
        disp('Dimension of space requested is currently not supported')
end

return