function[]=export_csv(mesh,filename,printflag)

c = clock;
basefilename = strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',filename,'_basedata','.csv']);
nodesfilename = strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',filename,'_nodedata','.csv']);
linesfilename = strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',filename,'_linedata','.csv']);
edgesfilename = strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',filename,'_edgedata','.csv']);
facesfilename = strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',filename,'_facedata','.csv']);
cellsfilename = strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3)),'_',filename,'_celldata','.csv']);

basefid = fopen(basefilename,'w');

dlmwrite(basefilename,mesh.D,'-append');
dlmwrite(basefilename,mesh.dim1min,'-append');
dlmwrite(basefilename,mesh.dim1max,'-append');
dlmwrite(basefilename,mesh.Ndim1,'-append');
dlmwrite(basefilename,mesh.dim2min,'-append');
dlmwrite(basefilename,mesh.dim2max,'-append');
dlmwrite(basefilename,mesh.Ndim2,'-append');
dlmwrite(basefilename,mesh.dim3min,'-append');
dlmwrite(basefilename,mesh.dim3max,'-append');
dlmwrite(basefilename,mesh.Ndim3,'-append');
dlmwrite(basefilename,mesh.Nlinesdim1,'-append');
dlmwrite(basefilename,mesh.Nlinesdim2,'-append');
dlmwrite(basefilename,mesh.Nlinesdim3,'-append');
dlmwrite(basefilename,mesh.totN,'-append');
dlmwrite(basefilename,mesh.totE,'-append');
dlmwrite(basefilename,mesh.totF,'-append');
dlmwrite(basefilename,mesh.totC,'-append');
dlmwrite(basefilename,mesh.edgeflag,'-append');
dlmwrite(basefilename,mesh.faceflag,'-append');
dlmwrite(basefilename,mesh.cellflag,'-append');

fclose(basefid);

switch mesh.D
%-----------------------------------------------------------------------------1D--------------------------------------------------------------------------------------------
    case 1
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('EXPORT MESH TO .csv FILE')
            disp(' ')
            disp('1D mesh')
            fprintf('Start point x0 = %d\n',mesh.dim1min)
            fprintf('End point xN = %6.2f\n',mesh.dim1max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        if mesh.edgeflag
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            edgefid = fopen(edgesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)' mesh.nodes.edges(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            for i=1:mesh.totE
                dlmwrite(filename,[mesh.edges.id(i) mesh.edges.centroid(:,i)' mesh.edges.nodes(:,i)'],'-append');
            end
            fclose(nodefid);
            fclose(linefid);
            fclose(edgefid);
        else
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            fclose(nodefid);
            fclose(linefid);
        end
%-----------------------------------------------------------------------------2D--------------------------------------------------------------------------------------------
    case 2
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('EXPORT MESH TO .csv FILE')
            disp(' ')
            disp('2D mesh')
            fprintf('Lower boundary along the first direction x0 = %d\n',mesh.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',mesh.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',mesh.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',mesh.dim2max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        if mesh.faceflag && mesh.edgeflag
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            edgefid = fopen(edgesfilename,'w');
            facefid = fopen(facesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)' mesh.nodes.edges(:,i)' mesh.nodes.faces(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            for i=1:mesh.totE
                dlmwrite(filename,[mesh.edges.id(i) mesh.edges.centroid(:,i)' mesh.edges.nodes(:,i)' mesh.edges.faces(:,i)'],'-append');
            end
            for i=1:mesh.totF
                dlmwrite(filename,[mesh.faces.id(i) mesh.faces.centroid(:,i)' mesh.faces.nodes(:,i)' mesh.faces.edges(:,i)'],'-append');
            end
            fclose(nodefid);
            fclose(linefid);
            fclose(edgefid);
            fclose(facefid);
        elseif mesh.edgeflag
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            edgefid = fopen(edgesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)' mesh.nodes.edges(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            for i=1:mesh.totE
                dlmwrite(filename,[mesh.edges.id(i) mesh.edges.centroid(:,i)' mesh.edges.nodes(:,i)'],'-append');
            end
            fclose(nodefid);
            fclose(linefid);
            fclose(edgefid);
        else
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            fclose(nodefid);
            fclose(linefid);
        end
%-----------------------------------------------------------------------------3D--------------------------------------------------------------------------------------------
    case 3
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('EXPORT MESH TO .csv FILE')
            disp(' ')
            disp('3D mesh')
            disp(' ')
            fprintf('Lower boundary along the first direction x0 = %d\n',mesh.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',mesh.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',mesh.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',mesh.dim2max)
            fprintf('Lower boundary along the second direction z0 = %d\n',mesh.dim3min)
            fprintf('Upper boundary along the second direction zN = %6.2f\n',mesh.dim3max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        if mesh.cellflag && mesh.faceflag && mesh.edgeflag
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            edgefid = fopen(edgesfilename,'w');
            facefid = fopen(facesfilename,'w');
            cellfid = fopen(cellsfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)' mesh.nodes.edges(:,i)' mesh.nodes.faces(:,i)' mesh.nodes.cells(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim3free,'-append');
            for i=1:mesh.totE
                dlmwrite(filename,[mesh.edges.id(i) mesh.edges.centroid(:,i)' mesh.edges.nodes(:,i)' mesh.edges.faces(:,i)' mesh.edges.cells(:,i)'],'-append');
            end
            for i=1:mesh.totF
                dlmwrite(filename,[mesh.faces.id(i) mesh.faces.centroid(:,i)' mesh.faces.nodes(:,i)' mesh.faces.edges(:,i)' mesh.faces.cells(:,i)'],'-append');
            end
            for i=1:mesh.totC
                dlmwrite(filename,[mesh.cells.id(i) mesh.cells.centroid(:,i)' mesh.cells.nodes(:,i)' mesh.cells.edges(:,i)' mesh.cells.faces(:,i)'],'-append');
            end
            fclose(nodefid);
            fclose(linefid);
            fclose(edgefid);
            fclose(facefid);
            fclose(cellfid);
        elseif mesh.faceflag && mesh.edgeflag
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            edgefid = fopen(edgesfilename,'w');
            facefid = fopen(facesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)' mesh.nodes.edges(:,i)' mesh.nodes.faces(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim3free,'-append');
            for i=1:mesh.totE
                dlmwrite(filename,[mesh.edges.id(i) mesh.edges.centroid(:,i)' mesh.edges.nodes(:,i)' mesh.edges.faces(:,i)'],'-append');
            end
            for i=1:mesh.totF
                dlmwrite(filename,[mesh.faces.id(i) mesh.faces.centroid(:,i)' mesh.faces.nodes(:,i)' mesh.faces.edges(:,i)'],'-append');
            end
            fclose(nodefid);
            fclose(linefid);
            fclose(edgefid);
            fclose(facefid);
        elseif mesh.edgeflag
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            edgefid = fopen(edgesfilename,'w');
            for i=1:mesh.totN
                dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)' mesh.nodes.edges(:,i)'],'-append');
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim3free,'-append');
            for i=1:mesh.totE
                dlmwrite(filename,[mesh.edges.id(i) mesh.edges.centroid(:,i)' mesh.edges.nodes(:,i)'],'-append');
            end
            fclose(nodefid);
            fclose(linefid);
            fclose(edgefid);
        else
            nodefid = fopen(nodesfilename,'w');
            linefid = fopen(linesfilename,'w');
            fclose(nodefid);
            fclose(linefid);
            for i=1:mesh.totN
               dlmwrite(nodesfilename,[mesh.nodes.id(i) mesh.nodes.meshcoordinates(:,i)'],'-append'); 
            end
            dlmwrite(linesfilename,mesh.coordlines.dim1free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim2free,'-append');
            dlmwrite(linesfilename,mesh.coordlines.dim3free,'-append');
            fclose(nodefid);
            fclose(linefid);
        end       
end

return