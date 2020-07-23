function[mesh]=load_csvmesh(filename,printflag)

basefilename = strcat([filename,'_basedata','.csv']);
nodesfilename = strcat([filename,'_nodedata','.csv']);
linesfilename = strcat([filename,'_linedata','.csv']);
edgesfilename = strcat([filename,'_edgedata','.csv']);
facesfilename = strcat([filename,'_facedata','.csv']);
cellsfilename = strcat([filename,'_celldata','.csv']);

R = 0;
C = 0;
R2 = 17 - 1;
C2 = 0;
RNG = [R C R2 C2];
basedata = csvread(basefilename,R,C,RNG);

mesh.D = basedata(1);
mesh.dim1min = basedata(2);
mesh.dim1max = basedata(3);
mesh.Ndim1 = basedata(4);
mesh.dim2min = basedata(5);
mesh.dim2max = basedata(6);
mesh.Ndim2 = basedata(7);
mesh.dim3min = basedata(8);
mesh.dim3max = basedata(9);
mesh.Ndim3 = basedata(10);
mesh.Nlinesdim1 = basedata(11);
mesh.Nlinesdim2 = basedata(12);
mesh.Nlinesdim3 = basedata(13);
mesh.totN = basedata(14);
mesh.totE = basedata(15);
mesh.totF = basedata(16);
mesh.totC = basedata(17);
mesh.edgeflag = basedata(18);
mesh.faceflag = basedata(19);
mesh.cellflag = basedata(20);

switch mesh.D
%-----------------------------------------------------------------------------1D--------------------------------------------------------------------------------------------
    case 1
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('LOAD MESH FROM .csv FILE FORMAT')
            disp(' ')
            disp('1D mesh')
            fprintf('Start point x0 = %d\n',mesh.dim1min)
            fprintf('End point xN = %6.2f\n',mesh.dim1max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        mesh.nodes.id = zeros(1,mesh.totN);
        mesh.nodes.meshcoordinates = zeros(3,mesh.totN);
        mesh.coordlines.dim1free = zeros(mesh.Nlinesdim1,3);
        R = 0;
        C = 0;
        R2 = mesh.Nlinesdim1 - 1;
        C2 = 3 - 1;
        RNG = [R C R2 C2];
        data = csvread(linesfilename,R,C,RNG);
        mesh.coordlines.dim1free = data;
        if mesh.edgeflag
            mesh.nodesperedge = 2;
            mesh.edgespernode = 2;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 6 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
                mesh.nodes.edges(1,i) = data(5);
                mesh.nodes.edges(2,i) = data(6);
            end
            for i=1:mesh.totE
                R = 0;
                C = 0;
                C2 = 6 - 1;
                RNG = [R C R C2];
                data = csvread(edgesfilename,R,C,RNG);
                mesh.edges.id(i) = data(1);
                mesh.edges.centroid(1,i) = data(2);
                mesh.edges.centroid(2,i) = data(3);
                mesh.edges.centroid(3,i) = data(4);
                mesh.edges.nodes(1,i) = data(5);
                mesh.edges.nodes(2,i) = data(6);
            end
        else
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 4 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
            end
        end
%-----------------------------------------------------------------------------2D--------------------------------------------------------------------------------------------
    case 2
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('LOAD MESH FROM .csv FILE FORMAT')
            disp(' ')
            disp('2D mesh')
            fprintf('Lower boundary along the first direction x0 = %d\n',mesh.dim1min)
            fprintf('Upper boundary along the first direction xN = %6.2f\n',mesh.dim1max)
            fprintf('Lower boundary along the second direction y0 = %d\n',mesh.dim2min)
            fprintf('Upper boundary along the second direction yN = %6.2f\n',mesh.dim2max)
            disp(' ')
            fprintf('Total number of nodes N = %6.2f\n',mesh.totN)
        end
        mesh.nodes.id = zeros(1,mesh.totN);
        mesh.nodes.meshcoordinates = zeros(3,mesh.totN);
        mesh.coordlines.dim1free = zeros(mesh.Nlinesdim1,3*mesh.Ndim2);
        mesh.coordlines.dim2free = zeros(mesh.Nlinesdim2,3*mesh.Ndim1);
        R = 0;
        C = 0;
        R2 = mesh.Nlinesdim1 - 1;
        C2 = 3*mesh.Ndim2 - 1;
        RNG = [R C R2 C2];
        data = csvread(linesfilename,R,C,RNG);
        mesh.coordlines.dim1free = data;
        R = mesh.Nlinesdim1+1 - 1;
        C = 0;
        R2 = mesh.Nlinesdim1+mesh.Nlinesdim2 - 1;
        C2 = 3*mesh.Ndim1 - 1;
        RNG = [R C R2 C2];
        data = csvread(linesfilename,R,C,RNG);
        mesh.coordlines.dim2free = data;
        if mesh.faceflag && mesh.edgeflag
            mesh.nodesperedge = 2;
            mesh.nodesperface = 4;
            mesh.edgespernode = 4;
            mesh.edgesperface = 4;
            mesh.facespernode = 4;
            mesh.facesperedge = 2;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.nodes.faces = zeros(mesh.facespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            mesh.edges.faces = zeros(mesh.facesperedge,mesh.totE);
            mesh.faces.id = zeros(1,mesh.totF);
            mesh.faces.centroid = zeros(3,mesh.totF);
            mesh.faces.nodes = zeros(mesh.nodesperface,mesh.totF);
            mesh.faces.edges = zeros(mesh.edgesperface,mesh.totF);
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 12 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
                mesh.nodes.edges(1,i) = data(5);
                mesh.nodes.edges(2,i) = data(6);
                mesh.nodes.edges(3,i) = data(7);
                mesh.nodes.edges(4,i) = data(8);
                mesh.nodes.faces(1,i) = data(9);
                mesh.nodes.faces(2,i) = data(10);
                mesh.nodes.faces(3,i) = data(11);
                mesh.nodes.faces(4,i) = data(12);
            end
            for i=1:mesh.totE
                R = 0;
                C = 0;
                C2 = 8 - 1;
                RNG = [R C R C2];
                data = csvread(edgesfilename,R,C,RNG);
                mesh.edges.id(i) = data(1);
                mesh.edges.centroid(1,i) = data(2);
                mesh.edges.centroid(2,i) = data(3);
                mesh.edges.centroid(3,i) = data(4);
                mesh.edges.nodes(1,i) = data(5);
                mesh.edges.nodes(2,i) = data(6);
                mesh.edges.faces(1,i) = data(7);
                mesh.edges.faces(2,i) = data(8);
            end
            for i=1:mesh.totF
                R = 0;
                C = 0;
                C2 = 12 - 1;
                RNG = [R C R C2];
                data = csvread(facesfilename,R,C,RNG);
                mesh.faces.id(i) = data(1);
                mesh.faces.centroid(1,i) = data(2);
                mesh.faces.centroid(2,i) = data(3);
                mesh.faces.centroid(3,i) = data(4);
                mesh.faces.nodes(1,i) = data(5);
                mesh.faces.nodes(2,i) = data(6);
                mesh.faces.nodes(3,i) = data(7);
                mesh.faces.nodes(4,i) = data(8);
                mesh.faces.edges(1,i) = data(9);
                mesh.faces.edges(2,i) = data(10);
                mesh.faces.edges(3,i) = data(11);
                mesh.faces.edges(4,i) = data(12);
            end
        elseif mesh.edgeflag
            mesh.nodesperedge = 2;
            mesh.edgespernode = 4;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 8 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
                mesh.nodes.edges(1,i) = data(5);
                mesh.nodes.edges(2,i) = data(6);
                mesh.nodes.edges(3,i) = data(7);
                mesh.nodes.edges(4,i) = data(8);
            end
            for i=1:mesh.totE
                R = 0;
                C = 0;
                C2 = 6 - 1;
                RNG = [R C R C2];
                data = csvread(edgesfilename,R,C,RNG);
                mesh.edges.id(i) = data(1);
                mesh.edges.centroid(1,i) = data(2);
                mesh.edges.centroid(2,i) = data(3);
                mesh.edges.centroid(3,i) = data(4);
                mesh.edges.nodes(1,i) = data(5);
                mesh.edges.nodes(2,i) = data(6);
            end
        else
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 4 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
            end
        end
%-----------------------------------------------------------------------------3D--------------------------------------------------------------------------------------------
    case 3
        c = clock;
        if printflag
            disp(strcat([num2str(c(1)),'_',num2str(c(2)),'_',num2str(c(3))]))
            disp(' ')
            disp('LOAD MESH FROM .csv FILE FORMAT')
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
        mesh.nodes.id = zeros(1,mesh.totN);
        mesh.nodes.meshcoordinates = zeros(3,mesh.totN);
        mesh.coordlines.dim1free = zeros(mesh.Nlinesdim1,3*mesh.Ndim2*mesh.Ndim3);
        mesh.coordlines.dim2free = zeros(mesh.Nlinesdim2,3*mesh.Ndim1*mesh.Ndim3);
        mesh.coordlines.dim3free = zeros(mesh.Nlinesdim3,3*mesh.Ndim1*mesh.Ndim2);
        R = 0;
        C = 0;
        R2 = mesh.Nlinesdim1 - 1;
        C2 = 3*mesh.Ndim2*mesh.Ndim3 - 1;
        RNG = [R C R2 C2];
        data = csvread(linesfilename,R,C,RNG);
        mesh.coordlines.dim1free = data;
        R = mesh.Nlinesdim1+1 - 1;
        C = 0;
        R2 = mesh.Nlinesdim1+mesh.Nlinesdim2 - 1;
        C2 = 3*mesh.Ndim1*mesh.Ndim3 - 1;
        RNG = [R C R2 C2];
        data = csvread(linesfilename,R,C,RNG);
        mesh.coordlines.dim2free = data;
        R = mesh.Nlinesdim1+mesh.Nlinesdim2+1 - 1;
        C = 0;
        R2 = mesh.Nlinesdim1+mesh.Nlinesdim2+mesh.Nlinesdim3 - 1;
        C2 = 3*mesh.Ndim1*mesh.Ndim2 - 1;
        RNG = [R C R2 C2];
        data = csvread(linesfilename,R,C,RNG);
        mesh.coordlines.dim3free = data;
        if mesh.cellflag && mesh.faceflag && mesh.edgeflag
            mesh.nodesperedge = 2;
            mesh.nodesperface = 4;
            mesh.nodespercell = 8;
            mesh.edgespernode = 6;
            mesh.edgesperface = 4;
            mesh.edgespercell = 12;
            mesh.facespernode = 12;
            mesh.facesperedge = 4;
            mesh.facespercell = 6;
            mesh.cellspernode = 8;
            mesh.cellsperedge = 4;
            mesh.cellsperface = 2;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.nodes.faces = zeros(mesh.facespernode,mesh.totN);
            mesh.nodes.cells = zeros(mesh.cellspernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            mesh.edges.faces = zeros(mesh.facesperedge,mesh.totE);
            mesh.edges.cells = zeros(mesh.cellsperedge,mesh.totE);
            mesh.faces.id = zeros(1,mesh.totF);
            mesh.faces.centroid = zeros(3,mesh.totF);
            mesh.faces.nodes = zeros(mesh.nodesperface,mesh.totF);
            mesh.faces.edges = zeros(mesh.edgesperface,mesh.totF);
            mesh.faces.cells = zeros(mesh.cellsperface,mesh.totF);
            mesh.cells.id = zeros(1,mesh.totC);
            mesh.cells.centroid = zeros(3,mesh.totC);
            mesh.cells.nodes = zeros(mesh.nodespercell,mesh.totC);
            mesh.cells.edges = zeros(mesh.edgespercell,mesh.totC);
            mesh.cells.faces = zeros(mesh.facespercell,mesh.totC);
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 30 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
                mesh.nodes.edges(1,i) = data(5);
                mesh.nodes.edges(2,i) = data(6);
                mesh.nodes.edges(3,i) = data(7);
                mesh.nodes.edges(4,i) = data(8);
                mesh.nodes.edges(5,i) = data(9);
                mesh.nodes.edges(6,i) = data(10);
                mesh.nodes.faces(1,i) = data(11);
                mesh.nodes.faces(2,i) = data(12);
                mesh.nodes.faces(3,i) = data(13);
                mesh.nodes.faces(4,i) = data(14);
                mesh.nodes.faces(5,i) = data(15);
                mesh.nodes.faces(6,i) = data(16);
                mesh.nodes.faces(7,i) = data(17);
                mesh.nodes.faces(8,i) = data(18);
                mesh.nodes.faces(9,i) = data(19);
                mesh.nodes.faces(10,i) = data(20);
                mesh.nodes.faces(11,i) = data(21);
                mesh.nodes.faces(12,i) = data(22);
                mesh.nodes.cells(1,i) = data(23);
                mesh.nodes.cells(2,i) = data(24);
                mesh.nodes.cells(3,i) = data(25);
                mesh.nodes.cells(4,i) = data(26);
                mesh.nodes.cells(5,i) = data(27);
                mesh.nodes.cells(6,i) = data(28);
                mesh.nodes.cells(7,i) = data(29);
                mesh.nodes.cells(8,i) = data(30);
            end
            for i=1:mesh.totE
                R = 0;
                C = 0;
                C2 = 14 - 1;
                RNG = [R C R C2];
                data = csvread(edgesfilename,R,C,RNG);
                mesh.edges.id(i) = data(1);
                mesh.edges.centroid(1,i) = data(2);
                mesh.edges.centroid(2,i) = data(3);
                mesh.edges.centroid(3,i) = data(4);
                mesh.edges.nodes(1,i) = data(5);
                mesh.edges.nodes(2,i) = data(6);
                mesh.edges.faces(1,i) = data(7);
                mesh.edges.faces(2,i) = data(8);
                mesh.edges.faces(3,i) = data(9);
                mesh.edges.faces(4,i) = data(10);
                mesh.edges.cells(1,i) = data(11);
                mesh.edges.cells(2,i) = data(12);
                mesh.edges.cells(3,i) = data(13);
                mesh.edges.cells(4,i) = data(14);
            end
            for i=1:mesh.totF
                R = 0;
                C = 0;
                C2 = 14 - 1;
                RNG = [R C R C2];
                data = csvread(facesfilename,R,C,RNG);
                mesh.faces.id(i) = data(1);
                mesh.faces.centroid(1,i) = data(2);
                mesh.faces.centroid(2,i) = data(3);
                mesh.faces.centroid(3,i) = data(4);
                mesh.faces.nodes(1,i) = data(5);
                mesh.faces.nodes(2,i) = data(6);
                mesh.faces.nodes(3,i) = data(7);
                mesh.faces.nodes(4,i) = data(8);
                mesh.faces.edges(1,i) = data(9);
                mesh.faces.edges(2,i) = data(10);
                mesh.faces.edges(3,i) = data(11);
                mesh.faces.edges(4,i) = data(12);
                mesh.faces.cells(1,i) = data(13);
                mesh.faces.cells(2,i) = data(14);
            end
            for i=1:mesh.totC
                R = 0;
                C = 0;
                C2 = 30 - 1;
                RNG = [R C R C2];
                data = csvread(cellsfilename,R,C,RNG);
                mesh.cells.id(i) = data(1);
                mesh.cells.centroid(1,i) = data(2);
                mesh.cells.centroid(2,i) = data(3);
                mesh.cells.centroid(3,i) = data(4);
                mesh.cells.nodes(1,i) = data(5);
                mesh.cells.nodes(2,i) = data(6);
                mesh.cells.nodes(3,i) = data(7);
                mesh.cells.nodes(4,i) = data(8);
                mesh.cells.nodes(5,i) = data(9);
                mesh.cells.nodes(6,i) = data(10);
                mesh.cells.nodes(7,i) = data(11);
                mesh.cells.nodes(8,i) = data(12);
                mesh.cells.edges(1,i) = data(13);
                mesh.cells.edges(2,i) = data(14);
                mesh.cells.edges(3,i) = data(15);
                mesh.cells.edges(4,i) = data(16);
                mesh.cells.edges(5,i) = data(17);
                mesh.cells.edges(6,i) = data(18);
                mesh.cells.edges(7,i) = data(19);
                mesh.cells.edges(8,i) = data(20);
                mesh.cells.edges(9,i) = data(21);
                mesh.cells.edges(10,i) = data(22);
                mesh.cells.edges(11,i) = data(23);
                mesh.cells.edges(12,i) = data(24);
                mesh.cells.faces(1,i) = data(25);
                mesh.cells.faces(2,i) = data(26);
                mesh.cells.faces(3,i) = data(27);
                mesh.cells.faces(4,i) = data(28);
                mesh.cells.faces(5,i) = data(29);
                mesh.cells.faces(6,i) = data(30);
            end
        elseif mesh.faceflag && mesh.edgeflag
            mesh.nodesperedge = 2;
            mesh.nodesperface = 4;
            mesh.edgespernode = 6;
            mesh.edgesperface = 4;
            mesh.facespernode = 12;
            mesh.facesperedge = 4;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.nodes.faces = zeros(mesh.facespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            mesh.edges.faces = zeros(mesh.facesperedge,mesh.totE);
            mesh.faces.id = zeros(1,mesh.totF);
            mesh.faces.centroid = zeros(3,mesh.totF);
            mesh.faces.nodes = zeros(mesh.nodesperface,mesh.totF);
            mesh.faces.edges = zeros(mesh.edgesperface,mesh.totF);
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 22 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
                mesh.nodes.edges(1,i) = data(5);
                mesh.nodes.edges(2,i) = data(6);
                mesh.nodes.edges(3,i) = data(7);
                mesh.nodes.edges(4,i) = data(8);
                mesh.nodes.edges(5,i) = data(9);
                mesh.nodes.edges(6,i) = data(10);
                mesh.nodes.faces(1,i) = data(11);
                mesh.nodes.faces(2,i) = data(12);
                mesh.nodes.faces(3,i) = data(13);
                mesh.nodes.faces(4,i) = data(14);
                mesh.nodes.faces(5,i) = data(15);
                mesh.nodes.faces(6,i) = data(16);
                mesh.nodes.faces(7,i) = data(17);
                mesh.nodes.faces(8,i) = data(18);
                mesh.nodes.faces(9,i) = data(19);
                mesh.nodes.faces(10,i) = data(20);
                mesh.nodes.faces(11,i) = data(21);
                mesh.nodes.faces(12,i) = data(22);
            end
            for i=1:mesh.totE
                R = 0;
                C = 0;
                C2 = 10 - 1;
                RNG = [R C R C2];
                data = csvread(edgesfilename,R,C,RNG);
                mesh.edges.id(i) = data(1);
                mesh.edges.centroid(1,i) = data(2);
                mesh.edges.centroid(2,i) = data(3);
                mesh.edges.centroid(3,i) = data(4);
                mesh.edges.nodes(1,i) = data(5);
                mesh.edges.nodes(2,i) = data(6);
                mesh.edges.faces(1,i) = data(7);
                mesh.edges.faces(2,i) = data(8);
                mesh.edges.faces(3,i) = data(9);
                mesh.edges.faces(4,i) = data(10);
            end
            for i=1:mesh.totF
                R = 0;
                C = 0;
                C2 = 12 - 1;
                RNG = [R C R C2];
                data = csvread(facesfilename,R,C,RNG);
                mesh.faces.id(i) = data(1);
                mesh.faces.centroid(1,i) = data(2);
                mesh.faces.centroid(2,i) = data(3);
                mesh.faces.centroid(3,i) = data(4);
                mesh.faces.nodes(1,i) = data(5);
                mesh.faces.nodes(2,i) = data(6);
                mesh.faces.nodes(3,i) = data(7);
                mesh.faces.nodes(4,i) = data(8);
                mesh.faces.edges(1,i) = data(9);
                mesh.faces.edges(2,i) = data(10);
                mesh.faces.edges(3,i) = data(11);
                mesh.faces.edges(4,i) = data(12);
            end
        elseif mesh.edgeflag
            mesh.nodesperedge = 2;
            mesh.edgespernode = 6;
            mesh.nodes.edges = zeros(mesh.edgespernode,mesh.totN);
            mesh.edges.id = zeros(1,mesh.totE);
            mesh.edges.centroid = zeros(3,mesh.totE);
            mesh.edges.nodes = zeros(mesh.nodesperedge,mesh.totE);
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 10 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
                mesh.nodes.edges(1,i) = data(5);
                mesh.nodes.edges(2,i) = data(6);
                mesh.nodes.edges(3,i) = data(7);
                mesh.nodes.edges(4,i) = data(8);
                mesh.nodes.edges(5,i) = data(9);
                mesh.nodes.edges(6,i) = data(10);
            end
            for i=1:mesh.totE
                R = 0;
                C = 0;
                C2 = 6 - 1;
                RNG = [R C R C2];
                data = csvread(edgesfilename,R,C,RNG);
                mesh.edges.id(i) = data(1);
                mesh.edges.centroid(1,i) = data(2);
                mesh.edges.centroid(2,i) = data(3);
                mesh.edges.centroid(3,i) = data(4);
                mesh.edges.nodes(1,i) = data(5);
                mesh.edges.nodes(2,i) = data(6);
            end
        else
            for i=1:mesh.totN
                R = 0;
                C = 0;
                C2 = 4 - 1;
                RNG = [R C R C2];
                data = csvread(nodesfilename,R,C,RNG);
                mesh.nodes.id(i) = data(1);
                mesh.nodes.meshcoordinates(1,i) = data(2);
                mesh.nodes.meshcoordinates(2,i) = data(3);
                mesh.nodes.meshcoordinates(3,i) = data(4);
            end
        end       
end

return