classdef Mesh3D < Geometry3D
% Parent class for 3D polygonal meshes.
%
%   Parent class for 3D polygonal meshes. 
%   Contains several static methods for creation of simples polygonal
%   meshes, and for reading from OFF file format.
%
%   Example
%     % Create 3D mesh representing an octahedron (6 vertices and 8 faces)
%     oct = Mesh3D.createOctahedron;
%     figure; draw(oct);
%
%     % Create 3D mesh representing an icosahedron (12 vertices and 20 faces)
%     ico = Mesh3D.createIcosahedron;
%     figure; draw(ico);
%
%     % Read triangular mesh in OFF format
%     mesh = Mesh3D.read_off('mushroom.off');
%     figure; draw(mesh, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
%
%   See also
%     TriMesh3D
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-02-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - BIA-BIBS.


%% Properties
properties
end % end properties


%% Constructor
methods (Access = protected)
    function obj = Mesh3D(varargin)
    % Constructor for Mesh3D class
    %   (will be called by subclasses)

    end

end % end constructors


%% Static factories
methods (Static)
    function mesh = createIcosahedron(varargin)
        % Create a new mesh representing an icosahedron.
        theta = 2*pi/5;
        l = 1/sin(theta/2)/2;
        z1 = sqrt(1-l*l);
        
        t1 = (0:2*pi/5:2*pi*(1-1/5))';
        x1 = l*cos(t1);
        y1 = l*sin(t1);
        
        t2 = t1 + 2*pi/10;
        x2 = l*cos(t2);
        y2 = l*sin(t2);
        
        h = sqrt(l*l-.5*.5);
        z2 = sqrt(3/4 - (l-h)*(l-h));
        
        nodes = [0 0 0;...
            [x1 y1 repmat(z1, [5 1])]; ...
            [x2 y2 repmat(z1+z2, [5 1])]; ...
            0 0 2*z1+z2];
        
        edges = [...
            1 2;1 3;1 4;1 5;1 6; ...
            2 3;3 4;4 5;5 6;6 2; ...
            2 7;7 3;3 8;8 4;4 9;9 5;5 10;10 6;6 11;11 2; ...
            7 8;8 9;9 10;10 11;11 7; ...
            7 12;8 12;9 12;10 12;11 12];
        
        % faces are ordered to have normals pointing outside of the mesh
        faces = [...
            1 3  2 ; 1 4  3 ; 1  5  4 ;  1  6  5 ;  1 2  6;...
            2 3  7 ; 3 4  8 ; 4  5  9 ;  5  6 10 ;  6 2 11;...
            7 3  8 ; 8 4  9 ; 9  5 10 ; 10  6 11 ; 11 2  7;...
            7 8 12 ; 8 9 12 ; 9 10 12 ; 10 11 12 ; 11 7 12];

        mesh = TriMesh3D(nodes, faces);
        mesh.Edges = edges;
    end
    
    function mesh = createOctahedron()
        % Create a new mesh representing an octahedron.
        nodes = [1 0 0;0 1 0;-1 0 0;0 -1 0;0 0 1;0 0 -1];
        edges = [1 2;1 4;1 5; 1 6;2 3;2 5;2 6;3 4;3 5;3 6;4 5;4 6];
        faces = [1 2 5;2 3 5;3 4 5;4 1 5;1 6 2;2 6 3;3 6 4;1 4 6];

        mesh = TriMesh3D(nodes, faces);
        mesh.Edges = edges;
    end
    
end % end methods


%% Serialization methods
methods

end

methods (Static)
    function mesh = read(fileName)
        % Read mesh for a file.
        
        [path, name, ext] = fileparts(fileName); %#ok<ASGLU>
        
        % check if off format
        if strcmpi(ext, '.off')
            mesh = read_off(fileName);
        else
            if exist('loadjson', 'file') == 0
                error('Requires the ''jsonlab'' library');
            end
            mesh = Geometry.read(fileName);
        end
    end
    
    function mesh = read_off(fileName)
        % Read a mesh stored in OFF format.
        %
        %   M = Mesh3D.read_off(FILENAME)
        %   Reads the data stored in file FILENAME and return the mesh into
        %   an instance of Mesh3D.
        %
        %   Example
        %     mesh = Mesh3D.read_off('mushroom.off');
        %     figure; 
        %     draw(mesh, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
        %     view([5 80]); light; lighting gouraud;
        
        % open file for reading
        f = fopen(fileName, 'r');
        if f == -1
            error('Geometry:Mesh3D:read_off:FileNotFound', ...
                ['Could not find file: ' fileName]);
        end
        
        % check format
        line = fgetl(f);   % -1 if eof
        if ~strcmp(line(1:3), 'OFF')
            error('Geometry:Mesh3D:read_off:FileFormatError', ...
                'Not a valid OFF file');
        end
        
        % number of faces and vertices
        line = fgetl(f);
        vals = sscanf(line, '%d %d');
        nVertices = vals(1);
        nFaces = vals(2);
        
        
        % Read vertex data
        [vertices, count] = fscanf(f, '%f ', [3 nVertices]);
        if count ~= nVertices * 3
            error('Geometry:Mesh3D:read_off:FileFormatError', ...
                ['Could not read all the ' num2str(nVertices) ' vertices']);
        end
        vertices = vertices';
        
        
        % Read Face data
        % First try to read faces as an homogeneous array. It if fails,
        % start from face offset and parse each face individually. In the
        % latter case, faces can have different number of vertices.
        
        % keep position of face info within file
        faceOffset = ftell(f);
        
        % read first face to assess number of vertices per face
        line = fgetl(f);
        if line == -1
            error('Geometry:Mesh3D:read_off:FileFormatError', ...
                'Unexpected end of file');
        end
        tokens = regexp(line, '\s+', 'split');
        face1 = str2double(tokens(2:end)) + 1;
        nv = length(face1);
        
        try
            % attenpt to read the remaining faces assuming they all have
            % the same number of vertices
            pattern = ['%d' repmat(' %d', 1, nv) '\n'];
            [faces, count] = fscanf(f, pattern, [(nv+1) (nFaces-1)]);
            if count ~= (nFaces-1) * (nv+1)
                error('Geometry:Mesh3D:read_off:FileFormatError', ...
                    'Could not read all the %d faces', nFaces);
            end
            
            % transpose, remove first column, use 1-indexing, and
            % concatenate with first face
            faces = [face1 ; faces(2:end,:)'+1];
            
            mesh = TriMesh3D(vertices, faces);
            
        catch
            % if attempt failed, switch to slower face-by-face parsing
            disp('read_off: Inhomogeneous number of vertices per face, switching to face-per-face parsing');
            
            fseek(f, faceOffset, 'bof');
            
            % allocate cell array
            faces = cell(1, nFaces);
            
            % iterate over faces
            for iFace = 1:nFaces
                % read next line
                line = fgetl(f);
                if line == -1
                    error('Geometry:Mesh3D:read_off:FileFormatError', ...
                        'Unexpected end of file');
                end
                
                % parse vertex indices for current face
                tokens = regexp(line, '\s+', 'split');
                faces{iFace} = str2double(tokens(2:end))' + 1;
            end
            
            disp('Non triangular meshes not yet implemented, returning stucture instead of class');
            mesh = struct('vertices', vertices, 'faces', faces);
        end
    end
end

end % end classdef

