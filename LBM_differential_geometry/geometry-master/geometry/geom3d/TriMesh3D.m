classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) TriMesh3D < Mesh3D
% Class for representing a 3D triangular mesh.
%
%   MESH = TriMesh3D(V, F)
%
%   Example
%   TriMesh3D
%
%   See also
%     Meshes3D
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-02-07,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.

properties
    % Coordinates of vertices, as a NV-by-3 array.
    Vertices;
    
    % Vertex indices for each edge, as a NE-by-2 array (optional).
    % Can be empty.
    Edges = [];
    
    % Vertex indices for each face, as a NF-by-3 array.
    Faces;
    
    % Mapping of faces associated to each edge (optional).
    % updated with method 'computeEdgeFaces'
    EdgeFaces = [];
end

%% Constructor
methods
    function obj = TriMesh3D(varargin)
        % Constructor for the TriMesh3D class.
        
        var1 = varargin{1};
        if isnumeric(var1)
            obj.Vertices = varargin{1};
            obj.Faces = varargin{2};
            
        elseif nargin == 1 && isa(var1, 'TriMesh3D')
            % Copy constructor from another TriMesh3D instance.
            obj.Vertices = var1.Vertices;
            obj.Edges = var1.Edges;
            obj.Faces = var1.Faces;
            obj.EdgeFaces = var1.EdgeFaces;

        elseif isstruct(var1)
            % Copy constructor from a structure.
            obj.Vertices = var1.vertices;
            obj.Faces = var1.faces;
        end
        
    end
end

%% Drawing functions
methods
    function h = drawFaceNormals(obj, varargin)
        % Draw the normal of each face.
        %
        % h = drawFaceNormals(mesh);
        pts = faceCentroids(obj);
        pos = pts.Coords;
        vn = faceNormals(obj);
        h = quiver3(pos(:, 1), pos(:, 2), pos(:, 3), ...
            vn(:, 1), vn(:, 2), vn(:, 3), 0, varargin{:});
    end
end

%% Global procesing of mesh
methods
    function res = smooth(obj, varargin)
        %SMOOTH Smooth a mesh.
        %
        % mesh2 = smooth(mesh, nIter);
        % Smoothes the mesh by replacing each vertex by the average of its
        % neighbors.
        
        % determine number of iterations
        nIter = 1;
        if ~isempty(varargin)
            nIter = varargin{1};
        end
        
        % compute adjacency matrix,
        % result is a Nv-by-Nv matrix with zeros on the diagonal
        adj = vertexAdjacencyMatrix(obj);
        
        % Add "self adjacencies"
        nv = size(obj.Vertices, 1);
        adj = adj + speye(nv);
        
        % weight each vertex by the number of its neighbors
        w = spdiags(full(sum(adj, 2).^(-1)), 0, nv, nv);
        adj = w * adj;
        
        % do averaging to smooth the field
        v2 = obj.Vertices;
        for k = 1:nIter
            v2 = adj * v2;
        end
        
        % return new TriMesh
        res = TriMesh3D(v2, obj.Faces);
    end
    
    function res = subdivide(obj, n)
        % Create a finer version of the mesh by subdividing each face.
        
        % compute the edge array
        computeEdges(obj);
        nEdges = size(obj.Edges, 1);
        
        % index of edges around each face
        faceEdgeIndices = meshFaceEdges(obj.Vertices, obj.Edges, obj.Faces);
        
        
        % Create new vertices on edges
        
        % several interpolated positions
        t = linspace(0, 1, n + 1)';
        coef2 = t(2:end-1);
        coef1 = 1 - t(2:end-1);
        
        % initialise the array of new vertices
        vertices2 = obj.Vertices;
        
        % keep an array containing index of new vertices for each original edge
        edgeNewVertexIndices = zeros(nEdges, n-1);
        
        % create new vertices on each edge
        for iEdge = 1:nEdges
            % extract each extremity as a point
            v1 = obj.Vertices(obj.Edges(iEdge, 1), :);
            v2 = obj.Vertices(obj.Edges(iEdge, 2), :);
            
            % compute new points
            newPoints = coef1 * v1 + coef2 * v2;
            
            % add new vertices, and keep their indices
            edgeNewVertexIndices(iEdge,:) = size(vertices2, 1) + (1:n-1);
            vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
        end
        
        
        % create array
        faces2 = zeros(0, 3);
        
        % iterate on faces of initial mesh
        nFaces = size(obj.Faces, 1);
        for iFace = 1:nFaces
            % compute index of each corner vertex
            face = obj.Faces(iFace, :);
            iv1 = face(1);
            iv2 = face(2);
            iv3 = face(3);
            
            % compute index of each edge
            faceEdges = faceEdgeIndices{iFace};
            ie1 = faceEdges(1);
            ie2 = faceEdges(2);
            ie3 = faceEdges(3);
            
            % indices of new vertices on edges
            edge1NewVertexIndices = edgeNewVertexIndices(ie1, :);
            edge2NewVertexIndices = edgeNewVertexIndices(ie2, :);
            edge3NewVertexIndices = edgeNewVertexIndices(ie3, :);
            
            % keep vertex 1 as reference for edges 1 and 3
            if obj.Edges(ie1, 1) ~= iv1
                edge1NewVertexIndices = edge1NewVertexIndices(end:-1:1);
            end
            if obj.Edges(ie3, 1) ~= iv1
                edge3NewVertexIndices = edge3NewVertexIndices(end:-1:1);
            end
            
            % create the first new face, on 'top' of the original face
            topVertexInds = [edge1NewVertexIndices(1) edge3NewVertexIndices(1)];
            newFace = [iv1 topVertexInds];
            faces2 = [faces2; newFace]; %#ok<AGROW>
            
            % iterate over middle strips
            for iStrip = 2:n-1
                % index of extreme vertices of current row
                ivr1 = edge1NewVertexIndices(iStrip);
                ivr2 = edge3NewVertexIndices(iStrip);
                
                % extreme vertices as points
                v1 = vertices2(ivr1, :);
                v2 = vertices2(ivr2, :);
                
                % create additional vertices within the bottom row of the strip
                t = linspace(0, 1, iStrip+1)';
                coef2 = t(2:end-1);
                coef1 = 1 - t(2:end-1);
                newPoints = coef1 * v1 + coef2 * v2;
                
                % compute indices of new vertices in result array
                newInds = size(vertices2, 1) + (1:iStrip-1);
                botVertexInds = [ivr1 newInds ivr2];
                
                % add new vertices
                vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
                
                % create top faces of current strip
                for k = 1:iStrip-1
                    newFace = [topVertexInds(k) botVertexInds(k+1) topVertexInds(k+1)];
                    faces2 = [faces2; newFace]; %#ok<AGROW>
                end
                
                % create bottom faces of current strip
                for k = 1:iStrip
                    newFace = [topVertexInds(k) botVertexInds(k) botVertexInds(k+1)];
                    faces2 = [faces2; newFace]; %#ok<AGROW>
                end
                
                % bottom vertices of current strip are top vertices of next strip
                topVertexInds = botVertexInds;
            end
            
            % for edge 2, keep vertex 2 of the current face as reference
            if obj.Edges(ie2, 1) ~= iv2
                edge2NewVertexIndices = edge2NewVertexIndices(end:-1:1);
            end
            
            % consider new vertices together with extremities
            botVertexInds = [iv2 edge2NewVertexIndices iv3];
            
            % create top faces for last strip
            for k = 1:n-1
                newFace = [topVertexInds(k) botVertexInds(k+1) topVertexInds(k+1)];
                faces2 = [faces2; newFace]; %#ok<AGROW>
            end
            
            % create bottom faces for last strip
            for k = 1:n
                newFace = [topVertexInds(k) botVertexInds(k) botVertexInds(k+1)];
                faces2 = [faces2; newFace]; %#ok<AGROW>
            end
        end

        % create the resulting data structure
        res = TriMesh3D(vertices2, faces2);
    end
end

%% Geometric information about mesh
methods
    function vol = volume(obj)
        % (signed) volume enclosed by this mesh.
        %
        % See Also
        %   surfaceArea

        % initialize an array of volume
        nFaces = size(obj.Faces, 1);
        vols = zeros(nFaces, 1);

        % Shift all vertices to the mesh centroid
        centroid = mean(obj.Vertices, 1);
        
        % compute volume of each tetraedron
        for iFace = 1:nFaces
            % consider the tetrahedron formed by face and mesh centroid
            tetra = obj.Vertices(obj.Faces(iFace, :), :);
            tetra = bsxfun(@minus, tetra, centroid);
            
            % volume of current tetrahedron
            vols(iFace) = det(tetra) / 6;
        end
        
        vol = sum(vols);
    end
    
    function area = surfaceArea(obj)
        % Surface area of this mesh, obtained by summing face areas.
        %
        % See Also
        %   volume
        
        % compute two direction vectors of each trinagular face, using the
        % first vertex of each face as origin
        v1 = obj.Vertices(obj.Faces(:, 2), :) - obj.Vertices(obj.Faces(:, 1), :);
        v2 = obj.Vertices(obj.Faces(:, 3), :) - obj.Vertices(obj.Faces(:, 1), :);
        
        % area of each triangle is half the cross product norm
        % see also crossProduct3d in MatGeom
        vn = zeros(size(v1));
        vn(:) = bsxfun(@times, v1(:,[2 3 1],:), v2(:,[3 1 2],:)) - ...
                bsxfun(@times, v2(:,[2 3 1],:), v1(:,[3 1 2],:));
        vn = sqrt(sum(vn .* vn, 2));
        
        % sum up and normalize
        area = sum(vn) / 2;
    end
    
%     function mb = meanBreadth(obj)
%         % Mean breadth of this mesh
%         % Mean breadth is proportionnal to the integral of mean curvature
%         %
%         % See Also
%         %   trimeshMeanBreadth
%         
%         mb = trimeshMeanBreadth(obj.Vertices, obj.Faces);
%     end
end


%% Vertex management methods
methods
    function nv = vertexNumber(obj)
        % Get the number of vertices in the mesh.
        nv = size(obj.Vertices, 1);
    end
    
    function verts = vertices(obj)
        % Return vertices in the mesh as a MultiPoint3D.
        verts = MultiPoint3D(obj.Vertices);
    end
    
    function adj = vertexAdjacencyMatrix(obj)
        % Get the adjacency matrix of mesh vertices.
        
        % forces faces to be floating point array, for sparse function
        faces = obj.Faces;
        if ~isfloat(faces)
            faces = double(obj.Faces);
        end
        
        % populate a sparse matrix
        adj = sparse(...
            [faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
            [faces(:,3); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1)], ...
            1.0);
        
        % remove double adjacencies
        adj = min(adj, 1);
        
        % ensure the size of the matrix is Nv-by-Nv
        % (this can happen if some vertices are not referenced)
        nv = size(obj.Vertices, 1);
        if size(adj, 1) < nv
            adj(nv, nv) = 0;
        end
    end
end

%% Edge management methods
methods
    function ne = edgeNumber(obj)
        % Get the number of edges in the mesh.
        
        % ne = edgeNumber(mesh)
        computeEdges(obj);
        ne = size(obj.Edges, 1);
    end
        
    function edgeList = edges(obj)
        % edgeList = edges(mesh);
        if isempty(obj.Edges)
            computeEdges(obj);
        end
        edgeList = obj.Edges;
    end
end

methods (Access = private)
    function computeEdges(obj)
        % Update the property Edges.
        
        % compute total number of edges
        % (3 edges per face)
        nFaces  = size(obj.Faces, 1);
        nEdges  = nFaces * 3;
        
        % create vertex indices for all edges (including duplicates)
        edges = zeros(nEdges, 2);
        for i = 1:nFaces
            f = obj.Faces(i, :);
            edges(((i-1)*3+1):i*3, :) = [f' f([2:end 1])'];
        end
        
        % remove duplicate edges, and sort the result
        obj.Edges = sortrows(unique(sort(edges, 2), 'rows'));
    end
    
    function edgeFaces = computeEdgeFaces(obj)
        % Update the property EdgeFaces.
        
        % ensure edge array is computed
        if isempty(obj.Edges)
            computeEdges(obj);
        end
        edges = obj.Edges;
        
        % allocate memory for result
        nEdges = size(obj.Edges, 1);
        obj.EdgeFaces = zeros(nEdges, 2);

        % iterate on faces
        nFaces = size(obj.Faces, 1);
        for iFace = 1:nFaces
            face = obj.Faces(iFace, :);
            
            % iterate on edges of current face
            for j = 1:length(face)
                % build edge: array of vertices
                j2 = mod(j, length(face)) + 1;
                
                % do not process edges with same vertices
                if face(j) == face(j2)
                    continue;
                end
                
                % vertex indices of current edge
                currentEdge = [face(j) face(j2)];
                
                % find index of current edge, assuming face is left-located
                b1 = ismember(obj.Edges, currentEdge, 'rows');
                indEdge = find(b1);
                if ~isempty(indEdge)
                    if obj.EdgeFaces(indEdge, 1) ~= 0
                        error('TriMesh3D:subdivide:IllegalTopology', ...
                            'Two faces were found on left side of edge %d ', indEdge);
                    end
                    
                    obj.EdgeFaces(indEdge, 1) = iFace;
                    continue;
                end
                
                % otherwise, assume the face is right-located
                b2 = ismember(edges, currentEdge([2 1]), 'rows');
                indEdge = find(b2);
                if ~isempty(indEdge)
                    if obj.EdgeFaces(indEdge, 2) ~= 0
                        error('TriMesh3D:subdivide:IllegalTopology', ...
                            'Two faces were found on left side of edge %d ', indEdge);
                    end
                    
                    obj.EdgeFaces(indEdge, 2) = iFace;
                    continue;
                end
                
                % If face was neither left nor right, error
                warning('TriMesh3D:subdivide:IllegalTopology', ...
                    'Edge %d of face %d was not found in edge array', ...
                    j, iFace);
                continue;
            end
        end
        
        edgeFaces = obj.EdgeFaces;
    end
end

%% Face management methods
methods
    function nf = faceNumber(obj)
        % Get the number of faces in the mesh.
        nf = size(obj.Faces, 1);
    end
    
    function normals = faceNormals(obj, inds)
        % Compute the normal vector to each face.
        %
        % vn = faceNormals(mesh);
        
        nf = size(obj.Faces, 1);
        if nargin == 1
            inds = 1:nf;
        end

        % compute vector of each edge
        v1 = obj.Vertices(obj.Faces(inds,2),:) - obj.Vertices(obj.Faces(inds,1),:);
        v2 = obj.Vertices(obj.Faces(inds,3),:) - obj.Vertices(obj.Faces(inds,1),:);

        % compute normals using cross product (vectors have same size)
        normals = cross(v1, v2, 2);
    end
    
    function pts = faceCentroids(obj, inds)
        % Compute the centroid of each face.
        %
        % pts = faceCentroids(mesh);
        nf = size(obj.Faces, 1);
        if nargin == 1
            inds = 1:nf;
        end
        pts = zeros(length(inds), 3);
        
        for i = 1:3
            pts = pts + obj.Vertices(obj.Faces(inds,i),:) / 3;
        end
        
        pts = MultiPoint3D(pts);
    end

    function poly = facePolygon(obj, ind)
        poly = obj.Vertices(obj.Faces(ind, :), :);
    end
end


%% Methods implementing Geometry3D
methods
    function box = boundingBox(obj)
        % Return the bounding box of this mesh.
        mini = min(obj.Vertices);
        maxi = max(obj.Vertices);
        box = Box3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current mesh, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});

        % add default drawing options
        options = {'FaceColor', [.7 .7 .7]};

        % extract optional drawing options
        if length(varargin) > 1 && ischar(varargin{1})
            options = [options varargin];
        end
        
        h = patch('Parent', ax, ...
            'vertices', obj.Vertices, 'faces', obj.Faces, ...
            options{:} );

        % optionnally add style processing
        if ~isempty(style)
            apply(style, hh);
        end
                
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transfo)
        % Apply a transform to this mesh.
        vt = transformPoint(transfo, obj.Vertices);
        res = TriMesh3D(vt, obj.Faces);
        res.Edges = obj.Edges;
        res.EdgeFaces = obj.EdgeFaces;
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this mesh.
        factor = varargin{1};
        res = TriMesh3D(obj.Vertices * factor, obj.Faces);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this mesh.
        shift = varargin{1};
        res = TriMesh3D(bsxfun(@plus, obj.Vertices, shift), obj.Faces);
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'TriMesh3D', ...
            'Vertices', obj.Vertices, ...
            'Faces', obj.Faces);
    end
end
methods (Static)
    function mesh = fromStruct(str)
        % Create a new instance from a structure.
        if ~(isfield(str, 'Vertices') && isfield(str, 'Faces'))
            error('Requires fields vertices and faces');
        end
        if size(str.Faces, 2) ~= 3
            error('Requires a triangular face array');
        end
        mesh = TriMesh3D(str);
    end
end

end