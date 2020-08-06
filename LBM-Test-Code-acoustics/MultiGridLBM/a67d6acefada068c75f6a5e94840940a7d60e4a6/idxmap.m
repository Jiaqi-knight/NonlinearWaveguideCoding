% ***************************************************************************************************

% Routines to flatten 2/3/4D array indices to a single index

% 4D array index
function [idx]=idxmap (varargin) 
% nargin 表示参数个数
if nargin == 7
    i = varargin{1};
    j = varargin{2};
    k = varargin{3};
    v = varargin{4};
    j_max = varargin{5};
    k_max = varargin{6};
    v_max = varargin{7};
    
    % Loop over vel then k then j then i
	idx = v + ((k-1)*v_max) + ((j-1)*v_max*k_max) + ((i-1)*v_max*k_max*j_max);

elseif nargin==5
    % 3D array index
    i = varargin{1};
    j = varargin{2};
    k = varargin{3};
    j_max = varargin{4};
    k_max = varargin{5};

    % Loop over k then j then i
	idx = k + ((j-1)*k_max) + ((i-1)*k_max*j_max);
    
elseif nargin==3
    % 2D array index
    i = varargin{1};
    j = varargin{2};
    j_max = varargin{3};
    
	% Loop over j then i
	idx = j + ((i-1)*j_max);
end
end
