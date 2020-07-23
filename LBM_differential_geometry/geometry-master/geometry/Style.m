classdef Style < handle
% Encapsulates information for drawing shapes.
%
%   Class Style
%   Contains the information for drawing a shape.
%
%   The different fields in Style are:
%     MarkerColor     = 'b';
%     MarkerStyle     = '+';
%     MarkerSize      = 6;
%     MarkerFillColor = 'none';
%     MarkerVisible   = false;
%     LineColor       = 'b';
%     LineWidth       = .5;
%     LineStyle       = '-';
%     LineVisible     = true;
%     FillColor       = 'y';
%     FillOpacity     = 1;
%     FillVisible     = false;
%
%
%   Example
%     % draw a polygon with default style
%     poly = [10 10;20 10;20 20;10 20];
%     figure; h1 = drawPolygon(poly, 'b');
%     axis equal; axis([0 50 0 50]);
%     % change style using Style class
%     style1 = Style('LineWidth', 2, 'LineColor', 'g');
%     apply(style1, h1)
%     % chage vertex style
%     hold on;
%     h2 = drawPolygon(poly);
%     style2 = Style('MarkerStyle', 's', 'MarkerColor', 'k', 'MarkerFillColor', 'w', 'MarkerVisible', true, 'LineVisible', false);
%     apply(style2, h2)
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - BIA-BIBS.



%% Properties
properties
    % global visibility of the shape
    Visible = true;
    
    % style for the markers / vertices
    MarkerColor     = 'b';
    MarkerStyle     = '+';
    MarkerSize      = 6;
    MarkerFillColor = 'none';
    MarkerVisible   = false;
    
    % style for the lines / curves / edges
    LineColor       = 'b';
    LineWidth       = .5;
    LineStyle       = '-';
    LineVisible     = true;
    
    % style for filling interior of polygons
    FillColor       = 'y';
    FillOpacity     = 1;
    FillVisible     = false;
    
    % style for polygonal surfaces
    FaceColor       = [.75 .75 .75];
    FaceOpacity     = 1;
    FaceVisible     = false;
    
end % end properties


%% Constructor
methods
    function obj = Style(varargin)
    % Constructor for Style class.

        if nargin == 0
            return;
        end
        
        var1 = varargin{1};
        if isa(var1, 'Style')
            copyFields(var1);
            return;
        end
        
        while length(varargin) >= 2
            name = varargin{1};
            value = varargin{2};
            
            if strcmpi(name, 'Visible')
                obj.Visible = value;
    
            elseif strcmpi(name, 'MarkerColor')
                obj.MarkerColor = value;
            elseif strcmpi(name, 'MarkerStyle')
                obj.MarkerStyle = value;
            elseif strcmpi(name, 'MarkerSize')
                obj.MarkerSize = value;
            elseif strcmpi(name, 'MarkerFillColor')
                obj.MarkerFillColor = value;
            elseif strcmpi(name, 'MarkerVisible')
                obj.MarkerVisible= value;
                
            elseif strcmpi(name, 'LineColor')
                obj.LineColor = value;
            elseif strcmpi(name, 'LineWidth')
                obj.LineWidth = value;
            elseif strcmpi(name, 'LineStyle')
                obj.LineStyle = value;
            elseif strcmpi(name, 'LineVisible')
                obj.LineVisible= value;
                
            elseif strcmpi(name, 'FillColor')
                obj.FillColor = value;
            elseif strcmpi(name, 'FillOpacity')
                obj.FillOpacity = value;
            elseif strcmpi(name, 'FillVisible')
                obj.FillVisible= value;
                
            elseif strcmpi(name, 'FaceColor')
                obj.FaceColor = value;
            elseif strcmpi(name, 'FaceOpacity')
                obj.FaceOpacity = value;
            elseif strcmpi(name, 'FaceVisible')
                obj.FaceVisible= value;
            end
            
            varargin(1:2) = [];
        end
        
        function copyFields(that)
            obj.Visible            = that.Visible;
            
            obj.MarkerColor        = that.MarkerColor;
            obj.MarkerStyle        = that.MarkerStyle;
            obj.MarkerSize         = that.MarkerSize;
            obj.MarkerFillColor    = that.MarkerFillColor;
            obj.MarkerVisible      = that.MarkerVisible;
            
            obj.LineColor          = that.LineColor;
            obj.LineWidth          = that.LineWidth;
            obj.LineStyle          = that.LineStyle;
            obj.LineVisible        = that.LineVisible;
            
            obj.FillColor          = that.FillColor;
            obj.FillOpacity        = that.FillOpacity;
            obj.FillVisible        = that.FillVisible;

            obj.FaceColor          = that.FaceColor;
            obj.FaceOpacity        = that.FaceOpacity;
            obj.FaceVisible        = that.FaceVisible;

        end
    end

end % end constructors


%% Methods
methods
    function apply(obj, h)
        % Apply the style to the given graphic handle(s).
        
        hType = get(h, 'Type');

        % setup marker style
        if obj.Visible && obj.MarkerVisible && ~strcmp(hType, 'patch')
            set(h, 'MarkerEdgeColor',   obj.MarkerColor);
            set(h, 'Marker',            obj.MarkerStyle);
            set(h, 'MarkerSize',        obj.MarkerSize);
            set(h, 'MarkerFaceColor',   obj.MarkerFillColor);            
            set(h, 'LineWidth',         obj.LineWidth);
        else
            set(h, 'Marker', 'none');
        end
        
        % setup line style
        if obj.Visible && obj.LineVisible && ~strcmp(hType, 'patch')
            set(h, 'LineStyle',         obj.LineStyle);
            set(h, 'Color',             obj.LineColor);
            set(h, 'LineWidth',         obj.LineWidth);
        else
            set(h, 'LineStyle', 'none');
        end
        
        % setup fill style
        if obj.Visible && obj.FillVisible && strcmp(hType, 'patch')
            set(h, 'FaceColor', obj.FillColor);
            set(h, 'FaceAlpha', obj.FillOpacity);
        end
        
        % setup face style
        if obj.Visible && obj.FaceVisible && strcmp(hType, 'patch')
            set(h, 'FaceColor', obj.FaceColor);
            set(h, 'FaceAlpha', obj.FaceOpacity);
        end
    end
    
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.

        % create empty struct
        str = struct();

        % global visibility
        if ~obj.Visible
            str.Visible = obj.Visible;
        end

        % update marker modifiers with values different from default
        if obj.MarkerVisible ~= false
            str.MarkerVisible = obj.MarkerVisible;
        end
        if obj.MarkerColor ~= 'b'
            str.MarkerColor = obj.MarkerColor;
        end
        if obj.MarkerStyle ~= '+'
            str.MarkerStyle = obj.MarkerStyle;
        end
        if obj.MarkerSize ~= 6
            str.MarkerSize = obj.MarkerSize;
        end
        if ~ischar(obj.MarkerFillColor) || ~strcmp(obj.MarkerFillColor, 'none')
            str.MarkerFillColor = obj.MarkerFillColor;
        end
        
        % update line modifiers with values different from default
        if obj.LineVisible ~= true
            str.LineVisible = obj.LineVisible;
        end
        if ~isSameColor(obj.LineColor, 'b')
            str.LineColor = obj.LineColor;
        end
        if obj.LineWidth ~= .5
            str.LineWidth = obj.LineWidth;
        end
        if ~strcmp(obj.LineStyle, '-')
            str.LineStyle = obj.LineStyle;
        end
        
        % update fill modifiers with values different from default
        if obj.FillVisible
            str.FillVisible = obj.FillVisible;
        end
        if ~isSameColor(obj.FillColor, 'y')
            str.FillColor = obj.FillColor;
        end
        if obj.FillOpacity ~= 1
            str.FillOpacity = obj.FillOpacity;
        end
        
        % update face modifiers with values different from default
        if ~obj.FaceVisible
            str.FaceVisible = obj.FaceVisible;
        end
        if ~isSameColor(obj.FaceColor, [.7 .7 .7])
            str.FaceColor = obj.FaceColor;
        end
        if obj.FaceOpacity ~= 1
            str.FaceOpacity = obj.FaceOpacity;
        end
        
        function b = isSameColor(color1, color2)
            if ischar(color1)
                color1 = colorFromName(color1);
            end
            if ischar(color2)
                color2 = colorFromName(color2);
            end
            b = all(color1 == color2);
        end
        
        function color = colorFromName(name)
            switch name(1)
                case 'k', color = [0 0 0];
                case 'r', color = [1 0 0];
                case 'g', color = [0 1 0];
                case 'b', color = [0 0 1];
                case 'y', color = [1 1 0];
                case 'm', color = [1 0 1];
                case 'c', color = [0 1 1];
                case 'w', color = [1 1 1];
                otherwise
                    error('unknown color name: %s', name);
            end
        end
        
    end
    
    function write(obj, fileName, varargin)
        % Write into a JSON file.
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
end

methods (Static)
    function style = fromStruct(str)
        % Create a new instance from a structure.
        
        % create default empty style
        style = Style();
        
        % global visibility
        if isfield(str, 'Visible')
            style.Visible = str.Visible;
        elseif isfield(str, 'visible')
            style.Visible = str.visible;
        end

        % parse marker style modifiers
        if isfield(str, 'MarkerVisible')
            style.MarkerVisible = str.MarkerVisible;
        elseif isfield(str, 'markerVisible')
            style.MarkerVisible = str.markerVisible;
        end
        if isfield(str, 'MarkerColor')
            style.MarkerColor = str.MarkerColor;
        elseif isfield(str, 'markerColor')
            style.MarkerColor = str.markerColor;
        end
        if isfield(str, 'MarkerStyle')
            style.MarkerStyle = str.MarkerStyle;
        elseif isfield(str, 'markerStyle')
            style.MarkerStyle = str.markerStyle;
        end
        if isfield(str, 'MarkerSize')
            style.MarkerSize = str.MarkerSize;
        elseif isfield(str, 'markerSize')
            style.MarkerSize = str.markerSize;
        end
        if isfield(str, 'MarkerFillColor')
            style.MarkerFillColor = str.MarkerFillColor;
        elseif isfield(str, 'markerFillColor')
            style.MarkerFillColor = str.markerFillColor;
        end
        
        
        % parse line style modifiers
        if isfield(str, 'LineVisible')
            style.LineVisible = str.LineVisible;
        elseif isfield(str, 'lineVisible')
            style.LineVisible = str.lineVisible;
        end
        if isfield(str, 'LineColor')
            style.LineColor = str.LineColor;
        elseif isfield(str, 'lineColor')
            style.LineColor = str.lineColor;
        end
        if isfield(str, 'LineWidth')
            style.LineWidth = str.LineWidth;
        elseif isfield(str, 'lineWidth')
            style.LineWidth = str.lineWidth;
        end
        if isfield(str, 'LineStyle')
            style.LineStyle = str.LineStyle;
        elseif isfield(str, 'lineStyle')
            style.LineStyle = str.lineStyle;
        end
         
        % parse fill style modifiers
        if isfield(str, 'FillVisible')
            style.FillVisible = str.FillVisible;
        elseif isfield(str, 'fillVisible')
            style.FillVisible = str.fillVisible;
        end
        if isfield(str, 'FillColor')
            style.FillColor = str.FillColor;
        elseif isfield(str, 'fillColor')
            style.FillColor = str.fillColor;
        end
        if isfield(str, 'FillOpacity')
            style.FillOpacity = str.FillOpacity;
        elseif isfield(str, 'fillOpacity')
            style.FillOpacity = str.fillOpacity;
        end
        
        % parse face style modifiers
        if isfield(str, 'FaceVisible')
            style.FaceVisible = str.FaceVisible;
        elseif isfield(str, 'faceVisible')
            style.FaceVisible = str.faceVisible;
        end
        if isfield(str, 'FaceColor')
            style.FaceColor = str.FaceColor;
        elseif isfield(str, 'faceColor')
            style.FaceColor = str.faceColor;
        end
        if isfield(str, 'FaceOpacity')
            style.FaceOpacity = str.FaceOpacity;
        elseif isfield(str, 'faceOpacity')
            style.FaceOpacity = str.faceOpacity;
        end
        
    end
    
    function style = read(fileName)
        % Read a style from a file in JSON format.
        style = Style.fromStruct(loadjson(fileName));
    end
end

end % end classdef
