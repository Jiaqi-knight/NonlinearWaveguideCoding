function figureFULL(varargin)
    if isempty(varargin)
        scrsz = 0.9*get(0,'ScreenSize'); % full screen looks better
    else
        scrsz = varargin{1};
    end
    figure('Position', scrsz, 'Units', 'normalized');
    axes('Position',[0 0 1 1], 'Units','normalized');
end