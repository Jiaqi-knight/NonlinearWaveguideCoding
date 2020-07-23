%cropfig
% Function which transforms a figure window to crop to the image contained
% within it.

function z = cropfig( fig, ax, method )

%Get original axis units
orig_ax_units = get(ax,'units');
orig_fig_units = get(fig,'units');

%Set units to pixels and crop figure to objects within it
set(ax,'units','pixels');
set(fig,'units','pixels');

%Get figure and axes bounding box dimensions in pixels
t = get(ax,'tightinset');
p = get(ax,'position');
po = get(ax,'outerposition');
fp = get(fig,'position');

if method==1 %Crop to the tight inset of the figure
    
    %Work out x and y shift to place bottom left hand corner axes + tight inset
    %in the bottom left hand corner of the figure
    dx = -p(1)+t(1);
    dy = -p(2)+t(2);
    
    %Crop figure to contents based upon tight inset
    po = [po(1)+dx,po(2)+dy,po(3),po(4)];
    set( ax,'outerposition',po);
    p = [p(1)+dx,p(2)+dy,p(3),p(4)];
    set( ax,'position',p);
    fp = [fp(1),fp(2),p(3)+t(1)+t(3),p(4)+t(2)+t(4)];
    set( fig,'position',fp);
    
else %Crop to the extent of the non-white pixels. Use this when axis labels are turned off.
    
    %Build in 2% border safety factor
    
    %Set background color of axes to be white
    background_color = get(fig,'color');
    set(fig,'color',[1 1 1]);
    
    %Obtain a bitmap image
    f=getframe(fig);
    
    %Determine mean colour matrix
    m = mean(f.cdata,3);
    dim = size(m);
    width_pix = dim(2);
    height_pix = dim(1);
    b = round(0.02*min([width_pix,height_pix]));
    
    %Determine maximum and minimum indices of pixels which are NOT white
    [row_ind,col_ind] = find(m~=255);
    min_row = min(row_ind);
    max_row = max(row_ind);
    min_col = min(col_ind);
    max_col = max(col_ind);
    
    %Shift axes so that bottom left hand corner is the bottom left hand corner
    %of the image
    p = [ p(1)-min_col+b,...
        p(2)-height_pix+max_row+b,...
        p(3),...
        p(4)] ;
    po = [ po(1)-min_col + b,...
        po(2)-height_pix+max_row + b,...
        po(3),...
        po(4)] ;
    set( ax,'outerposition',po);
    set( ax,'position',p,'clipping','on');
    
    %Modify figure to have the same aspect ratio as the image
    fp = [fp(1),fp(2),max_col-min_col + 2*b,max_row-min_row + 2*b];
    set( fig,'position',fp);
    
    %Restore figure background colour
    set(fig,'color',background_color);
end

% %Return units to original settings
set(fig,'units',orig_fig_units);
set(ax,'units',orig_ax_units);

%End of code