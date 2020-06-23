function edit_bespoke_colormap
try
    load('bespoke_colormap.mat')
    figure('units','normalized','position',[0,0,0.01,0.01],'toolbar','none','menubar','none');
    colormap(map);
    colormapeditor
    clear map
end