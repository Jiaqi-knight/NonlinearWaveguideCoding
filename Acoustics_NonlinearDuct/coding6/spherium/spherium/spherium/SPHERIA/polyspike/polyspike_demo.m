function polyspike_demo

%Number of azimuth spikes
M = 2;

%Number of elevation spikes
N = 5;

%Spikiness
k = 0;

%Number of parabola points between spikes
P = 40;

%Create polyspike plot
figure('color',[1 1 1],'name','polyspike test','renderer','opengl')
surface_handle = make_polyspike(N,M,k,P);
shading interp
lighting phong
camlight
axis vis3d
axis off
zoom(1.5)
print( gcf, 'polyspike_test', '-dpng', '-r300' )

%End of code