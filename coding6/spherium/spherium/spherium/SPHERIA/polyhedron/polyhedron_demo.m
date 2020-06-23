function polyhedron_demo

%Edges of a single cross section is 2*(N+1)
N = 5;

%Create polyspike plot
figure('color',[1 1 1],'name','polyhedron test','renderer','opengl')
surface_handle = make_polyhedron(N);
lighting phong
camlight
zoom(1.5)
print( gcf, 'polyhedron_test', '-dpng', '-r300' )

%End of code
