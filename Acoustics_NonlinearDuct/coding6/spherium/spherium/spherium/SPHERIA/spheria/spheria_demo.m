function spheria_demo

N = 300;
spherefunction = 'abs(cos(elev*azi)*sin(azi*elev))'; 

figure('name','spheria_demo','color',[1 1 1] )

subplot(1,2,1);
surface_handle1 = sphericalplot( 'Sphere', spherefunction,  N );
shading interp;
axis equal
axis vis3d
axis off
lighting phong
camlight
zoom(1.4)

subplot(1,2,2);
surface_handle2 = sphericalplot( 'Surface', spherefunction,  N );
shading interp;
axis equal
axis vis3d
axis off
lighting phong
camlight
zoom(1.4)

%End of code