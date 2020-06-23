%sphericalplot
% Demonstrates how to plot a function defined in elevation and azimuth
% co-ordinates. Two representations are offered. The first plots a
% sphere with a defined radius and colours this propotional to the value of
% the function evaluated at the radius defined. The second creates a
% surface with radius proportional to the value of the function (at a
% specified actual radius parameter) at each azimuth,elevation co-ordinate.
%
% Good functions!
% cos(elev)*atan(azi)                 - wierd shell like thingy
% cos(elev)*cos(elev)                 - gain of a dipole
% abs(cos(elev)*atan(elev))           - fruit bowl on inverted fruit bowl
% cos(elev)*atan(azi*elev)            - another strange shell-like thing
% abs(cos(elev*azi)*sin(azi*elev))    - bizzare pinecone

function surface_handle = sphericalplot( sphere_or_surface, spherefunction,  N )

%Define elevation and azimuth arrays
elev=linspace(-pi/2,pi/2,N);
azi=linspace(0,2*pi,N);
[azi,elev]=meshgrid(azi,elev);

%Define x,y,z coordinate matricies
x = cos(elev).*cos(azi);
y = cos(elev).*sin(azi);
z = sin(elev);

%Define sphere power function in terms of azimuth (azi) and elevation (elev). Add '.' to force
%element-wise matrix operations.
spherefunction = strrep( spherefunction,'*','.*' );
spherefunction = strrep( spherefunction,'/','./' );
spherefunction = strrep( spherefunction,'^','.^' );

%Evaluate spherefunction and genrate surface (or surface
%colouring) matrix
try
    eval( ['P=',spherefunction,';'] );
    P = P./max(max(abs(P)));
catch
    errordlg('Invalid spherefucntion!',...
        'Spherium Error','modal');
    return;
end

%Plot spherical surface
if strcmp( sphere_or_surface,'Surface' )
    %Plot 3D surface, radius is equal to the magnitude of the
    %spherefunction
    x = abs(P).*x;
    y = abs(P).*y;
    z = abs(P).*z;
end
surface_handle=surf(x,y,z,P);
axis equal
axis vis3d
axis tight
axis off
shading interp

%End of code.
