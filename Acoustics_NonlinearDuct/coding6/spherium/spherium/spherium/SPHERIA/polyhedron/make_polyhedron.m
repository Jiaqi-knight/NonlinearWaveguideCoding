%make_polyhedron
% Function which creates matrices of X,Y,Z vertices of a solid surface
% formed via a volume of revolution of a regular N-gon. The lengths of all
% vertices are set to be one unit. i.e. sqrt( x.^2 + y.^2 + z.^2 ) =
% ones(size(x)). A surface is then formed from these matrices.
%
% Note the regular N-gon will have 2(N+1) vertices as it will always have
% two vertices on the rotation axis, which is also the line of symmetry.
%
% LAST UPADTED by Andy French August 2012

function surface_handle = make_polyhedron(N)

%Number of vertices
M = 2*(N+1);

%Create normalized N-gon
elev = linspace( -pi, pi, M+1 );
azi = linspace( 0, 2*pi, M+1 );
x = zeros(M+1);
y = zeros(M+1);
z = zeros(M+1);

%Aximuth
for n=1:M+1
    %Elevation
    for m=1:M+1
        z(n,m) = cos( elev(m) );
        y(n,m) = sin( elev(m) )*cos( azi(n) );
        x(n,m) = sin( elev(m) )*sin( azi(n) );
    end
end

%Create surface
surface_handle = surf( x,y,z,rand(size(x)));
axis equal
axis vis3d
axis off
shading interp

%End of code