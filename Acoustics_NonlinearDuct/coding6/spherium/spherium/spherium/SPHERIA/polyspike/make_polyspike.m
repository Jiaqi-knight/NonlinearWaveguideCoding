%make_polyspike
% Function which creates a surface from matrices of X,Y,Z vertices of the
% solid formed via a volume of revolution of a regular N-gon. The lengths
% of all vertices are set to be one unit. i.e. sqrt( x.^2 + y.^2 + z.^2 ) =
% ones(size(x)).
%
% Note the regular N-gon will have 2(N+1) vertices as it will always have
% two vertices on the rotation axis, which is also the line of symmetry.
%
% k is a spikiness parameter. If abs(k)>0 then a parabola will be drawn
% between the vertices of the polygon.
%
% LAST UPADTED by Andy French August 2012

function surface_handle = make_polyspike(N,M,k,P)

%Azimuth part
theta = pi/M;
xx = linspace(-sin(theta),sin(theta),P);
yy = k*xx.^2 + cos(theta) - k*sin(theta)^2;
r = sqrt( xx.^2 + yy.^2 );
r_azi = repmat(r,[1,M]);

%Determine maximum k value to contain parbola between spikes
kmax_azi = cos(theta)/(1-cos(theta)^2);

%Elevation part
theta = 0.5*pi/N;
xx = linspace(-sin(theta),sin(theta),P);
yy = k*xx.^2 + cos(theta) - k*sin(theta)^2;
r = sqrt( xx.^2 + yy.^2 );
r_elev = repmat(r,[1,N]);

%Determine maximum k value to contain parbola between spikes
kmax_elev = cos(theta)/(1-cos(theta)^2);

%Combine range matrices
[r_elev,r_azi] = meshgrid(r_elev,r_azi);

%Create normalized N-gon
elev = linspace( 0, pi, N*P );
azi = linspace( 0, 2*pi, M*P );
[elev,azi] = meshgrid(elev,azi);
z = r_azi.*r_elev.*cos(elev);
y = r_azi.*r_elev.*sin(elev).*cos(azi);
x = r_azi.*r_elev.*sin(elev).*sin(azi);
r = sqrt( x.^2 + y.^2 + z.^2 );
rmax = max(max(r));
x = x/rmax;
y = y/rmax;
z = z/rmax;
C = sqrt(x.^2 + y.^2 + z.^2);

%Create surface and return a handle to it
surface_handle = surf(x,y,z,C);
shading interp;
axis equal
axis vis3d
axis off

%End of code