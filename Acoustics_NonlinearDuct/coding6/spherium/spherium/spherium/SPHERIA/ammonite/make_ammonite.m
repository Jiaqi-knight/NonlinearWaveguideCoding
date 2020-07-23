%make_ammonite
% Function which creates a 3D ammonite.

function surface_handle = make_ammonite( add_ridges, add_ridges_to_colour, add_bumps,...
    add_bumps_to_colour, spiral_type, spiral_turns,points_per_turn, cross_section_ratio,...
    bump_amplitude, ridge_frequency, spiral_bump_amplitude, spiral_bump_frequency, helicity )
points_per_turn=20

% Define spiral function. r is range, t is angle measured clockwise from y
% axis (in radians)
if strcmp( spiral_type, 'Archimedian' )
    spiral_func_str = 'r = (t/2*pi).^2';
elseif strcmp( spiral_type, 'Linear'  )
    spiral_func_str = 'r = t/(2*pi)';
else
    %Logarithmic spiral
    spiral_func_str = 'r = exp(t/(2*pi)) - 1';
end

%Define array of cross section angles alpha
alpha = linspace(0,2*pi,points_per_turn);

%Define array of spiral angles t
t = linspace(0,2*pi*spiral_turns,points_per_turn*spiral_turns);
scale=[0 :0.05: 1];
%Mesh these
[alpha,t]=meshgrid(alpha,t);
figure
for s=scale
surface_matrix = mesh(alpha,t,s*10*ones(size(t)));hold on
shading interp;
axis equal
axis vis3d
end
%Define cross section ellipse semi-major radii
R = 0.5*spiral(t+2*pi,spiral_func_str) - 0.5*spiral(t,spiral_func_str);
figure
for s=scale
R = 0.5*spiral(t+2*pi,spiral_func_str) - 0.5*spiral(t,spiral_func_str);
R=  s*R;  
%Define cross section ellipse radius
A = R.*( sin(alpha).^2 + cross_section_ratio*cross_section_ratio*cos(alpha).^2 ).^(-0.5);

%Add ridges
% if add_ridges==1
%     A = A.* ( bump_amplitude*cos(ridge_frequency*alpha) + 1 );
% end
% if add_ridges_to_colour==1
%     C = A.* ( bump_amplitude*cos(ridge_frequency*alpha) + 1 );
% else
%     C=A;
% end

%Add bumps along spiral
% if add_bumps==1
%     A = A.* ( spiral_bump_amplitude*cos(spiral_bump_frequency*t) + 1 );
% end
% if add_bumps_to_colour==1
%     C = C.* ( spiral_bump_amplitude*cos(spiral_bump_frequency*t) + 1 );
% end

%Define cross section ellipse Cartesian coordinates
a = A.*sin(alpha);
x = A.*cos(alpha)+helicity*t;
%Determine surface in x,y,z coordinates
k = a + 0.5*spiral(t+2*pi,spiral_func_str) + 0.5*spiral(t,spiral_func_str);
y = k.*sin(t);
z = k.*cos(t);
% R = sqrt( x.^2 + y.^2 + z.^2);
% Rmax = max(max(R));
% x = x./Rmax;
% y = y./Rmax;
% z = z./Rmax;
% C = C./Rmax;

%Create surface and return a handle to it
% figure
surface_handle = mesh(x,y,z);hold on


%surface_handle = surf(x,y,z,C);

shading interp;
axis equal
axis vis3d
% axis off
end
%%

%Function which evaluates spiral radius from angle t /radians.
function r = spiral(t,spiral_func_str)
eval([ spiral_func_str,';'] )

%End of code


