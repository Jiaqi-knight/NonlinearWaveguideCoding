%% ================================================================
% gdgproj:  For one of several parameterized surfaces, plots tangent-vector
% fields, unit normal field, all 4 metric coefficients, and scalar curvature.
% ----------------------------------------------------------------
% Usage: gdgproj {which_surface}
% where which_surface is one of the following names.
% 
% Name     Description                             Coordinate type   s    t
% ----     --------------------------------------  ---------------  --- -----
% flat     xy plane in R^3                         graph             x    y
% s2       two-sphere                              phi, theta       phi theta
% t2       torus                                   phi, theta       phi theta
% trough   parabola translation surface            graph             x    y
% pipe     circular cylinder                       y, theta          y  theta
% sor      surface of revolution of sin(y)/y       y, theta          y  theta
% prbd     paraboloid z = x^2 + y^2                graph             x    y
% uprbd    upside-down paraboloid                  graph             x    y
% hyp2top  top sheet of hyperboloid of two sheets  graph             x    y
% hyp1     hyperboloid of one sheet                z, theta          z  theta
% saddle   standard saddle                         graph             x    y
% monkey   monkey saddle                           graph             x    y
% ----------------------------------------------------------------
% To use this program:
% * Put this file and all the other *.m files from
%   http://math.arizona.edu/~kerl/gdg in a directory of your choosing.
% * Start Matlab in that directory.
% * Type 'help gdgproj' at the Matlab prompt.
% * Type 'gdgproj s2', 'gdgproj t2', etc. at the Matlab prompt.
% * Ooh and aah at the pretty pictures.
% That's about it.
% ----------------------------------------------------------------
% John Kerl
% Final project for Math 537A
% 2007-11-27
% ================================================================

%% ================================================================
function gdgproj(which_surface)

% Find out which surface to work with.
% flatsetup() is in flatsetup.m, etc.
if (strcmp(which_surface, 'flat'))
	[x,y,z,ds,dt] = flatsetup();
elseif (strcmp(which_surface, 's2'))
	[x,y,z,ds,dt] = s2setup();
elseif (strcmp(which_surface, 't2'))
	[x,y,z,ds,dt] = t2setup();

elseif (strcmp(which_surface, 'trough'))
	[x,y,z,ds,dt] = troughsetup();
elseif (strcmp(which_surface, 'pipe'))
	[x,y,z,ds,dt] = pipesetup();
elseif (strcmp(which_surface, 'sor'))
	[x,y,z,ds,dt] = sorsetup();

elseif (strcmp(which_surface, 'prbd'))
	[x,y,z,ds,dt] = prbdsetup();
elseif (strcmp(which_surface, 'uprbd'))
	[x,y,z,ds,dt] = uprbdsetup();
elseif (strcmp(which_surface, 'hyp2top'))
	[x,y,z,ds,dt] = hyp2topsetup();
elseif (strcmp(which_surface, 'hyp1'))
	[x,y,z,ds,dt] = hyp1setup();
elseif (strcmp(which_surface, 'saddle'))
	[x,y,z,ds,dt] = saddlesetup();
elseif (strcmp(which_surface, 'monkey'))
	[x,y,z,ds,dt] = monkeysetup();
elseif (strcmp(which_surface, 'wjq'))
	[x,y,z,ds,dt] =wjqsetup();
else
	error 'Unrecognized surface.  Please type ''help gdgproj'' for more information.'
end

% For clutter reduction, I will be writing .eps files to this subdirectory.
mkdir figures

%% ----------------------------------------------------------------
% For quiver plots we need to downsample -- else the arrows become too close
% together.  As my 8-year-old pointed out, they just look like lines until you
% zoom in.  I want the quiver plots to look like arrows, not flowlines, even
% when you haven't zoomed in yet.  Here I set the downsample ratio in one
% place.
%
% P.S. she thinks the surface plots are really cool.  Her favorite is the
% surface of revolution of sin(y)/y, followed closely by the hyperboloid of one
% sheet.  And she told her cousin that Matlab stands for Math Laboratory. :)
mxds=2;

%% ----------------------------------------------------------------
% The setup routines construct the x, y, and z meshgrids, and the ds & dt.
% Here we find the dimensions of those meshgrids.  The names 'ns' and 'nt' are
% short for 'number of mesh points for parameter s' and 'number of mesh points
% for parameter t'.
%
% The x, y, and z are ns x nt matrices of floating-point numbers.

% Here are the semantics of the ubiquitous meshgrid.  First, the example:
%
% >> [s,t]=meshgrid(1:4,5:9)
%
% s =
%
%      1     2     3     4
%      1     2     3     4
%      1     2     3     4
%      1     2     3     4
%      1     2     3     4
%
%
% t =
%
%      5     5     5     5
%      6     6     6     6
%      7     7     7     7
%      8     8     8     8
%      9     9     9     9
%
% Here is what just happened.
%
% * Matlab supports 1D linear meshes via statements such as:
%   o   1:4     expands to [1 2 3 4]
%   o   1:0.5:4 expands to [1.0 1.5 2.0 2.5 3.0 3.5 4.0]
%   o   etc.
%
% * The two arguments to meshgrid are 1D linear meshes:  call them s and t.
%
% * These 1D linear meshes have lengths: call them nt and ns.  Here, nt=4 and
%   ns=5.  (This seems backward but you'll see why in the next paragraph.)
%
% * There are two return values from meshgrid.  They are both matrices
%   of dimension ns x nt.  The rows of the first are ns copies of s;
%   the columns of the second are nt copies of t.  This seems like a bizarre
%   format but it's what Matlab surface-plot routines expect.
%
% * The setup routines establish meshgrids on s & t (x & y, z & theta, etc.).
%   They then compute x, y, and z via pointwise functions of s and t.  For
%   example, the standard saddle is:
%
%   x = s
%   y = t
%   z = s.^2 - t.^2

[ns,nt] = size(x);

%% ----------------------------------------------------------------
% Graph (a portion of) the surface.
% ...
% Actually, don't.  We'll be seeing many plots of the surface below, colored
% with interesting things such as the scalar curvature.  surf(x,y,z) with
% no additional arguments will color the surface using the height function,
% which isn't particularly interesting.

surf(x,y,z);
shading faceted;
rotate3d on;

% Create surfaceplot.eps
set(gcf, 'PaperPositionMode', 'Manual');
set(gcf, 'PaperPosition', [0 0 6.0 3.0])
print -depsc figures/surfaceplot.eps

%% ----------------------------------------------------------------
% Compute D/Ds's and D/Dt's.  (For the comments in this file, which are of
% course in plain ASCII rather than TeX, I use capital D for the
% partial-derivative symbol.)

% Matlab has a gradient routine.  Here is the relationship between gradients
% and tangent vectors:
%
% * x, y, and z each are scalar functions of the two arguments s and t.
%
% * The three gradients are
%   grad x = [ Dx/Ds, Dx/Dt ]
%   grad y = [ Dy/Ds, Dy/Dt ]
%   grad z = [ Dz/Ds, Dz/Dt ]
%
% * The tangent vectors we want are
%   D/Ds   = [ Dx/Ds, Dy/Ds, Dz/Ds ]
%   D/Dt   = [ Dx/Dt, Dy/Dt, Dz/Dt ].
%
% We'll actually do that reorganization below when we compute the g_ij's.

% Implementation detail:  the Matlab gradient uses centered differences except
% at endpoints; it uses left and right differences at left and right endpoints,
% respectively.  Type 'help gradient' at the Matlab prompt for more
% information.

[DxDs,DxDt] = gradient(x,ds,dt);
[DyDs,DyDt] = gradient(y,ds,dt);
[DzDs,DzDt] = gradient(z,ds,dt);

%% ----------------------------------------------------------------
% Plot D/Ds's and D/Dt's.
%
% Plot the surface with zeros for coloring to give the arrows something to sit
% on.  Else we have a bunch of arrows floating around in space & it's not as
% visually appealing.
%
% The 'figure' statement tells Matlab to open a new figure window.
%
% The 'hold on' statement tells Matlab to leave the surface plot showing and
% overlay the quiver plot, rather than having the quiver plot replace the
% surface plot.
%
% downsample is in downsample.m.  Its purpose was described
% above.

figure;
surf(x,y,z,zeros(ns,nt)); shading faceted; hold on;
quiver3(...
	downsample(x,    mxds),...
	downsample(y,    mxds),...
	downsample(z,    mxds),...
	downsample(DxDs, mxds),...
	downsample(DyDs, mxds),...
	downsample(DzDs, mxds));
title '\partial/\partial s';

% Create DDs.eps
set(gcf, 'PaperPositionMode', 'Manual');
set(gcf, 'PaperPosition', [0 0 3.0 3.0])
print -depsc figures/DDs.eps

figure;
surf(x,y,z,zeros(ns,nt)); shading faceted; hold on;
quiver3(...
	downsample(x,    mxds),...
	downsample(y,    mxds),...
	downsample(z,    mxds),...
	downsample(DxDt, mxds),...
	downsample(DyDt, mxds),...
	downsample(DzDt, mxds));
title '\partial / \partial t';

% Create DDt.eps
set(gcf, 'PaperPositionMode', 'Manual');
set(gcf, 'PaperPosition', [0 0 3.0 3.0])
print -depsc figures/DDt.eps

%% ----------------------------------------------------------------
% The metric g_{ij} is implicitly defined by the R^3 dot product.
% Compute g_ij's:
%   g11 = <DDs, DDs>   g12 = <DDs, DDt>
%   g21 = <DDt, DDs>   g22 = <DDt, DDt>

% Even though the g_ij's are symmetrical, so g_12 == g_21, compute and display
% them anyway.  I always find it satisfying when the numerics *do* the things
% the theory says they *should* do.

% Above I said "We'll actually do that reorganization below when we compute
% the g_ij's."  Here it is.
%
% Here is an example for 4x4 parameter mesh.  The Dx/Ds et al. are, for each
% point on the surface, just a scalar.  The collection of Dx/Ds over all points
% on the parameter mesh is a matrix (here, 4x4).  To find the D/Ds vector at a
% point (say, the 2,1 element of the parameter mesh) we need to collect the
% three numbers labeled with '*' below and put them in a vector, doing this for
% each point on the mesh.

% Dx/Ds o o o o     Dx/Dt o o o o
%       * o o o           o o o o
%       o o o o           o o o o
%       o o o o           o o o o

% Dy/Ds o o o o     Dy/Dt o o o o
%       * o o o           o o o o
%       o o o o           o o o o
%       o o o o           o o o o

% Dz/Ds o o o o     Dz/Dt o o o o
%       * o o o           o o o o
%       o o o o           o o o o
%       o o o o           o o o o

%       |
%       |
%       v

%     [ * ] D/Ds
%     [ * ]
%     [ * ]

% The triples DDs1, DDt1, N1, and N1hat are just single vectors at a point.
% DDs, DDt, and Nhat are vector fields.
%
% Storage details:  Matlab has always supported scalars, vectors, and matrices.
% Leaving questions of covariance and contravariance aside, one might think of
% those as rank-0, rank-1, and rank-2 arrays.  At some point in the development
% of the development of Matlab, they began to support higher-rank arrays.
% But some things are weird.  Consider the following asymmetry:
%
% >> A = zeros(5,6,7,8)   % Construct a zero-filled rank-4 array
%                         % with dimensions 5, 6, 7, 8.
% >> size(A)
% ans =
%      5   6   7   8
%
% >> size(A(1,1,:,:))
% ans =
%      1   1   7   8
%
% >> size(A(:,:,1,1))
% ans =
%      5   6
%
% I want to use :'s whenever possible in order to let Matlab do the looping for
% me.  I also want the indexing to be such that Matlab will give me a matrix
% (rank-2) array back:  e.g. the 5x6 rather than the silly 1x1x7x8.  This
% affects where I place my indices.

% Allocate space for single vectors.
DDs1 = [0 0 0]; DDt1 = [0 0 0];

% Allocate space for vector fields.
% Put the s and t indices first so that size(DDs(:,:,3)) is ns x nt.
DDs  = zeros(ns, nt, 3);
DDt  = zeros(ns, nt, 3);
Nhat = zeros(ns, nt, 3);

% Allocate space for metric coefficients.
% Put the s and t indices first so that size(g(:,:,1,1)) is ns x nt.
g    = zeros(ns,nt,2,2);
invg = zeros(ns,nt,2,2);

% Here is the one place in this file where I loop over s and t meshgrid
% indices.
for ss = 1:ns
	for tt = 1:nt

		DDs1(1) = DxDs(ss,tt); DDt1(1) = DxDt(ss,tt);
		DDs1(2) = DyDs(ss,tt); DDt1(2) = DyDt(ss,tt);
		DDs1(3) = DzDs(ss,tt); DDt1(3) = DzDt(ss,tt);

		% How to do the dot products:
		%
		% * Matlab vectors are row vectors when punctuated with commas or
		%   spaces, or column vectors when punctuated with semicolons.
		% * The dot product is implemented in Matlab using matrix
		%   multplication with transpose of one vector.
		% * Since I am using row vectors, I need to transpose the second
		%   vector.
		% * Matrix transpose, in turn, is done using the ' operator.

		% Metric coefficients.  Call them a, b, c, d for the moment; we'll
		% insert them in the g and ingv arrays a few lines below.
		a = DDs1 * DDs1'; b = DDs1 * DDt1';
		c = DDt1 * DDs1'; d = DDt1 * DDt1';

		g(ss,tt,1,1) = a; g(ss,tt,1,2) = b;
		g(ss,tt,2,1) = c; g(ss,tt,2,2) = d;

		% Inverse metric coefficients
		gdet = a*d - b*c;

		invg(ss,tt,1,1) =  d/gdet; invg(ss,tt,1,2) = -b/gdet;
		invg(ss,tt,2,1) = -c/gdet; invg(ss,tt,2,2) =  a/gdet;

		% Unit normals
		N1    = cross(DDs1, DDt1);
		N1hat = N1 / norm(N1);

		% Store the tangent basis and unit normals for later.
		DDs (ss,tt,:) = DDs1;
		DDt (ss,tt,:) = DDt1;
		Nhat(ss,tt,:) = N1hat;

	end
end

% Quiver-plot the unit normals.
%
% WARNING:  If z scaling is different than x & y scaling, the normals may not
% look perpendicular.  For the same reason, a parameterized plot of cos(t) and
% sin(t) may not look circular on your TI-83.

figure;
surf(x,y,z,zeros(ns,nt)); shading faceted; hold on;
quiver3(...
	downsample(x,           mxds),...
	downsample(y,           mxds),...
	downsample(z,           mxds),...
	downsample(Nhat(:,:,1), mxds),...
	downsample(Nhat(:,:,2), mxds),...
	downsample(Nhat(:,:,3), mxds));
title 'Nhat';

% Create Nhat.eps
set(gcf, 'PaperPositionMode', 'Manual');
set(gcf, 'PaperPosition', [0 0 3.0 3.0])
print -depsc figures/Nhat.eps

% Surface-plot the metric coefficients.
figure; surf(x,y,z,   g(:,:,1,1)); shading faceted; colorbar; title 'g_{11}';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/g11.eps
figure; surf(x,y,z,   g(:,:,1,2)); shading faceted; colorbar; title 'g_{12}';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/g12.eps
figure; surf(x,y,z,   g(:,:,2,1)); shading faceted; colorbar; title 'g_{21}';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/g21.eps
figure; surf(x,y,z,   g(:,:,2,2)); shading faceted; colorbar; title 'g_{22}';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/g22.eps

figure; surf(x,y,z,invg(:,:,1,1)); shading faceted; colorbar; title 'g_{11}^{-1}';
figure; surf(x,y,z,invg(:,:,1,2)); shading faceted; colorbar; title 'g_{12}^{-1}';
figure; surf(x,y,z,invg(:,:,2,1)); shading faceted; colorbar; title 'g_{21}^{-1}';
figure; surf(x,y,z,invg(:,:,2,2)); shading faceted; colorbar; title 'g_{22}^{-1}';

%% ----------------------------------------------------------------
% Compute derivatives of the metric coefficients, which are needed to compute
% the Christoffel symbols.  These are D/Di g_jk.

DD_g = zeros(ns,nt,2,2,2);

for j = 1:2
	for k = 1:2
		[DD_g(:,:,1,j,k), DD_g(:,:,2,j,k)] = gradient(g(:,:,j,k),ds,dt);
	end
end

%% ----------------------------------------------------------------
% Compute the Christoffel symbols:
%
% Gamma_ij^k = 1/2 g^km ( D/Di g_jm + D/DDj g_im - D/DDm g_ij )

Gamma = zeros(ns,nt,2,2,2);

for i = 1:2
	for j = 1:2
		for k = 1:2
			for m = 1:2
				Gamma(:,:,i,j,k) = Gamma(:,:,i,j,k) ...
					+ 0.5 * invg(:,:,k,m) .* ( ...
						DD_g(:,:,i,j,m) + ...
						DD_g(:,:,j,i,m) - ...
						DD_g(:,:,m,i,j)   ...
					);
			end
		end
	end
end

% Surface-plot the Christoffel symbols.

figure; surf(x,y,z, Gamma(:,:,1,1,1)); shading faceted; colorbar; title '\Gamma_{11}^1';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G111.eps

figure; surf(x,y,z, Gamma(:,:,1,2,1)); shading faceted; colorbar; title '\Gamma_{12}^1';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G121.eps

figure; surf(x,y,z, Gamma(:,:,2,1,1)); shading faceted; colorbar; title '\Gamma_{21}^1';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G211.eps

figure; surf(x,y,z, Gamma(:,:,2,2,1)); shading faceted; colorbar; title '\Gamma_{22}^1';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G221.eps

figure; surf(x,y,z, Gamma(:,:,1,1,2)); shading faceted; colorbar; title '\Gamma_{11}^2';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G112.eps

figure; surf(x,y,z, Gamma(:,:,1,2,2)); shading faceted; colorbar; title '\Gamma_{12}^2';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G122.eps

figure; surf(x,y,z, Gamma(:,:,2,1,2)); shading faceted; colorbar; title '\Gamma_{21}^2';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G212.eps

figure; surf(x,y,z, Gamma(:,:,2,2,2)); shading faceted; colorbar; title '\Gamma_{22}^2';
set(gcf, 'PaperPositionMode', 'Manual'); set(gcf, 'PaperPosition', [0 0 3.0 3.0]); print -depsc figures/G222.eps


%% ----------------------------------------------------------------
% Quiver-plot some covariant derivatives?

%% ----------------------------------------------------------------
% Compute D/Di Gamma_jk^l:  Needed for curvature.

DD_Gamma = zeros(ns,nt,2,2,2,2);

for j = 1:2
	for k = 1:2
		for l = 1:2
			[DD_Gamma(:,:,1,j,k,l), DD_Gamma(:,:,2,j,k,l)] = ...
				gradient(Gamma(:,:,j,k,l),ds,dt);
		end
	end
end

%% ----------------------------------------------------------------
% The 1,3 curvature tensor:
%
%	R_ijk^l =
%		  D/Di \Gamma_jk^l + \Gamma_ik^m \Gamma_jm^l
%		- D/Dj \Gamma_ik^l - \Gamma_jk^m \Gamma_im^l.

R13 = zeros(ns,nt,2,2,2,2);

for i = 1:2
	for j = 1:2
		for k = 1:2
			for l = 1:2
				for m = 1:2
					R13(:,:,i,j,k,l) = R13(:,:,i,j,k,l) ...
						+ DD_Gamma(:,:,i,j,k,l) ...
						- DD_Gamma(:,:,j,i,k,l) ...
						+ Gamma(:,:,i,k,m) .* Gamma(:,:,j,m,l) ...
						- Gamma(:,:,j,k,m) .* Gamma(:,:,i,m,l);
				end
			end
		end
	end
end

%% ----------------------------------------------------------------
% The 0,4 curvature tensor:
%
% Form R_ijkl from R_ijk^l: R_ijkl = R_ijk^m g_lm.
%
% Matlab doesn't support index contraction by a built-in function, but you
% certainly can do it.

R04 = zeros(ns,nt,2,2,2,2);
for i = 1:2
	for j = 1:2
		for k = 1:2
			for l = 1:2
				for m = 1:2
					R04(:,:,i,j,k,l) = R04(:,:,i,j,k,l) ...
						+ R13(:,:,i,j,k,m) .* g(:,:,l,m);
			 	end
			end
		end
	end
end

%% ----------------------------------------------------------------
% The 0,2 curvature tensor:
%
% R_ij = g^kl R_kijl.

R02 = zeros(ns,nt,2,2);
for i = 1:2
	for j = 1:2
		for k = 1:2
			for l = 1:2
				R02(:,:,i,j) = R02(:,:,i,j) ...
					+ invg(:,:,k,l) .* R04(:,:,k,i,j,l);
			end
		end
	end
end

%% ----------------------------------------------------------------
% Compute the scalar curvature:  S = g^ij R_ij
% (which in turn is g^kl g^ij R_kijl).

S = zeros(ns,nt);

for i = 1:2
	for j = 1:2
		S(:,:) = S(:,:) + invg(:,:,i,j) .* R02(:,:,i,j);
	end
end
figure; surf(x,y,z,S); shading faceted; colorbar; title 'Scalar curvature';

% Create scalarcurvature.eps
set(gcf, 'PaperPositionMode', 'Manual');
set(gcf, 'PaperPosition', [0 0 3.0 3.0])
print -depsc figures/scalarcurvature.eps

%% ================================================================
% To do:
% Sectional curvature:

%	                R(X,Y,Y,X)
%	K(X,Y) = --------------------------
%	         g(X,X) + g(Y,Y) - g(X,Y)^2

% Note that in dimension 2, R02 = 2K.  But compute & plot it anyway for
% verification.

%% ----------------------------------------------------------------
% To do:
% Second fundamental form

%% ----------------------------------------------------------------
% To do:
% * Draw some geodesics.  This probably should be done in the setup routines.
% * Compute D_t's.
% * Plot some Jacobi fields.
% * ...
% * ...
% * ...