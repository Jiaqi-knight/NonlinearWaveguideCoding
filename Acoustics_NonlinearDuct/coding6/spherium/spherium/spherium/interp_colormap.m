%interp_colormap
% Function which interpolates current colourmap to yield better graduated
% shading. N is number of possible colours.
%
% LAST UPDATED by Andy French 14 September 2008

function interp_colormap( N )

%Get current colourmap
map = colormap;

%Initialise new colormap
new_map = ones(N,3);

%Get size of current colormap and initalise red,green,blue vectors
dim = size(map);
R = ones(1,dim(1));
G = ones(1,dim(1));
B = ones(1,dim(1));
RR = ones(1,N);
GG = ones(1,N);
BB = ones(1,N);

%Populate these with current colormap
R(:) = map(:,1);
G(:) = map(:,2);
B(:) = map(:,3);

%Interpolate to yield new colour map
x = linspace( 1, dim(1), N );
RR = interp1( 1:dim(1), R, x );
GG = interp1( 1:dim(1), G, x );
BB = interp1( 1:dim(1), B, x );
new_map(:,1) = RR(:);
new_map(:,2) = GG(:);
new_map(:,3) = BB(:);

%Set colormap to be new map
colormap( new_map );

%End of code