%spharm
% Function which computes spherical harmonics, functions of azimuth and
% elevation. These are the angular solution to the Schrodinger Equation for
% the Hydrogen Atom.
%
% Syntax: Y = spharm( azi, elev, L, M )
%
% azi and elev are matrices of angles in radians.
% The orbital integers L and M characterize the order of
% the spherical harmonic.
%
% LAST UPDATED by Andy French 30/March/2011

function Y = spharm( azi, elev, L, M )

%Compute associated Legendre functions
y = legendre( L, cos( -elev + pi/2 ) );

%Initialise output and extract row M from y
dim = size(elev);
Y = NaN(dim);
if L==0
    %Check for vector or matrix input
    if isempty(find(dim==1, 1))
        %Matrix
        Y(:,:) = y( :, : );
    else
        %Vector, one dimension has length 1
        Y(:) = y(:);
    end
elseif ( L==round(L) ) && (L>0)
    %Check for vector or matrix input
    if isempty(find(dim==1, 1))
        %Matrix
        Y(:,:) = y( abs(M)+1, :, : );
    else
        %Vector, one dimension has length 1
        Y(:) = y( abs(M)+1, : );
    end
else
    return
end

%Consider negative values of M
if M<0
    MM = abs(M);
    Y = Y * ( (-1)^MM ) * factorial(L-MM)/ factorial(L+MM) ;
end

%Mutiply Y by normalization and phase (azi) term
Y = Y .* exp( 1i*M*azi ) * ( (-1)^M ) *...
    sqrt( (2*L+1)*factorial(L-M)/( 4*pi*factorial(L+M) ) );

%End of code

