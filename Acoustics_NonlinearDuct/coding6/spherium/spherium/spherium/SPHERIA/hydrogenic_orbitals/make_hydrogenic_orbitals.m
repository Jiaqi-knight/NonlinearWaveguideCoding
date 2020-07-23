%make_hydrogenic_orbitals
% Function which creates a matrix of radii (or size NxN) corresponding to
% matrices of azimuth and elevation corresponding to quantum orbitals of
% Hydrogenic atoms. The orbitals are the (angular) solutions
% ('wavefunctions') to the Schrodinger equation and the radius R is
% proportional to the square of the modulus of the wavefunctions.
%
% Reference: pp96 of Woan, "The Cambridge Handbook of Physics formulas"
%
% LAST UPDATED by Andy French 30/3/2011
%
% Energy levels in Hydrogenic type atoms are characterized by quantum
% numbers n = 1, 2, 3, 4....
% Orbital quantum numbers are:
% L = 0,1,2....(n-1)
% M = 0, +/-1, +/-2,...+/- L
%
% Orbitals are typically labelled in the form nXM e.g. 3P1
% X is characterized by L.
% L = 0, X = S   min(n) = 1   max(abs(M)) = 0
% L = 1, X = P   min(n) = 2   max(abs(M)) = 1
% L = 2, X = D   min(n) = 3   max(abs(M)) = 2
% L = 3, X = F   min(n) = 4   max(abs(M)) = 3
% L = 4, X = G   min(n) = 5   max(abs(M)) = 4
% L = 5, X = H   min(n) = 6   max(abs(M)) = 5
% L = 6, X = I   min(n) = 7   max(abs(M)) = 6
%
% By convention XM orbitals are the linear combination of +abs(M) and
% -abs(M) orbitals. For example 3P-1 will combine the (angular) orbitals
% Y(L=1,M=1) - Y(L=1,M=-1). Similarly 3P1 will correspond to Y(L=1,M=1) +
% Y(L=1,M=-1).

function surface_handle = make_hydrogenic_orbitals( L, M, N )

%Check M is valid
if isnumeric(M)==1
    if abs(M) <= L
        
        %Construct azi and elev matrices /radians
        azi = linspace(0,2*pi,N);
        elev = linspace( -pi/2,pi/2,N );
        [azi,elev] = meshgrid( azi,elev );
        
        
        
        %Construct orbitals
        if M==0
            R = spharm( azi, elev, L, 0 );
        elseif M<0
            M=abs(M);
            R = spharm( azi, elev, L, M ) - spharm( azi, elev, L, -M );
        else
            M=abs(M);
            R = spharm( azi, elev, L, M ) + spharm( azi, elev, L, -M );
        end
        R = R.*conj(R);
        
        %Normalize
        R = R/ max(max(R));
        
        %Define x,y,z coordinate matricies
        x = abs(R).*cos(elev).*cos(azi);
        y = abs(R).*cos(elev).*sin(azi);
        z = abs(R).*sin(elev);
    else
        errordlg('Invalid Hydrogenic orbital. Check magnetic quantum number',...
            'Spherium Error','modal');
        R = NaN(N,N);
    end
else
    errordlg('Invalid Hydrogenic orbital.',...
        'Spherium Error','modal');
    R = NaN(N,N);
    
end

%Plot spherical surface
surface_handle=surf(x,y,z,R);
axis equal
axis vis3d
axis tight
axis off
shading interp

%End of code