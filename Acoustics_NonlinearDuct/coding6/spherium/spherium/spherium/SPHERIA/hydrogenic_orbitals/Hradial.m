%Hradial
% Radial dependence of the wavefunction of quantum orbitals
% of Hydrogenic atoms. This assumes an atom consisting of a nucleus of
% charge -Ze, where e is the charge on the electron.
%
% LAST UPDATED by Andy French Dec 2011
%
% [w,E,r_mean,r2_mean,a0,a] = Hradial(r,n,L,Z,A)
%
% r is the radius from the nuclear centre in angstroms (10^-10 metres)
% N is the orbital quantum number
% L is the angular quantum number = 0,1,2....(n-1)
%
% Z is the atomic number (the total charge is Z x the charge on the electron)
% A is the atomic mass number. The atomic mass is A* mass of a carbon-12
% atom / 12 = A*u. u is the 'unified mass constant.'
%
% w is the wavefunction evaluated at all values of r. Note ww* is the
% probability density of the electron cloud.
% E is the orbital energy in eV
% r_mean is the mean radius of the electron cloud.
% r2_mean is the mean squared radius of the elctron cloud
% a0 is the Bohr radius /Angstroms
% a is the atomic radius in Angstroms using a 'Hydrogenic' approximation

function [w,E,r_mean,r2_mean,a0,a] = Hradial(r,n,L,Z,A)

%% Physical constants %%

%Permittivity of free space /Fm^-1
e0 = 8.854187817e-12 ;

%Planck's constant /Js
h = 6.6260755e-34 ;

%Charge on the electron /C
qe = 1.60217733e-19 ;

%Mass of the electron /kg
me = 9.1093897e-31 ;

%Unified mass constant (= 1/12 Mass of a Carbon-12 atom /kg)
u = 1.6605402e-27 ;

%

%Bohr radius /Angstrom's
a0 = (e0*h^2)/(pi*me*qe^2);
a0 = a0/1e-10;

%Reduced mass
mu = me*A*u/(me+A*u);

%Hydrogenic atom 'radius'
a = me*a0/(mu*Z);

%Total energy /eV
E = -mu*(qe^4)*(Z^2)/(8*(e0^2)*(h^2)*(n^2));
E = E/qe;

%Mean radius /Angstroms
r_mean = (a/2)*(3*n^2 -L*(L+1));

%Mean square radius /Angstrom^2
r2_mean = (a^2)*(n^2)*(1/2)*( 5*n^2 + 1 - 3*L*(L+1) );

%Compute scaled radius
x = 2*r/(a*n);

%Compute wavefunction
w1 = sqrt( factorial(n-L-1)/(2*n*factorial(n+L)) );
w2 = (2/(a*n))^(3/2);
w3 = (x.^L).*exp(-x/2);
w4 = laguerre(x,n,L);
w = w1.*w2.*w3.*w4;

%%

function test_Hradial

%Use carbon as an example
r = linspace(0, 10, 100);
Z = 6;
A = 12.011
n = 4;
L = 2;
[w,E,r_mean,r2_mean] = Hradial(r,n,L,Z,A);
figure('color',[1 1 1]);
plot( r, w.*conj(w) );
xlabel('radius /angstrom')
ylabel('Probability density')
title(['Hydrogenic atom: Z=',num2str(Z),...
    ', A=',num2str(A),', n=',num2str(n),', L=',num2str(L)])

%End of code