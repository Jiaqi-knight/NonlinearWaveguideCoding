function[data]=loadfromcsv(inpfilename,folder)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 2nd, 2014
%    Last update: July 9th, 2014
%
%    Description: 
%          Input: 
%         Output: 
%
%                 Ordering of flags:
%                 flags(1)  --> save fixed lattice (computational domain)
%                 flags(2)  --> save moving mesh (physical domain)
%                 flags(3)  --> save velocities (LSM)
%                 flags(4)  --> save masses (LSM)
%                 flags(5)  --> save elastic forces (LSM)
%                 flags(6)  --> save damping forces (LSM)
%                 flags(7)  --> save external forces (LSM)
%                 flags(8)  --> save total forces (LSM)
%                 flags(9)  --> save particle kinetic energy per dof (LSM)
%                 flags(10) --> save particle kinetic energy (LSM)
%                 flags(11) --> save total kinetic energy (LSM)
%                 flags(12) --> save particle potential energy (LSM)
%                 flags(13) --> save total potential energy (LSM)
%                 flags(14) --> save particle dissipation function (LSM)
%                 flags(15) --> save total dissipation function (LSM)
%                 flags(16) --> save covariant base vectors
%                 flags(17) --> save contravariant base vectors
%                 flags(18) --> save metric tensor
%                 flags(19) --> save metric tensor determinant
%                 flags(20) --> save g
%                 flags(21) --> save reciprocal metric coefficients
%                 flags(22) --> save first Christoffel symbols
%                 flags(23) --> save second Christoffel symbols
%                 flags(24) --> save Riemann tensor
%                 flags(25) --> save Ricci tensor
%                 flags(26) --> save Ricci curvature
%                 flags(27) --> save density (LBM)
%                 flags(28) --> save velocity in cartesian components (LBM)
%                 flags(29) --> save pressure (LBM)
%                 flags(30) --> save strain rate tensor (LBM)
%                 flags(31) --> save strain rate (LBM)
%                 flags(32) --> save stress tensor (LBM)
%                 flags(33) --> save viscosity (LBM)
%                 flags(34) --> save relaxation time (LBM)
%                 flags(35) --> save sound velocity (LBM)
%                 flags(36) --> save particle equilibrium populations (LBM)
%                 flags(37) --> save particle populations (LBM)

%%

filename = strcat(folder,'/',inpfilename,'.csv');

R1 = 0;
C1 = 0;
R2 = 0;
C2 = 9-1;
RNG = [R1 C1 R2 C2];
basicdata = csvread(filename,R1,C1,RNG);

data.N = basicdata(1,1);
data.Nx = basicdata(1,2);
data.Ny = basicdata(1,3);
data.Nz = basicdata(1,4);
data.dt = basicdata(1,5);
data.tmax = basicdata(1,6);
data.t = basicdata(1,7);
data.i = basicdata(1,8);
data.nvars = basicdata(1,9);

R1 = 1;
C1 = 0;
R2 = 1;
C2 = nvars-1;
RNG = [R1 C1 R2 C2];
dataindex = csvread(filename,R1,C1,RNG);

for i=1:nvars
    switch dataindex(1,i)
        case 1
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.lattice = csvread(filename,R1,C1,RNG);
        case 2
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.mesh = csvread(filename,R1,C1,RNG);
        case 3
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.vel = csvread(filename,R1,C1,RNG);
        case 4
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.m = csvread(filename,R1,C1,RNG);
        case 5
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Fel = csvread(filename,R1,C1,RNG);
        case 6
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Fdamp = csvread(filename,R1,C1,RNG);
        case 7
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Fext = csvread(filename,R1,C1,RNG);
        case 8
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Ftot = csvread(filename,R1,C1,RNG);
        case 9
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Ekinpdof = csvread(filename,R1,C1,RNG);
        case 10
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Ekinp = csvread(filename,R1,C1,RNG);
        case 11
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            Ekin = csvread(filename,R1,C1,RNG);
            data.Ekin = Ekin(1);
        case 12
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Epotp = csvread(filename,R1,C1,RNG);
        case 13
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            Epot = csvread(filename,R1,C1,RNG);
            data.Epot = Epot(1);
        case 14
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Dfuncp = csvread(filename,R1,C1,RNG);
        case 15
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            Dfunc = csvread(filename,R1,C1,RNG);
            data.Dfunc = Dfunc(1);
        case 16
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.covariantbase = csvread(filename,R1,C1,RNG);
        case 17
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.contravariantbase = csvread(filename,R1,C1,RNG);
        case 18
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.metriccoefficients = csvread(filename,R1,C1,RNG);
        case 19
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.J = csvread(filename,R1,C1,RNG);
        case 20
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.g = csvread(filename,R1,C1,RNG);
        case 21
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.reciprocalmetriccoefficients = csvread(filename,R1,C1,RNG);
        case 22
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.firstChristoffelsymbol = csvread(filename,R1,C1,RNG);
        case 23
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.secondChristoffelsymbol = csvread(filename,R1,C1,RNG);
        case 24
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Riemanntensor = csvread(filename,R1,C1,RNG);
        case 25
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.Riccitensor = csvread(filename,R1,C1,RNG);
        case 26
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.R = csvread(filename,R1,C1,RNG);
        case 27
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.rho = csvread(filename,R1,C1,RNG);
        case 28
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.u = csvread(filename,R1,C1,RNG);
        case 29
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.p = csvread(filename,R1,C1,RNG);
        case 30
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.S = csvread(filename,R1,C1,RNG);
        case 31
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.gammadot = csvread(filename,R1,C1,RNG);
        case 32
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.sigma = csvread(filename,R1,C1,RNG);
        case 33
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.mu = csvread(filename,R1,C1,RNG);
        case 34
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.tau = csvread(filename,R1,C1,RNG);
        case 35
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.cs = csvread(filename,R1,C1,RNG);
        case 36
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.feq = csvread(filename,R1,C1,RNG);
        case 37
            R1 = dataindex(2,i);
            C1 = dataindex(3,i);
            R2 = dataindex(4,i);
            C2 = dataindex(5,i);
            RNG = [R1 C1 R2 C2];
            data.f = csvread(filename,R1,C1,RNG);
    end
end

return