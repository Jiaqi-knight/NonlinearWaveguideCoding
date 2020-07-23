function[]=save3Dtocsv(baseoutfilename,folder,N,Nx,Ny,Nz,dt,tmax,t,i,flags,lattice,vel,m,Fel,Fdamp,Fext,Ekinpdof,Ekinp,Ekin,Epotp,Epot,Dfuncp,Dfunc,covariantbase,contravariantbase,metriccoefficients,g,sqrtg,reciprocalmetriccoefficients,firstChristoffelsymbol,secondChristoffelsymbol,Riemanntensor,Riccitensor,R,rho,u,p,S,gammadot,sigma,mu,tau,cs,Q,feq,f,defgrad,rightCauGretens,EulLagstrain,Cauchystress)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 2nd, 2014
%    Last update: July 31st, 2014
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
%                 flags(38) --> save deformation gradient (LSM)
%                 flags(39) --> save right Cauchy-Green tensor (LSM)
%                 flags(40) --> save Euler-Lagrange strain tensor (LSM)
%                 flags(41) --> save Cauchy stress tensor (LSM)


%%

if i~=0
    idigits = fix(abs(log10(abs(i))))+1;
else
    idigits = 1;
end
if itmax~=0
    itmaxdigits = fix(abs(log10(abs(itmax))))+1;
else
    itmaxdigits = 1;
end
if idigits<itmaxdigits
    itstring = num2str(i);
    diffdigits = itmaxdigits - idigits;
    for i=1:diffdigits
        itstring = strcat('0',itstring);
    end
end

filename = strcat(folder,'/',baseoutfilename,'_N',itstring,'.csv');

dataindex = [];
data = [];
C1 = -1;
C2 = -1;
for i=1:length(flags)
    if flags(i)
        R1 = 6;
        R2 = N+5;
        C1 = C2 + 1;
        switch i
            case 1
                C2 = C2 + 3;
                data = [data lattice(:,1:3)];
            case 2
                C2 = C2 + 3;
                data = [data lattice(:,7:9)];
            case 3
                C2 = C2 + 3;
                data = [data vel];
            case 4
                C2 = C2 + 1;
                data = [data m];
            case 5
                C2 = C2 + 3;
                data = [data Fel];
            case 6
                C2 = C2 + 3;
                data = [data Fdamp];
            case 7
                C2 = C2 + 3;
                data = [data Fext];
            case 8
                C2 = C2 + 3;
                data = [data Fel+Fdamp+Fext];
            case 9
                C2 = C2 + 3;
                data = [data Ekinpdof];
            case 10
                C2 = C2 + 1;
                data = [data Ekinp];
            case 11
                C2 = C2 + 1;
                data = [data Ekin*ones(N,1)];
            case 12
                C2 = C2 + 1;
                data = [data Epotp];
            case 13
                C2 = C2 + 1;
                data = [data Epot*ones(N,1)];
            case 14
                C2 = C2 + 1;
                data = [data Dfuncp];
            case 15
                C2 = C2 + 1;
                data = [data Dfunc*ones(N,1)];
            case 16
                C2 = C2 + 3;
                data = [data covariantbase];
            case 17
                C2 = C2 + 3;
                data = [data contravariantbase];
            case 18
                C2 = C2 + 6;
                data = [data metriccoefficients];
            case 19
                C2 = C2 + 1;
                data = [data g];
            case 20
                C2 = C2 + 1;
                data = [data sqrtg];
            case 21
                C2 = C2 + 6;
                data = [data reciprocalmetriccoefficients];
            case 22
                C2 = C2 + 27;
                data = [data firstChristoffelsymbol];
            case 23
                C2 = C2 + 27;
                data = [data secondChristoffelsymbol];
            case 24
                C2 = C2 + 81;
                data = [data Riemanntensor];
            case 25
                C2 = C2 + 6;
                data = [data Riccitensor];
            case 26
                C2 = C2 + 1;
                data = [data R];
            case 27
                C2 = C2 + 1;
                data = [data rho];
            case 28
                C2 = C2 + 3;
                data = [data u];
            case 29
                C2 = C2 + 1;
                data = [data p];
            case 30
                C2 = C2 + 6;
                data = [data S];
            case 31
                C2 = C2 + 1;
                data = [data gammadot];
            case 32
                C2 = C2 + 6;
                data = [data sigma];
            case 33
                C2 = C2 + 1;
                data = [data mu];
            case 34
                C2 = C2 + 1;
                data = [data tau];
            case 35
                C2 = C2 + 1;
                data = [data cs];
            case 36
                C2 = C2 + Q;
                data = [data feq];
            case 37
                C2 = C2 + Q;
                data = [data f];
            case 38
                C2 = C2 + 9;
                data = [data defgrad];
            case 39
                C2 = C2 + 9;
                data = [data rightCauGretens];
            case 40
                C2 = C2 + 9;
                data = [data EulLagstrain];
            case 41
                C2 = C2 + 9;
                data = [data Cauchystress];
        end
        dataindex = [dataindex [i;R1;C1;R2;C2]];
    end
end

nvars = size(dataindex,2);

basicdata = [N Nx Ny Nz dt tmax t i nvars];

fid = fopen(filename,'w');

dlmwrite(filename,basicdata,'-append');
dlmwrite(filename,dataindex,'-append');
dlmwrite(filename,data,'-append');

fclose(fid);

return