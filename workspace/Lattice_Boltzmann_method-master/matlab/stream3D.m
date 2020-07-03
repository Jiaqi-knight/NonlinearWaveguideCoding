function[f]=stream3D(N,Q,scheme,fold,lattice)

%%
%        Project: Duct Acoustic xx curvilinear coordination by LBM
%         Author: Jiaqi Wang
%    Institution: Shanghai Jiaotong University
%                 Sound and Vibration insitition
% Research group:
%        Version: 0.1
%  Creation date: July 3rd, 2020
%    Last update: J
%
%    Description: 圆柱-》Vertification code for g_ij, & christoffel symbol
%    R(u,n)=[sin(u)*cos(n), sin(u)*sin(n),cos(u)]
%    g=[dR/du*dR/du dR/du*dR/dn;
%       dR/dn*dR/du dR/dn*dR/dn;]
%    gamma_bc^a=1/2*g^ad*(g_bd,c+g_cd,c-g_bc,d)
%
%          Input:
%         Output:

%%

Nr=max(lattice(:,1))+1;
Nt=max(lattice(:,2))+1;
Nz=max(lattice(:,3))+1;

for k=1:N
    for kk=1:Q
        xx = lattice(k,1) + scheme(kk,1);
        yy = lattice(k,2) + scheme(kk,2);
        zz = lattice(k,3) + scheme(kk,3);
        if (xx > Nr-1 || xx < 0)%全周期性边界条件
            xx = mod(xx+Nr,Nr);
        end
        if (yy > Nt-1 || yy < 0)
            yy =mod(yy+Nt,Nt);
        end
        if (zz > Nz-1 || zz < 0)
            zz = mod(zz+Nz,Nz);
        end
        id=xx+1+(yy)*Nr+(zz)*Nr*Nt;
        f(id,kk)=fold(k,kk);
    end
end
end
