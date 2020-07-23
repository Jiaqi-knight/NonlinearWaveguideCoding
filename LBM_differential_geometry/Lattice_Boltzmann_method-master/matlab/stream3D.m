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
f=reshape(fold,[Nr,Nt,Nz,Q]);
f(:,       :,       :,      1)   =  f(  :,      :,       :,       1); %  0, 0, 0;
f(:,       :,      3:Nz,    2)   =  f(  :,      :,      2:Nz-1,   2); %  0, 0, 1;
f(:,       :,      1:Nz-2,  3)   =  f(  :,      :,      2:Nz-1,   3); %  0, 0,-1;
f(:,      3:Nt,     :,      4)   =  f(  :,     2:Nt-1,   :,       4); %  0, 1, 0;
f(:,      1:Nt-2,   :,      5)   =  f(  :,     2:Nt-1,   :,       5); %  0,-1, 0;
f(3:Nr,    :,       :,      6)   =  f( 2:Nr-1,  :,       :,       6); %  1, 0, 0;
f(1:Nr-2,  :,       :,      7)   =  f( 2:Nr-1,  :,       :,       7); % -1, 0, 0;
f(:,      3:Nt,    3:Nz,    8)   =  f(  :,     2:Nt-1,  2:Nz-1,   8); %  0, 1, 1;
f(:,      3:Nt,    1:Nz-2,  9)   =  f(  :,     2:Nt-1,  2:Nz-1,   9); %  0, 1,-1;
f(:,      1:Nt-2,  3:Nz,    10)  =  f(  :,     2:Nt-1,  2:Nz-1,   10);%  0,-1, 1;
f(:,      1:Nt-2,  1:Nz-2,  11)  =  f(  :,     2:Nt-1,  2:Nz-1,   11);%  0,-1,-1;
f(3:Nr,    :,      3:Nz,    12)  =  f( 2:Nr-1,  :,      2:Nz-1,   12);%  1, 0, 1;
f(3:Nr,    :,      1:Nz-2,  13)  =  f( 2:Nr-1,  :,      2:Nz-1,   13);%  1, 0,-1;
f(3:Nr,   3:Nt,     :,      14)  =  f( 2:Nr-1,  2:Nt-1,  :,       14);%  1, 1, 0;
f(3:Nr,   1:Nt-2,   :,      15)  =  f( 2:Nr-1,  2:Nt-1,  :,       15);%  1,-1, 0;
f(1:Nr-2,  :,      3:Nz,    16)  =  f( 2:Nr-1,   :,     2:Nz-1,   16);% -1, 0, 1;
f(1:Nr-2,  :,      1:Nz-2,  17)  =  f( 2:Nr-1,   :,     2:Nz-1,   17);% -1, 0,-1;
f(1:Nr-2, 3:Nt,     :,      18)  =  f( 2:Nr-1,  2:Nt-1,  :,       18);% -1, 1, 0;
f(1:Nr-2, 1:Nt-2,   :,      19)  =  f( 2:Nr-1,  2:Nt-1,  :,       19);% -1,-1, 0;
f(3:Nr,   3:Nt,    3:Nz,    20)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   20);%  1, 1, 1;
f(3:Nr,   3:Nt,    1:Nz-2,  21)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   21);%  1, 1,-1;
f(3:Nr,   1:Nt-2,  3:Nz,    22)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   22);%  1,-1, 1;
f(3:Nr,   1:Nt-2,  1:Nz-2,  23)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   23);%  1,-1,-1;
f(1:Nr-2, 3:Nt,    3:Nz,    24)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   24);% -1, 1, 1;
f(1:Nr-2, 3:Nt,    1:Nz-2,  25)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   25);% -1, 1,-1;
f(1:Nr-2, 1:Nt-2,  3:Nz,    26)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   26);% -1,-1, 1;
f(1:Nr-2, 1:Nt-2,  1:Nz-2,  27)  =  f( 2:Nr-1,  2:Nt-1, 2:Nz-1,   27);% -1,-1,-1;
f(:,       :,      7:Nz,    28)  =  f(  :,      :,      4:Nz-3,   28);%  0, 0, 3;
f(:,       :,      1:Nz-6,  29)  =  f(  :,      :,      4:Nz-3,   29);%  0, 0,-3;
f(:,      7:Nt,     :,      30)  =  f(  :,     4:Nt-3,   :,       30);%  0, 3, 0;
f(:,      1:Nt-6,   :,      31)  =  f(  :,     4:Nt-3,   :,       31);%  0,-3, 0;
f(7:Nr,    :,       :,      32)  =  f( 4:Nr-3,  :,       :,       32);%  3, 0, 0;
f(1:Nr-6,  :,       :,      33)  =  f( 4:Nr-3,  :,       :,       33);% -3, 0, 0;
f(7:Nr,   7:Nt,    7:Nz,    34)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   34);%  3, 3, 3;
f(7:Nr,   7:Nt,    1:Nz-6,  35)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   35);%  3, 3,-3;
f(7:Nr,   1:Nt-6,  7:Nz,    36)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   36);%  3,-3, 3;
f(7:Nr,   1:Nt-6,  1:Nz-6,  37)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   37);%  3,-3,-3;
f(1:Nr-6, 7:Nt,    7:Nz,    38)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   38);% -3, 3, 3;
f(1:Nr-6, 7:Nt,    1:Nz-6,  39)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   39);% -3, 3,-3;
f(1:Nr-6, 1:Nt-6,  7:Nz,    40)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   40);% -3,-3, 3;
f(1:Nr-6, 1:Nt-6,  1:Nz-6,  41)  =  f( 4:Nr-3, 4:Nt-3,  4:Nz-3,   41);% -3,-3,-3;

f=reshape(f,[Nr*Nt*Nz,Q]);



% for k=1:N
%     for kk=1:Q
%         xx = lattice(k,1) + scheme(kk,1);
%         yy = lattice(k,2) + scheme(kk,2);
%         zz = lattice(k,3) + scheme(kk,3);
%         if (xx > Nr-1 || xx < 0)%全周期性边界条件
%             xx = mod(xx+Nr,Nr);
%         end
%         if (yy > Nt-1 || yy < 0)
%             yy =mod(yy+Nt,Nt);
%         end
%         if (zz > Nz-1 || zz < 0)
%             zz = mod(zz+Nz,Nz);
%         end
%         id=xx+1+(yy)*Nr+(zz)*Nr*Nt;
%         f(id,kk)=fold(k,kk);
%     end
% end


end
