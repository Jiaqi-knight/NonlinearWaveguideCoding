function[Riemanntensor]=computeRiemanntensor3Dsurface(N,deltaq,secondChristoffelsymbol,metriccoefficients,firstdevneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 18th, 2014
%    Last update: July 18th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

g11 = metriccoefficients(:,1);
g22 = metriccoefficients(:,2);
g12 = metriccoefficients(:,3);

Gamma111 = secondChristoffelsymbol(:,1);
Gamma112 = secondChristoffelsymbol(:,2);
Gamma121 = secondChristoffelsymbol(:,3);
Gamma122 = secondChristoffelsymbol(:,4);
Gamma211 = secondChristoffelsymbol(:,5);
Gamma212 = secondChristoffelsymbol(:,6);
Gamma221 = secondChristoffelsymbol(:,7);
Gamma222 = secondChristoffelsymbol(:,8);

Gamma111d1 = zeros(N,1);
Gamma112d1 = zeros(N,1);
Gamma121d1 = zeros(N,1);
Gamma122d1 = zeros(N,1);
Gamma211d1 = zeros(N,1);
Gamma212d1 = zeros(N,1);
Gamma221d1 = zeros(N,1);
Gamma222d1 = zeros(N,1);

Gamma111d2 = zeros(N,1);
Gamma112d2 = zeros(N,1);
Gamma121d2 = zeros(N,1);
Gamma122d2 = zeros(N,1);
Gamma211d2 = zeros(N,1);
Gamma212d2 = zeros(N,1);
Gamma221d2 = zeros(N,1);
Gamma222d2 = zeros(N,1);

for i=1:N
    j1 = 1;
    switch firstdevneighbours(i,4*(j1-1)+1)
        case 1
            Gamma111d1(i,1) = 0.5.*(Gamma111(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma111(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
            Gamma112d1(i,1) = 0.5.*(Gamma112(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma112(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
            Gamma121d1(i,1) = 0.5.*(Gamma121(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma121(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
            Gamma122d1(i,1) = 0.5.*(Gamma122(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma122(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
            Gamma211d1(i,1) = 0.5.*(Gamma211(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma211(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
            Gamma212d1(i,1) = 0.5.*(Gamma212(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma212(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);           
            Gamma221d1(i,1) = 0.5.*(Gamma221(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma221(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
            Gamma222d1(i,1) = 0.5.*(Gamma222(firstdevneighbours(i,4*(j1-1)+3),1)-Gamma222(firstdevneighbours(i,4*(j1-1)+2),1))./deltaq(j1);
        case 2
            Gamma111d1(i,1) = (-1.5*Gamma111(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma111(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma111(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma112d1(i,1) = (-1.5*Gamma112(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma112(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma112(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma121d1(i,1) = (-1.5*Gamma121(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma121(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma121(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma122d1(i,1) = (-1.5*Gamma122(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma122(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma122(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma211d1(i,1) = (-1.5*Gamma211(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma211(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma211(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma212d1(i,1) = (-1.5*Gamma212(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma212(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma212(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma221d1(i,1) = (-1.5*Gamma221(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma221(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma221(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma222d1(i,1) = (-1.5*Gamma222(firstdevneighbours(i,4*(j1-1)+2),1)+2*Gamma222(firstdevneighbours(i,4*(j1-1)+3),1)-0.5*Gamma222(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
        case 3
            Gamma111d1(i,1) = (1.5*Gamma111(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma111(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma111(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma112d1(i,1) = (1.5*Gamma112(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma112(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma112(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma121d1(i,1) = (1.5*Gamma121(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma121(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma121(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma122d1(i,1) = (1.5*Gamma122(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma122(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma122(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma211d1(i,1) = (1.5*Gamma211(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma211(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma211(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma212d1(i,1) = (1.5*Gamma212(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma212(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma212(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma221d1(i,1) = (1.5*Gamma221(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma221(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma221(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
            Gamma222d1(i,1) = (1.5*Gamma222(firstdevneighbours(i,4*(j1-1)+2),1)-2*Gamma222(firstdevneighbours(i,4*(j1-1)+3),1)+0.5*Gamma222(firstdevneighbours(i,4*(j1-1)+4),1))./deltaq(j1);
    end
    j2 = 2;
    switch firstdevneighbours(i,4*(j2-1)+1)
        case 1
            Gamma111d2(i,1) = 0.5.*(Gamma111(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma111(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma112d2(i,1) = 0.5.*(Gamma112(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma112(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma121d2(i,1) = 0.5.*(Gamma121(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma121(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma122d2(i,1) = 0.5.*(Gamma122(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma122(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma211d2(i,1) = 0.5.*(Gamma211(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma211(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma212d2(i,1) = 0.5.*(Gamma212(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma212(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma221d2(i,1) = 0.5.*(Gamma221(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma221(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
            Gamma222d2(i,1) = 0.5.*(Gamma222(firstdevneighbours(i,4*(j2-1)+3),1)-Gamma222(firstdevneighbours(i,4*(j2-1)+2),1))./deltaq(j2);
        case 2
            Gamma111d2(i,1) = (-1.5*Gamma111(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma111(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma111(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma112d2(i,1) = (-1.5*Gamma112(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma112(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma112(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma121d2(i,1) = (-1.5*Gamma121(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma121(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma121(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma122d2(i,1) = (-1.5*Gamma122(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma122(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma122(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma211d2(i,1) = (-1.5*Gamma211(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma211(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma211(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma212d2(i,1) = (-1.5*Gamma212(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma212(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma212(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma221d2(i,1) = (-1.5*Gamma221(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma221(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma221(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma222d2(i,1) = (-1.5*Gamma222(firstdevneighbours(i,4*(j2-1)+2),1)+2*Gamma222(firstdevneighbours(i,4*(j2-1)+3),1)-0.5*Gamma222(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
        case 3
            Gamma111d2(i,1) = (1.5*Gamma111(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma111(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma111(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma112d2(i,1) = (1.5*Gamma112(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma112(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma112(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma121d2(i,1) = (1.5*Gamma121(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma121(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma121(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma122d2(i,1) = (1.5*Gamma122(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma122(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma122(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma211d2(i,1) = (1.5*Gamma211(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma211(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma211(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma212d2(i,1) = (1.5*Gamma212(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma212(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma212(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma221d2(i,1) = (1.5*Gamma221(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma221(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma221(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
            Gamma222d2(i,1) = (1.5*Gamma222(firstdevneighbours(i,4*(j2-1)+2),1)-2*Gamma222(firstdevneighbours(i,4*(j2-1)+3),1)+0.5*Gamma222(firstdevneighbours(i,4*(j2-1)+4),1))./deltaq(j2);
    end
end

% Rlijk = Gamma_lik/j - Gamma_lij/k + Gamma_ljs*Gamma_sik - Gamma_lks*Gamma_sij

R1111 = Gamma111d1 - Gamma111d1 + (Gamma111.*Gamma111 + Gamma112.*Gamma211) - (Gamma111.*Gamma111 + Gamma112.*Gamma211);
R1112 = Gamma112d1 - Gamma111d2 + (Gamma111.*Gamma112 + Gamma112.*Gamma212) - (Gamma121.*Gamma111 + Gamma122.*Gamma211);
R1121 = Gamma111d2 - Gamma112d1 + (Gamma121.*Gamma111 + Gamma122.*Gamma211) - (Gamma111.*Gamma112 + Gamma112.*Gamma212);
R1122 = Gamma112d2 - Gamma112d2 + (Gamma121.*Gamma112 + Gamma122.*Gamma212) - (Gamma121.*Gamma112 + Gamma122.*Gamma212);
R1211 = Gamma121d1 - Gamma121d1 + (Gamma111.*Gamma121 + Gamma112.*Gamma221) - (Gamma111.*Gamma121 + Gamma112.*Gamma221);
R1212 = Gamma122d1 - Gamma121d2 + (Gamma111.*Gamma122 + Gamma112.*Gamma222) - (Gamma121.*Gamma121 + Gamma122.*Gamma221);
R1221 = Gamma121d2 - Gamma122d1 + (Gamma121.*Gamma121 + Gamma122.*Gamma221) - (Gamma111.*Gamma122 + Gamma112.*Gamma222);
R1222 = Gamma122d2 - Gamma122d2 + (Gamma121.*Gamma122 + Gamma122.*Gamma222) - (Gamma121.*Gamma122 + Gamma122.*Gamma222);
R2111 = Gamma211d1 - Gamma211d1 + (Gamma211.*Gamma111 + Gamma212.*Gamma211) - (Gamma211.*Gamma111 + Gamma212.*Gamma211);
R2112 = Gamma212d1 - Gamma211d2 + (Gamma211.*Gamma112 + Gamma212.*Gamma212) - (Gamma221.*Gamma111 + Gamma222.*Gamma211);
R2121 = Gamma211d2 - Gamma212d1 + (Gamma221.*Gamma111 + Gamma222.*Gamma211) - (Gamma211.*Gamma112 + Gamma212.*Gamma212);
R2122 = Gamma212d2 - Gamma212d2 + (Gamma221.*Gamma112 + Gamma222.*Gamma212) - (Gamma221.*Gamma112 + Gamma222.*Gamma212);
R2211 = Gamma221d1 - Gamma221d1 + (Gamma211.*Gamma121 + Gamma212.*Gamma221) - (Gamma211.*Gamma121 + Gamma212.*Gamma221);
R2212 = Gamma222d1 - Gamma221d2 + (Gamma211.*Gamma122 + Gamma212.*Gamma222) - (Gamma221.*Gamma121 + Gamma222.*Gamma221);
R2221 = Gamma221d2 - Gamma222d1 + (Gamma221.*Gamma121 + Gamma222.*Gamma221) - (Gamma211.*Gamma122 + Gamma212.*Gamma222);
R2222 = Gamma222d2 - Gamma222d2 + (Gamma221.*Gamma122 + Gamma222.*Gamma222) - (Gamma221.*Gamma122 + Gamma222.*Gamma222);

% Riemanntensor_lijk = g_ls*R_sijk
Riemanntensor1111 = g11.*R1111 + g12.*R2111;
Riemanntensor1112 = g11.*R1112 + g12.*R2112;
Riemanntensor1121 = g11.*R1121 + g12.*R2121;
Riemanntensor1122 = g11.*R1122 + g12.*R2122;
Riemanntensor1211 = g11.*R1211 + g12.*R2211;
Riemanntensor1212 = g11.*R1212 + g12.*R2212;
Riemanntensor1221 = g11.*R1221 + g12.*R2221;
Riemanntensor1222 = g11.*R1222 + g12.*R2222;
Riemanntensor2111 = g21.*R1111 + g22.*R2111;
Riemanntensor2112 = g21.*R1112 + g22.*R2112;
Riemanntensor2121 = g21.*R1121 + g22.*R2121;
Riemanntensor2122 = g21.*R1122 + g22.*R2122;
Riemanntensor2211 = g21.*R1211 + g22.*R2211;
Riemanntensor2212 = g21.*R1212 + g22.*R2212;
Riemanntensor2221 = g21.*R1221 + g22.*R2221;
Riemanntensor2222 = g21.*R1222 + g22.*R2222;

Riemanntensor = [Riemanntensor1111 Riemanntensor1112 Riemanntensor1121 Riemanntensor1122  Riemanntensor1211 Riemanntensor1212 Riemanntensor1221 Riemanntensor1222 Riemanntensor2111 Riemanntensor2112 Riemanntensor2121 Riemanntensor2122 Riemanntensor2211 Riemanntensor2212 Riemanntensor2221 Riemanntensor2222];

return