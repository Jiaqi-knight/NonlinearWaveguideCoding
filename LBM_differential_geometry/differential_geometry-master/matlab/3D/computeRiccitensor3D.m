function[Riccitensor]=computeRiccitensor3D(Riemanntensor,reciprocalmetriccoefficients)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 29th, 2014
%    Last update: May 29th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

G11 = reciprocalmetriccoefficients(:,1);
G22 = reciprocalmetriccoefficients(:,2);
G33 = reciprocalmetriccoefficients(:,3);
G12 = reciprocalmetriccoefficients(:,4);
G13 = reciprocalmetriccoefficients(:,5);
G23 = reciprocalmetriccoefficients(:,6);

R1111 = Riemanntensor(:,1);
R1112 = Riemanntensor(:,2);
R1113 = Riemanntensor(:,3);
R1121 = Riemanntensor(:,4);
R1122 = Riemanntensor(:,5);
R1123 = Riemanntensor(:,6);
R1131 = Riemanntensor(:,7);
R1132 = Riemanntensor(:,8);
R1133 = Riemanntensor(:,9);
R1211 = Riemanntensor(:,10);
R1212 = Riemanntensor(:,11);
R1213 = Riemanntensor(:,12);
R1221 = Riemanntensor(:,13);
R1222 = Riemanntensor(:,14);
R1223 = Riemanntensor(:,15);
R1231 = Riemanntensor(:,16);
R1232 = Riemanntensor(:,17);
R1233 = Riemanntensor(:,18);
R1311 = Riemanntensor(:,19);
R1312 = Riemanntensor(:,20);
R1313 = Riemanntensor(:,21);
R1321 = Riemanntensor(:,22);
R1322 = Riemanntensor(:,23);
R1323 = Riemanntensor(:,24);
R1331 = Riemanntensor(:,25);
R1332 = Riemanntensor(:,26);
R1333 = Riemanntensor(:,27);
%R2111 = Riemanntensor(:,28);
%R2112 = Riemanntensor(:,29);
%R2113 = Riemanntensor(:,30);
R2121 = Riemanntensor(:,31);
R2122 = Riemanntensor(:,32);
R2123 = Riemanntensor(:,33);
R2131 = Riemanntensor(:,34);
R2132 = Riemanntensor(:,35);
R2133 = Riemanntensor(:,36);
%R2211 = Riemanntensor(:,37);
%R2212 = Riemanntensor(:,38);
%R2213 = Riemanntensor(:,39);
R2221 = Riemanntensor(:,40);
R2222 = Riemanntensor(:,41);
R2223 = Riemanntensor(:,42);
R2231 = Riemanntensor(:,43);
R2232 = Riemanntensor(:,44);
R2233 = Riemanntensor(:,45);
%R2311 = Riemanntensor(:,46);
%R2312 = Riemanntensor(:,47);
%R2313 = Riemanntensor(:,48);
R2321 = Riemanntensor(:,49);
R2322 = Riemanntensor(:,50);
R2323 = Riemanntensor(:,51);
R2331 = Riemanntensor(:,52);
R2332 = Riemanntensor(:,53);
R2333 = Riemanntensor(:,54);
%R3111 = Riemanntensor(:,55);
%R3112 = Riemanntensor(:,56);
%R3113 = Riemanntensor(:,57);
%R3121 = Riemanntensor(:,58);
%R3122 = Riemanntensor(:,59);
%R3123 = Riemanntensor(:,60);
R3131 = Riemanntensor(:,61);
R3132 = Riemanntensor(:,62);
R3133 = Riemanntensor(:,63);
%R3211 = Riemanntensor(:,64);
%R3212 = Riemanntensor(:,65);
%R3213 = Riemanntensor(:,66);
%R3221 = Riemanntensor(:,67);
%R3222 = Riemanntensor(:,68);
%R3223 = Riemanntensor(:,69);
R3231 = Riemanntensor(:,70);
R3232 = Riemanntensor(:,71);
R3233 = Riemanntensor(:,72);
%R3311 = Riemanntensor(:,73);
%R3312 = Riemanntensor(:,74);
%R3313 = Riemanntensor(:,75);
%R3321 = Riemanntensor(:,76);
%R3322 = Riemanntensor(:,77);
%R3323 = Riemanntensor(:,78);
R3331 = Riemanntensor(:,79);
R3332 = Riemanntensor(:,80);
R3333 = Riemanntensor(:,81);

% R_ij = R_ji
R11 = G11.*R1111 + G12.*R1112 + G13.*R1113 + G12.*R1211 + G22.*R1212 + G23.*R1213 + G13.*R1311 + G23.*R1312 + G33.*R1313;
R12 = G11.*R1121 + G12.*R1122 + G13.*R1123 + G12.*R1221 + G22.*R1222 + G23.*R1223 + G13.*R1321 + G23.*R1322 + G33.*R1323;
R13 = G11.*R1131 + G12.*R1132 + G13.*R1133 + G12.*R1231 + G22.*R1232 + G23.*R1233 + G13.*R1331 + G23.*R1332 + G33.*R1333;
R22 = G11.*R2121 + G12.*R2122 + G13.*R2123 + G12.*R2221 + G22.*R2222 + G23.*R2223 + G13.*R2321 + G23.*R2322 + G33.*R2323;
R23 = G11.*R2131 + G12.*R2132 + G13.*R2133 + G12.*R2231 + G22.*R2232 + G23.*R2233 + G13.*R2331 + G23.*R2332 + G33.*R2333;
R33 = G11.*R3131 + G12.*R3132 + G13.*R3133 + G12.*R3231 + G22.*R3232 + G23.*R3233 + G13.*R3331 + G23.*R3332 + G33.*R3333;

Riccitensor = [R11 R22 R33 R12 R13 R23];

return