function[circles]=intersectingcirclessection2D(x1,y1,R1,x2,y2,R2,R3,Neta1_1,Neta1_2,Neta1_3,Neta2,eta1min,eta2min,deltaeta1,deltaeta2,xmin,xmax,ymin,ymax,itmax,tol)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 16th, 2014
%    Last update: August 8th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

d = sqrt((x1-x2)^2+(y1-y2)^2);

delta = 0.5*d+0.5*(R2^2-R1^2)/d;
xa = x2 + delta;
%xb = xa;

ya = y2 + sqrt(R2^2-delta^2);
%yb = y2 - sqrt(R2^2-delta^2);

yc = ya;
err = 1;
it = 0;
tolNewton = 10^-14;
itmaxNewton = 100;
while err>=tolNewton && it<itmaxNewton
    f = sqrt((R1+R3)^2-yc^2) + sqrt((R2+R3)^2-yc^2) -d;
    df = -yc*(1/sqrt((R1+R3)^2-yc^2) + 1/sqrt((R2+R3)^2-yc^2));
    yc = yc - f/df;
    it = it + 1;
    err = abs(sqrt((R1+R3)^2-yc^2) + sqrt((R2+R3)^2-yc^2) -d);
end

yc = real(yc);
%yd = -yc;

xc = x2 + sqrt((R2+R3)^2-yc^2);
%xd = xc;

beta = asin(yc/(R1+R3));
gamma = asin(yc/(R2+R3));

%xq = x1 - R1*cos(beta);
%yq = y1 + R1*sin(beta);

%xp = x2 + R2*cos(gamma);
%yp = y2 + R2*sin(gamma);

base1 = sqrt(R2^2-ya^2);
base2 = sqrt(R1^2-ya^2);

%phi = atan(ya/base1);
%psi = atan(ya/base2);

Neta1 = Neta1_1 + Neta1_2 + Neta1_3;
eta1max = eta1min + Neta1*deltaeta1;
eta2max = eta2min + Neta2*deltaeta2;

N = Neta1*Neta2;

Nx = Neta1;
Ny = Neta2;

c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);

c1(1,1) = x2;
c1(1,2) = y2-R2;

c2(1,1) = x1;
c2(1,2) = y1-R1;

c3(1,1) = x1;
c3(1,2) = y1+R1;

c4(1,1) = x2;
c4(1,2) = y2+R2;

theta = (-0.5*pi:pi/(Neta2-1):0.5*pi)';
e4(:,1) = x1 + R1*cos(theta);
e4(:,2) = y1 + R1*sin(theta);

theta = (1.5*pi:-pi/(Neta2-1):0.5*pi)';
e2(:,1) = x2 + R2*cos(theta);
e2(:,2) = y2 + R2*sin(theta);

%gammacomp = 0.5*pi - gamma;
betacomp = 0.5*pi - beta;

theta1 = (0.5*pi:(gamma-0.5*pi)/(Neta1_1-1):gamma)';
theta2 = (pi+gamma:((1.5*pi+betacomp)-(pi+gamma))/(Neta1_2-1):1.5*pi+betacomp)';
theta3 = (pi-beta:-(0.5*pi-beta)/(Neta1_3-1):0.5*pi)';
e3(:,1) = [x2+R2.*cos(theta1); xc+R3.*cos(theta2); x1+R1.*cos(theta3)];
e3(:,2) = [y2+R2.*sin(theta1); yc+R3.*sin(theta2); y1+R1.*sin(theta3)];

e1(:,1) = e3(:,1);
e1(:,2) = -e3(:,2);

% figure();
% plot(x1,y1,'g*','LineWidth',2)
% hold on
% plot(x2,y2,'g*','LineWidth',2)
% hold on
% plot([x1;xc],[y1;yc],'b--','LineWidth',2)
% hold on
% plot([x2;xc],[y2;yc],'b--','LineWidth',2)
% hold on
% plot([x1;xq],[y1;yq],'b--','LineWidth',2)
% hold on
% plot([x2;xp],[y2;yp],'b--','LineWidth',2)
% hold on
% plot([x2;x1],[y2;y1],'b--','LineWidth',2)
% hold on
% plot([xa;xa],[yc+5;yd-5],'k--','LineWidth',2)
% hold on
% plot([x2-R2-5;x1+R1+5],[y2;y1],'k--','LineWidth',2)
% hold on
% plot(xa,ya,'g*','LineWidth',2)
% hold on
% plot(xb,yb,'g*','LineWidth',2)
% hold on
% plot(xc,yc,'g*','LineWidth',2)
% hold on
% plot(xd,yd,'g*','LineWidth',2)
% hold on
% plot(xp,yp,'g*','LineWidth',2)
% hold on
% plot(xq,yq,'g*','LineWidth',2)
% hold on
% plot(c1(1,1),c1(1,2),'k*','LineWidth',2)
% hold on
% plot(c2(1,1),c2(1,2),'k*','LineWidth',2)
% hold on
% plot(c3(1,1),c3(1,2),'k*','LineWidth',2)
% hold on
% plot(c4(1,1),c4(1,2),'k*','LineWidth',2)
% hold on
% plot(e1(:,1),e1(:,2),'r*')
% hold on
% plot(e2(:,1),e2(:,2),'r*')
% hold on
% plot(e3(:,1),e3(:,2),'r*')
% e3(:,1) = [x2+R2.*cos(theta1); xc+R3.*cos(theta2); x1+R1.*cos(theta3)];
% e3(:,2) = [y2+R2.*sin(theta1); yc+R3.*sin(theta2); y1+R1.*sin(theta3)];
% % plot(x2+R2.*cos(theta1),y2+R2.*sin(theta1),'c--')
% % hold on
% % plot(xc+R3.*cos(theta2),yc+R3.*sin(theta2),'m--')
% % hold on
% % plot(x1+R1.*cos(theta3),y1+R1.*sin(theta3),'c--')
% % hold on
% plot(e4(:,1),e4(:,2),'r*')
% hold on
% grid on
% xlabel('x')
% ylabel('y')
% title('Mesh in fluid domain (physical space)')
% axis equal

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

circles = transfiniteinterpolation2D(N,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

if itmax~=0
    [indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Nx,Ny);

    flagperiodicity = 0;
    periodicity = 0;
    flagintbounds = 0;
    indicesintbounds = 0;
    typeintbounds = 0;
    deltaq = [deltaeta1 deltaeta2];

    [structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(N,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

    circles = sparseellipticgridgen2D(Nx,N,circles,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0);
end

return
