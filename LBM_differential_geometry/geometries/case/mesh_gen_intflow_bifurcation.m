clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% DUCT WITH BIFURCATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta1min = 0;
Neta1 = 50;
deltaeta1 = 1;
eta1max = eta1min + (Neta1-1)*deltaeta1;

eta2min = 0;
Neta2 = 60;
deltaeta2 = 1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

eta3minsec1 = 0;
Neta3sec1 = 50;
deltaeta3sec1 = 1;
eta3maxsec1 = eta3minsec1 + (Neta3sec1-1)*deltaeta3sec1;

eta3minsec2 = eta3maxsec1+deltaeta3sec1;
Neta3sec2 = 50;
deltaeta3sec2 = 1;
eta3maxsec2 = eta3minsec2 + (Neta3sec2-1)*deltaeta3sec2;

eta3minsec3 = eta3maxsec2+deltaeta3sec2;
Neta3sec3 = 50;
deltaeta3sec3 = 1;
eta3maxsec3 = eta3minsec3 + (Neta3sec3-1)*deltaeta3sec3;

eta3minsec4 = eta3maxsec3+deltaeta3sec3;
Neta3sec4 = 50;
deltaeta3sec4 = 1;
eta3maxsec4 = eta3minsec4 + (Neta3sec4-1)*deltaeta3sec4;

xmin = 0;
xmax = 5;
ymin = -2.5;
ymax = 3;


x0 = 0;
y0 = 0;

r = 10;

itmax = 1;
tol = 10^-4;

s1min = 0;
s1max = 40;
Ns1 = Neta3sec1;

meanlinefunc1   = @(s) [0 0 s];
dmeanlinefunc1  = @(s) [0 0 1];
d2meanlinefunc1 = @(s) [1 0 0];

s2min = 0;
s2max = 10;
Ns2 = Neta3sec2;

meanlinefunc2   = @(s) [0 0 s];
dmeanlinefunc2  = @(s) [0 0 1];
d2meanlinefunc2 = @(s) [1 0 0];

eccentricityfunc = @(s,smin,smax,emin,emax) (emax-emin)*(s-smin)/(smax-smin) + emin;

s3min = 0;
s3max = 10;
Ns3 = Neta3sec3;

meanlinefunc3   = @(s) [0 0 s];
dmeanlinefunc3  = @(s) [0 0 1];
d2meanlinefunc3 = @(s) [1 0 0];

s4part1min = 0;
s4part1max = 25;
Ns4part1 = Neta3sec4;

meanlinefunc4part1   = @(s) [s 0 s];
dmeanlinefunc4part1  = @(s) [1 0 1];
d2meanlinefunc4part1 = @(s) [0 1 0];

s4part2min = 0;
s4part2max = 25;
Ns4part2 = Neta3sec4;

meanlinefunc4part2   = @(s) [-0.1*s 0 s];
dmeanlinefunc4part2  = @(s) [-0.1 0 1];
d2meanlinefunc4part2 = @(s) [0 -1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2D = Neta1*Neta2;
Neta3 =  Neta3sec1 + Neta3sec2 + Neta3sec3 + Neta3sec4;
Nx = Neta1;
Ny = Neta2;
Nz = Neta3sec1 + Neta3sec2 + Neta3sec3 + Neta3sec4;

flagperiodicity = 0;
periodicity = 0;

flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;

deltaq = [deltaeta1 deltaeta2];

lattice = generatelattice2D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);

% ----------> Section 1: single tube with circular-cross section

e1 = zeros(Neta1,6);
e2 = zeros(Neta2,6);
e3 = zeros(Neta1,6);
e4 = zeros(Neta2,6);
c1 = zeros(1,8);
c2 = zeros(1,8);
c3 = zeros(1,8);
c4 = zeros(1,8);

c1(1,1:2) = [x0-r/sqrt(2) y0-r/sqrt(2)];
c2(1,1:2) = [x0+r/sqrt(2) y0-r/sqrt(2)];
c3(1,1:2) = [x0+r/sqrt(2) y0+r/sqrt(2)];
c4(1,1:2) = [x0-r/sqrt(2) y0+r/sqrt(2)];

theta = (5*pi/4:0.5*pi/(Neta1-1):7*pi/4)';
e1(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
e1(:,2) = y0*ones(Neta1,1)+r.*sin(theta);

theta = (5*pi/4:-0.5*pi/(Neta2-1):3*pi/4)';
e2(:,1) = x0*ones(Neta2,1)+r.*cos(theta);
e2(:,2) = y0*ones(Neta2,1)+r.*sin(theta);

theta = (3*pi/4:-0.5*pi/(Neta1-1):0.25*pi)';
e3(:,1) = x0*ones(Neta1,1)+r.*cos(theta);
e3(:,2) = y0*ones(Neta1,1)+r.*sin(theta);

theta = (-0.25*pi:0.5*pi/(Neta2-1):0.25*pi)';
e4(:,1) = x0*ones(Neta2,1)+r.*cos(theta);
e4(:,2) = y0*ones(Neta2,1)+r.*sin(theta);

circle = transfiniteinterpolation2D(logfullfile,N2D,lattice(:,1:2),eta1min,eta1max,eta2min,eta2max,Neta1,Neta2,1,e1,e2,e3,e4,c1,c2,c3,c4);

[indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(logfullfile,Nx,Ny);

[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(logfullfile,N2D,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);

boundary = [indicesE1;indicesE2;indicesE3;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4];

if itmax~=0
    circle = sparseellipticgridgen2D(logfullfile,Nx,N2D,circle,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 
end

s1 = (s1min:(s1max-s1min)/(Ns1-1):s1max)';
meanline1 = zeros(Ns1,3);
T1 = zeros(Ns1,3);
N1 = zeros(Ns1,3);
B1 = zeros(Ns1,3);

for i = 1:Ns1
    meanline1(i,1:3) = meanlinefunc1(s1(i));
    T = dmeanlinefunc1(s1(i));
    N = d2meanlinefunc1(s1(i));
    Tmod = sqrt(sum(T.^2,2));
    Nmod = sqrt(sum(N.^2,2));
    if Tmod~=0
        T1(i,:) = T/Tmod;
    end
    if Nmod~=0
        N1(i,:) = N/Nmod;
    end
    vecprod = [T1(i,2)*N1(i,3)-T1(i,3)*N1(i,2) T1(i,3)*N1(i,1)-T1(i,1)*N1(i,3) T1(i,1)*N1(i,2)-T1(i,2)*N1(i,1)];
    Bmod = sqrt(sum(vecprod.^2,2));
    if Bmod~=0
        B1(i,:) = vecprod/Bmod;
    end
end

lattice3Dsec1 = zeros(N2D*Neta3sec1,9);

for i=1:Neta3sec1
    lattice3Dsec1((i-1)*N2D+1:i*N2D,1:2) = circle(:,1:2);
    lattice3Dsec1((i-1)*N2D+1:i*N2D,3) = eta3minsec1 + deltaeta3sec1*(i-1);
    lattice3Dsec1((i-1)*N2D+1:i*N2D,4:6) = [meanline1(i,1).*ones(N2D,1) meanline1(i,2).*ones(N2D,1) meanline1(i,3).*ones(N2D,1)] + [circle(:,5) circle(:,5) circle(:,5)].*[N1(i,1).*ones(N2D,1) N1(i,2).*ones(N2D,1) N1(i,3).*ones(N2D,1)] + [circle(:,6) circle(:,6) circle(:,6)].*[B1(i,1).*ones(N2D,1) B1(i,2).*ones(N2D,1) B1(i,3).*ones(N2D,1)];
    lattice3Dsec1((i-1)*N2D+1:i*N2D,7:9) = lattice3Dsec1((i-1)*N2D+1:i*N2D,4:6);
end

% ---------> Section 2: single tube with elliptical cross-section

x0sec2 = meanline1(end,1);
y0sec2 = meanline1(end,2);
z0sec2 = meanline1(end,3);

s2 = (s2min:(s2max-s2min)/(Ns2-1):s2max)';
meanline2 = zeros(Ns2,3);
T2 = zeros(Ns2,3);
N2 = zeros(Ns2,3);
B2 = zeros(Ns2,3);

eccentricity = zeros(Ns2,1);

for i = 1:Ns2
    meanline2(i,1:3) = meanlinefunc2(s2(i));
    eccentricity(i,1) = eccentricityfunc(s2(i),s2min,s2max,0,0.3);
    T = dmeanlinefunc2(s2(i));
    N = d2meanlinefunc2(s2(i));
    Tmod = sqrt(sum(T.^2,2));
    Nmod = sqrt(sum(N.^2,2));
    if Tmod~=0
        T2(i,:) = T/Tmod;
    end
    if Nmod~=0
        N2(i,:) = N/Nmod;
    end
    vecprod = [T2(i,2)*N2(i,3)-T2(i,3)*N2(i,2) T2(i,3)*N2(i,1)-T2(i,1)*N2(i,3) T2(i,1)*N2(i,2)-T2(i,2)*N2(i,1)];
    Bmod = sqrt(sum(vecprod.^2,2));
    if Bmod~=0
        B2(i,:) = vecprod/Bmod;
    end
end

lattice3Dsec2 = zeros(N2D*Neta3sec2,9);

b = r;

for i=1:Neta3sec2
    e = eccentricity(i,1);
    a = b/sqrt(1-e^2);
    ellipse = ellipticalsection2D(a,b,x0,y0,Neta1,Neta2,eta1min,eta1max,eta2min,eta2max,deltaeta1,deltaeta2,xmin,xmax,Nx,ymin,ymax,Ny);
    [indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4] = getindices2D(Nx,Ny);
    [structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours] = build_neighbourhoods2D(N2D,Nx,flagperiodicity,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesexternalE1,indicesexternalE2,indicesexternalE3,indicesexternalE4,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesC1,indicesC2,indicesC3,indicesC4,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4);
    boundary = [indicesE1;indicesE2;indicesE3;indicesE4;indicesC1;indicesC2;indicesC3;indicesC4];
    if itmax~=0
        ellipse = sparseellipticgridgen2D(Nx,N2D,ellipse,deltaq,0,0,indicesbulk,indicesE1,indicesE2,indicesE3,indicesE4,indicesC1,indicesC2,indicesC3,indicesC4,firstdevneighbours,itmax,tol,0); 
    end
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
%     plot(ellipse(boundary,5),ellipse(boundary,6),'.r')
%     hold on
%     for j=1:size(ellipse,1)
%         if ~any(boundary==j)
%             plot(ellipse(j,5),ellipse(j,6),'.b')
%             hold on
%         end
%     end
%     grid on
%     xlabel('x')
%     ylabel('y')
%     title('Mesh in fluid domain (physical space)')
%     axis equal
    lattice3Dsec2((i-1)*N2D+1:i*N2D,1:2) = ellipse(:,1:2);
    lattice3Dsec2((i-1)*N2D+1:i*N2D,3) = eta3minsec2 + deltaeta3sec2*(i-1);
    lattice3Dsec2((i-1)*N2D+1:i*N2D,4:6) = [x0sec2.*ones(N2D,1) y0sec2.*ones(N2D,1) z0sec2.*ones(N2D,1)] + [meanline2(i,1).*ones(N2D,1) meanline2(i,2).*ones(N2D,1) meanline2(i,3).*ones(N2D,1)] + [ellipse(:,5) ellipse(:,5) ellipse(:,5)].*[N2(i,1).*ones(N2D,1) N2(i,2).*ones(N2D,1) N2(i,3).*ones(N2D,1)] + [ellipse(:,6) ellipse(:,6) ellipse(:,6)].*[B2(i,1).*ones(N2D,1) B2(i,2).*ones(N2D,1) B2(i,3).*ones(N2D,1)];
    lattice3Dsec2((i-1)*N2D+1:i*N2D,7:9) = lattice3Dsec2((i-1)*N2D+1:i*N2D,4:6);
end

lattice3D = [lattice3Dsec1;lattice3Dsec2];

% ------------------> Section 3: tube splitting in two

x0sec3 = x0sec2 + meanline2(end,1);
y0sec3 = y0sec2 + meanline2(end,2);
z0sec3 = z0sec2 + meanline2(end,3);

s3 = (s3min:(s3max-s3min)/(Ns3-1):s3max)';
meanline3 = zeros(Ns3,3);
T3 = zeros(Ns3,3);
N3 = zeros(Ns3,3);
B3 = zeros(Ns3,3);

p = 0.95;

e = eccentricity(end,1);
c = b*e/sqrt(1-e^2);
R = b*sqrt(1-(c/(b/sqrt(1-e^2)))^2);
xs = (0.25*c:(p*R-0.25*c)/(Ns3-1):p*R)';
x1 = xs;
x2 = -xs;
y1 = zeros(length(x1),1);
y2 = zeros(length(x2),1);

Rjunc = 2;

for i = 1:Ns3
    meanline3(i,1:3) = meanlinefunc3(s3(i));
    T = dmeanlinefunc3(s3(i));
    N = d2meanlinefunc3(s3(i));
    Tmod = sqrt(sum(T.^2,2));
    Nmod = sqrt(sum(N.^2,2));
    if Tmod~=0
        T3(i,:) = T/Tmod;
    end
    if Nmod~=0
        N3(i,:) = N/Nmod;
    end
    vecprod = [T3(i,2)*N3(i,3)-T3(i,3)*N3(i,2) T3(i,3)*N3(i,1)-T3(i,1)*N3(i,3) T3(i,1)*N3(i,2)-T3(i,2)*N3(i,1)];
    Bmod = sqrt(sum(vecprod.^2,2));
    if Bmod~=0
        B3(i,:) = vecprod/Bmod;
    end
end

lattice3Dsec3 = zeros(N2D*Neta3sec3,9);

itmax = 1;
tol = 10^-4;

Neta1_1 = 16;
Neta1_2 = 18;
Neta1_3 = 16;

%writerObj = VideoWriter('bifurcation_section3.avi');
%open(writerObj);


for i=1:Neta3sec3
    circles = intersectingcirclessection2D(x1(i),y1(i),R,x2(i),y2(i),R,Rjunc,Neta1_1,Neta1_2,Neta1_3,Neta2,eta1min,eta2min,deltaeta1,deltaeta2,xmin,xmax,ymin,ymax,itmax,tol);
%    scrsz = get(0,'ScreenSize');
%     figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
%     grid on
%     xlabel('x')
%     ylabel('y')
%     title('Mesh in fluid domain (physical space)')
%     axis equal
%     plot(circles(boundary,5),circles(boundary,6),'.r')
%     hold on
%     for j=1:size(circles,1)
%         if ~any(boundary==j)
%             plot(circles(j,5),circles(j,6),'.b')
%             hold on
%         end
%     end
%     grid on
%     xlabel('x')
%     ylabel('y')
%     title('Mesh in fluid domain (physical space)')
%     axis equal
    %frame = getframe;
    %writeVideo(writerObj,frame);
    lattice3Dsec3((i-1)*N2D+1:i*N2D,1:2) = circles(:,1:2);
    lattice3Dsec3((i-1)*N2D+1:i*N2D,3) = eta3minsec3 + deltaeta3sec3*(i-1);
    lattice3Dsec3((i-1)*N2D+1:i*N2D,4:6) = [x0sec3.*ones(N2D,1) y0sec3.*ones(N2D,1) z0sec3.*ones(N2D,1)] + [meanline3(i,1).*ones(N2D,1) meanline3(i,2).*ones(N2D,1) meanline3(i,3).*ones(N2D,1)] + [circles(:,5) circles(:,5) circles(:,5)].*[N3(i,1).*ones(N2D,1) N3(i,2).*ones(N2D,1) N3(i,3).*ones(N2D,1)] + [circles(:,6) circles(:,6) circles(:,6)].*[B3(i,1).*ones(N2D,1) B3(i,2).*ones(N2D,1) B3(i,3).*ones(N2D,1)];
    lattice3Dsec3((i-1)*N2D+1:i*N2D,7:9) = lattice3Dsec3((i-1)*N2D+1:i*N2D,4:6);
end

%close(writerObj);

lattice3D = [lattice3D;lattice3Dsec3];

% ------------------> Section 4: separated tubes

r01sec4 = [x0sec3 y0sec3 z0sec3] + [meanline3(end,1) meanline3(end,2) meanline3(end,3)] + [1.025*R 0 0].*[N3(end,1) N3(end,2) N3(end,3)] + [1.025*R 0 0].*[B3(end,1) B3(end,2) B3(end,3)];
r02sec4 = [x0sec3 y0sec3 z0sec3] + [meanline3(end,1) meanline3(end,2) meanline3(end,3)] + [-1.025*R 0 0].*[N3(end,1) N3(end,2) N3(end,3)] + [-1.025*R 0 0].*[B3(end,1) B3(end,2) B3(end,3)];

s4part1 = (s4part1min:(s4part1max-s4part1min)/(Ns4part1-1):s4part1max)';
meanline4part1 = zeros(Ns4part1,3);
T4p1 = zeros(Ns4part1,3);
N4p1 = zeros(Ns4part1,3);
B4p1 = zeros(Ns4part1,3);

s4part2 = (s4part2min:(s4part2max-s4part2min)/(Ns4part2-1):s4part2max)';
meanline4part2 = zeros(Ns4part2,3);
T4p2 = zeros(Ns4part2,3);
N4p2 = zeros(Ns4part2,3);
B4p2 = zeros(Ns4part2,3);

%p = 0.95;
e = eccentricity(end,1);
c = b*e/sqrt(1-e^2);
R = b*sqrt(1-(c/(b/sqrt(1-e^2)))^2);


for i = 1:Ns4part1
    meanline4part1(i,1:3) = meanlinefunc4part1(s4part1(i));
    T = dmeanlinefunc4part1(s4part1(i));
    N = d2meanlinefunc4part1(s4part1(i));
    Tmod = sqrt(sum(T.^2,2));
    Nmod = sqrt(sum(N.^2,2));
    if Tmod~=0
        T4p1(i,:) = T/Tmod;
    end
    if Nmod~=0
        N4p1(i,:) = N/Nmod;
    end
    vecprod = [T4p1(i,2)*N4p1(i,3)-T4p1(i,3)*N4p1(i,2) T4p1(i,3)*N4p1(i,1)-T4p1(i,1)*N4p1(i,3) T4p1(i,1)*N4p1(i,2)-T4p1(i,2)*N4p1(i,1)];
    Bmod = sqrt(sum(vecprod.^2,2));
    if Bmod~=0
        B4p1(i,:) = vecprod/Bmod;
    end
end

for i = 1:Ns4part2
    meanline4part2(i,1:3) = meanlinefunc4part2(s4part2(i));
    T = dmeanlinefunc4part2(s4part2(i));
    N = d2meanlinefunc4part2(s4part2(i));
    Tmod = sqrt(sum(T.^2,2));
    Nmod = sqrt(sum(N.^2,2));
    if Tmod~=0
        T4p2(i,:) = T/Tmod;
    end
    if Nmod~=0
        N4p2(i,:) = N/Nmod;
    end
    vecprod = [T4p2(i,2)*N4p2(i,3)-T4p2(i,3)*N4p2(i,2) T4p2(i,3)*N4p2(i,1)-T4p2(i,1)*N4p2(i,3) T4p2(i,1)*N4p2(i,2)-T4p2(i,2)*N4p2(i,1)];
    Bmod = sqrt(sum(vecprod.^2,2));
    if Bmod~=0
        B4p2(i,:) = vecprod/Bmod;
    end
end

lattice3Dsec4 = zeros(N2D*Neta3sec4,9);

itmax = 1;
tol = 10^-4;

Neta1part1 = 25;
Neta1part2 = 25;

for i=1:Neta3sec4
    circle1 = circularsection2D(0,0,R,0,Neta1part1,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,ymin,ymax,1,10^-4);
    circle2 = circularsection2D(0,0,R,0,Neta1part2,deltaeta1,eta2min,Neta2,deltaeta2,xmin,xmax,ymin,ymax,1,10^-4);
    for j=1:Neta2
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2,1) = (eta1min:deltaeta1:eta1min+deltaeta1*(Neta1part2-1))';        
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+Neta1part2+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2+Neta1part1,1) = (eta1min+deltaeta1*Neta1part2:deltaeta1:eta1min+deltaeta1*(Neta1part2+Neta1part1-1))';
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+1:(i-1)*N2D+j*Neta1,2) = (eta2min + (j-1)*deltaeta2)*ones(Neta1,1);       
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2,4) = r02sec4(1,1)*ones(Neta1part2,1) + meanline4part2(i,1)*ones(Neta1part2,1) + circle2((j-1)*Neta1part2+1:j*Neta1part2,5).*N3(i,1) + circle2((j-1)*Neta1part2+1:j*Neta1part2,6).*B3(i,1);
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2,5) = r02sec4(1,2)*ones(Neta1part2,1) + meanline4part2(i,2)*ones(Neta1part2,1) + circle2((j-1)*Neta1part2+1:j*Neta1part2,5).*N3(i,2) + circle2((j-1)*Neta1part2+1:j*Neta1part2,6).*B3(i,2);
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2,6) = r02sec4(1,3)*ones(Neta1part2,1) + meanline4part2(i,3)*ones(Neta1part2,1) + circle2((j-1)*Neta1part2+1:j*Neta1part2,5).*N3(i,3) + circle2((j-1)*Neta1part2+1:j*Neta1part2,6).*B3(i,3);
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+Neta1part2+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2+Neta1part1,4) = r01sec4(1,1)*ones(Neta1part1,1) + meanline4part1(i,1)*ones(Neta1part1,1) + circle1((j-1)*Neta1part1+1:j*Neta1part1,5).*N3(i,1) + circle1((j-1)*Neta1part1+1:j*Neta1part1,6).*B3(i,1);
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+Neta1part2+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2+Neta1part1,5) = r01sec4(1,2)*ones(Neta1part1,1) + meanline4part1(i,2)*ones(Neta1part1,1) + circle1((j-1)*Neta1part1+1:j*Neta1part1,5).*N3(i,2) + circle1((j-1)*Neta1part1+1:j*Neta1part1,6).*B3(i,2);
        lattice3Dsec4((i-1)*N2D+(j-1)*Neta1+Neta1part2+1:(i-1)*N2D+(j-1)*Neta1+Neta1part2+Neta1part1,6) = r01sec4(1,3)*ones(Neta1part1,1) + meanline4part1(i,3)*ones(Neta1part1,1) + circle1((j-1)*Neta1part1+1:j*Neta1part1,5).*N3(i,3) + circle1((j-1)*Neta1part1+1:j*Neta1part1,6).*B3(i,3);
    end
    lattice3Dsec4((i-1)*N2D+1:i*N2D,3) = eta3minsec3 + deltaeta3sec3*(i-1);
    lattice3Dsec4((i-1)*N2D+1:i*N2D,7:9) = lattice3Dsec4((i-1)*N2D+1:i*N2D,4:6);
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
%     plot(lattice3Dsec4((i-1)*N2D+1:i*N2D,7),lattice3Dsec4((i-1)*N2D+1:i*N2D,8),'.b')
%     hold on
%     grid on
%     xlabel('x')
%     ylabel('y')
%     title('Mesh in fluid domain (physical space)')
%     axis equal
end

%close(writerObj);

scrsz = get(0,'ScreenSize');
figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(lattice3Dsec4((1-1)*N2D+1:1*N2D,7),lattice3Dsec4((1-1)*N2D+1:1*N2D,8),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

lattice3D = [lattice3D;lattice3Dsec4];

[indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8] = getindices3D(Nx,Ny,Nz);

boundary3D = [indicesF1;indicesF2;indicesF3;indicesF4;indicesF5;indicesF6;indicesE1;indicesE2;indicesE3;indicesE4;indicesE5;indicesE6;indicesE7;indicesE8;indicesE9;indicesE10;indicesE11;indicesE12;indicesC1;indicesC2;indicesC3;indicesC4;indicesC5;indicesC6;indicesC7;indicesC8];

for i=(Neta3sec1+Neta3sec2+Neta3sec3):Neta3
    for j=1:Neta2
        boundary3D = [boundary3D; Neta1part2 + (j-1)*Neta1 + (i-1)*Neta1*Neta2; Neta1part2+1 + (j-1)*Neta1 + (i-1)*Neta1*Neta2];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
f0 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot(circle(boundary,5),circle(boundary,6),'.r')
hold on
for i=1:size(circle,1)
   if ~any(boundary==i)
       plot(circle(i,5),circle(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f1 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f1,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot(circle(boundary,1),circle(boundary,2),'.r')
hold on
for i=1:size(circle,1)
   if ~any(boundary==i)
       plot(circle(i,1),circle(i,2),'.b')
       hold on
   end
end
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot(circle(boundary,5),circle(boundary,6),'.r')
hold on
for i=1:size(circle,1)
   if ~any(boundary==i)
       plot(circle(i,5),circle(i,6),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f2 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3Dsec1(:,7),lattice3Dsec1(:,8),lattice3Dsec1(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f3 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f3,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(lattice3Dsec1(:,1),lattice3Dsec1(:,2),lattice3Dsec1(:,3),'.b')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(lattice3Dsec1(:,7),lattice3Dsec1(:,8),lattice3Dsec1(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f4 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3Dsec2(:,7),lattice3Dsec2(:,8),lattice3Dsec2(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f5 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f5,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(lattice3Dsec2(:,1),lattice3Dsec2(:,2),lattice3Dsec2(:,3),'.b')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(lattice3Dsec2(:,7),lattice3Dsec2(:,8),lattice3Dsec2(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f6 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3D(:,7),lattice3D(:,8),lattice3D(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f7 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f7,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(lattice3D(:,1),lattice3D(:,2),lattice3D(:,3),'.b')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(lattice3D(:,7),lattice3D(:,8),lattice3D(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f8 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f8,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(lattice3Dsec3(:,1),lattice3Dsec3(:,2),lattice3Dsec3(:,3),'.b')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(lattice3Dsec3(:,7),lattice3Dsec3(:,8),lattice3Dsec3(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f9 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3Dsec3(:,7),lattice3Dsec3(:,8),lattice3Dsec3(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f10 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f10,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(lattice3Dsec4(:,1),lattice3Dsec4(:,2),lattice3Dsec4(:,3),'.b')
hold on
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(lattice3Dsec4(:,7),lattice3Dsec4(:,8),lattice3Dsec4(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f11 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3Dsec4(:,7),lattice3Dsec4(:,8),lattice3Dsec4(:,9),'.b')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f12 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f12,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);
plot3(lattice3D(boundary3D,1),lattice3D(boundary3D,2),lattice3D(boundary3D,3),'.r')
hold on
for i=1:size(lattice3D,1)
   if ~any(boundary3D==i)
       plot3(lattice3D(i,1),lattice3D(i,2),lattice3D(i,3),'.b')
       hold on
   end
end
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);
plot3(lattice3D(boundary3D,7),lattice3D(boundary3D,8),lattice3D(boundary3D,9),'.r')
hold on
for i=1:size(lattice3D,1)
   if ~any(boundary3D==i)
       plot3(lattice3D(i,7),lattice3D(i,8),lattice3D(i,9),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis([fcomp fphys],'equal')

scrsz = get(0,'ScreenSize');
f13 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3D(boundary3D,7),lattice3D(boundary3D,8),lattice3D(boundary3D,9),'.r')
hold on
for i=1:size(lattice3D,1)
   if ~any(boundary3D==i)
       plot3(lattice3D(i,7),lattice3D(i,8),lattice3D(i,9),'.b')
       hold on
   end
end
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f14 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3D(boundary3D,7),lattice3D(boundary3D,8),lattice3D(boundary3D,9),'.r')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')
axis equal

scrsz = get(0,'ScreenSize');
f15 = figure('Position',[1 1 ceil(0.5*scrsz(3)) ceil(0.5*scrsz(4))]);
plot3(lattice3D(boundary3D,1),lattice3D(boundary3D,2),lattice3D(boundary3D,3),'.r')
hold on
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in lattice domain (computational space)')
axis equal