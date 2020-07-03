%%
%        Project: Duct Acoustic in curvilinear coordination by LBM
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

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subfunction_path1=genpath('C:\Users\wjq\Desktop\NonlinearWaveguideCoding-master (1)\NonlinearWaveguideCoding-master\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\NonlinearWaveguideCoding-master (1)\NonlinearWaveguideCoding-master\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\NonlinearWaveguideCoding-master (1)\NonlinearWaveguideCoding-master\workspace\differential_geometry-master\matlab');
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% DUCT MESHES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta1min = 0;
Neta1 = 60;
deltaeta1 = 1;
eta1max = eta1min + (Neta1-1)*deltaeta1;

eta2min = 0;
Neta2 = 60;
deltaeta2 = 1;
eta2max = eta2min + (Neta2-1)*deltaeta2;

eta3min = 0;
Neta3 = 60;
deltaeta3 = 1;
eta3max = eta3min + (Neta3-1)*deltaeta3;


Nr = Neta1;
rmin = 1;
rmax = 2.5;
deltar = (rmax-rmin)/(Nr-1);

Nt = Neta2;
tmin = 0;
tmax = 2*pi;%0和2pi重复包含
deltat = (tmax-tmin)/(Nt);
tmax=tmax-deltat;
tmin:deltat:tmax;

Nz = Neta3;
zmin = 0;
zmax = 2;
deltaz = (zmax-zmin)/(Nz-1);
deltaq = [deltar deltat deltaz];

N = Nr*Nt*Nz;

lattice = generatelattice3D(eta1min,Neta1,deltaeta1,eta2min,Neta2,deltaeta2,eta3min,Neta3,deltaeta3,rmin,rmax,Nr,tmin,tmax,Nt,zmin,zmax,Nz);

%lattice第5，6,7列为transform后的坐标点
% Apply the parameterization to obtain R^3 coordinates
lattice(:,7) = lattice(:,4).*cos(lattice(:,5)); %x
lattice(:,8) = lattice(:,4).*sin(lattice(:,5)); %y
lattice(:,9) = lattice(:,6); %z

% lattice(:,7) = lattice(:,4).*cos(lattice(:,5)); %x
% lattice(:,8) = lattice(:,4).*sin(lattice(:,5)); %y
% lattice(:,9) = lattice(:,6); %z

[indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8]=getindices3D(Nr,Nt,Nz);
%% visualization
f1 = figure();
title('Computational and physical domains')
hold on
subplot(1,2,1);
subplot(1,2,2);
hsubfigs = get(f1,'Children');
fcomp = hsubfigs(2);
fphys = hsubfigs(1);
subplot(fcomp);%view(fcomp,[-42.3 -6.8]);
plot3(lattice(:,1),lattice(:,2),lattice(:,3),'b.')
hold on
plot3(lattice(indicesC1,1),lattice(indicesC1,2),lattice(indicesC1,3),'Marker','pentagram')
plot3(lattice(indicesC2,1),lattice(indicesC2,2),lattice(indicesC2,3),'Marker','pentagram')
plot3(lattice(indicesC3,1),lattice(indicesC3,2),lattice(indicesC3,3),'Marker','pentagram')
plot3(lattice(indicesC4,1),lattice(indicesC4,2),lattice(indicesC4,3),'Marker','pentagram')
plot3(lattice(indicesC5,1),lattice(indicesC5,2),lattice(indicesC5,3),'Marker','pentagram')
plot3(lattice(indicesC6,1),lattice(indicesC6,2),lattice(indicesC6,3),'Marker','pentagram')
plot3(lattice(indicesC7,1),lattice(indicesC7,2),lattice(indicesC7,3),'Marker','pentagram')
plot3(lattice(indicesC8,1),lattice(indicesC8,2),lattice(indicesC8,3),'Marker','pentagram')
text(Nr/2,-2,0,'E1');text(Nr+2,Nt/2,0,'E2');text(Nr/2,Nt+1,0,'E3');text(-4,Nt/2,0,'E4');
text(-2,-2,Nz/2,'E5');text(Nr,-1,Nz/2,'E6');text(Nr+1,Nt+1,Nz/2,'E7');text(-6,Nt,Nz/2,'E8');
text(Nr/2,-2,Nz,'E9');text(Nr+2,Nt/2,Nz,'E10');text(Nr/2,Nt+1,Nz,'E11');text(-4,Nt/2,Nz,'E12');
text(Nr/2,Nt/2,-1,'F1下','Color','red');text(Nr/2,-1,Nz/2,'F2前','Color','red');text(Nr+1,Nt/2,Nz/2,'F3右','Color','red');
text(Nr/2,Nt+1,Nz/2,'F4后','Color','red');text(-1,Nt/2,Nz/2,'F5左','Color','red');text(Nr/2,Nt/2,Nz+1,'F6上','Color','red');
text(lattice(indicesC1,1),lattice(indicesC1,2),lattice(indicesC1,3),'C1')
text(lattice(indicesC2,1),lattice(indicesC2,2),lattice(indicesC2,3),'C2')
text(lattice(indicesC3,1),lattice(indicesC3,2),lattice(indicesC3,3),'C3')
text(lattice(indicesC4,1),lattice(indicesC4,2),lattice(indicesC4,3),'C4')
text(lattice(indicesC5,1),lattice(indicesC5,2),lattice(indicesC5,3),'C5')
text(lattice(indicesC6,1),lattice(indicesC6,2),lattice(indicesC6,3),'C6')
text(lattice(indicesC7,1),lattice(indicesC7,2),lattice(indicesC7,3),'C7')
text(lattice(indicesC8,1),lattice(indicesC8,2),lattice(indicesC8,3),'C8')

for i=1:length(indicesE1)
    plot3(lattice(indicesE1(i),1),lattice(indicesE1(i),2),lattice(indicesE1(i),3),'or')
    hold on
end
for i=1:length(indicesF1)
    plot3(lattice(indicesF1(i),1),lattice(indicesF1(i),2),lattice(indicesF1(i),3),'*r')
    hold on
end
for i=1:length(indicesF2)
    plot3(lattice(indicesF2(i),1),lattice(indicesF2(i),2),lattice(indicesF2(i),3),'*k')
    hold on
end
for i=1:length(indicesF3)
    plot3(lattice(indicesF3(i),1),lattice(indicesF3(i),2),lattice(indicesF3(i),3),'*g')
    hold on
end
for i=1:length(indicesF4)
    plot3(lattice(indicesF4(i),1),lattice(indicesF4(i),2),lattice(indicesF4(i),3),'*y')
    hold on
end
grid on
xlabel('$\xi_{1}$','Interpreter','LaTex')
ylabel('$\xi_{2}$','Interpreter','LaTex')
zlabel('$\xi_{3}$','Interpreter','LaTex')
title('Mesh in lattice domain (computational space)')
subplot(fphys);%view(fphys,[-42.3 -6.8]);
plot3(lattice(:,7),lattice(:,8),lattice(:,9),'b.')
hold on
plot3(lattice(indicesC1,7),lattice(indicesC1,8),lattice(indicesC1,9),'Marker','pentagram')
plot3(lattice(indicesC2,7),lattice(indicesC2,8),lattice(indicesC2,9),'Marker','pentagram')
plot3(lattice(indicesC3,7),lattice(indicesC3,8),lattice(indicesC3,9),'Marker','pentagram')
plot3(lattice(indicesC4,7),lattice(indicesC4,8),lattice(indicesC4,9),'Marker','pentagram')
plot3(lattice(indicesC5,7),lattice(indicesC5,8),lattice(indicesC5,9),'Marker','pentagram')
plot3(lattice(indicesC6,7),lattice(indicesC6,8),lattice(indicesC6,9),'Marker','pentagram')
plot3(lattice(indicesC7,7),lattice(indicesC7,8),lattice(indicesC7,9),'Marker','pentagram')
plot3(lattice(indicesC8,7),lattice(indicesC8,8),lattice(indicesC8,9),'Marker','pentagram')
for i=1:length(indicesE1)
    plot3(lattice(indicesE1(i),7),lattice(indicesE1(i),8),lattice(indicesE1(i),9),'or')
    hold on
end
for i=1:length(indicesE2)
    plot3(lattice(indicesE2(i),7),lattice(indicesE2(i),8),lattice(indicesE2(i),9),'ok')
    hold on
end
for i=1:length(indicesE3)
    plot3(lattice(indicesE3(i),7),lattice(indicesE3(i),8),lattice(indicesE3(i),9),'xr')
    hold on
end
for i=1:length(indicesE4)
    plot3(lattice(indicesE4(i),7),lattice(indicesE4(i),8),lattice(indicesE4(i),9),'xk')
    hold on
end
for i=1:length(indicesE5)
    plot3(lattice(indicesE5(i),7),lattice(indicesE5(i),8),lattice(indicesE5(i),9),'xk')
    hold on
end
for i=1:length(indicesF1)
    plot3(lattice(indicesF1(i),7),lattice(indicesF1(i),8),lattice(indicesF1(i),9),'*r')
    hold on
end
for i=1:length(indicesF2)
    plot3(lattice(indicesF2(i),7),lattice(indicesF2(i),8),lattice(indicesF2(i),9),'*k')
    hold on
end
for i=1:length(indicesF3)
    plot3(lattice(indicesF3(i),7),lattice(indicesF3(i),8),lattice(indicesF3(i),9),'*g')
    hold on
end
for i=1:length(indicesF4)
    plot3(lattice(indicesF4(i),7),lattice(indicesF4(i),8),lattice(indicesF4(i),9),'oy')
    hold on
end
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title('Mesh in fluid domain (physical space)')

%% tensor calculation
%flagperiodicity = 1;
periodicity = 2;% 2  --> neighbour of F2 - F4
flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;
[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours]=build_neighbourhoods3D(N,Nr,Nt,Nz,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

%% covariantbase
xyz=reshape(lattice,[Nr,Nt,Nz,9]);
mxds=1;
covariantbase = computecovariantbase3D(N,deltaq,lattice,firstdevneighbours);
cobase=reshape(covariantbase,[Nr,Nt,Nz,9]);

%% metriccoefficients
[metriccoefficients,g,sqrtg] = computemetriccoefficients3D(covariantbase);
metric=reshape(metriccoefficients,[Nr,Nt,Nz,6]);
g_name={'g_{11}','g_{22}','g_{33}','g_{12}g_{21}','g_{13}g_{31}','g_{23}g_{32}'}

contravariantbase = computecontravariantbase3D(covariantbase,sqrtg);

[reciprocalmetriccoefficients,g,sqrtg] = computereciprocalmetriccoefficients3D(contravariantbase);

firstChristoffelsymbol = computefirstChristoffelsymbol3D(N,deltaq,covariantbase,firstdevneighbours);

secondChristoffelsymbol = computesecondChristoffelsymbol3D(N,deltaq,covariantbase,contravariantbase,firstdevneighbours);
Christoffelsymbol=reshape(secondChristoffelsymbol,[Nr,Nt,Nz,27]);





cs = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                       % Squared speed of sound cl^2

%% LBM
D=3;Q=41;
[scheme,cs,cssq,T0]=initializeELBMscheme(D,Q);
invcssq=1/cssq;
% ftrue = 1;  %potential on/off
%make_lattice(&fIn,&fOut,c,wi,nx,ny,nz,sd,lambda,T0);
%初始化rho和u
rho=ones(N,1);
u=zeros(N,3);
xx=Nr/2,yy=Nt/2,zz=Nz/2
id=xx+1+(yy)*Nr+(zz)*Nr*Nt;
rho(id)=0.1;
F=zeros(N,3);
[fin]=computefequilibrium3D(N,Q,rho,u,reciprocalmetriccoefficients,scheme,invcssq);
dt=1;
beta = 1.9;                                      % Relaxation frequency
% h=figure;
for solutime=1:50
    %eq(nx,ny,nz,&fIn,&fOut,&rho,&ux,&uy,&uz,c,wi,sd,omega,lambda,&ts,T0,dt,ca,ftrue);
    %stream(nx,ny,nz,&fIn,&fOut,c,sd);
    [feq]=computefequilibrium3D(N,Q,rho,u,reciprocalmetriccoefficients,scheme,invcssq);
    [fout]=stream3D(N,Q,scheme,fin,lattice);%propagation
    [rho]=computerho3D(fout);[u]=computevelocity3D(N,Q,rho,fout,scheme);
    [Fgeom]=computegeometricforcing3D(N,Q,secondChristoffelsymbol,scheme);
    [Fext]=computeexternalforcing3D(N,Q,contravariantbase,F);
    [Flambda]=computeforcingterm3D(N,Q,Fgeom,Fext,rho,u,reciprocalmetriccoefficients,scheme,invcssq);
    [fin]=collide3D(Q,fout,feq,Flambda,beta,dt);

    rrhhoo=reshape(rho,[Nr,Nt,Nz,1])-1;
    uu=reshape(u,[Nr,Nt,Nz,3]);
    
    %% tecplot output

     tsignal.cubes(solutime).zonename='mysurface zone';
     tsignal.cubes(solutime).x=xyz(:,:,:,7);    %size 3x3 
     tsignal.cubes(solutime).y=xyz(:,:,:,8);    %size 3x3
     tsignal.cubes(solutime).z=xyz(:,:,:,9);    %size 3x3
     tsignal.cubes(solutime).v(1,:,:,:)=rrhhoo;%可以为三维输入，p只有一维
     tsignal.cubes(solutime).v(2,:,:,:)=uu(:,:,:,1);
     tsignal.cubes(solutime).v(3,:,:,:)=uu(:,:,:,2);
     tsignal.cubes(solutime).v(4,:,:,:)=uu(:,:,:,3);
     tsignal.cubes(solutime).solutiontime=solutime;
     
end
 
    
    

 %% 整合数据，生成文件
%title=''; 
NAME = [date,'Duct'];  %存储文件夹   
NAME1 = ['Lattice_test','-',date];  %存储文件夹   
 
output_file_name=[NAME,'1.plt']; 
tsignal.Nvar=7;     
tsignal.varnames={'x','y','z','rho','ux','uy','uz'};
mat2tecplot(tsignal,output_file_name);
