%% 先修复管道为直线需要近似的bug
clc;clear;close all
subfunction_path1='./scr';
addpath(subfunction_path1);
%Surf_para=MATLAB4geomTurbo;

radius=20;
diameter=2*radius;
nx=2*diameter+10;
ny=2*diameter+10;
nz=4*diameter+60;
thickness_duct=2;
radius_intern=radius-2;

h_min=0;
h_max=100;

% Map definition
gridptsX = 0:1:nx-1;
gridptsY = 0:1:ny-1;
gridptsZ = 0:1:nz-1;


[mapX, mapY, mapZ] = meshgrid(gridptsX,gridptsY,gridptsZ);
N = numel(mapX);

s = 0:1:nz;
h=radius_intern*exp(linspace(0,0,length(s)));
kappa=(0)./h;tau=0.4./h;
a=kappa./(kappa.^2+tau.^2);
b=tau./(kappa.^2+tau.^2);
sw=sqrt(kappa.^2+tau.^2).*s;
mastX= a.*sin(sw);
mastY = a.*cos(sw);
mastZ = b.*sw;

if kappa~=0
    [t]=frenet(mastX(:),mastY(:),mastZ(:));
    [mastX,mastY,mastZ]=rotateGeo([mastX(:) mastY(:) mastZ(:)],t(1,:),'z',[nx/2,ny/2,0]);
    [surf_x,surf_y,surf_z,t,n]=tubeplot(mastX,mastY,mastZ,h,s,50);
else
    [t]=repmat( [0 0 1] , length(s) , 1 );
    mastX=mastX-mastX(1)+nx/2;
    mastY=mastY-mastY(1)+ny/2;
    mastZ=mastZ-mastZ(1)-0;
    [surf_x,surf_y,surf_z]=cylplot(mastX,mastY,mastZ,h,s,50);

end
% 




tic
%10核并行
p=10
parfor k=1:10
[k_in{k}]=inshape_tube(mapX((k-1)*(N/p)+1:(k)*N/p),mapY((k-1)*(N/p)+1:(k)*N/p),mapZ((k-1)*(N/p)+1:(k)*N/p),mastX(:),mastY(:),mastZ(:),t,h,h_min,h_max);
end

toc
kin=[];
for k=1:10
[kin]=[kin k_in{k}+(k-1)*(N/p)];
end


figure;view([-22 19.6]);
%plot3(mapX(:),mapY(:),mapZ(:),'.');hold on;

surf(surf_x,surf_y,surf_z,'LineStyle','none');
daspect([1,1,1]); camlight;
xlabel('x');ylabel('y');zlabel('z');
grid on
axis([0 nx,0 ny,0 nz])
plot3(mapX(kin), mapY(kin), mapZ(kin),'.');
axis equal

%%
% output to data file
Final3Ddata=zeros(size(mapX));
Final3Ddata(kin)=1;
[row,col,page] = ind2sub(size(mapX),cell2mat(k_in));
% figure
% plot3(row,col,page,'.')
% axis equal
Final3Ddata1=permute(Final3Ddata,[1,3,2]);
%Write output to 3D.dat with space delimiter.
dlmwrite('3D.dat',Final3Ddata1,'delimiter',' ');
% figure
% isosurface(Final3Ddata1)


function [newmastX,newmastY,newmastZ]=rotateGeo(M,t,base,destination) % 旋转Geo
% h:物体的坐标
% azel：旋转坐标轴的方向
% alpha：转的角度
% origin:旋转坐标轴的起点
origin=M(1,:);
if base=='x'
    axis=[1,0,0];
elseif base=='y'
    axis=[0,1,0];
elseif base=='z'
    axis=[0,0,1];
end
azel=cross(t,axis);
alph=acos(dot(t(1,:),axis)/(norm(t(1,:))*norm(axis)));
% find unit vector for axis of rotation
u = azel(:)/norm(azel);
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);y = u(2);z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
    x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
    x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
    [m1,n1]=size(M);
    B=repmat(origin,m1,1);
    M1 = (M-B)*rot+B;
    M2=((M1'))';%由于旋转接近180度，前尾缘倒置
    h1=M2-M(1,:)+destination;
    newmastX=h1(:,1);
    newmastY=h1(:,2);
    newmastZ=h1(:,3);

end




%%
function [k_in]=inshape_tube(x0,y0,z0,x,y,z,t,r,h_min,h_max)
for  k=1:length(x0)
    temp=sum(([x y z]-[x0(k) y0(k) z0(k)]).*t,2);
    k_zeros=find (temp(1:end-1).*temp(2:end)<=0);
    r0=sqrt((x(k_zeros)-x0(k)).^2+(y(k_zeros)-y0(k)).^2+(z(k_zeros)-z0(k)).^2);
    if isempty(r0)
        isin(k)=0;
    else
        isin(k)=max((r0.'-r(k_zeros))<h_max & (r0.'-r(k_zeros))>h_min);
    end
end
k_in=find(isin==1);
end




