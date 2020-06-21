
clc;clear;%close all

% Map definition
gridptsX = 0:1:299;
gridptsY = 0:1:299;
gridptsZ = 0:1:539;


[mapX, mapY, mapZ] = meshgrid(gridptsX,gridptsY,gridptsZ);
N = numel(mapX);

s = 0:1:300;
h=20*exp(linspace(0,1.1,length(s)));
kappa=(0.2)./h;tau=0.4./h;
a=kappa./(kappa.^2+tau.^2);
b=tau./(kappa.^2+tau.^2);
sw=sqrt(kappa.^2+tau.^2).*s;
mastX= a.*sin(sw);
mastY = a.*cos(sw);
mastZ = b.*sw;

figure;view([-22 19.6]);
[x,y,z,t,n]=tubeplot(mastX,mastY,mastZ,h,s,50);hold on;
% plot3(mastX,mastY,mastZ);plot3(mastX(1), mastY(1), mastZ(1),'r*');

% surf(x,y,z,'LineStyle','none');
% plot3(mapX(:), mapY(:), mapZ(:),'.');

axis equal
daspect([1,1,1]); camlight;
xlabel('x');ylabel('y');zlabel('z');
grid on
[newmastX,newmastY,newmastZ]=rotateGeo([mastX(:) mastY(:) mastZ(:)],t(1,:),'z',[250,250,0]);
plot3(newmastX,newmastY,newmastZ)
plot3(newmastX(1),newmastY(2),newmastZ(3),'b*');

[x,y,z,t,n]=tubeplot(newmastX,newmastY,newmastZ,h,s,50);hold on;
surf(x,y,z,'LineStyle','none');

tic
%10�˲���
p=10;
parfor k=1:10
[k_in{k}]=inshape_tube(mapX((k-1)*(N/p)+1:(k)*N/p),mapY((k-1)*(N/p)+1:(k)*N/p),mapZ((k-1)*(N/p)+1:(k)*N/p),newmastX(:),newmastY(:),newmastZ(:),t,h,-1,1);
end
toc
kin=[];
for k=1:10
[kin]=[kin k_in{k}+(k-1)*(N/p)];
end
figure
plot3(mapX(kin), mapY(kin), mapZ(kin),'.');

%% output to data file
Final3Ddata=zeros(size(mapX));
Final3Ddata(kin)=1;
[row,col,page] = ind2sub(size(mapX),cell2mat(k_in));
% figure
% plot3(row,col,page,'.')
% axis equal
Final3Ddata1=permute(Final3Ddata,[1,3,2]);
%Write output to 3D.dat with space delimiter.
dlmwrite('3D.dat',Final3Ddata1,'delimiter',' ');
figure
isosurface(Final3Ddata1)


function [newmastX,newmastY,newmastZ]=rotateGeo(M,t,base,destination) % ��תGeo
% h:���������
% azel����ת������ķ���
% alpha��ת�ĽǶ�
% origin:��ת����������
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
    M2=((M1'))';%������ת�ӽ�180�ȣ�ǰβԵ����
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


