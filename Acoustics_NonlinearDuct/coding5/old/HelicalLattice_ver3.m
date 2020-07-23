
clc;clear;
% Map definition
gridpts = linspace(-2, 2, 120);
[mapX, mapY, mapZ] = meshgrid(gridpts,gridpts,gridpts);
N = numel(mapX);

s = 0:0.01:4;
h=0.1*exp(linspace(0,1.5,length(s)));
kappa=(2/3)./h;tau=0.2./h;
a=kappa./(kappa.^2+tau.^2);
b=tau./(kappa.^2+tau.^2);
sw=sqrt(kappa.^2+tau.^2).*s;
mastX= a.*sin(sw);
mastY = a.*cos(sw);
mastZ = b.*sw;
figure;
[x,y,z,t]=tubeplot(mastX,mastY,mastZ,h,s,50);hold on;

tic
[k_in]=inshape_tube(mapX(:),mapY(:),mapZ(:),mastX(:),mastY(:),mastZ(:),t,h,-0.03,0.03);
toc

plot3(mapX(k_in), mapY(k_in), mapZ(k_in),'.');
axis equal
daspect([1,1,1]); camlight;


%% output to data file
Final3Ddata=zeros(size(mapX));
Final3Ddata(k_in)=1;
[row,col,page] = ind2sub(size(mapX),k_in)
figure
plot3(row,col,page,'.')
axis equal

%Write output to 3D.dat with space delimiter.
dlmwrite('3D.dat',Final3Ddata,'delimiter',' ');
isosurface(Final3Ddata)






%% 
    function [k_in]=inshape_tube(x0,y0,z0,x,y,z,t,r,h_min,h_max)
        for  k=1:length(x0)
            temp=sum(([x y z]-[x0(k) y0(k) z0(k)]).*t,2);
            k_zeros=find (temp(1:end-1).*temp(2:end)<=0);
            r0=sqrt((x(k_zeros)-x0(k)).^2+(y(k_zeros)-y0(k)).^2+(z(k_zeros)-z0(k)).^2);
            %disp('判断r0是否在r�?'),有一个为真即可isin
            if isempty(r0)
                isin(k)=0;
            else
                isin(k)=max((r0.'-r(k_zeros))<h_max & (r0.'-r(k_zeros))>h_min);
            end
        end
        k_in=find(isin==1);
    end


