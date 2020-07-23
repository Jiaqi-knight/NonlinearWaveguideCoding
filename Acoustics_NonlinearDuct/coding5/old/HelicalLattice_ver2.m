function HelicalLattice_ver2


% On the GPU, we define a number of variables including:
%
% * |map|: An N x 3 height field in a 10km x 10km grid (N = 10,000)
% * |masts|: An M x 3 array of antenna positions, at height |H|
% * |AntennaDirection|: A 3 x M array of vectors representing the
% orientation of each antenna.
clc;clear;
% Map definition
gridpts = linspace(-2, 2, 50);
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

% Antenna properties
M = numel(mastX);

% Set up indices into the data
mapIndex = gpuArray.colon(1,N)'; % N x 1 array of map indices
mastIndex = gpuArray.colon(1,M); % 1 x M array of mast indices
[RowIndex, ColIndex] = ndgrid(mapIndex, mastIndex);

% Put the map data on the GPU and concatenate the map positions into a
% single 3-column matrix containing all the coordinates [X, Y, Z].
map = gpuArray([mapX(:) mapY(:) mapZ(:)]);
masts = gpuArray([mastX(:) mastY(:) mastZ(:)]);
X = reshape(map, [], 1, 3);
A = reshape(masts, 1, [], 3);

%这是gpu擅长的事情
    function [d,dotProd] = powerKernel(mapIndex, mastIndex)
        
        % Implement norm and dot product calculations via loop
        dSq = 0;
        dotProd = 0;
        for coord = 1:3
            path = A(1,mastIndex,coord) - X(mapIndex,1,coord);
            dSq = dSq + path*path;
            dotProd = dotProd + path*t(mastIndex,coord);
        end
        d = sqrt(dSq);
        
    end


%%

tic 
[d,dotProd] = arrayfun(@powerKernel, mapIndex, mastIndex) ;
d2=dotProd(:,1:end-1).*dotProd(:,2:end);
D2=gather(d2);%gpu转换回cpu
D=gather(d);
%这是cpu擅长的事情
for k=1:N
    k_zeros=find (D2(k,:)<=0);
    r0=D(k,k_zeros);
    %disp('判断r0是否在r内'),有一个为真即可isin
    if isempty(r0)
        isin(k)=0;
    else
        isin(k)=(min(r0- h(k_zeros))<=0);
    end
end
k_in=find(isin==1);
toc

plot3(mapX(k_in), mapY(k_in), mapZ(k_in),'.');
axis equal
daspect([1,1,1]); camlight;
tic
[k_in]=inshape_tube(mapX(:),mapY(:),mapZ(:),mastX(:),mastY(:),mastZ(:),t,h);
toc
%% cpu
    function [k_in]=inshape_tube(x0,y0,z0,x,y,z,t,r)
        %x0=1;y0=1;z0=1;
        for  k=1:length(x0)
            temp=sum(([x y z]-[x0(k) y0(k) z0(k)]).*t,2);
            k_zeros=find (temp(1:end-1).*temp(2:end)<=0);
            r0=sqrt((x(k_zeros)-x0(k)).^2+(y(k_zeros)-y0(k)).^2+(z(k_zeros)-z0(k)).^2);
            %disp('判断r0是否在r内'),有一个为真即可isin
            if length(r0)==0
                isin(k)=0;
            else
                isin(k)=(min(r0.'-r(k_zeros))<=0);
            end
        end
        k_in=find(isin==1);
    end


end