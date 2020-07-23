function h=rotsurf(curve,direct,point,theta,f)
% rotsurf(curve,dirct,orgin,alpha,fun)用于绘制旋转曲面
%   curve=[x,y,z]为母线，其中x,y,z为列向量，分别代表母线的三维坐标
%   direct和origin分别代表旋转轴的方向和该旋转轴上的任意一点的坐标，这两个参数合起来确定了一条直线，即旋转轴，其中：
%       direct表示旋转轴的方向，有两种表示法[theta,phi]或[x0,y0,z0]，其中：
%           theta代表沿xoy平面从x轴正方向逆时针旋转的弧度，phi代表从xoy平面向z轴正方向旋转的弧度
%           [x0,y0,z0]代表方向向量
%           direct默认为[0 0 1]，即z轴方向
%       origin=[xo,yo,zo]为该旋转轴上的任意一点坐标，默认为[0 0 0]即原点
%   向量alpha为旋转的弧度，默认为0:pi/1:2*pi，采样点的范围和密度都可以手动控制
%   fun为绘图函数句柄，默认为@mesh
% h=rotsurf(...)
%  绘制曲面的同时返回该曲面的句柄h
%
%例1：绘制母线为x=0,y^2+z^2=1，旋转轴为x=1,z=-y-2的圆环
%t=linspace(-pi,pi,37).'; 
%y=sin(t);z=cos(t);x=y-y; 
%rotsurf([x,y,z],[0 -1 1],[1 -2 0],[],@surf) 
%xlabel('x');ylabel('y');zlabel('z');axis equal
%例2：绘制母线为z=x,y=1，旋转轴为z轴的单叶双曲面
% t=linspace(-pi,pi,37).'; 
% x=t;y=t-t+1;z=t;
% rotsurf([x y z]) 
% xlabel('x');ylabel('y');zlabel('z');axis equal

assert(nargin>=1 && nargin<=5,'参数个数错误！请看帮助！');
if nargin<5
    f=@mesh;
    if nargin<4
        theta=linspace(0,2*pi,37);
        if nargin<3
            point=[0,0,0];
            if nargin<2
                direct=[0,0,1];
            end
        end
    end
end
curve=squeeze(curve);
assert(ismatrix(curve),'参数1格式错误！请看帮助！');
if size(curve,2)~=3
    curve=curve.';
end
assert(size(curve,2)==3,'参数1格式错误！请看帮助！');
direct=squeeze(direct);
if isempty(direct)
    direct=[0,0,1];
end
assert(numel(direct)==2 || numel(direct)==3,'参数2格式错误！请看帮助！');
if numel(direct)==2
    direct=[cos(direct(2))*cos(direct(1)),cos(direct(2))*sin(direct(1)),sin(direct(2))];
end
direct=direct/norm(direct);
point=squeeze(point);
if isempty(point)
    point=[0 0 0];
end
assert(numel(point)==3,'参数3格式错误！请看帮助！');

theta=squeeze(theta);
if isempty(theta)
    theta=linspace(0,2*pi,37);
end
assert(length(theta)==numel(theta) ,'参数4格式错误！请看帮助！');

f=squeeze(f);
if isempty(f)
    f=@mesh;
end
if ischar(f)
    assert(length(f)==numel(f),'参数5格式错误！请看帮助！');
    f=str2func(f);
end
assert(numel(f)==1 && isa(f,'function_handle'),'参数5格式错误！请看帮助！');


x0=point(1);
y0=point(2);
z0=point(3);

x=curve(:,1);
y=curve(:,2);
z=curve(:,3);

nx=direct(1);
ny=direct(2);
nz=direct(3);

[X,~]=meshgrid(x,theta);
[Y,~]=meshgrid(y,theta);
[Z,T]=meshgrid(z,theta);

sint=sin(T);
cost=cos(T);

XX=(X-x0).*(nx^2*(1-cost)+cost)+(Y-y0).*(nx*ny*(1-cost)-nz*sint)+(Z-z0).*(nx*nz*(1-cost)+ny*sint)+x0;
YY=(X-x0).*(ny*nx*(1-cost)+nz*sint)+(Y-y0).*(ny^2*(1-cost)+cost)+(Z-z0).*(ny*nz*(1-cost)-nx*sint)+y0;
ZZ=(X-x0).*(nz*nx*(1-cost)-ny*sint)+(Y-y0).*(nz*ny*(1-cost)+nx*sint)+(Z-z0).*(nz^2*(1-cost)+cost)+z0;

hh=f(XX,YY,ZZ);
if nargout==1
    h=hh;
end
end
