% James-3.58
%h(s)
clc
clear
close all
s = 0:0.01:4;
h=0.1*exp(linspace(0,1.5,length(s)));
kappa=(2/3)./h;tau=0.2./h;
a=kappa./(kappa.^2+tau.^2);
b=tau./(kappa.^2+tau.^2);

sw=sqrt(kappa.^2+tau.^2).*s;


n=1;
for k=1:n
x(k,:) = a.*sin(sw+2*pi*k/n);
y(k,:) = a.*cos(sw+2*pi*k/n);
z(k,:) = b.*sw;
end
% figure(1)
%  
% xlabel('x'); ylabel('y'); title('Circula helix');
% axis equal

figure
for k=1:n
[X,Y,Z,x,y,z,t]=tubeplot(x(k,:),y(k,:),z(k,:),h,s,50);hold on;
end

[qx, qy ,qz] = meshgrid(-2:0.1:2, -2:0.1:2, -2:0.1:2);
[isin]=inshape_tube(qx(:),qy(:),qz(:),x,y,z,t,h);
k_in=find(isin==1);


%% 


%% 


plot3(x, y, z);
daspect([1,1,1]); camlight;
figure
plot3(qx(k_in), qy(k_in), qz(k_in),'.');
axis equal
%% gpu计算prefer
function [isin]=inshape_tube(x0,y0,z0,x,y,z,t,r)
%x0=1;y0=1;z0=1;  
for  k=1:length(x0)
temp=sum(([x y z]-[x0(k) y0(k) z0(k)]).*t,2);
[value,order]=min(abs(temp));
k_zeros=find (temp(1:end-1).*temp(2:end)<=0);
r0=sqrt((x(k_zeros)-x0(k)).^2+(y(k_zeros)-y0(k)).^2+(z(k_zeros)-z0(k)).^2);
%disp('判断r0是否在r内')
%有一个为真即可isin
%try
if length(r0)==0
    isin(k)=0;
else
    isin(k)=(min(r0.'-r(k_zeros))<=0);
end
%catch
%     k
% end

end
end

function [t,n,b,kappa,tau,theta_0]=frenet(x,y,z)

% FRENET Calculate the Frenet frame for a polygonal space curve
% [t,n,b]=frenet(x,y,z) returns the tangent unit vector, the normal
% and binormal of the space curve x,y,z. The curve may be a row or
% column vector, the frame vectors are each row vectors. 
%
% If two points coincide, the previous tangent and normal will be used.
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005

N=size(x,1);
if (N==1)
  x=x';
  y=y';
  z=z';
  N=size(x,1);
end

t=zeros(N,3);
b=zeros(N,3);
n=zeros(N,3);

q=[x y z];

for i=2:(N-1)
  t(i,:)=(q(i+1,:)-q(i-1,:));
  tl(i,:)=norm(t(i,:)); %ds=scalar
  if (tl(i,:)>0)
    t(i,:)=t(i,:)/tl(i,:);
  else
    t(i,:)=t(i-1,:);
  end
end

t(1,:)=q(2,:)-q(1,:);
tl(1,:)=norm(t(1,:));
t(1,:)=t(1,:)/tl(1,:);

t(N,:)=q(N,:)-q(N-1,:);
tl(N,:)=norm(t(N,:));
t(N,:)=t(N,:)/tl(N,:);

%t:nominalized tangent vector
for i=2:(N-1)
  n(i,:)=(t(i+1,:)-t(i-1,:));
  nl(i,:)=norm(n(i,:));
  if (nl(i,:)>0)
    n(i,:)=n(i,:)/nl(i,:);
  else
    n(i,:)=n(i-1,:);
  end
end
%n:norm(dt/ds)
n(1,:)=t(2,:)-t(1,:);
nl(1,:)=norm(n(1,:));
n(1,:)=n(1,:)/nl(1,:);

n(N,:)=t(N,:)-t(N-1,:);
nl(N,:)=norm(n(N,:));
n(N,:)=n(N,:)/nl(N,:);

for i=1:N
  b(i,:)=cross(t(i,:),n(i,:)); %normlized
  %b=txn
end;

for i=2:(N-1)
  bl(i,:)=norm((b(i+1,:)-b(i-1,:)));
end
bl(1,:)=norm(b(2,:)-b(1,:));
bl(N,:)=norm(b(N,:)-b(N-1,:));



kappa=nl./tl;
tau=bl./tl;
tll=tl/2;tll(1)=tl(1);tll(end)=tl(end);
theta_0=cumsum(tau.*tll);
end


function [varargout]=tubeplot(x,y,z,varargin)  

% TUBEPLOT - plots a tube r along the space curve x,y,z.
%
% tubeplot(x,y,z) plots the basic tube with radius 1
% tubeplot(x,y,z,r) plots the basic tube with variable radius r (either a vector or a value)
% tubeplot(x,y,z,r,v) plots the basic tube with coloring dependent on the values in the vector v
% tubeplot(x,y,z,r,v,s) plots the tube with s tangential subdivisions (default is 6)
%
% [X,Y,Z]=tubeplot(x,y,z) returns [Nx3] matrices suitable for mesh or surf
%
% Note that the tube may pinch at points where the normal and binormal 
% misbehaves. It is suitable for general space curves, not ones that 
% contain straight sections. Normally the tube is calculated using the
% Frenet frame, making the tube minimally twisted except at inflexion points.
%
% To deal with this problem there is an alternative frame:
% tubeplot(x,y,z,r,v,s,vec) calculates the tube by setting the normal to
% the cross product of the tangent and the vector vec. If it is chosen so 
% that it is always far from the tangent vector the frame will not twist unduly
%
% Example:
%
%  t=0:(2*pi/100):(2*pi);
%  x=cos(t*2).*(2+sin(t*3)*.3);
%  y=sin(t*2).*(2+sin(t*3)*.3);
%  z=cos(t*3)*.3;
%  tubeplot(x,y,z,0.14*sin(t*5)+.29,t,10)
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005


  subdivs = 6;

  N=size(x,1);
  if (N==1)
    x=x';
    y=y';
    z=z';
    N=size(x,1);
  end

  if (nargin == 3)
    r=x*0+1;
  else
    r=varargin{1};
    if (size(r,1)==1 & size(r,2)==1)
      r=r*ones(N,1);
    end
  end
  if (nargin > 5)
    subdivs=varargin{3}+1;
  end
  if (nargin > 6)
    vec=varargin{4};
    [t,n,b]=frame(x,y,z,vec);
  else
    [t,n,b,kappa,tau,theta_0]=frenet(x,y,z);
  end


  X=zeros(N,subdivs);
  Y=zeros(N,subdivs);
  Z=zeros(N,subdivs);
  theta=linspace(2*pi/subdivs,2*pi,subdivs)+1/2*pi;
  %theta_0=0*theta_0;
  for i=1:N
    X(i,:)=x(i) + r(i)*(n(i,1)*cos(theta-theta_0(i)) + b(i,1)*sin(theta-theta_0(i)));
    Y(i,:)=y(i) + r(i)*(n(i,2)*cos(theta-theta_0(i)) + b(i,2)*sin(theta-theta_0(i)));
    Z(i,:)=z(i) + r(i)*(n(i,3)*cos(theta-theta_0(i)) + b(i,3)*sin(theta-theta_0(i)));
  end

  %if (nargout==0)
    if (nargin > 4)
      V=varargin{2};
      if (size(V,1)==1)
	V=V';
      end
      V=V*ones(1,subdivs);
      surf(X,Y,Z,V,'LineStyle','none');

    else
      surf(X,Y,Z,'LineStyle','none');

    end
  %else
    varargout(1) = {X}; 
    varargout(2) = {Y}; 
    varargout(3) = {Z}; 
    varargout(4) = {x}; 
    varargout(5) = {y}; 
    varargout(6) = {z}; 
    varargout(7) = {t}; 

  end

