clear;clc;close all
cs1=1/sqrt(3);
h1=plot([0,1,2],[cs1,cs1,cs1],'b')
legend('exact');
hold on;
for j=0.1:0.1:2;
omega=j;
 XMAX=1000;w=[1/6,2/3,1/6]; 
c=[-1,0,1];
rho=ones(XMAX,1);u=zeros(XMAX,1);fneq=zeros(XMAX,3);k=2*pi/50;x=1:50;
rho(x)=1+0.01*(1-cos(k*x));a1=max(rho);
u(x)=0.01/sqrt(3)*(1-cos(k*x));
feq=rho*w+3*(rho*w).*(u*c)+9/2*(rho*w).*(u*c).^2-3/2*(rho.*(u.^2))*w;
f=feq;%接下来是为了迭代初值
  for i=1:300
     f=f-omega*(f-feq);
      f(:,1)=circshift(f(:,1),-1);f(:,3)=circshift(f(:,3),1);
   end
%  fneq(x,:)=k/omega*0.01*sqrt(3)*sin(k*x')*(w.*(3*c.^2-1));
%  f=feq+fneq;
%.....f~=feq
% plot(rho);figure;
for i = 1:1000
    feq=rho*w+3.*rho*w.*(u*c)+9/2.*rho*w.*(u.^2*c.^2)-3/2.*rho*w.*(u.^2*ones(1,3));
    f=f-omega*(f-feq);%collision
    f(:,1)=circshift(f(:,1),-1);f(:,3)=circshift(f(:,3),1);%streaming&&Periodic
    rho=sum(f,2);
    u=(f(:,3)-f(:,1))./rho;
%     plot(rho) 
%     pause(0.01);
   
end
l=find(rho==max(rho));
ix=l-2:l+2;
ixx=linspace(l-2,l+2);
iyy=spline(ix,rho(ix),ixx);
% hold on;

%plot(ixx,iyy,'r')
m=find(iyy==max(iyy));
t=ixx(m);
cs=(t-25)/1000;
% v=(log(a1)-log(max(iyy))/((t-25)*(k^2)));
h2=plot(omega,cs,'b*') 
ylim([0,1]);
hold on;
end
legend([h1,h2],'exact', 'test');
xlabel('omega')
ylabel('cs')
title('test of sound speed cs')