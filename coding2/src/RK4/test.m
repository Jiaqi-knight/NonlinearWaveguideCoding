clc;
clear;
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);
% Theda=30*pi/180;r=0.075;c0=0.0404;
% c1_0=1/(4/3*pi*r^3)*c0;H=(7-tan(Theda))*cos(Theda);
[T,Y]=ode45(@myode,[0 20],[0.25 0],options);
% [T,Y]=ode45(@myode,[0 20],[c1_0 H]);
% [T,Y]=ode45(@myode,[0 20],[0 1],options);
figure(1)
plot(T,Y(:,1),'r');
hold on
plot(T,Y(:,2),'b');