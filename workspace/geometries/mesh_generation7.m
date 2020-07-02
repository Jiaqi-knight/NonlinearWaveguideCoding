clear all
close all
clc
subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
% subfunction_path4=genpath('C:\Users\wjq\Desktop\workspace\geometry-master\geometry')
addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
% addpath(subfunction_path4);

formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];
%%

N = 50;

x0 = 0;
y0 = 0;

r = 10;

b = r;

e = (0:0.01:0.3)';

colors = {'b.','g.','r.','c.','m.','y.','k.'};
numcolor = length(colors);
colorcount = 1;

f0 = figure;
grid on
hold on
xlabel('x')
ylabel('y')
title('Ellipse as a function of eccentricity')
for i=1:length(e)
    a = b/sqrt(1-e(i)^2);
    ellipse = ellipse2D(a,b,0,x0,y0,N);
    plot(ellipse(:,1),ellipse(:,2),colors{colorcount})
    hold on
    if colorcount<numcolor
        colorcount = colorcount + 1;
    else
        colorcount = 1;
    end
end
axis equal

c = e(end)*b/sqrt(1-e(end)^2);

center1 = [0.5*c 0;0.25*c 0;c 0];
center2 = [-0.5*c 0;-0.25*c 0;-c 0];
R1 = b/sqrt(1-e(end)^2) - c;
R2 = b/sqrt(1-e(end)^2) - c;
R3 = b*sqrt(1-(c/(b/sqrt(1-e(end)^2)))^2);
R4 = b*sqrt(1-(c/(b/sqrt(1-e(end)^2)))^2);

theta = (0:2*pi/100:2*pi)';

f1 = figure;
plot(r*cos(theta),r*sin(theta),'c*')
hold on
%plot(center1(1)+R1*cos(theta),center1(2)+R1*sin(theta),'b.')
%hold on
%plot(center2(1)+R2*cos(theta),center2(2)+R2*sin(theta),'b.')
%hold on
for i=1:size(center1,1)
    plot(center1(i,1)+R3*cos(theta),center1(i,2)+R3*sin(theta),colors{i})
    hold on
    plot(center2(i,1)+R4*cos(theta),center2(i,2)+R4*sin(theta),colors{i})
    hold on
    plot(center1(i,1),center1(i,2),'kx')
    hold on
    plot(center2(i,1),center2(i,2),'kx')
    hold on
end
a = b/sqrt(1-e(end)^2);
ellipse = ellipse2D(a,b,0,x0,y0,N);
plot(ellipse(:,1),ellipse(:,2),'m.')
hold on
grid on
xlabel('x')
ylabel('y')
title('Ellipse and circles')
axis equal