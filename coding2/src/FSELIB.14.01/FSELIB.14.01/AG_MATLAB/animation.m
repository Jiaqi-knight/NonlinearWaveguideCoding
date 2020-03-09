close all
clear all

N = 32;
step = 2*pi/N;
rad1 = 1.0;
rad2 = 0.5;

for i=1:N+1
  phi = (i-1.0)*step;
  x(i) = rad1*cos(phi);
  y(i) = rad2*sin(phi);
end

h = plot(x,y,'EraseMode','xor')
axis([0 5 0 5])

tic
for i=1:100
  x = x+0.10;
  y = y+0.10;
  set(h,'YData',y,'XData',x)
  drawnow 
  pause(0.1)
end
toc

