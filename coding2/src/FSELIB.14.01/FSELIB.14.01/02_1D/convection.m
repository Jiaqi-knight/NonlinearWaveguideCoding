close all
clear all

%==========================
% FSELIB
%
% animation of the solution
% of the convection equation
% using different methods
%===========================

%---
% parameters
%---

n = 128;
U = 1.0;
conv = 0.2;
L = 1.0;
nstep = 20000;

method = 2; % FTFS
method = 3; % FTCS
method = 1; % FTBS

method

%---
% prepare
%---

Dx = L/n;
Dt = conv*Dx/abs(U);

%-------------------------------
% absisas and initial condition
%-------------------------------
                                                                                
for i=1:n+2
  x(i) = (i-1.0)*Dx;
  y(i) = 0.5D0*cos(2.0*pi*x(i)/L);
end

time = 0;

% animation:

h = plot(x,y,'-o','EraseMode','xor');
axis([0 1 -0.6 0.6])
xlabel('x','fontsize',14);
ylabel('f','fontsize',14);
set(gca,'fontsize',14)
                                                                                
%--------------
% time stepping
%--------------

  for istep=1:nstep

%---
% save the old values
%---

  for i=1:n+1
    ysv(i) = y(i);
  end

  ysv(n+2) = ysv(2);

  for i=2:n+1
    if(method==1)
     y(i) = ysv(i)-Dt*U*(ysv(i)-ysv(i-1))/(x(i)-x(i-1));  % FTBS
    elseif(method==2)
     y(i) = ysv(i)-Dt*U*(ysv(i+1)-ysv(i))/(x(i+1)-x(i));  % FTFS
    elseif(method==3)
     y(i) = ysv(i)-Dt*U*(ysv(i+1)-ysv(i-1))/(x(i+1)-x(i-1));  % FTCS
    end
  end
  y(1) = y(n+1);

%---
% animation
%---

  set(h,'XData',x)
  set(h,'YData',y)
  drawnow

%---
  end  % of time stepping
%---
