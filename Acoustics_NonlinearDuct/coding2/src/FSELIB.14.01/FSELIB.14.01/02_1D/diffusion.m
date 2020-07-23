close all
clear all

%=================================
% FSELIB
%
% Animation of the solution of the
% unsteady diffusion equation
%=================================

n = 32;
alpha=0.501;
L=1.0;
nstep=20000;

%---
% prepare
%---

Dx = L/n;
Dt = alpha*Dx^2;

%-------------------------------
% absisas and initial condition
%-------------------------------
                                                                                
for i=1:n+1
  x(i) = (i-1.0)*Dx;
  y(i) = 0.0;
end

time = 0;

% animation:

h = plot(x,y,'-o','EraseMode','xor')
axis([0 1 0 1])
xlabel('x','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)

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

  for i=2:n
    y(i) = alpha*ysv(i-1)+ (1.0-2.0*alpha)*ysv(i)...
         + alpha*ysv(i+1);
  end

%---
% boundary conditions
%---
                                                                                
  y(1) = (25.0*exp(-3.0*time) ...
       +  80.0*(1-exp(-2.0*time))-25.0)/60.0;
  y(n+1) = 0.0;
  time = time+Dt;

  if(y(1)>0.99) break; end

%---
% animation
%---

  set(h,'XData',x)
  set(h,'YData',y)
  drawnow

%---
  end  % of time stepping
%---
