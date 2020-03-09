clear all
close all

%==========
% FSELIB
%
% plot the Bessel function: J_{-1/3}
%==========

L = 20.0;
Np = 64;
Dx = (L-0.001)/Np;

for i=1:Np+1
 x(i) = (i-1.0)*Dx+0.001;
 y(i) = besselj(-1/3,x(i));
 y1(i) = x(i)^(1/3)*y(i)
end

%---
% plot
%---

figure(1)
hold on;
plot(x,y,'k')
plot(x,y1,'k--')
xlabel('x','fontsize',14);
ylabel('J_{-1/3}         x^{1/3}J_{-1/3}','fontsize',14);
set(gca,'fontsize',14)
axis([0 L -1 1])
plot([0,L],[0,0],'k:')
box on

