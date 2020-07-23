clear all
close all

%=======
% FSLELIB
%
% plot the exact and approximate dispersion relation
%=======

ne = 32;
pi2 = 2*pi;

x     = zeros(ne,1);
s     = zeros(ne,1);
slump = zeros(ne,1);

for n=1:ne
   x(n) = n/ne;kh = pi2*x(n);
   s(n) = sin(kh) * 3/(2+cos(kh));
   slump(n) = sin(kh);
end

figure(1)
hold on
plot(x, s,'ko-');           % consistent
plot(x, slump,'kx-');  % lumped
xe=[0 0.5]; ye=[0 0.5*pi2]; plot(xe,ye,'k--'); % exact
axis([0 1 -2 2]);
xlabel('n/N_E','fontsize',15);
ylabel('s/N_E','fontsize',15);
set(gca,'fontsize',15)
box on

