close all
clear all

%==========================================
% FSELIB
%
% Code sdl
%
% Steady one-dimensional diffusion
% with linear elements
%==========================================

%---
% prepare
%---

figure(1)
hold on
xlabel('x','fontsize',14);
ylabel('f','fontsize',14);
set(gca,'fontsize',14)
box on

%-----------
% input data
%-----------

L=1.0; k=1.0;

icase = 3;
icase = 1;

if(icase==1)
 q0=-1.0; fL=0.0;
 ne=32; ratio=1.0;
 ne=16; ratio=5.0;
end

if(icase==2)
 q0=-1; fL=exp(1.0);  % CASE B
 ne=512; ratio=1.0;
end

if(icase==3)
 q0=-100/26^2;    % CASE R and R1
 fL=1/26.0;
 ne=10; ratio=1.0;
 nplot=64;
 for i=1:nplot
  xex(i) = (i-1)/nplot;
  xhat = 2*xex(i)-1;
  yex(i) = 1.0/(1+25*xhat^2);
 end
 plot(xex,yex,'-k','LineWidth',2)
end

%----------------
% grid generation
%----------------

xe = elm_line1(0,L,ne,ratio);

%-------------------
% specify the source
%-------------------

for i=1:ne+1
 if(icase==1)
  s(i) = 0.0;
  s(i) = 1.0;
  s(i) = 10.0*exp(-5.0*xe(i)^2/L^2);
 elseif(icase==2)
  s(i) = -exp(xe(i));   % CASE B
 elseif(icase==3)
  xhat = 2*xe(i)-1;
  s(i) = 1.0/(1+25*xhat^2);
  s(i) = 200*(1-75*xhat^2)/(1+25*xhat^2)^3 ;
 end

end

%-----------------
% compact assembly
%-----------------

[at,bt,ct,b] = sdl_sys (ne,xe,q0,fL,k,s);

%--------------
% linear solver
%--------------

f = thomas (ne,at,bt,ct,b);

f(ne+1) = fL;

%-----
% plot
%-----

%plot(xe, f,'-ko');
plot(xe, f,'-k+');

%-----
% done
%-----
