clear all
close all

%==========================================
% FSELIB
%
% Code hlml
%
% One-dimensional Helmholtz
% equation with linear elements
%
% f'' + alpha f = 0
%
% ne: number of elements
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

L = 1.0; k = 1.0;
ne = 32; ratio = 2.0;

alpha = 87.4;
q0 = -1.0; fL=0.2;

alpha = -54.2;
q0 = 1.0; fL=0.5;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

%-----------------
% system assembly
%-----------------

[at,bt,ct,b] = hlml_sys (ne,xe,q0,fL,k,alpha);

%--------------
% linear solver
%--------------

f = thomas (ne,at,bt,ct,b);

f(ne+1) = fL;

%-----
% plot
%-----

plot(xe, f,'-ko');

%-----
% done
%-----

%=================
% exact solution
%=================

sa = sqrt(abs(alpha));
npl = 64;
Dx = L/npl;

for i=1:npl+1
  xex(i) = (i-1.0)*Dx;
end

%---
if(alpha>0)
%---
  c1 = -q0/(k*sa);
  c2 = (fL -c1*sin(sa*L))/cos(sa*L); 
  for i=1:npl+1
    yex(i) = c1*sin(sa*xex(i))+c2*cos(sa*xx(i));
  end
%---
elseif(alpha<0)
%---
  den = exp(sa*L)+exp(-sa*L);
  c1 = (fL-q0/(k*sa)*exp(-sa*L))/den;
  c2 = (fL+q0/(k*sa)*exp( sa*L))/den;
  for i=1:npl+1
    yex(i) = c1*exp(sa*xex(i))+c2*exp(-sa*xex(i));
  end
%---
end
%---

plot(xex,yex,'r:')

