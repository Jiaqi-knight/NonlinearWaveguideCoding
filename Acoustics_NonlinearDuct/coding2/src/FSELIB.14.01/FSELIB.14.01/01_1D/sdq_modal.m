clear all
close all

%=================================
% FSELIB
%
% CODE sdq_modal
%
% steady one-dimensional diffusion
% with quadratic modal expansion
%=================================

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

L=1.0;
k=1.0; q0 = -1.0; fL=0.0;

ne=4; ratio=5.0;
ne=1; ratio=5.0;
ne=32; ratio=2.0;
ne=16; ratio=5.0;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

%------------------------
% number of global modes
%-----------------------

ng = 2*ne+1;

%--------------------------
% generate the global nodes
%--------------------------

Ic = 0;  % counter

for i=1:ne
  Ic = Ic+1; xg(Ic) = xe(i);
  Ic = Ic+1; xg(Ic) = 0.5*(xe(i)+xe(i+1));
end

xg(ng) = xe(ne+1);

%-------------------
% specify the source
%-------------------

for i=1:ng
 s(i) = 10.0*exp(-5.0*xg(i)^2/L^2);
end

%---------------------------------------------
% compute the modal coefficients of the source
% by collocation at the element mid-point
%---------------------------------------------

for i=1:ng
 smodal(i) = s(i);
end

for i=1:ne
 j=2*i;
 smodal(j) = (-s(j-1)+2.0*s(j)-s(j+1))/2.0;
end

%-----------------------------------------------------
% compact element assembly of the pentadiagonal system
%-----------------------------------------------------

[ap,bp,cp,dp,ep,rp] = sdq_modal_sys (ne,xe,q0,fL,k,smodal);

%------
% formulate and solve a tridiaginal system
%------

for i=1:ne
 j=2*i-1;
 at(i) = ap(j); bt(i) = cp(j); ct(i) = ep(j); rt(i) = rp(j);
end

sol = thomas (ne,at,bt,ct,rt);

for i=1:ne
 j = 2*i-1;
 f(j) = sol(i);
end

f(ng) = fL;

%------
% solve for the even-numbered unknowns
% to recover the mid-point values
%------

for i=1:ne
 j=2*i;
 f(j) = rp(j)/ap(j);
 f(j) = (f(j-1)+2.0*f(j)+f(j+1))/2.0;
end

%-----
% plot
%-----

plot(xg, f,'-k+');
for i=1:ne+1
 plot([xe(i),xe(i)], [0, 2.5],'k:');
end

%-----
% done
%-----
