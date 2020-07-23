clear all
close all

%============================================
% Code hlm_lob
%
% Helholtz equation with spectral elements
% corresponding to the Lobatto nodes
%===========================================

%---
% prepare
%---

figure(1)
hold on
xlabel('x','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)
box on

%-----------
% input data
%-----------

L=1.0; k=1.0;
fL=0; q0=-pi; helm=pi^2; 

ne=2; ratio=1.0;     % number of elements, stretch ratio

for i=1:ne
 np(i)=6;   % element polynomial order
end

%-----
% trap
%-----

for i=1:ne
   if(np(i)>7)
   disp(' -->');
   disp(' max polynomial order is 6');
   error(' hlm_lob: Sorry this high order not yet implemented');
   end
end

%------------------------
% element node generation 
%------------------------

[xe,xen,xg,c,ng] = discr_lob (0,L,ne,ratio,np);

%-----------------
% element assembly
%-----------------

[gdm,b] = hlm_lob_sys (ne,xe,np,ng,c,q0,fL,k,helm);

%--------------
% linear solver
%--------------

gdm(:,ng) = [];  % remove the last (ng) column
gdm(ng,:) = [];  % remove the last (ng) row
b(:,ng) = [];    % remove the last (ng) element

f = b/gdm';       % solve the linear system

f = [f fL];      % add the value at the right end

%-----
% plot
%-----

plot(xg, f,'k-+')

%ye = zeros(ne+1,1);
%plot(xe,ye,'ko');
%for i=1:ne
% plot([xe(i),xe(i)],[0,2.5],'k:');
%end

%-----
% done
%-----
