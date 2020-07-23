
clear all
close all

%==========
% Graphs of the Lobatto polynomials
%==========

nplot = 2*64;

for i=1:nplot+1
   xi(i) = -1.0+2.0*(i-1)/nplot;
end

figure(1)
hold on
xlabel('\xi','fontsize',15);
ylabel('Lo','fontsize',15);
set(gca,'fontsize',15)
box on
 
iplot = 1;

%----
for degree=0:5 
%----
 
for i=1:nplot+1
  x=xi(i);
  if(degree==0)
   Le(i) = 1.0;
  elseif(degree==1)
   Le(i) = 3*x;
  elseif(degree==2)
   Le(i) = 1.5*(5*x^2-1);
  elseif(degree==3)
   Le(i) = 2.5*(7*x^2-3)*x;
  elseif(degree==4)
   Le(i) = 15/8*(21*x^4-14*x^2+1);
  elseif(degree==5)
   Le(i) = 0.125*(693*x^4-630*x^2+105)*x;
  end
end

if(iplot==1)
   plot(xi,Le,'k-','linewidth',1);
elseif(iplot==2)
   plot(xi,Le,'k--','linewidth',1);
elseif(iplot==3)
   plot(xi,Le,'k-.','linewidth',1);
elseif(iplot==4)
   plot(xi,Le,'k:','linewidth',1);
elseif(iplot==5)
   plot(xi,Le,'k-','linewidth',2);
else
   plot(xi,Le,'k--','linewidth',2);
   iplot = 0;
end

iplot = iplot+1;

%----
end % of degree
%----

axis([-1 1 -10 10.0])
 
%---
% done
%---
