
clear all
close all

%==========
% Graphs of the Legendre polynomials
%==========

nplot = 2*64;

for i=1:nplot+1
   xi(i) = -1.0+2.0*(i-1)/nplot;
end

figure(1)
hold on
xlabel('\xi','fontsize',15);
ylabel('L','fontsize',15);
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
   Le(i) = x;
  elseif(degree==2)
   Le(i) = 0.5*(3*x^2-1);
  elseif(degree==3)
   Le(i) = 0.5*(5*x^2-3)*x;
  elseif(degree==4)
   Le(i) = 0.125*(35*x^4-30*x^2+3);
  elseif(degree==5)
   Le(i) = 0.125*(63*x^4-70*x^2+15)*x;
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

axis([-1 1 -1 1.1])
 
%---
% done
%---
