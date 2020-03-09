clear all
close all

%==========
% Lagrange interpolating functions psi
%==========

%-----
% data:
%-----
 
m = input('Will use m+1 nodes please enter m: ');

ibase = 1;  % evenly spaced
ibase = 2;  % lobatto

%---------------
% evenly spaced:
%---------------

if(ibase==1)

 for i=1:m+1
   xi(i) = -1.0+2.0*(i-1)/m;
 end

end

%--------
% lobatto:
%--------

if(ibase==2)

 xi(1) = -1.0;
 if(m>1)
   [zL, wL]= lobatto(m-1);
   for i=2:m
    xi(i) = zL(i-1);
   end
 end
xi(m+1)=1.0;

axis([-1 1 -0.2 1])
axis([-1 1 -0.25 1])

end

%----------
% Chebyshev
%----------

if(ibase==3)

 for i=1:m+1
  xi(i) = cos((i-1)*pi/m);
 end

end

%----------------
% function values
%----------------

for i=1:m+1
  y(i) = 1.0/(1.0+xi(i)^2);
end
 
%-----------------
% prepare a graph:
%-----------------

figure(1)
hold on
xlabel('\xi','fontsize',15);
ylabel('\psi','fontsize',15);
set(gca,'fontsize',15)
box on
 
a = xi(1);
b = xi(m+1);

nplot = 2*64;
step = (b-a)/nplot;

iplot = 1;

%----
for i=1:m+1 % will plot m+1 polynomials
%----
 
for l=1:nplot+1
  xint(l)=a+step*(l-1);
  yint(l) = 1.0;
  for j=1:m+1
    if(j~=i)
       yint(l) = yint(l)*(xint(l)-xi(j))/(xi(i)-xi(j));
    end
  end
end
 
if(iplot==1)
   plot(xint,yint,'k-','linewidth',1);
elseif(iplot==2)
   plot(xint,yint,'k--','linewidth',1);
elseif(iplot==3)
   plot(xint,yint,'k-.','linewidth',1);
else
   plot(xint,yint,'k:','linewidth',1);
   iplot = 0;
end

iplot = iplot+1;

%----
end % of i
%----

% axis([-1 1 -10.5 11.0])
 
% done
