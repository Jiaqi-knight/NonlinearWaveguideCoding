
clear all
close all

%======================================
% Lobatto nodes on the triangle
%
% display the Lobatto node distribution
% on the triangle
%======================================

option=0; % biased
option=1; % symmetrical

ichoose = 2; % chebyshev
ichoose = 0; % uniform
ichoose = 1; % lobatto

idot = 0;
idot = 1;

m = 3;
m = 7;

%-------------
% uniform grid
%-------------

if(ichoose==0)

 for i=1:m+1
   v(i) = (i-1.0)/m;
 end

end

%-------------
% lobatto grid
%-------------

if(ichoose==1)

 [t, W] = lobatto(m-1);    % W will not be needed

 v(1) = 0.0;
 for i=2:m
  v(i) = 0.5*(1.0+t(i-1));
 end
 v(m+1) = 1.0;

end

%--------------
% chebyshev grid
%---------------

if(ichoose==2)

 for i=1:m+1
  c = cos(pi/(m+1) * (i-0.5));
  v(m+2-i)= 0.5*(1.0+c);
 end

end

%----
% horizontal, vertical, and diagonal lines
%----

for i=2:m
 xv1(i) = v(i);
 yv1(i) = 0.0;
 xv2(i) = v(i);
 yv2(i) = 1.0-xv2(i);
 xo1(i) = yv1(i); 
 yo1(i) = xv1(i); 
 xo2(i) = yv2(i); 
 yo2(i) = xv2(i); 
 xd1(i) = v(i);
 yd1(i) = 0.0;
 xd2(i) = 0.0;
 yd2(i) = v(i);
end

%-----
% plot
%-----

figure(1)
hold on
axis equal
axis([0 1 0 1 0 1]);
axis off
set(gca,'fontsize',14)
xlabel('\xi','fontsize',14);
ylabel('\eta','fontsize',14);
box on

for i=1:m+1
   for j=1:m+2-i
       k = m+3-i-j;
          if(option==0)
            xx = v(i);                          % biased
            yy = v(j);
          elseif(option==1)
            xx = (1.0+2.0*v(i)-v(j)-v(k))/3.0;   % symmetric
            yy = (1.0+2.0*v(j)-v(i)-v(k))/3.0;
           end
     if(ichoose==0) plot(xx,yy,'ks','Markersize',10); end
     if(ichoose==1) plot(xx,yy,'ko','Markersize',10); end
     if(ichoose==2) plot(xx,yy,'kx'); end
   end
end

%------------------
% plot the triangle
%------------------

xx(1)=0.0;  yy(1)=0.0;
xx(2)=1.0;  yy(2)=0.0;
xx(3)=0.0;  yy(3)=1.0;
xx(4)=0.0;  yy(4)=0.0;

plot(xx,yy,'k-');

if(idot==1)
 for i=2:m
  plot([xv1(i),xv2(i)],[yv1(i),yv2(i)],'k:')
  plot([xo1(i),xo2(i)],[yo1(i),yo2(i)],'k:')
  plot([xd1(i),xd2(i)],[yd1(i),yd2(i)],'k:')
 end
end

%-------------------------------
% map to an equilateral triangle
%-------------------------------

figure(2)
hold on
axis equal
axis([0 1 0 1 0 1]);
set(gca,'fontsize',14)
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
axis off
box on

for i=1:m+1
   for j=1:m+2-i
       k = m+3-i-j;
        if(option==0)
         xx = v(i);      % biased
         yy = v(j);
        elseif(option==1)
         xx = (1.0+2.0*v(i)-v(j)-v(k))/3.0;    % symmetric
         yy = (1.0+2.0*v(j)-v(i)-v(k))/3.0;
       end
      xx = xx+0.5*yy; 
      yy = sqrt(3)*yy/2.0;
     if(ichoose==0) plot(xx,yy,'ks','Markersize',10); end
     if(ichoose==1) plot(xx,yy,'ko','Markersize',10); end
     if(ichoose==2) plot(xx,yy,'kx'); end
   end
end

%------------------------------
% plot the equilateral triangle
%------------------------------

xxp(1)=0.0;   yyp(1)=0.0;
xxp(2)=1.0;   yyp(2)=0.0;
xxp(3)=0.5;   yyp(3)=sqrt(3)/2.0;
xxp(4)=0.0;   yyp(4)=0.0;

%xxp = xxp - 0.50;yyp = yyp - 1.0/3.0;

plot(xxp,yyp,'k-');

for i=2:m
 xvh1(i) = xv1(i)+0.5*yv1(i);
 yvh1(i) = 0.5*sqrt(3)*yv1(i);
 xvh2(i) = xv2(i)+0.5*yv2(i);
 yvh2(i) = 0.5*sqrt(3)*yv2(i);
 xoh1(i) = xo1(i)+0.5*yo1(i);
 yoh1(i) = 0.5*sqrt(3)*yo1(i);
 xoh2(i) = xo2(i)+0.5*yo2(i);
 yoh2(i) = 0.5*sqrt(3)*yo2(i);
 xdh1(i) = xd1(i)+0.5*yd1(i);
 ydh1(i) = 0.5*sqrt(3)*yd1(i);
 xdh2(i) = xd2(i)+0.5*yd2(i);
 ydh2(i) = 0.5*sqrt(3)*yd2(i);
end

if(idot==1)
 for i=2:m
  plot([xvh1(i),xvh2(i)],[yvh1(i),yvh2(i)],'k:')
  plot([xoh1(i),xoh2(i)],[yoh1(i),yoh2(i)],'k:')
  plot([xdh1(i),xdh2(i)],[ydh1(i),ydh2(i)],'k:')
 end
end


%-----
% done
%-----
