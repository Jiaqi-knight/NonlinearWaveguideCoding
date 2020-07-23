function plot_t4 (ne,ng,p,c,f)

%================================================
% Color mapped visualization of a function f
% in a domain discretized into 4-node tetrahedra
%================================================

%----------------
% plotting option
%--------------

option = 4;
option = 1;
option = 2; % wireframe

%--------------------------------------------------
% compute the maximim and minimum of the function f
%--------------------------------------------------

fmin =  100.0;
fmax = -100.0;

for i=1:ng
 if(f(i) > fmax) fmax = f(i); end
 if(f(i) < fmin) fmin = f(i); end
end

range = 1.2*(fmax-fmin);
shift = fmin;

%------------------------------------------
% shift the color index in the range (0, 1)
% and plot 4-node elements
%------------------------------------------

for l=1:ne

 j=c(l,1); xp(1)=p(j,1); yp(1)=p(j,2); zp(1)=p(j,3); 
           cp(1)=(f(j)-shift)/range;
 j=c(l,2); xp(2)=p(j,1); yp(2)=p(j,2); zp(2)=p(j,3);
           cp(2)=(f(j)-shift)/range;
 j=c(l,3); xp(3)=p(j,1); yp(3)=p(j,2); zp(3)=p(j,3);
           cp(3)=(f(j)-shift)/range;
 j=c(l,4); xp(4)=p(j,1); yp(4)=p(j,2); zp(4)=p(j,3);
           cp(4)=(f(j)-shift)/range;

 xmean = 0.25*(xp(1)+xp(2)+xp(3)+xp(4));
 ymean = 0.25*(yp(1)+yp(2)+yp(3)+yp(4));
 zmean = 0.25*(zp(1)+zp(2)+zp(3)+zp(4));

% if(xmean>0 & ymean>0)
% if(xmean>0)
    volume = tetra4_vis (xp,yp,zp,cp,option);
% end

end

%-----
% done
%-----

return;
