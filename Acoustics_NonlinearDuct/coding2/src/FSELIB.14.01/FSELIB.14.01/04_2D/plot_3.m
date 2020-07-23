function plot_3 (ne,ng,p,c,f)

%-=============================================
% color mapped visualization of a function f
% in a domain discretized into 3-node triangles
%==============================================

%----------------------
% compute the maximim and minimum of the function f
%----------------------

fmin =  100.0;
fmax = -100.0;

for i=1:ng
 if(f(i) > fmax) fmax = f(i); end
 if(f(i) < fmin) fmin = f(i); end
end

range = 1.2*(fmax-fmin);
shift = fmin;

%---
% shift the color index in the range (0, 1)
% and plot 4-point patches
%---

for l=1:ne
 j=c(l,1); xp(1)=p(j,1); yp(1)=p(j,2); cp(1)=(f(j)-shift)/range;
 j=c(l,2); xp(2)=p(j,1); yp(2)=p(j,2); cp(2)=(f(j)-shift)/range;
 j=c(l,3); xp(3)=p(j,1); yp(3)=p(j,2); cp(3)=(f(j)-shift)/range;
 j=c(l,1); xp(4)=p(j,1); yp(4)=p(j,2); cp(4)=(f(j)-shift)/range;
 patch(xp,yp,cp);
end

%-----
% done
%-----

return;
