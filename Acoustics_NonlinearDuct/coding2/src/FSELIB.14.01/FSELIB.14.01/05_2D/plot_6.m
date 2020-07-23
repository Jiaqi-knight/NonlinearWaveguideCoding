function plot6 (ne,ng,p,c,f)

%===========================================
% Color mapped visualization of a function f
% in a domain discretized into 6-node triangles
%===========================================

%---
% compute the maximim and minimum of the function f
%---

fmax = -100.0; fmin =  100.0;

for i=1:ng
 if(f(i) > fmax) fmax = f(i); end
 if(f(i) < fmin) fmin = f(i); end
end

range = 1.2*(fmax-fmin);shift = fmin;
if(abs(range) < 0.0001) range = 1.0; end

%---
% shift the color index in the range (0, 1)
% and graph 7-point patches
%---

for l=1:ne

 j=c(l,1); xp(1)=p(j,1); yp(1)=p(j,2); zp(1)=f(j);
  cp(1)=(f(j)-shift)/range;
 j=c(l,4); xp(2)=p(j,1); yp(2)=p(j,2); zp(2)=f(j);
  cp(2)=(f(j)-shift)/range;
 j=c(l,2); xp(3)=p(j,1); yp(3)=p(j,2); zp(3)=f(j);
  cp(3)=(f(j)-shift)/range;
 j=c(l,5); xp(4)=p(j,1); yp(4)=p(j,2); zp(4)=f(j);
  cp(4)=(f(j)-shift)/range;
 j=c(l,3); xp(5)=p(j,1); yp(5)=p(j,2); zp(5)=f(j);
  cp(5)=(f(j)-shift)/range;
 j=c(l,6); xp(6)=p(j,1); yp(6)=p(j,2); zp(6)=f(j);
  cp(6)=(f(j)-shift)/range;
 j=c(l,1); xp(7)=p(j,1); yp(7)=p(j,2); zp(7)=f(j);
  cp(7)=(f(j)-shift)/range;

 zp(1) = 0.0;
 zp(2) = 0.0;
 zp(3) = 0.0;
 zp(4) = 0.0;
 zp(5) = 0.0;
 zp(6) = 0.0;
 zp(7) = 0.0;
% disp (cp)
% patch(xp,yp,cp);
 patch(xp,yp,zp,cp);
end

%-----
% done
%-----

return;
