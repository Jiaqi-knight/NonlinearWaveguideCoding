function yint = lagrange_interp (N,x,y,xint)

%-----------------------
% lagrange interpolation
%-----------------------

yint = 0.0;

for i=1:N+1
  lpoly = 1.0;
  for j=1:N+1
    if j~=i
       lpoly = lpoly*(xint-x(j))/(x(i)-x(j));
    end
  end
  yint = yint + lpoly*y(i);
end

%-----
% done
%-----

return
