%=======================
% Evaluation of the Lebesque function
% and Lebesque constant 
%=======================

N = 7   % polynomial order
M = 64  % evaluation f the Leb function

%---
% evenly spaced points
%---

for i=1:N+1
 xi(i) = -1.0+(i-1)*2.0/N;
end

%---
% compute the leb function
% at M+1 evaluation points
%---

for i=1:M+1
 xev(i) = -1.0+(i-1)*2.0/M;
end

%---
% evaluate leb 
%---

for l=1:M+1

  leb(l) = 0.0;

  for i=1:N+1  % run over nodes
    lpoly = 1.0;
     for j=1:N+1
       if j~=i
          lpoly = lpoly*(xev(l)-xi(j))/(xi(i)-xi(j));
       end
     end
    leb(l) = leb(l) + abs(lpoly);
  end

end

%---
% plot the Lebesque function
%---

plot(xev, leb,'-o');
xlabel('\xi','fontsize',15);
ylabel('Lebesque function','fontsize',15);
set(gca,'fontsize',15)

%----
% find the Lebesque constant
%----

max(leb)
