function [t,n,b,kappa,tau,theta_0]=frenet(x,y,z)

% FRENET Calculate the Frenet frame for a polygonal space curve
% [t,n,b]=frenet(x,y,z) returns the tangent unit vector, the normal
% and binormal of the space curve x,y,z. The curve may be a row or
% column vector, the frame vectors are each row vectors. 
%
% If two points coincide, the previous tangent and normal will be used.
%
% Written by Anders Sandberg, asa@nada.kth.se, 2005

N=size(x,1);
if (N==1)
  x=x';
  y=y';
  z=z';
  N=size(x,1);
end

t=zeros(N,3);
b=zeros(N,3);
n=zeros(N,3);

q=[x y z];

for i=2:(N-1)
  t(i,:)=(q(i+1,:)-q(i-1,:));
  tl(i,:)=norm(t(i,:)); %ds=scalar
  if (tl(i,:)>0)
    t(i,:)=t(i,:)/tl(i,:);
  else
    t(i,:)=t(i-1,:);
  end
end

t(1,:)=q(2,:)-q(1,:);
tl(1,:)=norm(t(1,:));
t(1,:)=t(1,:)/tl(1,:);

t(N,:)=q(N,:)-q(N-1,:);
tl(N,:)=norm(t(N,:));
t(N,:)=t(N,:)/tl(N,:);

%t:nominalized tangent vector
for i=2:(N-1)
  n(i,:)=(t(i+1,:)-t(i-1,:));
  nl(i,:)=norm(n(i,:));
  if (nl(i,:)>0)
    n(i,:)=n(i,:)/nl(i,:);
  else
    n(i,:)=n(i-1,:);
  end
end
%n:norm(dt/ds)
n(1,:)=t(2,:)-t(1,:);
nl(1,:)=norm(n(1,:));
n(1,:)=n(1,:)/nl(1,:);

n(N,:)=t(N,:)-t(N-1,:);
nl(N,:)=norm(n(N,:));
n(N,:)=n(N,:)/nl(N,:);

for i=1:N
  b(i,:)=cross(t(i,:),n(i,:)); %normlized
  %b=txn
end;

for i=2:(N-1)
  bl(i,:)=norm((b(i+1,:)-b(i-1,:)));
end
bl(1,:)=norm(b(2,:)-b(1,:));
bl(N,:)=norm(b(N,:)-b(N-1,:));



kappa=nl./tl;
tau=bl./tl;
tll=tl/2;tll(1)=tl(1);tll(end)=tl(end);
theta_0=cumsum(tau.*tll);






