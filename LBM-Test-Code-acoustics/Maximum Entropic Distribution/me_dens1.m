function [lambda,p,entr]=me_dens1(mu,x,lambda0)
%ME_DENS1
% [LAMBDA,P,ENTR]=ME_DENS1(MU,X,LAMBDA0)
% This program calculates the Lagrange Multipliers of the ME
% probability density functions p(x) from the knowledge of the
% N contstraints in the form:
% E{fin(x)}=MU(n) n=0:N with fi0(x)=1, MU(0)=1.
%
% MU is a table containing the constraints MU(n),n=1:N.
% X is a table defining the range of the variation of x.
% LAMBDA0 is a table containing the first estimate of the LAMBDAs.
% (This argument is optional.)
% LAMBDA is a table containing the resulting Lagrange parameters.
% P is a table containing the resulting pdf p(x).
% ENTR is a table containing the entropy values at each
% iteration.
%
% Author: A. Mohammad-Djafari
% Date : 10-01-1991
%
mu=mu(:); mu=[1;mu]; % add mu(0)=1
x=x(:); lx=length(x); % x axis
xmin=x(1); xmax=x(lx); dx=x(2)-x(1);
%
if(nargin == 2) % initialize LAMBDA
lambda=zeros(size(mu)); % This produces a uniform
lambda(1)=log(xmax-xmin); % distribution.
else
lambda=lambda0(:);
end
N=length(lambda);
%
fin=fin1_x(x); % fin1_x(x) is an external
% % function which provides fin(x).
iter=0;
while 1 % start iterations
iter=iter+1;
disp('---------------'); disp(['iter=',num2str(iter)]);
%
p=exp(-(fin*lambda)); % Calculate p(x)-test,eq-3
plot(x,p); % plot it
%
G=zeros(N,1); % Calculate Gn
for n=1:N
G(n)=dx*sum(fin(:,n).*p);
end
%
entr(iter)=lambda'*G(1:N); % Calculate the entropy value
disp(['Entropy=',num2str(entr(iter))])
%
gnk=zeros(N,N); % Calculate gnk
gnk(1,:)=-G'; gnk(:,1)=-G; % first line and first column
for i=2:N % lower triangle part of the
for j=2:i % matrix G
gnk(i,j)=-dx*sum(fin(:,j).*fin(:,i).*p);
end
end
for i=2:N % uper triangle part of the
for j=i+1:N % matrix G
gnk(i,j)=gnk(j,i);
end
end
%
v=mu-G; % Calculate v eq-5-down
delta=gnk\v; % Calculate delta %eq-7
lambda=lambda+delta; % Calculate lambda%eq-5-down
eps=1e-6; % Stopping rules
if(abs(delta./lambda)<eps), break, end
if(iter>2)
if(abs((entr(iter)-entr(iter-1))/entr(iter))<eps),break, end
end
end
%
p=exp(-(fin*lambda)); % Calculate the final p(x)
plot(x,p); % plot it
entr=entr(:);
disp('----- END -------')
end
