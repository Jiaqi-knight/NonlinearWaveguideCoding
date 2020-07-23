function [lambda,p,entr]=me_dens2(mu,x,lambda0)
%ME_DENS2
% [LAMBDA,P,ENTR]=ME_DENS2(MU,X,LAMBDA0)
% This program calculates the Lagrange Multipliers of the ME
% probability density functions p(x) from the knowledge of the
% N moment contstraints in the form:
% E{x^n}=mu(n) n=0:N with mu(0)=1.
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
M=2*N-1; % Calcul de fin(x)=x.^n
fin=zeros(length(x),M); %
fin(:,1)=ones(size(x)); % fi0(x)=1
for n=2:M
fin(:,n)=x.*fin(:,n-1);
end
%
iter=0;
while 1 % start iterations
iter=iter+1;
disp('---------------'); disp(['iter=',num2str(iter)]);
%
p=exp(-(fin(:,1:N)*lambda)); % Calculate p(x)
plot(x,p); % plot it
%
G=zeros(M,1); % Calculate Gn
for n=1:M
G(n)=dx*sum(fin(:,n).*p);
end
%
entr(iter)=lambda'*G(1:N); % Calculate the entropy value
disp(['Entropy=',num2str(entr(iter))])
%
gnk=zeros(N,N); % Calculate gnk
for i=1:N % Matrix G is a Hankel matrix
gnk(:,i)=-G(i:N+i-1);
end
%
v=mu-G(1:N); % Calculate v
delta=gnk\v; % Calculate delta
lambda=lambda+delta; % Calculate lambda
eps=1e-6; % Stopping rules
if(abs(delta./lambda)<eps), break, end
if(iter>2)
if(abs((entr(iter)-entr(iter-1))/entr(iter))<eps),break, end
end
end
%
p=exp(-(fin(:,1:N)*lambda)); % Calculate the final p(x)
plot(x,p); % plot it
entr=entr(:);
disp('----- END -------')
end
