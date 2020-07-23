function [lambda,p,entr]=me_dens3(mu,x,lambda0)
%ME_DENS3
% [LAMBDA,P,ENTR]=ME_DENS3(MU,X,LAMBDA0)
% This program calculates the Lagrange Multipliers of the ME
% probability density functions p(x) from the knowledge of the
% Fourier moments values :
% E{exp[-j n w0 x]}=mu(n) n=0:N with mu(0)=1.
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
mu=mu(:);mu=[1;mu]; % add mu(0)=1
x=x(:); lx=length(x); % x axis
xmin=x(1); xmax=x(lx); dx=x(2)-x(1);
if(nargin == 2) % initialize LAMBDA
lambda=zeros(size(mu)); % This produces a uniform
lambda(1)=log(xmax-xmin); % distribution.
else
lambda=lambda0(:);
end
N=length(lambda);
%
M=2*N-1; % Calculate fin(x)=exp[-jnw0x]
fin=fin3_x(x,M); % fin3_x(x) is an external
% % function which provides fin(x).
iter=0;
while 1 % start iterations
iter=iter+1;
disp('---------------'); disp(['iter=',num2str(iter)]);
%
% Calculate p(x)
p=exp(-real(fin(:,1:N))*real(lambda)+imag(fin(:,1:N))*imag(lambda));
plot(x,p); % plot it
%
G=zeros(M,1); % Calculate Gn
for n=1:M
G(n)=dx*sum(fin(:,n).*p);
end
%plot([real(G(1:N)),real(mu),imag(G(1:N)),imag(mu)])
%
entr(iter)=lambda'*G(1:N); % Calculate the entropy
disp(['Entropy=',num2str(entr(iter))])
%
gnk=zeros(N,N); % Calculate gnk
for n=1:N % Matrix gnk is a Hermitian
for k=1:n % Toeplitz matrix.
gnk(n,k)=-G(n-k+1); % Lower triangle part
end
end
for n=1:N
for k=n+1:N
gnk(n,k)=-conj(G(k-n+1)); % Upper triangle part
end
end
%
v=mu-G(1:N); % Calculate v
delta=gnk\v; % Calculate delta
lambda=lambda+delta; % Calculate lambda
eps=1e-3; % Stopping rules
if(abs(delta)./abs(lambda)<eps), break, end
if(iter>2)
if(abs((entr(iter)-entr(iter-1))/entr(iter))<eps),break, end
end
end
% Calculate p(x)
p=exp(-real(fin(:,1:N))*real(lambda)+imag(fin(:,1:N))*imag(lambda));
plot(x,p); % plot it
entr=entr(:);
disp('----- END -------')
end


