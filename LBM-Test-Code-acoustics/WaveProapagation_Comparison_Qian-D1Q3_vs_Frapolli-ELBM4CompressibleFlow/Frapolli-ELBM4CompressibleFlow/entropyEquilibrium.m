function [ feq] = entropyEquilibrium(N,D,Q,T,rho,w,c,u)
%%
% This program calculates the Lagrange Multipliers of the ME
% probability density functions feq(x) from the knowledge of the
% N contstraints in the form:
% {rho,rho*u,2*rho*E}=sum_{i=1}^{n_Q}{1,v_i,v_i^2}f_i(rho,u,T) ->Gn-eq-3.3
% f_i^eq=rho*W_i*exp[-sum \lambda_n*\phi_n(c_i)]  ->eq-3
% H=-\intfeq(c_i)ln{feq(c_i)dx ->eq-1
%\phi_n(c_i)={1,c_i,c_i^2} 
%In discret form 
%Compare with:
% feq=rho*w+3*(rho*w).*(u*c)+9/2*(rho*w).*(u*c).^2-3/2*(rho.*(u.^2))*w;

w=repmat(w,N,1);
% xmin=-1; xmax=1; dx=0.01; % define the x axis
% x=[xmin:dx:xmax]';
E=1/2*D.*T+1/2*u.^2;
mu=[rho rho.*u 2*rho.*E]; % define the mu values
N_mu=3;
% mu=mu(:); mu=[1;mu]; % add mu(0)=1
% x=x(:); lx=length(x); % x axis
% xmin=x(1); xmax=x(lx); dx=x(2)-x(1);
%
% if(nargin == 2) % initialize LAMBDA
% lambda=zeros(size(mu)); % This produces a uniform
% lambda(1)=log(xmax-xmin); % distribution.
% else
% lambda=lambda0(:);
% end
% N=length(lambda);
%
lambda=zeros(size(mu)); % This produces a uniform
lambda(:,:)=0.1;% distribution.%可能根据平均密度决定好的初始值

M=2*N_mu-1; % Calcul de fin(x)=x.^n
fin=zeros(N,M,Q); %
fin(:,1,:)=1; % fi0(x)=1
for k=1:Q
    fin(:,2,k)=c(k);
    fin(:,3,k)=c(k)^2;
    fin(:,4,k)=c(k)^3;
    fin(:,5,k)=c(k)^4;
end
%
% T=(G(:,3)-G(:,2).^2./G(:,1))./(D*G(:,1));
% if T>(1+sqrt(2/5)) | T<(1-sqrt(2/5))
%     return;
% end
% w(:,1)=(36-49*T+42*T.^2-15*T.^3)/36;
% w(:,2)=(12-13*T+5*T.^2).*T/16;
% w(:,3)=w(:,2);
% w(:,4)=(-3+10*T-5*T.^2).*T/40;
% w(:,5)=w(:,4);
% w(:,6)=(4-15*T+15*T.^2).*T/720;
% w(:,7)=w(:,6);


iter=0;
while 1 % start iterations
iter=iter+1;
disp('---------------'); disp(['iter=',num2str(iter)]);

%eq-3.24
f1=exp(permute(multiprod(permute(lambda,[1,3,2]),fin(:,1:N_mu,:),[2,3]),[1,3,2])); % Calculate p(x)
feq=(repmat(rho,1,Q).*w).*f1; 

G=multiprod(fin,feq,[2,3]);

%
entr(:,iter)=multiprod(permute(G(:,1:N_mu),[1,3,2]),lambda,[2,3]); % Calculate the entropy value
% disp(['Entropy=',num2str(entr(iter))])
%
gnk=zeros(N,N_mu,N_mu); % Calculate gnk
for i=1:N_mu % Matrix G is a Hankel matrix
gnk(:,:,i)=-G(:,i:N_mu+i-1);
end
%
v=mu-G(:,1:N_mu); % Calculate v
for k=1:N
delta(k,:)=permute(gnk(k,:,:),[2,3,1])\v(k,:).'; % Calculate delta
end
lambda=lambda+delta; % Calculate lambda
eps=1e-6; % Stopping rules
if(max(abs(delta./lambda))<eps)
    break
end
if(iter>2)
if(max(abs((entr(:,iter)-entr(:,iter-1))./entr(:,iter)))<eps)
    break;
end
end
end
%
f1=exp(-permute(multiprod(permute(lambda,[1,3,2]),fin(:,1:N_mu,:),[2,3]),[1,3,2])); % Calculate p(x)
feq=(repmat(rho,1,Q).*w).*f1; 

% plot(x,feq); % plot it
% entr=entr(:);
disp('---------------'); disp(['iter=',num2str(iter)]);



end