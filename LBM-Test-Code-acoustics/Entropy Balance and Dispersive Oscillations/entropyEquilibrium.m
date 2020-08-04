function [feq,T] = entropyEquilibrium(N,D,Q,rho,w,c,u,T)
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
Er=rho*D.*T+rho.*u.^2;
mu=[rho rho.*u Er]; % define the mu values
N_mu=3;
lambda=zeros(size(mu)); % This produces a uniform
lambda(:,1)=0.1;% distribution.%可能根据平均密度决定好的初始值
M=2*N_mu-1; % Calcul de fin(x)=x.^n
fin=zeros(N,M,Q); %
for k=1:Q
    fin(:,1,:)=1; 
    fin(:,2,k)=c(k);
    fin(:,3,k)=c(k)^2;
    fin(:,4,k)=c(k)^3;
    fin(:,5,k)=c(k)^4;
end
iter=0;
while 1 % start iterations
iter=iter+1;
% disp('---------------'); disp(['iter=',num2str(iter)]);
f1=exp(-permute(multiprod(permute(lambda,[1,3,2]),fin(:,1:N_mu,:),[2,3]),[1,3,2])); % Calculate p(x)
feq=(repmat(rho,1,Q).*w).*f1; 
G=multiprod(fin,feq,[2,3]);
%
% entr(:,iter)=multiprod(permute(G(:,1:N_mu),[1,3,2]),lambda,[2,3]); % Calculate the entropy value
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
eps=1e-15; % Stopping rules
if(max(max(abs(delta)))<eps)
    break
end
% if(iter>2)
% if(max(abs((entr(:,iter)-entr(:,iter-1))./entr(:,iter)))<eps)
%     break;
% end
% end
%
f1=exp(-permute(multiprod(permute(lambda,[1,3,2]),fin(:,1:N_mu,:),[2,3]),[1,3,2])); % Calculate p(x)
feq=(repmat(rho,1,Q).*w).*f1; 

% plot(x,feq); % plot it
% entr=entr(:);

end
disp('---------------'); disp(['iter=',num2str(iter)]);
T=(G(:,3)-G(:,2).^2./G(:,1))./(D*G(:,1));
