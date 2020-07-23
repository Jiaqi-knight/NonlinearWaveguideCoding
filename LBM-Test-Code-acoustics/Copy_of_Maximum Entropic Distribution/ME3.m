%ME3
% This scripts shows how to use the function ME_DENS3
% in the case of the trigonometric moments.
clear;clf
xmin=-5; xmax=5; dx=0.5; % define the x axis
x=[xmin:dx:xmax]';lx=length(x); 
p=(1/sqrt(2*pi))*exp(-.5*(x.*x));% Gaussian distribution
plot(x,p);title('p(x)')
%
M=3;fin=fin3_x(x,M); % Calculate fin(x)
%
mu=zeros(M,1); % Calculate mun
for n=1:M
mu(n)=dx*sum(fin(:,n).*p);
end
%
w0=2*pi/(xmax-xmin);w=w0*[0:M-1]'; % Define the w axis
%
mu=mu(2:M); % Attention : mu(0) is added
% in ME_DENS3
[lambda,p,entr]=me_dens3(mu,x);
disp([mu;lambda;entr(length(entr))]')
