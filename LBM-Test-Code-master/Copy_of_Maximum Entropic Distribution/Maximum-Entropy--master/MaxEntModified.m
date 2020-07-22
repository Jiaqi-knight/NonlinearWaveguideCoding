% This code is developed by Foad Karimian.
% It is using Newton's method to solve non-linear MaxEnt equation.
% It is developed based on following paper.
% Mohammad-Djafari, A. (1992). A Matlab program to calculate the maximum entropy distributions. 
% In Maximum Entropy and Bayesian Methods (pp. 221-233). Springer, Dordrecht.

% To verify the code, following data can be used based on example provided
% in paper: if x=[-1:0.01:1], mu=[0 0.3 0 0.15], then independent of
% initial value for lambda, final lambda values should be lambda=[0.9392 0
% -3.3414 0 4.6875].
clear all;

xmin=-2;                %X lower bound
xmax=2;                 %X upper bound
dx=0.01;                %Integral percision
mu=[0 0.3 0 0.15];      %constrains
lambda=[1 2 -1 -1];       %Initial lambda values

x=[xmin:dx:xmax];       %Xs
mu=mu(:);
mu=[1;mu]
x=x(:);
lambda=lambda(:);
lambda=[1 ;lambda];     %add lambda0

M=length(mu);
v=zeros(M,1);
gnk=zeros(M);
phlam=zeros(length(x),M);
p=zeros(length(x),1);

phi=ones(length(x),M);       %function to generate moments (mean, variance,...)
phi(:,2)=phi(:,2).*x;        %first column is x


for i=3:M
    phi(:,i)=phi(:,i-1).*x; %phi contains x^i in "i" column
end

iter=0;
while 1
 
    iter=iter+1;
    disp("---------------"); 
    disp(["iter=",num2str(iter)]);
    
    phlam=phi*lambda;            %multipication of lambdas to x, generates sumation within exponential equatoin in G and f(x) and ...
    p=exp(-phlam);              %generates exponential function
    
    G=ones(M,1);
    
    for i=1:M
        G(i,1)=sum(phi(:,i).*p.*dx);
    end
    
    v=mu-G;
    
    for k=1:M
        for i=1:M
            gnk(k,i)=-sum(phi(:,i).*phi(:,k).*p.*dx);
        end
    end
    
    entropy(iter)=-sum(p.*log(p).*dx);
    
    delta=gnk^-1*v;
    lambda=lambda+delta;
    eps=1e-6;
    
    if(abs(delta./lambda)<10000*eps),break,end
         if(iter>2)
          if(abs((entropy(iter)-entropy(iter-1))/entropy(iter))<eps),break,end
         end
    plot(x,p)
    
    
end
disp(["lambda0=",num2str(lambda(1));
      "lambda1=",num2str(lambda(2));
      "lambda2=",num2str(lambda(3));
      "lambda3=",num2str(lambda(4));
      "lambda4=",num2str(lambda(5))]);
 Entropy=entropy(iter)
