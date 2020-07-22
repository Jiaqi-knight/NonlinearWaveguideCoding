% Applying Nelder-Mead method to find Lagrange multipliers to generate
% distribution with maximum entropy regarding constraints. 
% This code is developed by Foad Karimian.
% To verify the code, following data can be used based on example provided
% in Mohammad-Djafari, A. (1992). A Matlab program to calculate the maximum entropy 
% distributions: if x=[-1:0.01:1], mu=[0 0.3 0 0.15], then independent of
% initial value for lambda, final lambda values should be lambda=[0.9392 0
% -3.3414 0 4.6875].

clear all;
% clc;

xmin=-5;
xmax=5;
dx=0.001;
x=[xmin:dx:xmax];
mu=[0 1 -0.3491 2.5835];

mu=mu(:);                    %import mu and make a vector
x=x(:);                      %make a vector of x
M=length(mu);                %determines summation over indicies
phi=ones(length(x),M);       %function to generate moments (mean, variance,...)
phi(:,1)=phi(:,1).*x;        %first column is x

for i=2:M
    phi(:,i)=phi(:,i-1).*x;  %generate x^i
end

phmu=zeros(length(x),M);

for i=1:M
    phmu(:,i)=phi(:,i)-mu(i);       %generates x^i - mu(i)
end

l0=zeros(M,1);

Q = @(l) sum(exp(-phmu*l).*dx);

options = optimset('Display','iter','PlotFcns',@optimplotfval);
lambda = fminsearch(Q,l0,options);  %minimizing potential function to find Lagrangian multiplyers

q = sum(exp(-phmu*lambda).*dx);     %calculate potential value

pdf=exp(-phmu*lambda)./q;             %generate distribution  

lambda0 = log(q.*exp(-lambda.'*mu));    %find lambda0, normalizing factor

lambda=[lambda0;lambda];

cdf=zeros(length(x),1);
cdf(1)=pdf(1).*dx;

for i=2:length(x)
    cdf(i)=pdf(i).*dx+cdf(i-1);     %generate cdf of MaxEnt
end

ent=0;
entropy=0;
for i=1:length(pdf)
    if pdf(i) == 0
        ent = ent;
    else
        ent = pdf(i).*log(pdf(i)).*dx; 
        entropy = entropy + ent;
    end
end


disp('Lagrange multipliers (lambda0, lambda1, ...) are:'); disp(num2str(lambda));

disp('Entropy: '); disp(num2str(-entropy));

figure(1);          %draw pdf MaxEnt distribution
plot(x,pdf)
title('MaxEnt PDF')
xlabel('Data')
ylabel('PDF')


figure(2);            %draw pdf MaxEnt distribution
plot(x,cdf)
title('MaxEnt CDF')
xlabel('Data')
ylabel('CDF')

movegui(figure(1),'center')      %place figure to center on screen
movegui(figure(2),'east')      %place figure to right on screen