% This code is developed by Foad Karimian
% This part of code calculates MLE estimations of parameters for Lognormal
% and 3-paramter Weibull distributions

clear all;
clc;

%Actual data in form of Log(cycles) is entered in variable "data"
data=[4.13944	4.15100	4.16708	4.17609	4.18327	4.18563	4.20412	4.20863	4.21181	4.22050	4.22461	4.23277	4.24460	4.25154	4.26007	4.26423	4.29144	4.29776	4.30270	4.31260	4.32672	4.33630	4.33630	4.37024	4.38322	4.41162	4.41192	4.41664	4.42488	4.43993];

Lognorm=mle(data);      %MLE estimation of mean and std
Lognorm=Lognorm.';      %make a vector of dist. parameters


xmin=-5.*Lognorm(2)+Lognorm(1);     %set x interval with respect to std, extending from 4*std before mean to 4*std after mean
xmax=5.*Lognorm(2)+Lognorm(1);
dx=0.001;
x=[xmin:dx:xmax];
x=x(:);

normp = normpdf(x,Lognorm(1),Lognorm(2));    %generate pdf of lognormal distribution
normc = normcdf(x,Lognorm(1),Lognorm(2));    %generate cdf of lognormal distribution

figure(3);                          %draw pdf lognormal distribution
plot(x,normp)
title('Lognormarl PDF')
xlabel('Log(N)')
ylabel('PDF')

figure(4);                          %draw cdf lognormal distribution
plot(x,normc)
title('Lognormarl CDF')
xlabel('Log(N)')
ylabel('CDF')

movegui(figure(3),'northwest')      %place lognormal pdf on top-left on screen 
movegui(figure(4),'north')          %place lognormal cdf on top-center on screen 

custpdf = @(x,a,b,c) (x>c).*(b/a).*(((x-c)/a).^(b-1)).*exp(-((x-c)/a).^b);  %PDF for 3-parameter Weibull
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'FunValCheck','off');
weibull= mle(data,'pdf',custpdf,'start',[0.2 2.4 3],'Options',opt,...       %MLE estimation of parameters
    'LowerBound',[0 0 -Inf],'UpperBound',[Inf Inf min(data)]);

weibull=weibull.';                  %make a vector of parameters

weipdf=custpdf(x,weibull(1),weibull(2),weibull(3));    %generate PDF of 3-parameter weibull distribtion

weicdf=zeros(length(x),1);
weicdf(1)=weipdf(1).*dx;

for i=2:length(x)
    weicdf(i)=weipdf(i).*dx+ weicdf(i-1);                   %generate CDF of 3-parameter weibull distribtion
end

figure(5)           %draw pdf 3-parameter weibull distribution
plot(x,weipdf)
title('3 Parameter Weibull PDF')
xlabel('Log(N)')
ylabel('PDF')

figure(6)           %draw cdf 3-parameter weibull distribution
plot(x,weicdf)
title('3 Parameter Weibull CDF')
xlabel('Log(N)')
ylabel('CDF')

movegui(figure(5),'west')      %place figure on left-center on screen
movegui(figure(6),'center')      %place figure on center on screen


% Applying Nelder-Mead method to find Lagrange Multipliers to generate distribution with highst
% entropy with constraint. 
% This code is developed by Foad Karimian.
% To verify the code, following data can be used based on example provided
% in Mohammad-Djafari, A. (1992). A Matlab program to calculate the maximum entropy 
% distributions: if x=[-1:0.01:1], mu=[0 0.3 0 0.15], then independent of
% initial value for lambda, final lambda values should be lambda=[0.9392 0
% -3.3414 0 4.6875].

data = data(:);
S = std(data);
Mean = mean(data);
data_st=(data-Mean)/S;

Mst=mean(data_st);
Sst=std(data_st);
sk=skewness(data_st);
kr=kurtosis(data_st);
mu=[Mst Sst sk kr];

xmin=-5;
xmax=5;
dxst=(xmax-xmin)/(length(x)+1);
x=[xmin+dxst:dxst:xmax-dxst];

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
    
lambda0 = log(sum(exp(-phi*lambda).*dx));   %find lambda0, normalizing factor

cdf=zeros(length(x),1);
cdf(1)=pdf(1).*dx;

for i=2:length(x)
    cdf(i)=pdf(i).*dx+cdf(i-1);     %generate cdf of MaxEnt
end


figure(1);          %draw pdf MaxEnt distribution
plot(x,pdf)
title('MaxEnt PDF')
xlabel('Standardized Log(N)')
ylabel('PDF')


figure(2);            %draw pdf MaxEnt distribution
plot(x,cdf)
title('MaxEnt CDF')
xlabel('Standardized Log(N)')
ylabel('CDF')

movegui(figure(1),'southwest')      %place figure on lower-left on screen
movegui(figure(2),'south')      %place figure on lower-center on screen

% Following part of code calculates MxEnt pdf and cdf for unstandardized
% data

x=(x.*std(data))+mean(data);        %Changing X interval to match interval for Lognormal and Weibull

phmun=zeros(length(x),M);           %generating x^i-mu(i) for new x's and mu's
mun=[Mean;S;sk;kr];                 %new mu is actual mean and std of data

for i=1:M
    phmun(:,i)=phi(:,i)-mun(i);       %generates x^i - mu(i) for actual data (unstandardiz)
end


qin = sum(exp(-(phmun*lambda)).*dx);     %calculate potential function value

pdfin=exp(-(phmun*lambda))./qin;         %generate PDF  

cdfin=zeros(length(x),1);                %generate CDF 
cdfin(1)=pdfin(1).*dx;

for i=2:length(x)
    cdfin(i)=pdfin(i).*dx+cdfin(i-1);     %generate cdf of MaxEnt
end

cdf_data=zeros(length(data),1);

for i=1:length(data)
    cdf_data(i)=(i-0.3)/(length(data)+0.4);
end

figure(7)                                   %Lognormal, weibull and MaxEnt PDF in one figure
plot(x,normp,':',x,weipdf,'--',x,pdfin)
title('Lognormal vs 3-p Weibull vs MaxEnt PDF')
legend('Lognormal','3 Paramter Weibull','MaxEnt')
xlabel('Log(N)')
ylabel('PDF')

figure(8)                                   %Lognormal, weibull and MaxEnt CDF in one figure
plot(x,normc,':',x,weicdf,'--',x,cdfin,data,cdf_data,'s')
title('Lognormal vs 3-p Weibull vs MaxEnt CDF')
legend({'Lognormal','3 Parameter Weibull','MaxEnt','Log(cycles)'},'Location','northwest')
xlabel('Log(N)')
ylabel('CDF')


movegui(figure(7),'northeast')      %place figure on upper-right on screen
movegui(figure(8),'east')      %place figure on right on screen

disp('------------------------- MLE parameters for Lognormal and 3-Parameter Weibull distributions -------------------------');

disp('MLE parameters for Lognormal distribution are:');
disp(['mean= ',num2str(Lognorm(1)),'      ','std= ',num2str(Lognorm(2))]);      %display Lognormal Distribution parameters
     
 
disp('MLE parameters for 3-Parameter Weibull distributions are:');
disp(['alpha= ',num2str(weibull(1)),'      ','beta= ',num2str(weibull(2)),'      ','Log(n0)= ',num2str(weibull(3))]);    %display 3-Parameter Weibull Distribution parameters
 
disp('------------------------- Lagrange multipliers for MaxEnt distribution -------------------------');
 
disp('Lagrange multipliers (lambda0, lambda1, ...) are: ');       %display Lagrangian Multipliers parameters
disp(num2str(lambda0));
disp(num2str(lambda));

% This part of the code find r.m.s. of all distributions from data points

d_lognorm=zeros(length(data),1);        % for lognormal distribution

for i=1:length(data)
    d_lognorm(i,1)=(i-0.3)/(length(data)+0.4) - interp1(x,normc,data(i)); 
end

rms_lognorm=(sum(d_lognorm.^2)/length(data))^0.5;

d_wei=zeros(length(data),1);        % for weibull distribution

for i=1:length(data)
    d_wei(i,1)=(i-0.3)/(length(data)+0.4) - interp1(x,weicdf,data(i)); 
end

rms_wei=(sum(d_wei.^2)/length(data))^0.5;

d_maxent=zeros(length(data),1);        % for maxent distribution

for i=1:length(data)
    d_maxent(i,1)=(i-0.3)/(length(data)+0.4) - interp1(x,cdfin,data(i)); 
end

rms_maxent=(sum(d_maxent.^2)/length(data))^0.5;

% Following part of code calculates Likelihood values for distriutions (models) given data and
% AIC criteria for model selection and
% BIC criteria for model selection

lnorm = interp1(x,normp,data(1));       % lnorm calculates likelihood value for normal dist. Initial value is set
for i=2:length(data)
    lnorm = lnorm * interp1(x,normp,data(i));  % Likelihood value is generated by multiplying normal pdf value for each data point
end


lwei = interp1(x,weipdf,data(1));       % lwei calculates likelihood value for weibull dist. Initial value is set
for i=2:length(data)
    lwei = lwei * interp1(x,weipdf,data(i));  % Likelihood value is generated by multiplying weibull pdf value for each data point
end


lmaxent = interp1(x,pdfin,data(1));       % lnorm calculates likelihood value for MaxEnt dist. Initial value is set
for i=2:length(data)
    lmaxent = lmaxent * interp1(x,pdfin,data(i));  % Likelihood value is generated by multiplying MaxEnt pdf value for each data point
end

[aic,bic]=aicbic([log(lnorm);log(lwei);log(lmaxent)],[2;3;length(lambda)],length(data));       %calculates AIC and BIC for dist.s (models)

amin=min(aic);     % Find AIC min
deltaaic=aic-amin; % Calculate delta_i=AICi-AICmin
rlaic=exp(-0.5.*deltaaic);          % Compute relative likelihood of models

bmin=min(bic);      % Find BIC min
deltabic=bic-bmin;  % Calculate delta_i=BICi-BICmin

% Following part of code selects best model based on Maximum Likelihood ratio or Bayesian model
% comparison when prior probabilities for both models are equal

lmin=min([lnorm;lwei;lmaxent]);     % Find minimum value of likelihood
k=[lnorm;lwei;lmaxent]./lmin;       % divide each likelihood by minimum likelihood

% Following part of code calculates covariance matrix for distributions

lcov = mlecov(Lognorm,data,'pdf',@normpdf);         % Covariance matrix for Lognormal distribution
wcov = mlecov(weibull,data,'pdf',custpdf);         % Covariance matrix for 3-p Weibull distribution

% Following part of code calculates Hessian matrix of MaxEnt distribution.
% It should be a positive definite matrix (symetric and positive value
% pivot) which results in positive eigen values and positive determinant
% that shows lagrangian multiliers are minimum

g=zeros(length(mu));

for i=1:length(mu)
    for j=1:length(mu)
        g(i,j)=qin*sum(phmun(:,i).*phmun(:,j).*pdfin.*dx);
    end
end

dg=det(g);

% This Part of code calculates Entropy of each distributions

entnorm = sum(normp.*log(normp).*dx);
entmaxent = sum(pdfin.*log(pdfin).*dx);
ent = 0;
entwei = 0;
for i=1:length(weipdf)
    if weipdf(i) == 0
        ent = ent;
    else
        ent = weipdf(i).*log(weipdf(i)).*dx; 
        entwei = entwei + ent;
    end
end



disp('------------------------- Covariance Matrix for distribution parameters -------------------------')
disp('Covaraince matrix for Lognormal distribution is (mean std):')
disp(num2str(lcov));
disp('Covaraince matrix for 3-Parameter Weibull distribution is (alpha beta log(n0)):')
disp(num2str(wcov));

disp('------------------------- Hessian Matrix for MaxEnt distribution -------------------------')
disp('Hessian matrix of MaxEnt distribution is (lambda1, lambda2,...): ');
disp(num2str(g));
disp(['Determinant of Hessian matrix is: ',num2str(dg)])

disp('------------------------- Entropy of distributions -------------------------')
disp(['Lognormal: ',num2str(entnorm)])
disp(['3-Parameter Weibull: ',num2str(entwei)])
disp(['MaxEnt: ',num2str(entmaxent)])

disp('------------------------- Tail fit and Root Mean Square for distributions -------------------------')

disp('Tail fit for Lognormal distribution are (for two first data points):');
disp([num2str(d_lognorm(1)),'      ',num2str(d_lognorm(2))]);
disp('R.M.S. for Lognormal distribution is:');
disp(num2str(rms_lognorm));

disp('Tail fit for 3 Parameter Weibull distribution are (for two first data points):');
disp([num2str(d_wei(1)),'      ',num2str(d_wei(2))]);
disp('R.M.S. for 3 Parameter Weibull distribution is:');
disp(num2str(rms_wei));

disp('Tail fit for MaxEnt distribution are (for two first data points):');
disp([num2str(d_maxent(1)),'      ',num2str(d_maxent(2))]);
disp('R.M.S. for MaxEnt distribution is:');
disp(num2str(rms_maxent));

disp('------------------------- Likelihood values for distributions -------------------------')

disp('Likelihood value for dist.s are:');
disp(['Lognormal= ',num2str(lnorm),'      ','3-Parameter Weibull= ',num2str(lwei),'      ','MaxEnt= ',num2str(lmaxent)]);

disp('------------------------- Model selection based on Akaike Information Criterion (AIC) -------------------------')

disp('AIC values for models are:');
disp(['Lognormal= ',num2str(aic(1)),'      ','3-Parameter Weibull= ',num2str(aic(2)),'      ','MaxEnt= ',num2str(aic(3))]);

if amin==aic(1)
    disp('Best model based on AIC is Lognormal distribution');
    disp(['3-Parameter Weibull is as ',num2str(rlaic(2)),' times probable as Lognormal to minimize information loss']);
    disp(['MaxEnt is as ',num2str(rlaic(3)),' times probable as Lognormal to minimize information loss']);
end
if amin==aic(2)
    disp('Best model based on AIC is 3-Parameter Weibull distribution');
    disp(['Lognormal is as ',num2str(rlaic(1)),' times probable as 3-Parameter Weibull to minimize information loss']);
    disp(['MaxEnt is as ',num2str(rlaic(3)),' times probable as 3-Parameter Weibull to minimize information loss']);
end
if amin==aic(3)
    disp('Best model based on AIC is MaxEnt distribution');
    disp(['Lognormal is as ',num2str(rlaic(1)),' times probable as MaxEnt to minimize information loss']);
    disp(['3-Parameter Weibull is as ',num2str(rlaic(2)),' times probable as MaxEnt to minimize information loss']);
end

disp('------------------------- Model selection based on Bayesian Information Criterion (BIC) -------------------------')

disp('BIC values for models are:');
disp(['Lognormal= ',num2str(bic(1)),'      ','3-Parameter Weibull= ',num2str(bic(2)),'      ','MaxEnt= ',num2str(bic(3))]);

if bmin==bic(1)
    disp('Best model based on BIC is Lognormal distribution');
    disp(['BIC difference between 3-Parameter Weibull and Lognormal is ',num2str(deltabic(2))]);
    disp(['BIC difference between MaxEnt and Lognormal is ',num2str(deltabic(3))]);
end
if bmin==bic(2)
    disp('Best model based on BIC is 3-Parameter Weibull distribution');
    disp(['BIC difference between Lognormal and 3-Parameter Weibull is ',num2str(deltabic(1))]);
    disp(['BIC difference between MaxEnt and 3-Parameter Weibull is ',num2str(deltabic(3))]);
end
if bmin==bic(3)
    disp('Best model based on BIC is MaxEnt distribution');
    disp(['BIC difference between Lognormal and MaxEnt is ',num2str(deltabic(1))]);
    disp(['BIC difference between 3-Parameter Weibull and MaxEnt is ',num2str(deltabic(2))]);
end

disp('------------------------- Model selection based on Maximum Likelihood (Bayesian Model Comparison) criteria -------------------------')

if lmin==lnorm
    disp(['3-Parameter Weibull is prefered to Lognormal with K value of ',num2str(k(2))]);
    disp(['MaxEnt is prefered to Lognormal with K value of ',num2str(k(3))]);
end
if lmin==lwei
    disp(['Lognormal is preferred to 3-Parameter Weibull with K value of ',num2str(k(1))]);
    disp(['MaxEnt is preferred to 3-Parameter Weibull with K value of ',num2str(k(3))]);
end
if lmin==lmaxent
    disp(['Lognormal is prefered to MaxEnt with K value of ',num2str(k(1))]);
    disp(['3-Parameter Weibull is prefered to MaxEnt with K value of ',num2str(k(2))]);
end