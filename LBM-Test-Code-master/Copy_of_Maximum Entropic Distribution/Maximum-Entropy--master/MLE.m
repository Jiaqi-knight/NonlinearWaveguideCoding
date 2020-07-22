% This code is developed by Foad Karimian
% This part of code calculates MLE estimations of parameters for Lognormal
% and 3-paramter Weibull distributions

%Actual data in form of Log(cycles) is entered in variable "data"
data=[4.139438274	4.151001908	4.16708092	4.176091259	4.183269844	4.185627136	4.204119983	4.208629438	4.21181443	4.220500346	4.224610747	4.232767474	4.244598701	4.251540889	4.260071388	4.264227348	4.291435455	4.297760511	4.302698819	4.312600439	4.326724913	4.336299596	4.336299596	4.370235437	4.383222742	4.411619706	4.411922596	4.416640507	4.424881637	4.439932396];

Lognorm=mle(data);      %MLE estimation of mean and std
Lognorm=Lognorm.';      %make a vector of dist. parameters
disp(["mean:",num2str(Lognorm(1));      %display parameters
     "std:",num2str(Lognorm(2))]);

xmin=-5.*Lognorm(2)+Lognorm(1);     %set x interval with respect to std, extending from 4*std before mean to 4*std after mean
xmax=5.*Lognorm(2)+Lognorm(1);
x=[xmin:0.01:xmax];
dx=0.01
normp = normpdf(x,Lognorm(1),Lognorm(2));    %generate pdf of lognormal distribution
normc = normcdf(x,Lognorm(1),Lognorm(2));    %generate cdf of lognormal distribution

figure(3);                          %draw pdf lognormal distribution
plot(x,normp)
title('Lognormarl PDF')
xlabel('data')
ylabel('PDF')

figure(4);                          %draw cdf lognormal distribution
plot(x,normc)
title('Lognormarl CDF')
xlabel('data')
ylabel('CDF')

movegui(figure(3),'northwest')      %place lognormal pdf on top-left on screen 
movegui(figure(4),'north')      %place lognormal cdf on top-center on screen 

custpdf = @(x,a,b,c) (x>c).*(b/a).*(((x-c)/a).^(b-1)).*exp(-((x-c)/a).^b);  %PDF for 3-parameter Weibull
opt = statset('MaxIter',1e5,'MaxFunEvals',1e5,'FunValCheck','off');
weibull= mle(data,'pdf',custpdf,'start',[0.2 2.4 4],'Options',opt,...       %MLE estimation of parameters
    'LowerBound',[0 0 -Inf],'UpperBound',[Inf Inf min(data)]);

weibull=weibull.';                  %make a vector of parameters

disp(["alpha:",num2str(weibull(1)); %display parameters
     "beta:",num2str(weibull(2));
     "Log(n0)",num2str(weibull(3))]);
weipdf=custpdf(x,weibull(1),weibull(2),weibull(3));    %generate PDF of 3-parameter weibull distribtion

weicdf=zeros(length(x),1);
weicdf(1)=weipdf(1).*dx;

for i=2:length(x)
    weicdf(i)=weipdf(i).*dx+ weicdf(i-1);                   %generate CDF of 3-parameter weibull distribtion
end

figure(5)           %draw pdf 3-parameter weibull distribution
plot(x,weipdf)
title('3 Parameter Weibull PDF')
xlabel('data')
ylabel('PDF')

figure(6)           %draw cdf 3-parameter weibull distribution
plot(x,weicdf)
title('3 Parameter Weibull CDF')
xlabel('data')
ylabel('CDF')

movegui(figure(5),'west')      %place figure on left-center on screen
movegui(figure(6),'center')      %place figure on center on screen

figure(7)
plot(x,normp,x,weipdf)
title('Lognormal vs 3-p Weibull PDF')
legend('Lognormal','3 Paramter Weibull')

figure(8)
plot(x,normc,x,weicdf)
title('Lognormal vs 3-p Weibull CDF')
legend('Lognormal','3 Paramter Weibull')

movegui(figure(7),'northeast')      %place figure on upper-right on screen
movegui(figure(8),'east')      %place figure on right on screen
