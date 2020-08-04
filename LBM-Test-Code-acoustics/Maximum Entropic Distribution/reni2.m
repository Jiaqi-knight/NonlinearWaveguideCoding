function [lambda,p,beta,conditions] = reni2(mu,xprime,alpha)
%UNTITLED Summary of this function goes here
% Detailed explanation goes here
%p has the values of PDF
Eps=1e-4;
xmin=xprime(1); %in order to avoid conflict we defined xprime instead of x
xmax=xprime(length(xprime));
dx=xprime(2)-xprime(1);
lambda=zeros(size(mu));
n=length(lambda);
phi2=fin_3(xprime,n);
r=1;
[row,column]=size(xprime);
finalx=xprime;
phin=zeros(n,column);
muprime=zeros(n,column);
lambdaprime=zeros(n,column);
phi1=zeros(n,1);
not=0;
y=0;
for myi=1:length(xprime) %Calculation of muprime matrix which is n*r
    muprime(:,myi)=mu';
end %End of muprime calculation
for ex=xmin:dx:xmax %start of loop for calculating Phi-en
    for i=1:n %caculation of matrix phi1 in order to use in sigma
        phi1(i,1)=phi2{i,1}(ex);
    end %end of phi1 calculation
    phin(:,r)=phi1(:); %Calculation of Phi-en matrix which is r*n
    r=r+1;
end %end of loop for calculating Phi-en
while 1 %first of while loop
    for myi=1:length(xprime) %Calculation of lambdaprime matrix which is n*r
        lambdaprime(:,myi)=lambda';
    end %End of lambdaprime calculation
    if (n==1)
        lmphi=((1-alpha)/alpha)*(lambdaprime.*(muprime-phin));
        %Calculation of sum of ((1-alpha)/alpha)*lambda.*(mu-phi) in which each column stands for single x
    else
        lmphi=((1-alpha)/alpha)*sum(lambdaprime.*(muprime-phin));
    end
    beta=dx*sum((1-lmphi).^(1/(alpha-1))); %calculation of beta
    p=1/beta*((1-lmphi).^(1/(alpha-1))); %Calculation of P(x)
    a=(1-lmphi).^((1/(alpha-1))-1);
    g=zeros(1,n);
    for myi=1:n %calculation of g(1,i)
        g(1,myi)=dx*sum(phin(myi,:).*p);
    end %for end
    if(isnan(g)|isinf(g)|(~isreal(g))) %If any element in g was not a number or was infinite or was complex the
        % calculation is broken.
        not=1;
        break
    end
    gmk=zeros(n,n);
    for m=1:n
        for k=1:n
            first=(((-(1-lmphi)/beta)*sum((muprime(k,:)-phin(k,:)).*a)+(muprime(k,:)-phin(k,:))).*a)/(alpha*beta);
            gmk=dx*sum(first.*phin(m,:));
        end
    end %end of forming g(i,j)
    v=(mu-g)';
    delta=gmk\v;
    lambda=lambda+delta';
    if(abs(delta./(lambda'))<eps)
        for myi=1:length(finalx) %Calculation of lambdaprime matrix which is n*r
            lambdaprime(:,myi)=lambda';
        end %End of lambdaprime calculation
        if (n==1)
            lmphi=((1-alpha)/alpha)*(lambdaprime.*(muprime-phin));
            %Calculation of sum of ((1-alpha)/alpha)*lambda.*(mu-phi) in which each column stands for a single x
        else
            lmphi=((1-alpha)/alpha)*sum(lambdaprime.*(muprime-phin));
        end %Calculation of sum of ((1-alpha)/alpha)*lambda.*(mu-phi) in which each column stands for a single x
        beta=dx*sum((1-lmphi).^(1/(alpha-1))); %calculation of beta
        p=1/beta*((1-lmphi).^(1/(alpha-1))); %Calculation of P(x)
        break,
    end %End of delta if
end%while end
if(not==0)
    p_integral=dx*sum(p);
    conditions=zeros(1,n);
    for myj=1:n
        conditions(1,myj)=dx*sum(phin(myj,:).*p);
    end
    conditions=[p_integral,conditions]; %Conditions shows how correct our calculations has been
    if(size(p)==size(finalx))
        plot(finalx,p)
        disp('End of the program');
    else
        finalx(end)=[];
        plot(finalx,p)
        disp('End of the program');
    end
else
    disp('No answer exists')
end
end % function end