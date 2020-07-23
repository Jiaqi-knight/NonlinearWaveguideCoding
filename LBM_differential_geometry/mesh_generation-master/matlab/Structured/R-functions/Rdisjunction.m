function[x]=Rdisjunction(x1,x2,alpha,m)

x = (1./(1+alpha)).*(x1+x2+sqrt(x1.^2+x2.^2-2*alpha.*x1.*x2)).*(x1.^2+x2.^2).^(0.5*m);

return