clc
clear
close all
subfunction_path1='.\chebfun-master'
addpath(genpath(subfunction_path1));

m=0:10;
n=4;
for k=1:length(m)
temp1=roots(diff(chebfun(@(t) besselj(m(k),t),[0,300]))+0.00001);
temp2= besselj(m(k),temp1);
jmn(:,k)=temp1(1:n);
BesselValue(:,k)=temp2(1:n);
end



%Orthogonal vertified-Rienstra:
Nmn=sqrt(2)./(BesselValue.*sqrt(1-m.^2./jmn.^2));
for km=1:length(m)
    for kn=1:n
      temp(kn,km)=sum(chebfun(@(r) besselj(m(km),jmn(kn,km)*r)  *  besselj(m(km),jmn(2,km)*r)*r  *  Nmn(2,km) * Nmn(kn,km),[0,1]));
    end
end

%Orthogonal vertified-James:
h=1;
Cmn=sqrt((-1).^m)./(sqrt(pi)*h*BesselValue.*sqrt(1-m.^2./jmn.^2));
for km=1:length(m)
    for kn=1:n
      temp3(kn,km)=sum(chebfun(@(theta) exp(sqrt(-1)*m(km)*theta)*exp(sqrt(-1)*(-m(1))*theta),[0,2*pi]));
      temp4(kn,km)=sum(chebfun(@(r) besselj(m(km),jmn(kn,km)*r)*besselj(-m(km),jmn(2,km)*r)*r*Cmn(2,km)*Cmn(kn,km),[0,1]));
    end
end
temp5=temp3.*temp4;

