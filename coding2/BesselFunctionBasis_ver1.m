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

Nmn=sqrt(2)./(BesselValue.*sqrt(1-m.^2./jmn.^2));

%Orthogonal vertified:
for km=1:length(m)
    for kn=1:n
      temp(kn,km)=sum(chebfun(@(r) besselj(m(km),jmn(kn,km)*r)^2*r*Nmn(kn,km)^2,[0,1]));
    end
end
