%-----------------------------------------------------------------%
%this Code is to prove the orthogonal property of Bessel functions
%provided by Jiaqi, email:Jiaqi_Wang@sjtu.edu.cn
%2020-03-06
%-----------------------------------------------------------------%

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
jmn_pm=[fliplr(jmn) jmn(:,2:end)];
BesselValue_pm=[fliplr(repmat((-1).^m,n,1).*BesselValue) BesselValue(:,2:end)];



%Orthogonal vertified-Rienstra:
Nmn=sqrt(2)./(BesselValue.*sqrt(1-m.^2./jmn.^2));
for km=1:length(m)
    for kv=1:n
      temp(kv,km)=sum(chebfun(@(r) besselj(m(km),jmn(kv,km)*r)  *  besselj(m(km),jmn(2,km)*r)*r  *  Nmn(2,km) * Nmn(kv,km),[0,1]));
    end
end

%Orthogonal vertified-James:
h=1;

m2=-10:10;
Cmn=(i).^m2./(sqrt(pi)*h*BesselValue_pm.*sqrt(1-m2.^2./jmn_pm.^2));

% for km=1:length(m2)
%     for ku=1:
%         for kn=1:length(m2)
%             for kv=1:n
%       temp3(kv,km)=sum(chebfun(@(theta) exp(sqrt(-1)*m2(km)*theta)*exp(sqrt(-1)*(-m2(km))*theta),[0,2*pi]));
%       temp4(kv,km)=sum(chebfun(@(r) besselj(m2(km),jmn(kv,km)*r)*besselj(-m2(km),jmn(2,km)*r)*r*Cmn(2,km)*Cmn(kv,km),[0,1]));
%       temp5=temp3.*temp4;
%             end
%         end
%     end
% end

for km1=1:length(m2)
        for km2=1:length(m2)
        temp3(km1,km2)=sum(chebfun(@(theta) exp(sqrt(-1)*m2(km1)*theta)*exp(sqrt(-1)*(m2(km2))*theta),[0,2*pi]));
        end
end
figure(1);image(real(temp3),'CDataMapping','scaled');

%same m, diff n
%proved that: m=2n1+1 * m=2n2 =0, m=2n1 * m=2n2 !=0, except for n1=n2
for km=1:length(m2)
for kn1=1:n
    for kn2=1:n
        temp4{km}(kn1,kn2)=sum(chebfun(@(r) besselj(m2(km),jmn_pm(kn1,km)*r)*besselj(m2(km),jmn_pm(kn2,km)*r)*r*Cmn(kn1,km)*Cmn(kn2,km),[0,1]));
    end
end
%figure(km);image(real(temp4{km}),'CDataMapping','scaled');
end

%same n, diff m;
%non-orthogonal!!
for kn=1:n
for km1=1:length(m2)
    for km2=1:length(m2)
        temp5{kn}(km1,km2)=sum(chebfun(@(r) besselj(m2(km1),jmn_pm(kn,km1)*r)*besselj(m2(km2),jmn_pm(kn,km2)*r)*r*Cmn(kn,km1)*Cmn(kn,km2),[0,1]));
    end
end
figure;image(real(temp5{kn}),'CDataMapping','scaled');
end



