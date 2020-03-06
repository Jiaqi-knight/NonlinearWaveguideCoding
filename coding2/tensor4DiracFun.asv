%-----------------------------------------------------------------%
%this Code is to prove the tensor calculataion with Dirac's funcrion
%reconstrcut: \delta_{m,-n}\delta_{u,v}=[I] for easy our model
%provided by Jiaqi, email:Jiaqi_Wang@sjtu.edu.cn
%2020-03-06
%-----------------------------------------------------------------%

clc
clear
close all
subfunction_path1='.\tensor_toolbox-master'
addpath(genpath(subfunction_path1));
subfunction_path2='.\chebfun-master'
addpath(genpath(subfunction_path2));

m=-2:2;u=1:3;
delta_mn=eye(length(m));
delta_m0n=fliplr(delta_mn);
delta_uv=eye(length(u));
delta_u0v=fliplr(delta_uv);

%kron_fliplr the second line in order to form I matrix
delta_mn_u0v=kron_fliplr(delta_mn,delta_u0v);
delta_m0n_uv=kron_fliplr(delta_m0n,delta_uv); 



%we began to consrtruct 4D tensor
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
h=1;
m2=-10:10;
Cmn=(i).^m2./(sqrt(pi)*h*BesselValue_pm.*sqrt(1-m2.^2./jmn_pm.^2));


for km1=1:length(m2)
for km2=1:length(m2)
for kn1=1:n
    for kn2=1:n
        %cost much time
        temp4(km1,km2,kn1,kn2)=sum(chebfun(@(r) besselj(m2(km1),jmn_pm(kn1,km1)*r)*besselj(m2(km2),jmn_pm(kn2,km2)*r)*r*Cmn(kn1,km1)*Cmn(kn2,km2),[0,1]));
    end
end
end
end



figure;image(real(reshape(temp4(2,2,:,:),4,4)),'CDataMapping','scaled');


