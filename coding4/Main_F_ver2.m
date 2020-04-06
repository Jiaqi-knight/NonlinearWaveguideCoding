%Main-Function
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
%-----------------------------------------------------------------%
%For Boudnary, only end point is need to solved!
clc
clear
close all
opengl('save','software')
subfunction_path1='.\subfunction1'
addpath(genpath(subfunction_path1));
subfunction_path2='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding2\src\chebfun-master';
addpath(genpath(subfunction_path2));
subfunction_path3='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding2\src\tensor_toolbox-master';
addpath(genpath(subfunction_path3));
subfunction_path4='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding4\subfunction1\Multiprod_2009';
addpath(genpath(subfunction_path4));
save_directory='./data';
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('ÎÄ¼þ¼Ð´æÔÚ£¡');
end
% %% #######Geometry########%
ns=51;%odd
Geo.s =linspace(0,1,ns);
or_n=1:2:ns;or_r=2:2:ns;
%Geo.h=0.1*exp(linspace(0,1.5,length(Geo.s)));
Geo.h=1*ones(size(Geo.s));

Geo.tau=0.0./Geo.h;
%Geo.tau=linspace(0,4,500)
%Geo.tau =logspace(0,0.6,50)-1;
%Geo.tau =0.0*ones(size(Geo.s));
%Geo.tau =logspace(0,0.1,500)-1;

Geo.kappa=(0.0001)./Geo.h;
%Geo.kappa=logspace(0,0.6,500)-1

Geo.sw=sqrt(Geo.kappa.^2+Geo.tau.^2).*Geo.s;
Geo.x = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*sin(Geo.sw+0);Geo.y = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*cos(Geo.sw+0);Geo.z = Geo.tau./(Geo.kappa.^2+Geo.tau.^2).*Geo.sw;
Geo.theta_0=cumsum(Geo.tau.*[0 diff(Geo.s)]);
Geo.h_diff=[(Geo.h(2)-Geo.h(1))/(Geo.s(2)-Geo.s(1)) diff(Geo.h)./diff(Geo.s)];
h=figure
tubeplot(Geo.x,Geo.y,Geo.z,Geo.h,Geo.s,50);hold on;plot3(Geo.x, Geo.y, Geo.z);daspect([1,1,1]); camlight;
saveas(h,[save_directory,'\','Geo','.png'])


% %% #######Wave########%
Geo.m=[-1:1];
Geo.n=1;
Wave.a=[1 2 3]; %P^{a}=P^{a*},U^{a}=U^{a*}
Wave.b=[-3 -2 -1 1 2 3];
Wave.a_b=Wave.a-Wave.b.';
Wave.k=0.95*1.8412/Geo.h(end)+0.00001*i; %+-;
Wave.gamma=1.4;

%%
Geo_b.m=Geo.m;
Geo_b.n=Geo.n;
Geo_b.s=Geo.s(or_n);                  %[Geo.s(1) Geo.s(end)];
Geo_b.h=Geo.h(or_n);                  %[Geo.h(1) Geo.h(end)];
Geo_b.kappa=Geo.kappa(or_n);          %[Geo.kappa(1) Geo.kappa(end)];
Geo_b.tau=Geo.tau(or_n);              %[Geo.tau(1) Geo.tau(end)];
Geo_b.h_diff=Geo.h_diff(or_n)         %[Geo.h_diff(1) Geo.h_diff(end)]
Geo_b.theta_0=Geo.theta_0(or_n);

RK_b.m=Geo.m;
RK_b.n=Geo.n;
RK_b.s=Geo.s(or_r);                  %[Geo.s(1) Geo.s(end)];
RK_b.h=Geo.h(or_r);                  %[Geo.h(1) Geo.h(end)];
RK_b.kappa=Geo.kappa(or_r);          %[Geo.kappa(1) Geo.kappa(end)];
RK_b.tau=Geo.tau(or_r);              %[Geo.tau(1) Geo.tau(end)];
RK_b.h_diff=Geo.h_diff(or_r)         %[Geo.h_diff(1) Geo.h_diff(end)]
RK_b.theta_0=Geo.theta_0(or_r);

%% #######Name########%
Name.op_m=  {'1','r','r2'};
Name.op=    {'ab','pr_ab','a_pr_b','ps_ab','a_ps_b'};
Name.X2=    {'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}','X_{\alpha\{\beta\}}'};
Name.X3=    {'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}','X_{\alpha\{\beta\}\gamma}'};
Name.T=     {'ab','pt_ab','a_pt_b','ps_ab','a_ps_b','ab_cos','ab_sin','pt_ab_cos','pt_ab_sin','ps_ab_cos','ps_ab_sin'};
Name.T2=    {'\Theta_{\alpha\beta}','\Theta_{(\alpha)\beta}','\Theta_{\alpha(\beta)}','\Theta_{\{\alpha\}\beta}','\Theta_{\alpha\{\beta\}}','\Theta_{\alpha\beta}[cos\phi]',...
    '\Theta_{\alpha\beta}[sin\phi]','\Theta_{(\alpha)\beta}[cos\phi]','\Theta_{(\alpha)\beta}[sin\phi]','\Theta_{\{\alpha\}\beta}[cos\phi]','\Theta_{\{\alpha\}\beta}[sin\phi]'};
Name.T3=    {'\Theta_{\alpha\beta\gamma}','\Theta_{(\alpha)\beta\gamma}','\Theta_{\alpha(\beta)\gamma}','\Theta_{\{\alpha\}\beta\gamma}','\Theta_{\alpha\{\beta\}\gamma}','\Theta_{\alpha\beta\gamma}[cos\phi]',...
    '\Theta_{\alpha\beta\gamma}[sin\phi]','\Theta_{(\alpha)\beta\gamma}[cos\phi]','\Theta_{(\alpha)\beta\gamma}[sin\phi]','\Theta_{\{\alpha\}\beta\gamma}[cos\phi]','\Theta_{\{\alpha\}\beta\gamma}[sin\phi]'};
%% #######Function########%
[Psi]= PFun(Geo_b)
[ePsi]= PFun(RK_b)
[Out]= EFun(Geo_b,Psi,Wave,1)
[OutRK]= EFun(RK_b,ePsi,Wave,0)
[Ad]= RKFUN(Geo_b,RK_b,Out,OutRK,Wave,or_n,or_r) 
%%

%Initial pressure modes

%suppose:
%a=1
P0=zeros(length(Geo_b.m)*Geo_b.n,1,1,length(Wave.a));

P0(2,1,1,1)=1

L=multiprod(Out.N_a,Ad.Y_a,[1,2])+Out.H_a;
L1=L(:,:,1:end-1,:)
ds=diff(Geo_b.s)
for ka=1:length(Wave.a)
for ks=1:length(ds)
    L_inv(:,:,ks,ka)=inv(L1(:,:,ks,ka));
end
end

exp1=exp(bsxfun(@times,L1,reshape(ds.'*ones(size(Wave.a)),1,1,length(ds),length(Wave.a))))  

exp2=multiprod(L_inv,(Fun2_a.I2(:,:,1:end-1,:)-exp1),[1,2])


for ka=1:length(Wave.a)
for  ks=1:length(ds)

    P0(:,:,ks+1,ka)=exp1(:,:,ks,ka)*P0(:,:,ks,ka)-exp2(:,:,ks,ka)*P0(:,:,ks,ka);

end
end
Ad.YY_ab


Wave.a_b

Wave.b