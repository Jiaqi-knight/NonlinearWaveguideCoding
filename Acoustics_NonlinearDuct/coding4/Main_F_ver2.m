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
ns=101;%odd
Geo.s =linspace(0,5,ns);
or_n=1:2:ns;or_r=2:2:ns;
%Geo.h=0.1*exp(linspace(0,1.5,length(Geo.s)));
Geo.h=1*ones(size(Geo.s));

Geo.tau=0.0./Geo.h;
%Geo.tau=linspace(0,4,500)
%Geo.tau =logspace(0,0.6,ns)-1;
%Geo.tau =0.0*ones(size(Geo.s));
%Geo.tau =logspace(0,0.1,500)-1;

Geo.kappa=(0)./Geo.h;
%Geo.kappa=logspace(0,0.6,500)-1

Geo.sw=sqrt(Geo.kappa.^2+Geo.tau.^2).*Geo.s;
Geo.x = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*sin(Geo.sw+0);Geo.y = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*cos(Geo.sw+0);Geo.z = Geo.tau./(Geo.kappa.^2+Geo.tau.^2).*Geo.sw;
Geo.theta_0=cumsum(Geo.tau.*[0 diff(Geo.s)]);
Geo.h_diff=[(Geo.h(2)-Geo.h(1))/(Geo.s(2)-Geo.s(1)) diff(Geo.h)./diff(Geo.s)];
h=figure
tubeplot(Geo.x,Geo.y,Geo.z,Geo.h,Geo.s,50);hold on;plot3(Geo.x, Geo.y, Geo.z);daspect([1,1,1]); camlight;
saveas(h,[save_directory,'\','Geo','.png'])


% %% #######Wave########%
Geo.m=[-6:6];
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
beta=1+(Wave.gamma-1)/2;
M=0.1;
Wave.k
sigma=Geo_b.s./(beta*M*Wave.k);

L=multiprod(Out.N_a,Ad.Y_a,[1,2])+Out.H_a%-permute(bsxfun(@times,eye(length(Geo.m)*Geo.n),reshape(beta*M*Wave.k.*Wave.a.^2./max(Wave.a),1,1,length(Wave.a))),[1,2,4,3]);
L1=L(:,:,1:end-1,:);
ds=diff(Geo_b.s);
exp1=bsxfun(@times,L1,reshape(ds.'*ones(size(Wave.a)),1,1,length(ds),length(Wave.a)));%expm

for ka=1:length(Wave.a)
    for ks=1:length(ds)
        L_inv(:,:,ks,ka)=inv(L1(:,:,ks,ka));
        Exp1(:,:,ks,ka)=expm(exp1(:,:,ks,ka));
    end
end
Exp2=multiprod(L_inv,(Out.I2(:,:,1:end-1,:)-Exp1),[1,2]);

%Initial pressure modes
or_ab=max(Wave.a)-min(Wave.b);
P_list=bsxfun(@times,zeros(length(Geo_b.m)*Geo_b.n,or_ab*2+1),reshape(ones(size(Geo_b.s)),1,1,length(Geo_b.s)));
kka=1;
%p_intial=[0 0 0.5 0 0 0];
P_list(7,or_ab+kka+1,1)=M;
P_list(:,or_ab-kka+1,1)=conj(P_list(:,or_ab+kka+1,1));
tic
for  ks=1:length(ds)
    
    P_a_b=permute(reshape(P_list(:,Wave.a_b+or_ab+1,ks),length(Geo_b.m)*Geo_b.n,length(Wave.b),length(Wave.a)),[1,4,5,6,2,3]);
    P_b=permute(bsxfun(@times,P_list(:,Wave.b+or_ab+1,ks),reshape(ones(size(Wave.a)),1,1,length(Wave.a))),[4,1,5,6,2,3]);
    
    term1=multiprod(permute(bsxfun(@times,Out.N_a(:,:,ks,:),reshape(ones(size(Wave.b)),1,1,1,1,length(Wave.b))),[1,2,6,3,5,4]),...
        multiprod(multiprod(Ad.YY_ab(:,:,:,ks,:,:),P_a_b,[1,2]),P_b,[2,3]),[1,2]);
    term2=multiprod(multiprod(Out.C(:,:,:,ks,:,:),multiprod(permute(Ad.Y_a_b(:,:,ks,:,:),[1,2,6,3,4,5]),P_a_b,[1,2]),[1,2]),P_b,[2,3]);
    term3=multiprod(multiprod(Out.D(:,:,:,ks,:,:),multiprod(permute(Ad.Y_a_b(:,:,ks,:,:),[1,2,6,3,4,5]),P_a_b,[1,2]),[1,2]),...
        multiprod(permute(bsxfun(@times,Ad.Y_b(:,:,ks,:),reshape(ones(size(Wave.a)),1,1,1,1,length(Wave.a))),[6,1,2,3,4,5]),P_b,[2,3]),[2,3]);
    
    temp_p=permute(multiprod(Exp1(:,:,ks,:),permute(P_list(:,Wave.a+or_ab+1,ks),[1,3,4,2]),[1,2])...
        -multiprod(Exp2(:,:,ks,:),permute(sum(term1+term2+term3,5),[1,2,4,6,5,3]),[1,2]),[1,4,2,3]);
    %temp_p=permute(multiprod(Exp1(:,:,ks,:),permute(P_list(:,Wave.a+or_ab+1,ks),[1,3,4,2]),[1,2]),[1,4,2,3]);
    P_list(:,or_ab+1-max(Wave.a):or_ab+max(Wave.a)+1,ks+1)=[fliplr(conj(temp_p)) zeros(length(Geo_b.m)*Geo_b.n,1) temp_p];
    
end
toc


figure;
plot(abs((permute(P_list(:,8,:),[1,3,2]))).')
