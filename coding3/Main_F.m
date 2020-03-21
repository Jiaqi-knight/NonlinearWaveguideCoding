%Main-Function
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
%-----------------------------------------------------------------%
clc
clear
close all
%opengl('save','software')
load('Database_X2.mat');load('Database_X3.mat');
subfunction_path1='.\subfunction1'
addpath(genpath(subfunction_path1));
subfunction_path2='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding2\src\chebfun-master';
addpath(genpath(subfunction_path2));
subfunction_path3='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding2\src\tensor_toolbox-master';
addpath(genpath(subfunction_path3));
subfunction_path4='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding3\subfunction1\Multiprod_2009';
addpath(genpath(subfunction_path4));

%% #######Geometry########%
Geo.s =logspace(0,1,50);
Geo.h=0.1*exp(linspace(0,1.5,length(Geo.s)));
%Geo.h=0.1*ones(size(Geo.h));

Geo.kappa=(2/3)./Geo.h;Geo.tau=0.2./Geo.h;
Geo.sw=sqrt(Geo.kappa.^2+Geo.tau.^2).*Geo.s;
Geo.x = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*sin(Geo.sw+0);Geo.y = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*cos(Geo.sw+0);Geo.z = Geo.tau./(Geo.kappa.^2+Geo.tau.^2).*Geo.sw;
Geo.theta_0=cumsum(Geo.tau.*[0 diff(Geo.s)]);
%tubeplot(x,y,z,h,s,50);hold on;plot3(x, y, z);daspect([1,1,1]); camlight;

%% #######Wave########%
Geo.m=-1:1;
Geo.n=1;
a=[1:2]; %P^{a}=P^{a*},U^{a}=U^{a*}
b=[-3:3];
k=4/norm(Geo.h);
%% #######Name########%
Name.op_m={'1','r','r2'};
Name.op={'ab','pr_ab','a_pr_b','ps_ab','a_ps_b'};
Name.X2={'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}','X_{\alpha\{\beta\}}'};
Name.X3={'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}','X_{\alpha\{\beta\}\gamma}'};
Name.T={'ab','pt_ab','a_pt_b','ps_ab','a_ps_b','ab_cos','ab_sin','pt_ab_cos','pt_ab_sin','ps_ab_cos','ps_ab_sin'};
Name.T2={'\Theta_{\alpha\beta}','\Theta_{(\alpha)\beta}','\Theta_{\alpha(\beta)}','\Theta_{\{\alpha\}\beta}','\Theta_{\alpha\{\beta\}}','\Theta_{\alpha\beta}[cos\phi]',...
    '\Theta_{\alpha\beta}[sin\phi]','\Theta_{(\alpha)\beta}[cos\phi]','\Theta_{(\alpha)\beta}[sin\phi]','\Theta_{\{\alpha\}\beta}[cos\phi]','\Theta_{\{\alpha\}\beta}[sin\phi]'};
Name.T3={'\Theta_{\alpha\beta\gamma}','\Theta_{(\alpha)\beta\gamma}','\Theta_{\alpha(\beta)\gamma}','\Theta_{\{\alpha\}\beta\gamma}','\Theta_{\alpha\{\beta\}\gamma}','\Theta_{\alpha\beta\gamma}[cos\phi]',...
    '\Theta_{\alpha\beta\gamma}[sin\phi]','\Theta_{(\alpha)\beta\gamma}[cos\phi]','\Theta_{(\alpha)\beta\gamma}[sin\phi]','\Theta_{\{\alpha\}\beta\gamma}[cos\phi]','\Theta_{\{\alpha\}\beta\gamma}[sin\phi]'};

%% #######Function########%

[Base]=BaseJ(Geo.m,Geo.n,Geo.h);
tic
%2D-bsxfun(@times, 3D(\alpha*\beta*s), 1*1*s*a)-->4D[\alpha*\beta*s*a]
Fun.I2=bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ab'),reshape(ones(size(a)),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun.V=bsxfun(@times,X(Base,Geo,2,'a_pr_b','r').*Theta(Geo,2,'ab'),reshape(1./(sqrt(-1)*k*a),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun.W=bsxfun(@times,X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab'),reshape(-1./(sqrt(-1)*k*a),1,1,1,length(a))); %(James-3.27,Jiaqi-77)
Fun.N=bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ab'),reshape((sqrt(-1)*k*a),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','r2').*Theta(Geo,2,'ab_cos'),reshape((sqrt(-1)*k*Geo.kappa.'.*a),1,1,length(Geo.kappa),length(a)));%(James-3.35b)
Fun.G=-bsxfun(@times,X(Base,Geo,2,'ps_ab','r').*Theta(Geo,2,'ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ps_ab'),reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun.H=-bsxfun(@times,X(Base,Geo,2,'a_ps_b','r').*Theta(Geo,2,'ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'a_ps_b'),reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun.M_2_1=(bsxfun(@times,X(Base,Geo,2,'pr_ab','r').*Theta(Geo,2,'ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'pr_ab','r2').*Theta(Geo,2,'ab_cos'),reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a))));
Fun.M_3_1=(bsxfun(@times,X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'pt_ab_cos'),reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a))));
Fun.M=-Fun.N-mtimesx(Fun.M_2_1,Fun.V)-mtimesx(Fun.M_3_1,Fun.W);%(James-3.35a)
%3D-bsxfun(@times, 4D(\alpha*\beta*\gamma*s), 1*1*1*s*b*a)-->6D[\alpha*\beta*gamma*s*b*a]
Fun.I3=bsxfun(@times,X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab'),...
    reshape(bsxfun(@times,ones(size(Geo.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun.A_1=bsxfun(@times,X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab'),...
    reshape(bsxfun(@times,ones(size(Geo.s.'))*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)))...
    -bsxfun(@times,X(Base,Geo,3,'ab','r2').*Theta(Geo,3,'ab_cos'),...
    reshape(bsxfun(@times,Geo.kappa.'*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun.A_2=permute(bsxfun(@times,X(Base,Geo,2,'pr_ab','r').*Theta(Geo,2,'ab'),...
    reshape(ones(size(b.'))*ones(size(a)),1,1,1,length(b),length(a))),[6,1,2,3,4,5]);
for i=1:length(Geo.s);for j=1:length(a);Fun.N_inv(:,:,i,j)=Fun.N(:,:,i,j)^(-1);end;end;
Fun.N_inv=permute(bsxfun(@times,Fun.N_inv,reshape(b,1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun.A_2_1_1= specialFun(Geo.s,Geo.h,Geo.kappa,Geo.m,Geo.n,Base.Cmn1,Base.jmn_pm,'hh`^2/[1-\kappa*h*cos\psi]');
Fun.A_2_1_2= specialFun(Geo.s,Geo.h,Geo.kappa,Geo.m,Geo.n,Base.Cmn1,Base.jmn_pm,'h(1-\kappa*h*cos\psi)');
Fun.A_2_1_3=bsxfun(@times,X(Base,Geo,2,'ab','1').*Theta(Geo,2,'ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ab_cos'),reshape((Geo.kappa.'.*ones(size(a))),1,1,length(Geo.kappa),length(a)));
Fun.A_2_1_4=Fun.M_2_1;
Fun.A_2_1=permute(bsxfun(@times,Fun.A_2_1_1+Fun.A_2_1_2+Fun.A_2_1_3+Fun.A_2_1_4,reshape(b,1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun.A_2_2=bsxfun(@times,X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab_cos'),...
    reshape(bsxfun(@times,Geo.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun.A_3=bsxfun(@times,X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','r').*Theta(Geo,2,'pt_ab_cos'),reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a)));
Fun.A_3=permute(bsxfun(@times,Fun.A_3,reshape(b,1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun.A_3_1_2=bsxfun(@times,X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab'),reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab_cos'),reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a)));
Fun.A_3_1_2=permute(bsxfun(@times,Fun.A_3_1_2,reshape(b,1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun.A_3_2=bsxfun(@times,X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab_sin'),...
    reshape(bsxfun(@times,Geo.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun.A = -Fun.A_1-multiprod(multiprod(Fun.A_2,Fun.N_inv,[2,3]),-multiprod(Fun.I3,Fun.A_2_1,[2,3]),[2,3])...
       -multiprod(multiprod(Fun.A_3,Fun.N_inv,[2,3]),(multiprod(Fun.I3,Fun.A_3_1_2,[2,3])+Fun.A_3_2),[2,3]);%(James-3.35e)

   
   
   
size(Fun.A_2)

multiprod(Fun.I3,Fun.A_2_1,[2,3],[2,3]);

size(Fun.I3)
size(Fun.A_2_1)
size(R)
A=rand(1,3,4,2,3)
B=rand(1,2,3,3)
multiprod(bsxfun(@times,rand(3,4,2),reshape(ones(1,1),1,1,1,1)),bsxfun(@times,B,reshape(ones(1,1),1,1,1,1)),[2,3],[2,3])



size(ans)



mtimesx(permute(T,[2,1,3,4]),permute(Fun.I3,[3,2,1,4,5,6]))
% permute(Fun.I3,)
size(T)
% mtimesx(,)
% length(size(-Fun.I3))
% 
% size(Fun.A_2)
% size(Fun.N_inv)
%Fun.A_2=
% Fun.B
% Fun.C
% Fun.D
% Fun.E
toc

ones(1,2,3,4)+ones(1,2)