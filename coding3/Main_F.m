%Main-Fun_action
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
gamma=1.4;
Geo.kappa=(2/3)./Geo.h;Geo.tau=0.2./Geo.h;
Geo.sw=sqrt(Geo.kappa.^2+Geo.tau.^2).*Geo.s;
Geo.x = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*sin(Geo.sw+0);Geo.y = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*cos(Geo.sw+0);Geo.z = Geo.tau./(Geo.kappa.^2+Geo.tau.^2).*Geo.sw;
Geo.theta_0=cumsum(Geo.tau.*[0 diff(Geo.s)]);
%tubeplot(x,y,z,h,s,50);hold on;plot3(x, y, z);daspect([1,1,1]); camlight;

%% #######Wave########%
Geo.m=-1:1;
Geo.n=1;
a=[1 2]; %P^{a}=P^{a*},U^{a}=U^{a*}
b=[-3 -2 -1 1 2 3];
a_b=a-b.';
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
Psi.ab_r=X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ab');
Psi.a_pr_b_r =X(Base,Geo,2,'a_pr_b','r').*Theta(Geo,2,'ab');
Psi.ab_pt_ab =X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab');
Psi.ab_r2_cos =X(Base,Geo,2,'ab','r2').*Theta(Geo,2,'ab_cos');
Psi.ps_ab_r_1 =X(Base,Geo,2,'ps_ab','r').*Theta(Geo,2,'ab');
Psi.ps_ab_r_2 =X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ps_ab');
Psi.a_ps_b_r_1 =X(Base,Geo,2,'a_ps_b','r').*Theta(Geo,2,'ab');
Psi.a_ps_b_r_2 =X(Base,Geo,2,'ab','r').*Theta(Geo,2,'a_ps_b');
Psi.pr_ab_r=X(Base,Geo,2,'pr_ab','r').*Theta(Geo,2,'ab');
Psi.pr_ab_r2_cos=X(Base,Geo,2,'pr_ab','r2').*Theta(Geo,2,'ab_cos');
Psi.pt_ab=X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab');
Psi.pt_ab_r_cos=X(Base,Geo,2,'ab','r').*Theta(Geo,2,'pt_ab_cos');
Psi.ab_s1= specialFun(Geo.s,Geo.h,Geo.kappa,Geo.m,Geo.n,Base.Cmn1,Base.jmn_pm,'hh`^2/[1-\kappa*h*cos\psi]');
Psi.ab_s2= specialFun(Geo.s,Geo.h,Geo.kappa,Geo.m,Geo.n,Base.Cmn1,Base.jmn_pm,'h(1-\kappa*h*cos\psi)');
Psi.ab=X(Base,Geo,2,'ab','1').*Theta(Geo,2,'ab');
Psi.ab_r_cos=X(Base,Geo,2,'ab','r').*Theta(Geo,2,'ab_cos');
Psi.pt_ab_cos=X(Base,Geo,2,'ab','1').*Theta(Geo,2,'pt_ab_cos');

Psi.abc=X(Base,Geo,3,'ab','1').*Theta(Geo,3,'ab');
Psi.abc_r=X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab');
Psi.abc_r2_cos=X(Base,Geo,3,'ab','r2').*Theta(Geo,3,'ab_cos');
Psi.abc_r_cos= X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab_cos');
Psi.abc_r_sin=X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab_sin');
Psi.pr_abc_r=X(Base,Geo,3,'pr_ab','r').*Theta(Geo,3,'ab');
Psi.pr_abc_r2_cos=X(Base,Geo,3,'pr_ab','r2').*Theta(Geo,3,'ab_cos');
Psi.pt_abc=X(Base,Geo,3,'ab','1').*Theta(Geo,3,'pt_ab');
Psi.pt_abc_r_cos=X(Base,Geo,3,'ab','r').*Theta(Geo,3,'pt_ab_cos');
Psi.abc=X(Base,Geo,3,'ab','1').*Theta(Geo,3,'ab');
Psi.abc_r_cos=X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab_cos');
Psi.abc_r_sin=X(Base,Geo,3,'ab','r').*Theta(Geo,3,'ab_sin');
Psi.ps_abc_ps_r=X(Base,Geo,3,'ps(ab)','r').*Theta(Geo,3,'ab');  %(Jiaqi-153)
Psi.ps_abc_r=X(Base,Geo,3,'ps_ab','r').*Theta(Geo,3,'ab');
Psi.ps_abc=X(Base,Geo,3,'ps_ab','1').*Theta(Geo,3,'ab');

tic
%2D-bsxfun(@times, 3D(\alpha*\beta*s), 1*1*s*a)-->4D[\alpha*\beta*s*a]
Fun2_a.I2=bsxfun(@times,Psi.ab_r,reshape(ones(size(a)),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a.V=bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_b.V=bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*b),1,1,1,length(b)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a_b.V=bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a_b),1,1,1,length(b),length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a.W=-bsxfun(@times,Psi.ab_pt_ab,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a))); %(James-3.27,Jiaqi-77)
Fun2_b.W=-bsxfun(@times,Psi.ab_pt_ab,reshape(1./(sqrt(-1)*k*b),1,1,1,length(b))); %(James-3.27,Jiaqi-77)
Fun2_a_b.W=-bsxfun(@times,Psi.ab_pt_ab,reshape(1./(sqrt(-1)*k*a_b),1,1,1,length(b),length(a))); %(James-3.27,Jiaqi-77)
Fun2_a.N=bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo.kappa.'.*a),1,1,length(Geo.kappa),length(a)));%(James-3.35b)
Fun2_a_b.N=bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a_b),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo.kappa.'.*permute(a_b,[3,1,2])),1,1,length(Geo.kappa),length(b),length(a)));%(James-3.35b)
for i=1:length(Geo.s);for j=1:length(a);Fun2_a.N_inv(:,:,i,j)=Fun2_a.N(:,:,i,j)^(-1);end;end;
Fun2_a.G=-bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_b.G=-bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(b)),1,1,1,length(b)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(b)),1,1,1,length(b)));
Fun2_a_b.G=-bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a_b)),1,1,1,length(b),length(a)));
Fun2_a.H=-bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_b.H=-bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(b)),1,1,1,length(b)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(b)),1,1,1,length(b)));%(James-3.35d)
Fun2_a_b.H=-bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a_b)),1,1,1,length(b),length(a)));%(James-3.35d)
Fun2_a.M_2_1=bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a)));
Fun2_a_b.M_2_1=bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo.kappa.'.*permute(ones(size(a_b)),[3,1,2]),1,1,length(Geo.kappa),length(b),length(a)));
Fun2_a.M_3_1=bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a)));
Fun2_a_b.M_3_1=bsxfun(@times,Psi.pt_ab,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo.kappa.'.*permute(ones(size(a_b)),[3,1,2]),1,1,length(Geo.kappa),length(b),length(a)));


Fun2_a.M  =-Fun2_a.N-multiprod(Fun2_a.M_2_1,Fun2_a.V,[1,2])-multiprod(Fun2_a.M_3_1,Fun2_a.W,[1,2]);%(James-3.35a)
Fun2_a_b.M=-Fun2_a_b.N-multiprod(Fun2_a_b.M_2_1,Fun2_a_b.V,[1,2])-multiprod(Fun2_a_b.M_3_1,Fun2_a_b.W,[1,2]);%(James-3.35a)



Fun2_a.A_2_1_1= bsxfun(@times,Psi.ab_s1,reshape(ones(size(a)),1,1,1,length(a)));
Fun2_a.A_2_1_2= bsxfun(@times,Psi.ab_s2,reshape(ones(size(a)),1,1,1,length(a)));
Fun2_a.A_2_1_3= bsxfun(@times,Psi.ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ab_r_cos,reshape((Geo.kappa.'.*ones(size(a))),1,1,length(Geo.kappa),length(a)));
Fun2_a.A_2_1_4=Fun2_a.M_2_1;
Fun2_a.A_3_1_2=bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pt_ab_cos,reshape(Geo.kappa.'.*ones(size(a)),1,1,length(Geo.kappa),length(a)));


%3D-bsxfun(@times, 4D(\alpha*\beta*\gamma*s), 1*1*1*s*b*a)-->6D[\alpha*\beta*gamma*s*b*a]
Fun3_ab.I3_1=bsxfun(@times,Psi.abc,reshape(bsxfun(@times,ones(size(Geo.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.I3_r=bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));

Fun3_ab.A_1=bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo.s.'))*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)))...
    -bsxfun(@times,Psi.abc_r2_cos,reshape(bsxfun(@times,Geo.kappa.'*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.A_2=  permute(bsxfun(@times,Fun2_a.M_2_1,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.N_inv=permute(bsxfun(@times,Fun2_a.N_inv,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_2_1=permute(bsxfun(@times,Fun2_a.A_2_1_1+Fun2_a.A_2_1_2+Fun2_a.A_2_1_3+Fun2_a.A_2_1_4,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_2_2=bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.A_3=  permute(bsxfun(@times,Fun2_a.M_3_1,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_3_1_2=permute(bsxfun(@times,Fun2_a.A_3_1_2,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_3_2=bsxfun(@times,Psi.abc_r_sin,reshape(bsxfun(@times,Geo.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.A = -Fun3_ab.A_1...
    -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-multiprod(Fun3_ab.I3_r,Fun3_ab.A_2_1,[2,3])-Fun3_ab.A_2_2),[1,2])...
    -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(multiprod(Fun3_ab.I3_r,Fun3_ab.A_3_1_2,[2,3])+Fun3_ab.A_3_2),[1,2]);%(James-3.35e)
Fun3_ab.B_1=bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo.s.'))*((gamma-1)/2*sqrt(-1)*k*ones(size(b))),reshape(a,1,1,[])),1,1,1,length(Geo.s),length(b),length(a)))...
    -bsxfun(@times,Psi.abc_r2_cos,reshape(bsxfun(@times,Geo.kappa.'*((gamma-1)/2*sqrt(-1)*k*ones(size(b))),reshape(a,1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.B_2=Fun3_ab.A_1;
Fun3_a_b.V= permute(Fun2_a_b.V,[1,2,6,3,4,5]);
Fun3_b.V= permute(bsxfun(@times,Fun2_b.V,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_ab.B_3=multiprod(multiprod(Fun3_ab.B_2,Fun3_a_b.V,[1,2]),Fun3_b.V,[2,3]);
Fun3_a_b.W= permute(Fun2_a_b.W,[1,2,6,3,4,5]);
Fun3_b.W= permute(bsxfun(@times,Fun2_b.W,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.G= permute(Fun2_a_b.G,[1,2,6,3,4,5]);
Fun3_b.G= permute(bsxfun(@times,Fun2_b.G,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.H= permute(Fun2_a_b.H,[1,2,6,3,4,5]);
Fun3_b.H= permute(bsxfun(@times,Fun2_b.H,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_ab.B_4=multiprod(multiprod(Fun3_ab.B_2,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
Fun3_a_b.M= permute(Fun2_a_b.M,[1,2,6,3,4,5]);
Fun3_a_b.I=permute(bsxfun(@times,Fun2_a.I2,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_b.I=permute(bsxfun(@times,Fun2_a.I2,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);

Fun3_ab.B_5_1=multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_2=multiprod(multiprod(Fun3_ab.B_1/((gamma-1)/2),Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_31=bsxfun(@times,Psi.pr_abc_r,reshape(bsxfun(@times,ones(size(Geo.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)))...
    -bsxfun(@times,Psi.pr_abc_r2_cos,reshape(bsxfun(@times,Geo.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.B_5_3=multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.V,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_41=bsxfun(@times,Psi.pt_abc,reshape(bsxfun(@times,ones(size(Geo.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)))...
    -bsxfun(@times,Psi.pt_abc_r_cos,reshape(bsxfun(@times,Geo.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.B_5_4=multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.W,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_51=bsxfun(@times,Psi.abc,reshape(bsxfun(@times,ones(size(Geo.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)))...
    -bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.B_5_5=multiprod(multiprod(Fun3_ab.B_5_51,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);

Fun3_ab.B_6_1=multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_2=multiprod(multiprod(Fun3_ab.B_1/((gamma-1)/2),Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_3=multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.V,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_4=multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_5=multiprod(multiprod(Fun3_ab.B_5_51,Fun3_a_b.W,[1,2]),Fun3_b.V,[2,3]);

Fun3_ab.B = -Fun3_ab.B_1-Fun3_ab.B_2-Fun3_ab.B_3-Fun3_ab.B_4...
    -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(Fun3_ab.B_5_1+Fun3_ab.B_5_2+Fun3_ab.B_5_3+Fun3_ab.B_5_4+Fun3_ab.B_5_5),[1,2])...
    -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(Fun3_ab.B_6_1+Fun3_ab.B_6_2+Fun3_ab.B_6_3+Fun3_ab.B_6_4+Fun3_ab.B_6_5),[1,2]);

Fun3_ab.C_1=multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.I,[1,2]),Fun3_b.M,[2,3]);
Fun3_ab.C_2=multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.I,[2,3]);
Fun3_ab.C_3=Fun3_ab.B_1/((gamma-1)/2);
Fun3_ab.C_4=multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.C_5=multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_1=bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.C_6=multiprod(multiprod(Fun3_ab.B_6_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_7_1=bsxfun(@times,Psi.abc_r_sin,reshape(bsxfun(@times,Geo.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo.s),length(b),length(a)));
Fun3_ab.C_7=multiprod(multiprod(Fun3_ab.B_7_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.C=Fun3_ab.C_1+Fun3_ab.C_2+Fun3_ab.C_3+Fun3_ab.C_4+Fun3_ab.C_5+Fun3_ab.C_6-Fun3_ab.C_7; %(James-3.35g)

Fun3_ab.D_1=bsxfun(@times,Psi.ps_abc_ps_r, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
Fun3_ab.D_2=multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.G,[1,2]),Fun3_b.I,[2,3]);
Fun3_ab.D_3=multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.I,[1,2]),Fun3_b.G,[2,3]);
Fun3_ab.D_4=bsxfun(@times,Psi.ps_abc_r, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
Fun3_ab.D=-Fun3_ab.D_1+Fun3_ab.D_2+Fun3_ab.D_3+Fun3_ab.D_4;

Fun3_ab.E_1_1=multiprod(multiprod(Fun3_ab.D_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.E_1_2=multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.G,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.E_1_3=multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.I,[1,2]),multiprod(Fun3_b.G,Fun3_b.V,[2,3]),[2,3]);
Fun3_ab.E_1_4_1=bsxfun(@times,Psi.ps_abc, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
Fun3_ab.E_1_4=multiprod(multiprod(Fun3_ab.E_1_4_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);

Fun3_ab.E_2_1=multiprod(multiprod(Fun3_ab.D_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.E_2_2=multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.G,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.E_2_3=multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.I,[1,2]),multiprod(Fun3_b.H,Fun3_b.W,[2,3]),[2,3]);
Fun3_ab.E_2_4=multiprod(multiprod(Fun3_ab.E_1_4_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);

Fun3_ab.E =-multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-Fun3_ab.E_1_1+Fun3_ab.E_1_2+Fun3_ab.E_1_3+Fun3_ab.E_1_4),[1,2])...
           -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-Fun3_ab.E_2_1+Fun3_ab.E_2_2-Fun3_ab.E_2_3+Fun3_ab.E_2_4),[1,2]);


%Warning:a-b!=0, be NaN, need to be deleted.


toc