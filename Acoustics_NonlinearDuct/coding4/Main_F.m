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
    disp('文件夹存在！');
end
% %% #######Geometry########%
Geo.s =linspace(0,10,50);
%Geo.h=0.1*exp(linspace(0,1.5,length(Geo.s)));
Geo.h=1*ones(size(Geo.s));

Geo.tau=0.0./Geo.h;
%Geo.tau=linspace(0,4,500)
%Geo.tau =logspace(0,0.6,50)-1;
%Geo.tau =0.0*ones(size(Geo.s));
%Geo.tau =logspace(0,0.1,500)-1;

Geo.kappa=(0.000)./Geo.h;
%Geo.kappa=logspace(0,0.6,500)-1

Geo.sw=sqrt(Geo.kappa.^2+Geo.tau.^2).*Geo.s;
Geo.x = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*sin(Geo.sw+0);Geo.y = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*cos(Geo.sw+0);Geo.z = Geo.tau./(Geo.kappa.^2+Geo.tau.^2).*Geo.sw;
Geo.theta_0=cumsum(Geo.tau.*[0 diff(Geo.s)]);
Geo.h_diff=[(Geo.h(2)-Geo.h(1))/(Geo.s(2)-Geo.s(1)) diff(Geo.h)./diff(Geo.s)];
gamma=1.4;
h=figure
tubeplot(Geo.x,Geo.y,Geo.z,Geo.h,Geo.s,50);hold on;plot3(Geo.x, Geo.y, Geo.z);daspect([1,1,1]); camlight;
saveas(h,[save_directory,'\','Geo','.png'])


% %% #######Wave########%
Geo.m=[-3:3];
Geo.n=3;
n_matrix=length(Geo.m)*Geo.n;
% a=[1 2 3]; %P^{a}=P^{a*},U^{a}=U^{a*}
% b=[-3 -2 -1 1 2 3];
% a_b=a-b.';
% k=0.95*1.8412/Geo.h(end)+0.00001*i; %+-;
% 
% %%
% Geo_b.m=Geo.m;
% Geo_b.n=Geo.n;
% Geo_b.s=Geo.s;                  %[Geo.s(1) Geo.s(end)];
% Geo_b.h=Geo.h;                  %[Geo.h(1) Geo.h(end)];
% Geo_b.kappa=Geo.kappa;          %[Geo.kappa(1) Geo.kappa(end)];
% Geo_b.tau=Geo.tau;              %[Geo.tau(1) Geo.tau(end)];
% Geo_b.h_diff=Geo.h_diff         %[Geo.h_diff(1) Geo.h_diff(end)]
% Geo_b.theta_0=Geo.theta_0;
% %% #######Name########%
% Name.op_m=  {'1','r','r2'};
% Name.op=    {'ab','pr_ab','a_pr_b','ps_ab','a_ps_b'};
% Name.X2=    {'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}','X_{\alpha\{\beta\}}'};
% Name.X3=    {'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}','X_{\alpha\{\beta\}\gamma}'};
% Name.T=     {'ab','pt_ab','a_pt_b','ps_ab','a_ps_b','ab_cos','ab_sin','pt_ab_cos','pt_ab_sin','ps_ab_cos','ps_ab_sin'};
% Name.T2=    {'\Theta_{\alpha\beta}','\Theta_{(\alpha)\beta}','\Theta_{\alpha(\beta)}','\Theta_{\{\alpha\}\beta}','\Theta_{\alpha\{\beta\}}','\Theta_{\alpha\beta}[cos\phi]',...
%     '\Theta_{\alpha\beta}[sin\phi]','\Theta_{(\alpha)\beta}[cos\phi]','\Theta_{(\alpha)\beta}[sin\phi]','\Theta_{\{\alpha\}\beta}[cos\phi]','\Theta_{\{\alpha\}\beta}[sin\phi]'};
% Name.T3=    {'\Theta_{\alpha\beta\gamma}','\Theta_{(\alpha)\beta\gamma}','\Theta_{\alpha(\beta)\gamma}','\Theta_{\{\alpha\}\beta\gamma}','\Theta_{\alpha\{\beta\}\gamma}','\Theta_{\alpha\beta\gamma}[cos\phi]',...
%     '\Theta_{\alpha\beta\gamma}[sin\phi]','\Theta_{(\alpha)\beta\gamma}[cos\phi]','\Theta_{(\alpha)\beta\gamma}[sin\phi]','\Theta_{\{\alpha\}\beta\gamma}[cos\phi]','\Theta_{\{\alpha\}\beta\gamma}[sin\phi]'};
% 
% %% #######Function########%
% %% different h can be similarly deal with
% tic
% X2.ab=         X(Geo_b,2,'ab','1');        %1/h
% X2.ab_r=       X(Geo_b,2,'ab','r');        %1
% X2.ab_r2 =     X(Geo_b,2,'ab','r2');       %h
% X2.pr_ab_r=    X(Geo_b,2,'pr_ab','r');     %1/h
% X2.pr_ab_r2=   X(Geo_b,2,'pr_ab','r2');    %1
% X2.a_pr_b_r =  X(Geo_b,2,'a_pr_b','r');    %1
% X2.a_ps_b_r =  X(Geo_b,2,'a_ps_b','r');    %h'/h
% X2.ps_ab_r =   X(Geo_b,2,'ps_ab','r');     %h'/h
% toc;tic
% X3.abc=        X(Geo_b,3,'ab','1');        %1/h^2
% X3.abc_r=      X(Geo_b,3,'ab','r');        %1/h
% X3.abc_r2=     X(Geo_b,3,'ab','r2');       %1
% X3.ps_abc_r=   X(Geo_b,3,'ps_ab','r');     %h'/h^2
% X3.ps_abc=     X(Geo_b,3,'ps_ab','1');     %h'/h^3
% X3.pr_abc_r=   X(Geo_b,3,'pr_ab','r');     %1/h^2
% X3.pr_abc_r2=  X(Geo_b,3,'pr_ab','r2');    %1/h
% X3.ps_abc_ps_r=X(Geo_b,3,'ps(ab)','r');    %h'/h^2
% toc
% %%
% T2.ab=         Theta(Geo_b,2,'ab');
% T2.ps_ab=      Theta(Geo_b,2,'ps_ab');
% T2.pt_ab=      Theta(Geo_b,2,'pt_ab');
% T2.ab_cos=     Theta(Geo_b,2,'ab_cos');
% T2.a_ps_b=     Theta(Geo_b,2,'a_ps_b');
% T2.pt_ab_cos=  Theta(Geo_b,2,'pt_ab_cos');
% 
% T3.abc=        Theta(Geo_b,3,'ab');
% T3.pt_abc=     Theta(Geo_b,3,'pt_ab');
% T3.abc_cos=    Theta(Geo_b,3,'ab_cos');
% T3.abc_sin=    Theta(Geo_b,3,'ab_sin');
% T3.pt_abc_cos= Theta(Geo_b,3,'pt_ab_cos');
% %%
% Psi.ab_r=           X2.ab_r         .*T2.ab;
% Psi.a_pr_b_r =      X2.a_pr_b_r     .*T2.ab;
% Psi.ab_r2_cos =     X2.ab_r2        .*T2.ab_cos;
% Psi.ps_ab_r_1 =     X2.ps_ab_r      .*T2.ab;
% Psi.ps_ab_r_2 =     X2.ab_r         .*T2.ps_ab;
% Psi.a_ps_b_r_1 =    X2.a_ps_b_r     .*T2.ab;
% Psi.a_ps_b_r_2 =    X2.ab_r         .*T2.a_ps_b;
% Psi.pr_ab_r=        X2.pr_ab_r      .*T2.ab;
% Psi.pr_ab_r2_cos=   X2.pr_ab_r2     .*T2.ab_cos;
% Psi.pt_ab=          X2.ab           .*T2.pt_ab;
% Psi.pt_ab_r_cos=    X2.ab_r         .*T2.pt_ab_cos;
% Psi.ab=             X2.ab           .*T2.ab;
% Psi.ab_r_cos=       X2.ab_r         .*T2.ab_cos;
% Psi.pt_ab_cos=      X2.ab           .*T2.pt_ab_cos;
% Psi.ab_s1= specialFun(Geo_b.theta_0,Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'hh`^2/[1-\kappa*h*cos\psi]');
% Psi.ab_s2= specialFun(Geo_b.theta_0,Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h(1-\kappa*h*cos\psi)');
% 
% 
% Psi.abc=                X3.abc          .*T3.abc;
% Psi.abc_r=              X3.abc_r        .*T3.abc;
% Psi.abc_r2_cos=         X3.abc_r2       .*T3.abc_cos;
% Psi.abc_r_cos=          X3.abc_r        .*T3.abc_cos;
% Psi.abc_r_sin=          X3.abc_r        .*T3.abc_sin;
% Psi.pr_abc_r=           X3.pr_abc_r     .*T3.abc;
% Psi.pr_abc_r2_cos=      X3.pr_abc_r2    .*T3.abc_cos;
% Psi.pt_abc=             X3.abc          .*T3.pt_abc;
% Psi.pt_abc_r_cos=       X3.abc_r        .*T3.pt_abc_cos;
% Psi.abc=                X3.abc          .*T3.abc;
% Psi.abc_r_cos=          X3.abc_r        .*T3.abc_cos;
% Psi.abc_r_sin=          X3.abc_r        .*T3.abc_sin;
% Psi.ps_abc_ps_r=        X3.ps_abc_ps_r  .*T3.abc;  %(Jiaqi-153)
% Psi.ps_abc_r=           X3.ps_abc_r     .*T3.abc;
% Psi.ps_abc=             X3.ps_abc       .*T3.abc;
% 
load('Data.mat')
% %
% tic
% %2D-bsxfun(@times, 3D(\alpha*\beta*s), 1*1*s*a)-->4D[\alpha*\beta*s*a]
% Fun2_a.I2=      bsxfun(@times,Psi.ab_r,reshape(ones(size(a)),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
% Fun2_a.V=       bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
% Fun2_b.V=       bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*b),1,1,1,length(b)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
% Fun2_a_b.V=     bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a_b),1,1,1,length(b),length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
% Fun2_a.W=       -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a))); %(James-3.27,Jiaqi-77)
% Fun2_b.W=       -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*b),1,1,1,length(b))); %(James-3.27,Jiaqi-77)
% Fun2_a_b.W=     -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*a_b),1,1,1,length(b),length(a))); %(James-3.27,Jiaqi-77)
% Fun2_a.N=       bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*a),1,1,length(Geo_b.kappa),length(a)));%(James-3.35b)
% Fun2_a_b.N=     bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a_b),1,1,1,length(b),length(a)))...
%     -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*permute(a_b,[3,1,2])),1,1,length(Geo_b.kappa),length(b),length(a)));%(James-3.35b)
% Fun2_b.N=       bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*b),1,1,1,length(b)))...
%     -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*b),1,1,length(Geo_b.kappa),length(b )));%(James-3.35b)
% for i=1:length(Geo_b.s);for j=1:length(a);Fun2_a.N_inv(:,:,i,j)=Fun2_a.N(:,:,i,j)^(-1);end;end;
% for i=1:length(Geo_b.s);for j=1:length(b); for k=1:length(a);if a(k)-b(j)==0;Fun2_a_b.N_inv(:,:,i,j,k)=NaN+NaN*sqrt(-1);else;Fun2_a_b.N_inv(:,:,i,j,k)=Fun2_a_b.N(:,:,i,j,k)^(-1);end;end;end;end;
% for i=1:length(Geo_b.s);for j=1:length(b);Fun2_b.N_inv(:,:,i,j)=Fun2_b.N(:,:,i,j)^(-1);end;end;
% 
% Fun2_a.G=       -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
% Fun2_b.G=       -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(b)),1,1,1,length(b)))...
%     -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(b)),1,1,1,length(b)));
% Fun2_a_b.G=     -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
%     -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a_b)),1,1,1,length(b),length(a)));
% Fun2_a.H=       -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
% Fun2_b.H=       -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(b)),1,1,1,length(b)))...
%     -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(b)),1,1,1,length(b)));%(James-3.35d)
% Fun2_a_b.H=     -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
%     -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a_b)),1,1,1,length(b),length(a)));%(James-3.35d)
% Fun2_a.M_2_1=   bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a)),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));
% Fun2_a_b.M_2_1= bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
%     -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*permute(ones(size(a_b)),[3,1,2]),1,1,length(Geo_b.kappa),length(b),length(a)));
% Fun2_b.M_2_1=   bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(b)),1,1,1,length(b)))...
%     -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*ones(size(b)),1,1,length(Geo_b.kappa),length(b)));
% 
% Fun2_a.M_3_1=   bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));
% Fun2_a_b.M_3_1= bsxfun(@times,Psi.pt_ab,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
%     -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*permute(ones(size(a_b)),[3,1,2]),1,1,length(Geo_b.kappa),length(b),length(a)));
% Fun2_b.M_3_1=   bsxfun(@times,Psi.pt_ab,reshape(ones(size(b)),1,1,1,length(b)))...
%     -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*ones(size(b)),1,1,length(Geo_b.kappa),length(b)));
% 
% Fun2_a.M  =     -Fun2_a.N-multiprod(Fun2_a.M_2_1,Fun2_a.V,[1,2])-multiprod(Fun2_a.M_3_1,Fun2_a.W,[1,2]);%(James-3.35a)
% Fun2_a_b.M=     -Fun2_a_b.N-multiprod(Fun2_a_b.M_2_1,Fun2_a_b.V,[1,2])-multiprod(Fun2_a_b.M_3_1,Fun2_a_b.W,[1,2]);%(James-3.35a)
% Fun2_b.M=       -Fun2_b.N-multiprod(Fun2_b.M_2_1,Fun2_b.V,[1,2])-multiprod(Fun2_b.M_3_1,Fun2_b.W,[1,2]);%(James-3.35a)
% 
% Fun2_a.A_2_1_1= bsxfun(@times,Psi.ab_s1,reshape(ones(size(a)),1,1,1,length(a)));
% Fun2_a.A_2_1_2= bsxfun(@times,Psi.ab_s2,reshape(ones(size(a)),1,1,1,length(a)));
% Fun2_a.A_2_1_3= bsxfun(@times,Psi.ab,reshape(ones(size(a)),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.ab_r_cos,reshape((Geo_b.kappa.'.*ones(size(a))),1,1,length(Geo_b.kappa),length(a)));
% Fun2_a.A_2_1_4= Fun2_a.M_2_1;
% Fun2_a.A_3_1_2= bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
%     -bsxfun(@times,Psi.pt_ab_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));
% 
% 
% %3D-bsxfun(@times, 4D(\alpha*\beta*\gamma*s), 1*1*1*s*b*a)-->6D[\alpha*\beta*gamma*s*b*a]
% Fun3_ab.I3_1=   bsxfun(@times,Psi.abc,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.I3_r=   bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% 
% Fun3_ab.A_1=    bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
%     -bsxfun(@times,Psi.abc_r2_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.A_2=    permute(bsxfun(@times,Fun2_a.M_2_1,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
% Fun3_ab.N_inv=  permute(bsxfun(@times,Fun2_a.N_inv,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
% Fun3_ab.A_2_1=  permute(bsxfun(@times,Fun2_a.A_2_1_1+Fun2_a.A_2_1_2+Fun2_a.A_2_1_3+Fun2_a.A_2_1_4,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
% Fun3_ab.A_2_2=  bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.A_3=    permute(bsxfun(@times,Fun2_a.M_3_1,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
% Fun3_ab.A_3_1_2=permute(bsxfun(@times,Fun2_a.A_3_1_2,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
% Fun3_ab.A_3_2=  bsxfun(@times,Psi.abc_r_sin,reshape(bsxfun(@times,Geo_b.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.A =     -Fun3_ab.A_1...
%     -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-multiprod(Fun3_ab.I3_r,Fun3_ab.A_2_1,[2,3])-Fun3_ab.A_2_2),[1,2])...
%     -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(multiprod(Fun3_ab.I3_r,Fun3_ab.A_3_1_2,[2,3])+Fun3_ab.A_3_2),[1,2]);%(James-3.35e)
% Fun3_ab.B_1=    bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*((gamma-1)/2*sqrt(-1)*k*ones(size(b))),reshape(a,1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
%     -bsxfun(@times,Psi.abc_r2_cos,reshape(bsxfun(@times,Geo_b.kappa.'*((gamma-1)/2*sqrt(-1)*k*ones(size(b))),reshape(a,1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.B_2=    Fun3_ab.A_1;
% Fun3_a_b.V=     permute(Fun2_a_b.V,[1,2,6,3,4,5]);
% Fun3_b.V=       permute(bsxfun(@times,Fun2_b.V,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
% Fun3_ab.B_3=    multiprod(multiprod(Fun3_ab.B_2,Fun3_a_b.V,[1,2]),Fun3_b.V,[2,3]);
% Fun3_a_b.W=     permute(Fun2_a_b.W,[1,2,6,3,4,5]);
% Fun3_b.W=       permute(bsxfun(@times,Fun2_b.W,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
% Fun3_a_b.G=     permute(Fun2_a_b.G,[1,2,6,3,4,5]);
% Fun3_b.G=       permute(bsxfun(@times,Fun2_b.G,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
% Fun3_a_b.H=     permute(Fun2_a_b.H,[1,2,6,3,4,5]);
% Fun3_b.H=       permute(bsxfun(@times,Fun2_b.H,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
% Fun3_ab.B_4=    multiprod(multiprod(Fun3_ab.B_2,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
% Fun3_a_b.M=     permute(Fun2_a_b.M,[1,2,6,3,4,5]);
% Fun3_b.M=       permute(bsxfun(@times,Fun2_b.M,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
% 
% Fun3_a_b.I=     permute(bsxfun(@times,Fun2_a.I2,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
% Fun3_b.I=       permute(bsxfun(@times,Fun2_a.I2,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
% 
% Fun3_ab.B_5_1=  multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.B_5_2=  multiprod(multiprod(Fun3_ab.B_1/((gamma-1)/2),Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.B_5_31= bsxfun(@times,Psi.pr_abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
%     -bsxfun(@times,Psi.pr_abc_r2_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.B_5_3=  multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.V,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.B_5_41= bsxfun(@times,Psi.pt_abc,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
%     -bsxfun(@times,Psi.pt_abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.B_5_4=  multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.W,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.B_5_51= bsxfun(@times,Psi.abc,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
%     -bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.B_5_5=  multiprod(multiprod(Fun3_ab.B_5_51,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
% 
% Fun3_ab.B_6_1=  multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.B_6_2=  multiprod(multiprod(Fun3_ab.B_1/((gamma-1)/2),Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.B_6_3=  multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.V,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.B_6_4=  multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.B_6_5=  multiprod(multiprod(Fun3_ab.B_5_51,Fun3_a_b.W,[1,2]),Fun3_b.V,[2,3]);
% 
% Fun3_ab.B =     -Fun3_ab.B_1-Fun3_ab.B_2-Fun3_ab.B_3-Fun3_ab.B_4...
%     -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(Fun3_ab.B_5_1+Fun3_ab.B_5_2+Fun3_ab.B_5_3+Fun3_ab.B_5_4+Fun3_ab.B_5_5),[1,2])...
%     -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(Fun3_ab.B_6_1+Fun3_ab.B_6_2+Fun3_ab.B_6_3+Fun3_ab.B_6_4+Fun3_ab.B_6_5),[1,2]);
% 
% Fun3_ab.C_1=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.I,[1,2]),Fun3_b.M,[2,3]);
% Fun3_ab.C_2=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.I,[2,3]);
% Fun3_ab.C_3=    Fun3_ab.B_1/((gamma-1)/2);
% Fun3_ab.C_4=    multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.C_5=    multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.B_6_1=  bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.C_6=    multiprod(multiprod(Fun3_ab.B_6_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.B_7_1=  bsxfun(@times,Psi.abc_r_sin,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
% Fun3_ab.C_7=    multiprod(multiprod(Fun3_ab.B_7_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.C=      Fun3_ab.C_1+Fun3_ab.C_2+Fun3_ab.C_3+Fun3_ab.C_4+Fun3_ab.C_5+Fun3_ab.C_6-Fun3_ab.C_7; %(James-3.35g)
% 
% Fun3_ab.D_1=    bsxfun(@times,Psi.ps_abc_ps_r, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
% Fun3_ab.D_2=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.G,[1,2]),Fun3_b.I,[2,3]);
% Fun3_ab.D_3=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.I,[1,2]),Fun3_b.G,[2,3]);
% Fun3_ab.D_4=    bsxfun(@times,Psi.ps_abc_r, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
% Fun3_ab.D=      -Fun3_ab.D_1+Fun3_ab.D_2+Fun3_ab.D_3+Fun3_ab.D_4;
% 
% Fun3_ab.E_1_1=  multiprod(multiprod(Fun3_ab.D_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.E_1_2=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.G,[1,2]),Fun3_b.V,[2,3]);
% Fun3_ab.E_1_3=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.I,[1,2]),multiprod(Fun3_b.G,Fun3_b.V,[2,3]),[2,3]);
% Fun3_ab.E_1_4_1=bsxfun(@times,Psi.ps_abc, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
% Fun3_ab.E_1_4=  multiprod(multiprod(Fun3_ab.E_1_4_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
% 
% Fun3_ab.E_2_1=  multiprod(multiprod(Fun3_ab.D_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.E_2_2=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.G,[1,2]),Fun3_b.W,[2,3]);
% Fun3_ab.E_2_3=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.I,[1,2]),multiprod(Fun3_b.H,Fun3_b.W,[2,3]),[2,3]);
% Fun3_ab.E_2_4=  multiprod(multiprod(Fun3_ab.E_1_4_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
% 
% Fun3_ab.E =     -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-Fun3_ab.E_1_1+Fun3_ab.E_1_2+Fun3_ab.E_1_3+Fun3_ab.E_1_4),[1,2])...
%     -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-Fun3_ab.E_2_1+Fun3_ab.E_2_2-Fun3_ab.E_2_3+Fun3_ab.E_2_4),[1,2]);


%Warning:a-b!=0, be NaN, need to be deleted.
%toc

%% #######Boudnary########%
Fun2_a.NM=              multiprod(Fun2_a.N,Fun2_a.M,[1,2]);
Fun2_a_b.NM=            multiprod(Fun2_a_b.N,Fun2_a_b.M,[1,2]);
Fun2_b.NM=              multiprod(Fun2_b.N,Fun2_b.M,[1,2]);
Fun2_a.L=               [-Fun2_a.G  -Fun2_a.M;...
                          Fun2_a.N   Fun2_a.H];
Fun2_b.L=               [-Fun2_b.G  -Fun2_b.M;...
                          Fun2_b.N   Fun2_b.H];
Fun2_a_b.L=               [-Fun2_a_b.G  -Fun2_a_b.M;...
                            Fun2_a_b.N   Fun2_a_b.H];
Fun2_a.L2=              multiprod(Fun2_a.L,Fun2_a.L,[1,2]);

%Case1: Torsion Free Outlet
%a

for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        Fun2_a.Y(:,:,kh,ka) =                                 sqrt(-1)*Fun2_a.N_inv(:,:,kh,ka)*sqrtm(Fun2_a.NM(:,:,kh,ka));
        [Fun2_a.Vb(:,:,kh,ka),Fun2_a.Lambda(:,:,kh,ka)]=      eigs(sqrt(-1)*sqrtm(Fun2_a.NM(:,:,kh,ka)),n_matrix);
        Fun2_a.YN(:,:,kh,ka)=                                 Fun2_a.Y(:,:,kh,ka)*Fun2_a.N(:,:,kh,ka);
        [Fun2_a.Wb(:,:,kh,ka),Fun2_a.Lambda1(:,:,kh,ka)]=     eigs(Fun2_a.YN(:,:,kh,ka),n_matrix);        
        Fun2_a.Vb_inv(:,:,kh,ka)=                             inv(Fun2_a.Vb(:,:,kh,ka));
        Fun2_a.Wb_inv(:,:,kh,ka)=                             inv(Fun2_a.Wb(:,:,kh,ka));
    end
end
%b
for kh=1:length(Geo_b.h)
    for kb=1:length(b)
        Fun2_b.Y(:,:,kh,kb) = sqrt(-1)*Fun2_b.N_inv(:,:,kh,kb)*sqrtm(Fun2_b.NM(:,:,kh,kb));
        [Fun2_b.Vb(:,:,kh,kb),Fun2_b.Lambda(:,:,kh,kb)]=eigs(sqrt(-1)*sqrtm(Fun2_b.NM(:,:,kh,kb)),n_matrix);
        Fun2_b.YN(:,:,kh,kb)=Fun2_b.Y(:,:,kh,kb)*Fun2_b.N(:,:,kh,kb);
        [Fun2_b.Wb(:,:,kh,kb),temp]=eigs(Fun2_b.YN(:,:,kh,kb),n_matrix);
        Fun2_b.Vb_inv(:,:,kh,kb)=inv(Fun2_b.Vb(:,:,kh,kb));
        Fun2_b.Wb_inv(:,:,kh,kb)=inv(Fun2_b.Wb(:,:,kh,kb));
    end
end
%a_b
temp1=(NaN+NaN*sqrt(-1))*ones(size(Fun2_a_b.N(:,:,1,1,1)));
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        for kb=1:length(b)
            if a(ka)-b(kb)==0
                Fun2_a_b.Y(:,:,kh,kb,ka)=         temp1;
                Fun2_a_b.Vb(:,:,kh,kb,ka)=        temp1;
                Fun2_a_b.Lambda(:,:,kh,kb,ka)=    temp1;
                Fun2_a_b.Wb(:,:,kh,kb,ka)=        temp1;
                Fun2_a_b.Vb_inv(:,:,kh,kb,ka)=    temp1;
                Fun2_a_b.Wb_inv(:,:,kh,kb,ka)=    temp1;
            else
                Fun2_a_b.Y(:,:,kh,kb,ka) = sqrt(-1)*Fun2_a_b.N_inv(:,:,kh,kb,ka)*sqrtm(Fun2_a_b.NM(:,:,kh,kb,ka));
                [Fun2_a_b.Vb(:,:,kh,kb,ka),Fun2_a_b.Lambda(:,:,kh,kb,ka)]=eigs(sqrt(-1)*sqrtm(Fun2_a_b.NM(:,:,kh,kb,ka)),n_matrix);
                Fun2_a_b.YN(:,:,kh,kb,ka)=Fun2_a_b.Y(:,:,kh,kb,ka)*Fun2_a_b.N(:,:,kh,kb,ka);
                [Fun2_a_b.Wb(:,:,kh,kb,ka),temp]=eigs(Fun2_a_b.YN(:,:,kh,kb,ka),n_matrix);
                Fun2_a_b.Vb_inv(:,:,kh,kb,ka)=inv(Fun2_a_b.Vb(:,:,kh,kb,ka));
                Fun2_a_b.Wb_inv(:,:,kh,kb,ka)=inv(Fun2_a_b.Wb(:,:,kh,kb,ka));
            end
        end
    end
end
%2D->3D
Fun2_a.Y_minus=        -Fun2_a.Y;
Fun3_a_b.Vb=           permute(Fun2_a_b.Vb,[1,2,6,3,4,5]);
Fun3_b.Vb=             permute(bsxfun(@times,Fun2_b.Vb,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.Vb_inv=       permute(Fun2_a_b.Vb_inv,[1,2,6,3,4,5]);
Fun3_b.Vb_inv=         permute(bsxfun(@times,Fun2_b.Vb_inv,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.YV=           permute(multiprod(Fun2_a_b.Y,Fun2_a_b.Vb,[1,2]),[1,2,6,3,4,5]);
Fun3_b.YV=             permute(bsxfun(@times,multiprod(Fun2_b.Y,Fun2_b.Vb,[1,2]),reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a.Wb=             permute(bsxfun(@times,Fun2_a.Wb,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Wb_inv=         permute(bsxfun(@times,Fun2_a.Wb_inv,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Y=              permute(bsxfun(@times,Fun2_a.Y,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_ab.II=            ones(size(Fun3_ab.A));
Fun3_a.Lambda=         permute(bsxfun(@times,Fun2_a.Lambda,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a_b.Lambda=       permute(Fun2_a_b.Lambda,[1,2,6,3,4,5]);
Fun3_b.Lambda=         permute(bsxfun(@times,Fun2_b.Lambda,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3.Lambda_alpha_a=   multiprod(Fun3_a.Lambda,multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.I,[2,3]),[1,2]);
Fun3.Lambda_beta_a_b=  multiprod(multiprod(Fun3_ab.II,Fun3_a_b.Lambda,[1,2]),Fun3_b.I,[2,3]);
Fun3.Lambda_gamma_b=   multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.Lambda,[2,3]);
Fun3_ab.Yp=            (multiprod(Fun3_a.Wb_inv,multiprod(multiprod(Fun3_ab.A,Fun3_a_b.YV,[1,2]),Fun3_b.YV,[2,3]),[1,2])...
                       +multiprod(Fun3_a.Wb_inv,multiprod(multiprod(Fun3_ab.B,Fun3_a_b.Vb,[1,2]),Fun3_b.Vb,[2,3]),[1,2])...
                       -multiprod(multiprod(Fun3_a.Wb_inv,Fun3_a.Y,[1,2]),multiprod(multiprod(Fun3_ab.C,Fun3_a_b.YV,[1,2]),Fun3_b.Vb,[2,3]),[1,2]))...
                        ./(Fun3.Lambda_alpha_a+Fun3.Lambda_beta_a_b+Fun3.Lambda_gamma_b);
Fun3_ab.YY=            multiprod(Fun3_a.Wb,multiprod(multiprod(Fun3_ab.Yp,Fun3_a_b.Vb_inv,[1,2]),Fun3_b.Vb_inv,[2,3]),[1,2]);
Fun3_ab.YY_minus=      -Fun3_ab.YY;

%Warning:a-b!=0, be NaN, need to be deleted.



%Case2: Torsion Helical Duct
%a
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        [Fun2_a.Vb_H(:,:,kh,ka),Fun2_a.Lambda_H(:,:,kh,ka)]= eigs(( Fun2_a.L(:,:,kh,ka)),n_matrix*2); %sqrt(-1)*
        eigenvalues_a=diag(Fun2_a.Lambda_H(:,:,kh,ka));
        order_a=find(real(eigenvalues_a)<0);
        %Fun2_a.Vb_HH(:,:,kh,ka)=Fun2_a.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order_a,kh,ka);
        Fun2_a.Vb_HH(:,:,kh,ka)=Fun2_a.Vb(:,:,kh,ka);%test!
        Fun2_a.Vb_inv_HH(:,:,kh,ka)=inv(Fun2_a.Vb_HH(:,:,kh,ka));
        %Fun2_a.Lambda_HH(:,:,kh,ka)=diag(eigenvalues_a(order_a));
        Fun2_a.Lambda_HH(:,:,kh,ka)=Fun2_a.Lambda(:,:,kh,ka);
    end
end
%b
for kh=1:length(Geo_b.h)
    for kb=1:length(b)
        [Fun2_b.Vb_H(:,:,kh,kb),Fun2_b.Lambda_H(:,:,kh,kb)]= eigs(( Fun2_b.L(:,:,kh,kb)),n_matrix*2); %sqrt(-1)*
        eigenvalues_b=diag(Fun2_b.Lambda_H(:,:,kh,kb));
        order_b=find(real(eigenvalues_b)<0);
        %Fun2_b.Vb_HH(:,:,kh,kb)=Fun2_b.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order_b,kh,kb);
        Fun2_b.Vb_HH(:,:,kh,kb)=Fun2_b.Vb(:,:,kh,kb);%test!
        Fun2_b.Vb_inv_HH(:,:,kh,kb)=inv(Fun2_b.Vb_HH(:,:,kh,kb));
        %Fun2_b.Lambda_HH(:,:,kh,kb)=diag(eigenvalues_b(order_b)); %test!
        Fun2_b.Lambda_HH(:,:,kh,kb)=Fun2_b.Lambda(:,:,kh,kb);
    end
end
%a_b
temp1=(NaN+NaN*sqrt(-1))*ones(size(Fun2_a_b.N(:,:,1,1,1)));
temp2=(NaN+NaN*sqrt(-1))*ones(size(Fun2_a_b.L(:,:,1,1,1)));

for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        for kb=1:length(b)
            if a(ka)-b(kb)==0
                Fun2_a_b.Lambda_H(:,:,kh,kb,ka)=  temp2;
                Fun2_a_b.Vb_H(:,:,kh,kb,ka)=      temp2;
                Fun2_a_b.Vb_HH(:,:,kh,kb,ka)=     temp1;
                Fun2_a_b.Vb_inv_HH(:,:,kh,kb,ka)= temp1;
                Fun2_a_b.Vb_inv(:,:,kh,kb,ka)=    temp1;
                Fun2_a_b.Lambda_HH(:,:,kh,kb,ka)= temp1;
            else
                [Fun2_a_b.Vb_H(:,:,kh,kb,ka),Fun2_a_b.Lambda_H(:,:,kh,kb,ka)]= eigs((Fun2_a_b.L(:,:,kh,kb,ka)),n_matrix*2); %sqrt(-1)*
                eigenvalues_a_b=diag(Fun2_a_b.Lambda_H(:,:,kh,kb,ka));
                order_a_b=find(real(eigenvalues_a_b)<0);
                %Fun2_a_b.Vb_HH(:,:,kh,kb,ka)=Fun2_a_b.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order_a_b,kh,kb,ka);
                Fun2_a_b.Vb_HH(:,:,kh,kb,ka)=Fun2_a_b.Vb(:,:,kh,kb,ka);%test!
                Fun2_a_b.Vb_inv_HH(:,:,kh,kb,ka)=inv(Fun2_a_b.Vb_HH(:,:,kh,kb,ka));
                %Fun2_a_b.Lambda_HH(:,:,kh,kb,ka)=diag(eigenvalues_a_b(order_a_b));
                 Fun2_a_b.Lambda_HH(:,:,kh,kb,ka)=Fun2_a_b.Lambda(:,:,kh,kb,ka);

            end
        end
    end
end
Fun2_a.Y_H=    multiprod(Fun2_a.N_inv,multiprod(multiprod(Fun2_a.Vb_HH,Fun2_a.Lambda_HH,[1,2]),Fun2_a.Vb_inv_HH,[1,2])-Fun2_a.H,[1,2]);
Fun2_b.Y_H=    multiprod(Fun2_b.N_inv,multiprod(multiprod(Fun2_b.Vb_HH,Fun2_b.Lambda_HH,[1,2]),Fun2_b.Vb_inv_HH,[1,2])-Fun2_b.H,[1,2]);
Fun2_a_b.Y_H=  multiprod(Fun2_a_b.N_inv,multiprod(multiprod(Fun2_a_b.Vb_HH,Fun2_a_b.Lambda_HH,[1,2]),Fun2_a_b.Vb_inv_HH,[1,2])-Fun2_a_b.H,[1,2]);
max(max(Fun2_a.Vb_HH(:,:,2,1)-Fun2_a_b.Vb_HH(:,:,2,6,2)))

max(max(Fun2_a.NM(:,:,2,1)-Fun2_a_b.NM(:,:,2,6,2)))


max(max(Fun2_a.Lambda(:,:,2,3)-Fun2_b.Lambda(:,:,2,6)))



Fun2_a.YNG=multiprod(Fun2_a.Y_H,Fun2_a.N,[1,2])+Fun2_a.G;     %James-3.52
Fun2_b.YNG=multiprod(Fun2_b.Y_H,Fun2_b.N,[1,2])+Fun2_b.G;
Fun2_a_b.YNG=multiprod(Fun2_a_b.Y_H,Fun2_a_b.N,[1,2])+Fun2_a_b.G;
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        [Fun2_a.Wb_H(:,:,kh,ka),Fun2_a.Lambda1_H(:,:,kh,ka)]=eigs(Fun2_a.YN(:,:,kh,ka),n_matrix);%!Fun2_a.YNG
         Fun2_a.Wb_inv_H(:,:,kh,ka)=inv(Fun2_a.Wb_H(:,:,kh,ka));
    end
end
for kh=1:length(Geo_b.h)
    for kb=1:length(b)
        [Fun2_b.Wb_H(:,:,kh,kb),Fun2_b.Lambda1_H(:,:,kh,kb)]=eigs(Fun2_b.YN(:,:,kh,kb),n_matrix);%!Fun2_a.YNG
         Fun2_b.Wb_inv_H(:,:,kh,kb)=inv(Fun2_b.Wb_H(:,:,kh,kb));
    end
end

for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        for kb=1:length(b)
            if a(ka)-b(kb)==0
                Fun2_a_b.Wb_H(:,:,kh,kb,ka)=temp1;
                Fun2_a_b.Lambda1_H(:,:,kh,kb,ka)=temp1;               
            else
                [Fun2_a_b.Wb_H(:,:,kh,kb,ka),Fun2_a_b.Lambda1_H(:,:,kh,kb,ka)]=eigs(Fun2_a_b.YN(:,:,kh,kb,ka),n_matrix);%!Fun2_a.YNG
                 Fun2_a_b.Wb_inv_H(:,:,kh,kb,ka)=inv(Fun2_a_b.Wb_H(:,:,kh,kb,ka));
            end
        end
    end
end

%2D->3D
Fun3_a_b.Vb_H=            permute(Fun2_a_b.Vb_HH,[1,2,6,3,4,5]);

Fun3_b.Vb_H=              permute(bsxfun(@times,Fun2_b.Vb_HH,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.Vb_inv_H=        permute(Fun2_a_b.Vb_inv_HH,[1,2,6,3,4,5]);
Fun3_b.Vb_inv_H=          permute(bsxfun(@times,Fun2_b.Vb_inv_HH,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);

Fun3_a_b.YV_H=            permute(multiprod(Fun2_a_b.Y_H,Fun2_a_b.Vb_HH,[1,2]),[1,2,6,3,4,5]);%!!!!!!!!!均为0，bug
Fun3_b.YV_H=              permute(bsxfun(@times,multiprod(Fun2_b.Y_H,Fun2_b.Vb_HH,[1,2]),reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a.Wb_H=              permute(bsxfun(@times,Fun2_a.Wb_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Wb_inv_H=          permute(bsxfun(@times,Fun2_a.Wb_inv_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Y_H=               permute(bsxfun(@times,Fun2_a.Y_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Lambda_HH=         permute(bsxfun(@times,Fun2_a.Lambda_HH,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a_b.Lambda_HH=       permute(Fun2_a_b.Lambda_HH,[1,2,6,3,4,5]);
Fun3_b.Lambda_HH=         permute(bsxfun(@times,Fun2_b.Lambda_HH,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a.Lambda1_H=         permute(bsxfun(@times,Fun2_a.Lambda1_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);

Fun3.Lambda_H_alpha_a=    multiprod(Fun3_a.Lambda1_H,multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.I,[2,3]),[1,2]);
Fun3.Lambda_H_beta_a_b=   multiprod(multiprod(Fun3_ab.II,Fun3_a_b.Lambda_HH,[1,2]),Fun3_b.I,[2,3]);
Fun3.Lambda_H_gamma_b=    multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.Lambda_HH,[2,3]);


Fun3_ab.Yp_H=            (multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.A,Fun3_a_b.YV_H,[1,2]),Fun3_b.YV_H,[2,3]),[1,2])...
                         +multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.B,Fun3_a_b.Vb_H,[1,2]),Fun3_b.Vb_H,[2,3]),[1,2])...
                         -multiprod(multiprod(Fun3_a.Wb_inv_H,Fun3_a.Y_H,[1,2]),multiprod(multiprod(Fun3_ab.C,Fun3_a_b.YV_H,[1,2]),Fun3_b.Vb_H,[2,3]),[1,2])...
                         -multiprod(multiprod(Fun3_a.Wb_inv_H,Fun3_a.Y_H,[1,2]),multiprod(multiprod(Fun3_ab.D,Fun3_a_b.YV_H,[1,2]),Fun3_b.YV_H,[2,3]),[1,2])...
                         +multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.E,Fun3_a_b.YV_H,[1,2]),Fun3_b.Vb_H,[2,3]),[1,2]))...
                         ./(Fun3.Lambda_H_alpha_a+Fun3.Lambda_H_beta_a_b+Fun3.Lambda_H_gamma_b);
Fun3_ab.YY_H=            multiprod(Fun3_a.Wb_H,multiprod(multiprod(Fun3_ab.Yp_H,Fun3_a_b.Vb_inv_H,[1,2]),Fun3_b.Vb_inv_H,[2,3]),[1,2]);

%Warning:a-b!=0, be NaN, need to be deleted.

%% Vertification
%Case: for staright duct, should be 0
max(max(max(max(Fun2_a.Y_H-Fun2_a.Y))))
max(max(max(max(Fun2_b.Y_H-Fun2_b.Y))))
max(max(max(max(max(Fun2_a_b.Y_H-Fun2_a_b.Y)))))
Fun3_a_b.YV_H-Fun3_a_b.YV_H

t1=multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.A,Fun3_a_b.YV_H,[1,2]),Fun3_b.YV_H,[2,3]),[1,2]);
t2=multiprod(Fun3_a.Wb_inv,multiprod(multiprod(Fun3_ab.A,Fun3_a_b.YV_H,[1,2]),Fun3_b.YV,[2,3]),[1,2]);
for kh=2
for k=1:3
    for k1=1:6
max(max(max(max(max(max(Fun2_a_b.Y_H(:,:,kh,k1,k)-Fun2_a_b.Y(:,:,kh,k1,k)))))))
end
end
end


%问题出在Fun3_a_b.YV_H，Fun3_a_b.YV
Fun3_a_b.YV_H=            permute(multiprod(Fun2_a_b.Y_H,Fun2_a_b.Vb_HH,[1,2]),[1,2,6,3,4,5]);%!!!!!!!!!均为0，bug

Fun3_a_b.YV=              permute(multiprod(Fun2_a_b.Y,Fun2_a_b.Vb,[1,2]),[1,2,6,3,4,5]);


t11=multiprod(Fun2_a_b.Y_H,Fun2_a_b.Vb_HH,[1,2]);
t12=multiprod(Fun2_a_b.Y,  Fun2_a_b.Vb,[1,2]);




%Fun2_a_b.Y_H,Fun2_a_b.Y
Fun2_a_b.NM-Fun2_a_b.NM_H



Fun2_a_b.Y_H-Fun2_a_b.Y

Fun2_a_b.Vb_HH-Fun2_a_b.Vb

size(Fun3_a.Wb_inv_H)
size(Fun3_a.Wb_inv)

max(max(Fun3_a_b.YV_H-Fun3_a_b.YV))




max(max(max(max(max(max(max(Fun3_ab.YY_H)))))))
max(max(max(max(max(max(max(Fun3_ab.YY)))))))   %bug!!

max(max(max(max(max(Fun2_a.G)))))
max(max(max(max(Fun2_a.Lambda1_H-Fun2_a.Lambda1))))

max(max(max(max(Fun2_a.YN-Fun2_a.YNG))))

max(max(max(max(Fun2_a.YN))))
max(max(max(max(Fun2_a.YNG))))

max(max(max(max(Fun2_a.Vb_HH-Fun2_a.Vb))))

max(max(max(max(Fun2_a.Vb_HH))))
max(max(max(max(Fun2_a.Vb))))

max(max(max(max(max(max(Fun3.Lambda_H_alpha_a-Fun3.Lambda_alpha_a))))))
max(max(max(max(max(max(Fun3.Lambda_beta_a_b-Fun3.Lambda_H_beta_a_b))))))
max(max(max(max(max(max(Fun3.Lambda_gamma_b-Fun3.Lambda_H_gamma_b))))))


Fun3.Lambda_H_beta_a_b=   multiprod(multiprod(Fun3_ab.II,Fun3_a_b.Lambda_HH,[1,2]),Fun3_b.I,[2,3]);
Fun3.Lambda_beta_a_b=     multiprod(multiprod(Fun3_ab.II,Fun3_a_b.Lambda,[1,2]),Fun3_b.I,[2,3]);

max(max(max(max(max(max(Fun3_a_b.Lambda_HH-Fun3_a_b.Lambda))))))


max(max(max(max(max(Fun2_a.Wb_inv_H-Fun2_a.Wb_inv)))))

ts_max(Fun2_a.Lambda1_H,Fun2_a.Lambda)



Fun3_ab.E
Fun3_ab.D

%% RK4,a,ab,b

 
 


%s只计算一半
or=1:2:50;ds=Geo.s(2)-Geo.s(1);
%a
Y_a_0=Fun2_a.Y_H(:,:,end,:);
RK.Y_a(:,:,1,:)=Y_a_0;
for k=1:length(or)-1
k_1= - multiprod(multiprod(Y_a_0,Fun2_a.N(:,:,end-or(k)+1,:),[1,2]),Y_a_0,[1,2])- Fun2_a.M(:,:,end-or(k)+1,:) - multiprod(Y_a_0,Fun2_a.H(:,:,end-or(k)+1,:),[1,2]) - multiprod(Fun2_a.G(:,:,end-or(k)+1,:),Y_a_0,[1,2]);
Y_1= Y_a_0-k_1*ds;
max(max(max(k_1)))
k_2= - multiprod(multiprod(Y_1,Fun2_a.N(:,:,end-or(k),:),[1,2]),Y_1,[1,2])-Fun2_a.M(:,:,end-or(k),:) - multiprod(Y_1,Fun2_a.H(:,:,end-or(k),:),[1,2]) - multiprod(Fun2_a.G(:,:,end-or(k),:),Y_1,[1,2]);
Y_2= Y_a_0-k_2*ds;
max(max(max(k_2)))
k_3= - multiprod(multiprod(Y_2,Fun2_a.N(:,:,end-or(k),:),[1,2]),Y_2,[1,2])-Fun2_a.M(:,:,end-or(k),:) - multiprod(Y_2,Fun2_a.H(:,:,end-or(k),:),[1,2]) - multiprod(Fun2_a.G(:,:,end-or(k),:),Y_2,[1,2]);
Y_3= Y_a_0-k_3*ds*2;
max(max(max(k_3)))
k_4= - multiprod(multiprod(Y_3,Fun2_a.N(:,:,end-or(k)-1,:),[1,2]),Y_3,[1,2])-Fun2_a.M(:,:,end-or(k)-1,:) - multiprod(Y_3,Fun2_a.H(:,:,end-or(k)-1,:),[1,2]) - multiprod(Fun2_a.G(:,:,end-or(k)-1,:),Y_3,[1,2]);
max(max(max(k_4)))
Y_a_0= Y_a_0 - (1/6)*(k_1+2*k_2+2*k_3+k_4)*ds*2;
RK.Y_a(:,:,k+1,:)=Y_a_0;
Y_a_0
end
%b
Y_b_0=Fun2_b.Y_H(:,:,end,:);
RK.Y_b(:,:,1,:)=Y_b_0;
for k=1:length(or)-1
k_1= - multiprod(multiprod(Y_b_0,Fun2_b.N(:,:,end-or(k)+1,:),[1,2]),Y_b_0,[1,2])- Fun2_b.M(:,:,end-or(k)+1,:) - multiprod(Y_b_0,Fun2_b.H(:,:,end-or(k)+1,:),[1,2]) - multiprod(Fun2_b.G(:,:,end-or(k)+1,:),Y_b_0,[1,2]);
Y_1= Y_b_0-k_1*ds;
k_2= - multiprod(multiprod(Y_1,Fun2_b.N(:,:,end-or(k),:),[1,2]),Y_1,[1,2])-Fun2_b.M(:,:,end-or(k),:) - multiprod(Y_1,Fun2_b.H(:,:,end-or(k),:),[1,2]) - multiprod(Fun2_b.G(:,:,end-or(k),:),Y_1,[1,2]);
Y_2= Y_b_0-k_2*ds;
k_3= - multiprod(multiprod(Y_2,Fun2_b.N(:,:,end-or(k),:),[1,2]),Y_2,[1,2])-Fun2_b.M(:,:,end-or(k),:) - multiprod(Y_2,Fun2_b.H(:,:,end-or(k),:),[1,2]) - multiprod(Fun2_b.G(:,:,end-or(k),:),Y_2,[1,2]);
Y_3= Y_b_0-k_3*ds*2;
k_4= - multiprod(multiprod(Y_3,Fun2_b.N(:,:,end-or(k)-1,:),[1,2]),Y_3,[1,2])-Fun2_b.M(:,:,end-or(k)-1,:) - multiprod(Y_3,Fun2_b.H(:,:,end-or(k)-1,:),[1,2]) - multiprod(Fun2_b.G(:,:,end-or(k)-1,:),Y_3,[1,2]);
Y_b_0= Y_b_0 - (1/6)*(k_1+2*k_2+2*k_3+k_4)*ds*2;
RK.Y_b(:,:,k+1,:)=Y_b_0;    
end
%a_b
Y_a_b_0=Fun2_a_b.Y_H(:,:,end,:,:);
RK.Y_a_b(:,:,1,:,:)=Y_a_b_0;
for k=1:length(or)-1
k_1= - multiprod(multiprod(Y_a_b_0,Fun2_a_b.N(:,:,end-or(k)+1,:,:),[1,2]),Y_a_b_0,[1,2])- Fun2_a_b.M(:,:,end-or(k)+1,:,:) - multiprod(Y_a_b_0,Fun2_a_b.H(:,:,end-or(k)+1,:,:),[1,2]) - multiprod(Fun2_a_b.G(:,:,end-or(k)+1,:,:),Y_a_b_0,[1,2]);
Y_1= Y_a_b_0-k_1*ds;
k_2= - multiprod(multiprod(Y_1,Fun2_a_b.N(:,:,end-or(k),:,:),[1,2]),Y_1,[1,2])-Fun2_a_b.M(:,:,end-or(k),:,:) - multiprod(Y_1,Fun2_a_b.H(:,:,end-or(k),:,:),[1,2]) - multiprod(Fun2_a_b.G(:,:,end-or(k),:,:),Y_1,[1,2]);
Y_2= Y_a_b_0-k_2*ds;
k_3= - multiprod(multiprod(Y_2,Fun2_a_b.N(:,:,end-or(k),:,:),[1,2]),Y_2,[1,2])-Fun2_a_b.M(:,:,end-or(k),:,:) - multiprod(Y_2,Fun2_a_b.H(:,:,end-or(k),:,:),[1,2]) - multiprod(Fun2_a_b.G(:,:,end-or(k),:,:),Y_2,[1,2]);
Y_3= Y_a_b_0-k_3*ds*2;
k_4= - multiprod(multiprod(Y_3,Fun2_a_b.N(:,:,end-or(k)-1,:,:),[1,2]),Y_3,[1,2])-Fun2_a_b.M(:,:,end-or(k)-1,:,:) - multiprod(Y_3,Fun2_a_b.H(:,:,end-or(k)-1,:,:),[1,2]) - multiprod(Fun2_a_b.G(:,:,end-or(k)-1,:,:),Y_3,[1,2]);
Y_a_b_0= Y_a_b_0 - (1/6)*(k_1+2*k_2+2*k_3+k_4)*ds*2;
RK.Y_a_b(:,:,k+1,:,:)=Y_a_b_0;    
end



% RK.Y_a(:,:,end,1)%-RK.Y_b(:,:,end,3)
% size(RK.Y_a)