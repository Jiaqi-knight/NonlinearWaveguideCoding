%Main-Function
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
%-----------------------------------------------------------------%
clc
%clear
close all
%opengl('save','software')
load('Database_X2.mat');load('Database_X3.mat');
subfunction_path1='.\subfunction1'
addpath(genpath(subfunction_path1));
subfunction_path2='.\chebfun-master'
addpath(genpath(subfunction_path2));

%% #######Geometry########%
Geo.s =logspace(0,1,50);
Geo.h=0.1*exp(linspace(0,1.5,length(s)));
Geo.kappa=(2/3)./h;Geo.tau=0.2./h;
Geo.sw=sqrt(kappa.^2+tau.^2).*s;
Geo.x = kappa./(kappa.^2+tau.^2).*sin(sw+0);Geo.y = kappa./(kappa.^2+tau.^2).*cos(sw+0);Geo.z = tau./(kappa.^2+tau.^2).*sw;
Geo.theta_0=cumsum(tau.*[0 diff(s)]);
%tubeplot(x,y,z,h,s,50);hold on;plot3(x, y, z);daspect([1,1,1]); camlight;
Geo.m=-5:5;
Geo.n=4;
Geo.a=[1:5]; %P^{a}=P^{a*},U^{a}=U^{a*}

%% #######Name########%
Name.op_m={'1','r','r2'};
Name.op={'ab','pr_ab','ps_ab','a_ps_b'};
Name.X2={'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}','X_{\alpha\{\beta\}}'};
Name.X3={'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}','X_{\alpha\{\beta\}\gamma}'};
Name.T={'ab','pt_ab','a_pt_b','ps_ab','a_ps_b','ab_cos','ab_sin','pt_ab_cos','pt_ab_sin','ps_ab_cos','ps_ab_sin'};
Name.T2={'\Theta_{\alpha\beta}','\Theta_{(\alpha)\beta}','\Theta_{\alpha(\beta)}','\Theta_{\{\alpha\}\beta}','\Theta_{\alpha\{\beta\}}','\Theta_{\alpha\beta}[cos\phi]',...
    '\Theta_{\alpha\beta}[sin\phi]','\Theta_{(\alpha)\beta}[cos\phi]','\Theta_{(\alpha)\beta}[sin\phi]','\Theta_{\{\alpha\}\beta}[cos\phi]','\Theta_{\{\alpha\}\beta}[sin\phi]'};
Name.T3={'\Theta_{\alpha\beta\gamma}','\Theta_{(\alpha)\beta\gamma}','\Theta_{\alpha(\beta)\gamma}','\Theta_{\{\alpha\}\beta\gamma}','\Theta_{\alpha\{\beta\}\gamma}','\Theta_{\alpha\beta\gamma}[cos\phi]',...
    '\Theta_{\alpha\beta\gamma}[sin\phi]','\Theta_{(\alpha)\beta\gamma}[cos\phi]','\Theta_{(\alpha)\beta\gamma}[sin\phi]','\Theta_{\{\alpha\}\beta\gamma}[cos\phi]','\Theta_{\{\alpha\}\beta\gamma}[sin\phi]'};

%% #######Function########%

[Base]=BaseJ(m,n,h);
[X2]= X(Base,Geo,dimention,op,op_m);


X(Base,m,n,2,'a_ps_b','r')
tic
Fun.V=bsxfun(@times,X2{4,2}.*Theta(Geo,2,'ab'),reshape(1./(sqrt(-1)*a.'*kappa),1,1,length(a),length(kappa)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun.W=bsxfun(@times,X2{1,1}.'.*Theta(Geo,2,'pt_ab'),reshape(1./(sqrt(-1)*a.'*kappa),1,1,length(a),length(kappa))); %(James-3.27,Jiaqi-77)
Fun.N=bsxfun(@times,X2{1,2}.*Theta(Geo,2,'ab'),reshape((sqrt(-1)*a.'*kappa),1,1,length(a),length(kappa)))...
    -bsxfun(@times,X2{1,3}.*Theta(Geo,2,'ab_cos'),reshape((sqrt(-1)*a.'*kappa.^2),1,1,length(a),length(kappa)));%(James-3.35b)
Fun.G=-bsxfun(@times,X2{1,2}.*Theta(Geo,2,'ab'),reshape((sqrt(-1)*a.'*kappa),1,1,length(a),length(kappa)))
% Fun.H
% Fun.A
% Fun.B
% Fun.C
% Fun.D
% Fun.E
toc