%deltaT(m,n,dimention,1,multiply_factor)
%X Function for \mathacal{X} functions.
%m:circumferiential mode
%n:radial mode
%dimension:
%   2D           -  \X_{\alpha\beta}
%   3D           -  \X_{\alpha\beta\gamma}
%X Options(op):
%   ab           -  \X_{\alpha\beta}[rr]
%   pr_ab        -  \X_{[\alpha]\beta}[rr]
%   ps_ab        -  \X_{\{\alpha\}\beta}[rr]
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
%-----------------------------------------------------------------%

clc
clear
close all
% subfunction_path1='.\subfunction1'
% addpath(genpath(subfunction_path1));
subfunction_path2='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\NonlinearWaveguideCoding\coding2\src\chebfun-master'
addpath(genpath(subfunction_path2));

%% #######Geometry########%
s =logspace(0,1,50);
h=0.1*exp(linspace(0,1.5,length(s)));

kappa=(2/3)./h;tau=0.2./h;
sw=sqrt(kappa.^2+tau.^2).*s;
x = kappa./(kappa.^2+tau.^2).*sin(sw+0);y = kappa./(kappa.^2+tau.^2).*cos(sw+0);z = tau./(kappa.^2+tau.^2).*sw;
theta_0=cumsum(tau.*[0 diff(s)]);

%%
% %kron_fliplr the second line in order to form I matrix
% delta_mn_u0v=kron_fliplr(delta_mn,delta_u0v);
% delta_m0n_uv=kron_fliplr(delta_m0n,delta_uv);

m=-1:1;
n=1;
%h=0.1;


%%
m0=0:m(end);
for k=1:length(m0)
    temp1=roots(diff(chebfun(@(t) besselj(m0(k),t),[0,300]))+0.00001);
    temp2= besselj(m0(k),temp1);
    jmn(:,k)=temp1(1:n);
    BesselValue(:,k)=temp2(1:n);
end
jmn_pm=[fliplr(jmn) jmn(:,2:end)];
BesselValue_pm=[fliplr(repmat((-1).^m0,n,1).*BesselValue) BesselValue(:,2:end)];
Cmn1=bsxfun(@times,(sqrt(-1)).^m./(sqrt(pi)*BesselValue_pm.*sqrt(1-m.^2./jmn_pm.^2)),reshape(1./h,1,1,[]));%£¨Jiaqi-130£©4*11*50


%% X-2D
op_m={'1','r','r2'};
op={'ab','pr_ab','ps_ab','a_ps_b'};
op_name2={'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}','X_{\alpha\{\beta\}}'};
op_name3={'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}','X_{\alpha\{\beta\}\gamma}'};







tic
figure
for k=1:length(op)
    for kk=1:length(op_m)
        [X0]= deltaJ(s,h,Cmn1,jmn_pm,m,n,2,0,op{k},op_m{kk});
        [X1]= deltaJ(s,h,Cmn1,jmn_pm,m,n,2,1,op{k},op_m{kk});
        [X_1]= deltaJ(s,h,Cmn1,jmn_pm,m,n,2,-1,op{k},op_m{kk});
        X2{k,kk}=X0+X1+X_1;
        h1=subplot(4,12,(k-1)*length(op_m)+kk);image(real(X2{k,kk}(:,:,1)),'CDataMapping','scaled');grid on;axis equal;axis off;title([op_name2{k},'[',op_m{kk},']'])
        set(h1,'CLim',[-(2*pi) (2*pi)],'DataAspectRatio',[1 1 1]);
        h2=subplot(4,12,12+(k-1)*length(op_m)+kk);image(imag(X2{k,kk}(:,:,1)),'CDataMapping','scaled');grid on;axis equal;axis off;title([op_name2{k},'[',op_m{kk},']'])
        set(h1,'CLim',[-(2*pi) (2*pi)],'DataAspectRatio',[1 1 1]);
        colormap(jet);
    end
end
toc
tic

for k=1:length(op)
    for kk=1:length(op_m)
        [X0]= deltaJ(s,h,Cmn1,jmn_pm,m,n,3,0,op{k},op_m{kk});
        [X1]= deltaJ(s,h,Cmn1,jmn_pm,m,n,3,1,op{k},op_m{kk});
        [X_1]= deltaJ(s,h,Cmn1,jmn_pm,m,n,3,-1,op{k},op_m{kk});
        X3{k,kk}=X0+X1+X_1;
        h01=subplot(4,12,24+(k-1)*length(op_m)+kk);
        h1 = slice(real(X3{k,kk}(:,:,:,1)), [], [], 1:size(X3{k,kk}(:,:,:,1),3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
        set(h01,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
        title([op_name3{k},'[',op_m{kk},']'])
        h02=subplot(4,12,36+(k-1)*length(op_m)+kk);
        h1 = slice(imag(X3{k,kk}(:,:,:,1)), [], [], 1:size(X3{k,kk}(:,:,:,1),3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
        set(h02,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
        title([op_name3{k},'[',op_m{kk},']'])
        colorbar('Position',[0.05 0.35 0.02 0.35]);
    end
end

toc

% figure;
% h01=subplot(1,2,1)
% h1 = slice(real(X), [], [], 1:size(X,3));
% set(h1, 'EdgeColor','none')%, 'FaceColor','interp'
% alpha(.1);axis equal;
% %set(h01,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
%
% h02=subplot(1,2,2)
% h2= slice(imag(X), [], [], 1:size(X,3));
% set(h2, 'EdgeColor','none')%, 'FaceColor','interp'
% alpha(.05);axis equal;colormap(jet);
% %set(h02,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
% colorbar('Position',[0.05 0.35 0.02 0.35]);


