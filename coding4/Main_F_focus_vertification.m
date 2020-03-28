% Main-Function
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
% -----------------------------------------------------------------%
% For Boudnary, only end point is need to solved!
clc
clear
%close all
opengl('save','software')
%load('Database_X2.mat');load('Database_X3.mat');
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
%% #######Geometry########%
%Geo.s =logspace(0,1,10);
Geo.s =linspace(0,10,500);
%Geo.h=0.1*exp(linspace(0,1.5,length(Geo.s)));
Geo.h=1*ones(size(Geo.s));

%Geo.tau=1/2./Geo.h;
%Geo.tau=linspace(0,4,500)
%Geo.tau =logspace(0,0.6,500)-1;
Geo.tau =0.0*ones(size(Geo.s));
%Geo.tau =logspace(0,0.1,500)-1;

%Geo.kappa=(2/3)./Geo.h;
Geo.kappa=logspace(0,0.6,500)-1

Geo.sw=sqrt(Geo.kappa.^2+Geo.tau.^2).*Geo.s;
Geo.x = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*sin(Geo.sw+0);Geo.y = Geo.kappa./(Geo.kappa.^2+Geo.tau.^2).*cos(Geo.sw+0);Geo.z = Geo.tau./(Geo.kappa.^2+Geo.tau.^2).*Geo.sw;
Geo.theta_0=cumsum(Geo.tau.*[0 diff(Geo.s)]);
Geo.h_diff=[(Geo.h(2)-Geo.h(1))/(Geo.s(2)-Geo.s(1)) diff(Geo.h)./diff(Geo.s)];
gamma=1.4;
h=figure
tubeplot(Geo.x,Geo.y,Geo.z,Geo.h,Geo.s,50);hold on;plot3(Geo.x, Geo.y, Geo.z);daspect([1,1,1]); camlight;
saveas(h,[save_directory,'\','Geo','.png'])
%% #######Wave########%
Geo.m=[-3:3];
Geo.n=3;
a=[1 2 3]; %P^{a}=P^{a*},U^{a}=U^{a*}
b=[-3 -2 -1 1 2 3];
a_b=a-b.';
k=0.95*1.8412/Geo.h(end);

%%
Geo_b.m=Geo.m;
Geo_b.n=Geo.n;
Geo_b.s=Geo.s;                  %[Geo.s(1) Geo.s(end)];
Geo_b.h=Geo.h;                  %[Geo.h(1) Geo.h(end)];
Geo_b.kappa=Geo.kappa;          %[Geo.kappa(1) Geo.kappa(end)];
Geo_b.tau=Geo.tau;              %[Geo.tau(1) Geo.tau(end)];
Geo_b.h_diff=Geo.h_diff         %[Geo.h_diff(1) Geo.h_diff(end)]

%% #######Function########%
%% different h can be similarly deal with
tic
X2.ab=         X(Geo_b,2,'ab','1');        %1/h
X2.ab_r=       X(Geo_b,2,'ab','r');        %1
X2.ab_r2 =     X(Geo_b,2,'ab','r2');       %h
X2.pr_ab_r=    X(Geo_b,2,'pr_ab','r');     %1/h
X2.pr_ab_r2=   X(Geo_b,2,'pr_ab','r2');    %1
X2.a_pr_b_r =  X(Geo_b,2,'a_pr_b','r');    %1
X2.a_ps_b_r =  X(Geo_b,2,'a_ps_b','r');    %h'/h
X2.ps_ab_r =   X(Geo_b,2,'ps_ab','r');     %h'/h
toc;

%%
T2.ab=         Theta(Geo_b,2,'ab');
T2.ps_ab=      Theta(Geo_b,2,'ps_ab');
T2.pt_ab=      Theta(Geo_b,2,'pt_ab');
T2.ab_cos=     Theta(Geo_b,2,'ab_cos');
T2.a_ps_b=     Theta(Geo_b,2,'a_ps_b');
T2.pt_ab_cos=  Theta(Geo_b,2,'pt_ab_cos');

%%
Psi.ab_r=           X2.ab_r         .*T2.ab;
Psi.a_pr_b_r =      X2.a_pr_b_r     .*T2.ab;
Psi.ab_r2_cos =     X2.ab_r2        .*T2.ab_cos;
Psi.ps_ab_r_1 =     X2.ps_ab_r      .*T2.ab;
Psi.ps_ab_r_2 =     X2.ab_r         .*T2.ps_ab;
Psi.a_ps_b_r_1 =    X2.a_ps_b_r     .*T2.ab;
Psi.a_ps_b_r_2 =    X2.ab_r         .*T2.a_ps_b;
Psi.pr_ab_r=        X2.pr_ab_r      .*T2.ab;
Psi.pr_ab_r2_cos=   X2.pr_ab_r2     .*T2.ab_cos;
Psi.pt_ab=          X2.ab           .*T2.pt_ab;
Psi.pt_ab_r_cos=    X2.ab_r         .*T2.pt_ab_cos;
Psi.ab=             X2.ab           .*T2.ab;
Psi.ab_r_cos=       X2.ab_r         .*T2.ab_cos;
Psi.pt_ab_cos=      X2.ab           .*T2.pt_ab_cos;

%%
tic
%2D-bsxfun(@times, 3D(\alpha*\beta*s), 1*1*s*a)-->4D[\alpha*\beta*s*a]
Fun2_a.V=       bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a.W=       -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a))); %(James-3.27,Jiaqi-77)
Fun2_a.N=       bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*a),1,1,length(Geo_b.kappa),length(a)));%(James-3.35b)
Fun2_a.G=       -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_a.H=       -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_a.M_2_1=   bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));

Fun2_a.M_3_1=   bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));

Fun2_a.M  =     -Fun2_a.N-multiprod(Fun2_a.M_2_1,Fun2_a.V,[1,2])-multiprod(Fun2_a.M_3_1,Fun2_a.W,[1,2]);%(James-3.35a)


Fun2_a.L=               [-Fun2_a.G  -Fun2_a.M;...
    Fun2_a.N   Fun2_a.H];
Fun2_a.L2=              multiprod(Fun2_a.L,Fun2_a.L,[1,2]);
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        [Fun2_a.Vb_H(:,:,kh,ka),Fun2_a.Lambda_H(:,:,kh,ka)]=    eig(( Fun2_a.L(:,:,kh,ka))); %sqrt(-1)*
    end
end
toc



close all
fieldnames = {'LL2.gif'};
k=1;
figure(3);    set(gcf,'outerposition',get(0,'screensize'));%最大化
suptitle('Spectrum of L \kappa=[0,2],h=1,k=1.75,\tau=0,m=[-5:5],n=[1:3]');
for kh=1
    Lambda_H1=diag(Fun2_a.Lambda_H(:,:,kh,1));
    Lambda_H2=diag(Fun2_a.Lambda_H(:,:,kh,2));
    Lambda_H3=diag(Fun2_a.Lambda_H(:,:,kh,3));
    
    subplot(4,6,[1 2 7 8]);grid on
    plot(real(Lambda_H1),imag(Lambda_H1),'MarkerSize',1,'Marker','x','LineWidth',1,'LineStyle','none','Color',[1,0,0]);hold on
    %ylim([-50 50]);xlim([-20 20]);
    title('a=1');ylabel('\lambda')
    
    subplot(4,6,[3 4 9 10]);grid on
    plot(real(Lambda_H2),imag(Lambda_H2),'MarkerSize',1,'Marker','x','LineWidth',1,'LineStyle','none','Color',[1,0,0]);hold on
    %ylim([-50 50]);xlim([-20 20]);
    
    subplot(4,6,[5 6 11 12]);grid on
    plot(real(Lambda_H3),imag(Lambda_H3),'MarkerSize',1,'Marker','x','LineWidth',1,'LineStyle','none','Color',[1,0,0]);hold on
    %ylim([-50 50]); xlim([-20 20]);
    title('a=3')
    
    
    subplot(4,6,[13 14]);grid on;
    plot(real(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
        'Color',[1 0 0]);hold on;ylabel('real')
    %ylim([-10 10]);
    subplot(4,6,[19 20]);grid on
    plot(imag(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
        'Color',[1 0 0]);hold on;ylabel('imag')
    %ylim([-10 10]);
    subplot(4,6,[15 16]);grid on
    plot(real(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
        'Color',[1 0 0]);hold on;ylabel('real')
    %ylim([-10 10]);
    subplot(4,6,[21 22]);grid on
    plot(imag(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
        'Color',[1 0 0]);hold on;ylabel('imag')
    %ylim([-10 10]);
    subplot(4,6,[17 18]);grid on
    plot(real(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
        'Color',[1 0 0]);hold on;ylabel('real')
    %ylim([-10 10]);
    subplot(4,6,[23 24]);grid on
    plot(imag(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
        'Color',[1 0 0]);hold on;ylabel('imag')
    
    pause(1)
end



frame = getframe(gcf); % 获取整个窗口内容的图像
im=frame2im(frame);
[I{k},map{k}]=rgb2ind(im,256);
imwrite(I{k},map{k},fieldnames{1},'gif','Loopcount',Inf,'DelayTime',1);

waitingTime=1-log(1:length(Geo_b.h))/log(length(Geo_b.h));


for kh=2:length(Geo_b.h)
    Lambda_H1=diag(Fun2_a.Lambda_H(:,:,kh,1));
    Lambda_H2=diag(Fun2_a.Lambda_H(:,:,kh,2));
    Lambda_H3=diag(Fun2_a.Lambda_H(:,:,kh,3));
    k=k+1;
    subplot(4,6,[1 2 7 8]);grid on
    plot(real(Lambda_H1),imag(Lambda_H1),'MarkerSize',0.2,'Marker','x','LineWidth',0.15,'LineStyle','none','Color',[0,0,0]);hold on
    %ylim([-50 50]);xlim([-20 20]);
    
    subplot(4,6,[3 4 9 10]);grid on
    plot(real(Lambda_H2),imag(Lambda_H2),'MarkerSize',0.2,'Marker','x','LineWidth',0.15,'LineStyle','none','Color',[0,0,0]);hold on
    %ylim([-50 50]);xlim([-20 20]);
    
    subplot(4,6,[5 6 11 12]);grid on
    plot(real(Lambda_H3),imag(Lambda_H3),'MarkerSize',0.2,'Marker','x','LineWidth',0.15,'LineStyle','none','Color',[0,0,0]);hold on
    %ylim([-50 50]); xlim([-20 20]);
    
    
    subplot(4,6,[13 14]);grid on
    plot(real(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
        'Color',[0 0 0]);hold on;ylabel('real')
    %ylim([-10 10]);
    subplot(4,6,[19 20]);grid on
    plot(imag(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
        'Color',[0 0 0]);hold on;ylabel('imag')
    %ylim([-10 10]);
    subplot(4,6,[15 16]);grid on
    plot(real(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
        'Color',[0 0 0]);hold on;ylabel('real')
    %ylim([-10 10]);
    subplot(4,6,[21 22]);grid on
    plot(imag(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
        'Color',[0 0 0]);hold on;ylabel('imag')
    %ylim([-10 10]);
    subplot(4,6,[17 18]);grid on
    plot(real(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
        'Color',[0 0 0]);hold on;ylabel('real')
    %ylim([-10 10]);
    subplot(4,6,[23 24]);grid on
    plot(imag(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
        'Color',[0 0 0]);hold on;ylabel('imag')
    
    
    frame = getframe(gcf);% 获取整个窗口内容的图像
    im=frame2im(frame);
    [I{k},map{k}]=rgb2ind(im,256);
    %追加模式
    imwrite(I{k},map{k},fieldnames{1},'gif','WriteMode','append','DelayTime',waitingTime(k));
    
    %pause(0.1)
end

% for kh=2:length(Geo_b.h)
% subplot(1,3,1);grid on
% Lambda_H=diag(Fun2_a.Lambda_H(:,:,kh,1));
% plot(real(Lambda_H),imag(Lambda_H),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none');hold on
% ylim([-50 50]);xlim([-20 20]);
%
% subplot(1,3,2);grid on
% Lambda_H=diag(Fun2_a.Lambda_H(:,:,kh,2));
% plot(real(Lambda_H),imag(Lambda_H),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none');hold on
% ylim([-50 50]);xlim([-20 20]);
%
% subplot(1,3,3);grid on
% Lambda_H=diag(Fun2_a.Lambda_H(:,:,kh,3));
% plot(real(Lambda_H),imag(Lambda_H),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none');hold on
% ylim([-50 50]); xlim([-20 20]);
% pause(0.1)
% end



