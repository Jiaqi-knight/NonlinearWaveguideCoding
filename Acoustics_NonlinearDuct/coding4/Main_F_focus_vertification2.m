% Main-Function
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
% -----------------------------------------------------------------%
% For Boudnary, only end point is need to solved!
%only for straight duct

clc
clear
%close all
opengl('save','software')
%load('Database_X2.mat');load('Database_X3.mat');
subfunction_path4='C:\Users\Jiaqi-knight\Documents\GitHub\NonlinearWaveguideCoding\coding4\subfunction1';
addpath(genpath(subfunction_path4));

save_directory='./data';
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
%% #######Geometry########%
Geo.s =linspace(0,10,1);
Geo.h=1*ones(size(Geo.s));
Geo.tau=0./Geo.h;
Geo.kappa=(0)./Geo.h;
%Geo.h_diff=[(Geo.h(2)-Geo.h(1))/(Geo.s(2)-Geo.s(1)) diff(Geo.h)./diff(Geo.s)];
Geo.h_diff=[0]
%% #######Wave########%
Geo.m=[-23:23];
Geo.n=1;
a=[1]; %P^{a}=P^{a*},U^{a}=U^{a*}
%b=[-3 -2 -1 1 2 3];
%a_b=a-b.';
%                                                                                                                                                              
Wave.k=3*1.8411/Geo.h(end)-0.0000000001*i; %+-;
                                                                              
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
%X2.a_ps_b_r =  X(Geo_b,2,'a_ps_b','r');    %h'/h
%X2.ps_ab_r =   X(Geo_b,2,'ps_ab','r');     %h'/h
toc;

%%
T2.ab=         Theta(Geo_b,2,'ab');
%T2.ps_ab=      Theta(Geo_b,2,'ps_ab');
T2.pt_ab=      Theta(Geo_b,2,'pt_ab');
T2.ab_cos=     Theta(Geo_b,2,'ab_cos');
%T2.a_ps_b=     Theta(Geo_b,2,'a_ps_b');
T2.pt_ab_cos=  Theta(Geo_b,2,'pt_ab_cos');
%%
 Psi.ab_r=           X2.ab_r         .*T2.ab;
 Psi.a_pr_b_r =      X2.a_pr_b_r     .*T2.ab;
 Psi.ab_r2_cos =     X2.ab_r2        .*T2.ab_cos;
% Psi.ps_ab_r_1 =     X2.ps_ab_r      .*T2.ab;
% Psi.ps_ab_r_2 =     X2.ab_r         .*T2.ps_ab;
% Psi.a_ps_b_r_1 =    X2.a_ps_b_r     .*T2.ab;
% Psi.a_ps_b_r_2 =    X2.ab_r         .*T2.a_ps_b;
 Psi.pr_ab_r=        X2.pr_ab_r      .*T2.ab;
 Psi.pr_ab_r2_cos=   X2.pr_ab_r2     .*T2.ab_cos;
 Psi.pt_ab=          X2.ab           .*T2.pt_ab;
 Psi.pt_ab_r_cos=    X2.ab_r         .*T2.pt_ab_cos;
 Psi.ab=             X2.ab           .*T2.ab;
 Psi.pt_ab_r=        X2.ab_r         .*T2.pt_ab;

 Psi.ab_r_cos=       X2.ab_r         .*T2.ab_cos;
 Psi.a_pr_b_r1=specialFun(zeros(size(Geo_b.h)),Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h')-Psi.ab- Psi.pr_ab_r
 (max(Psi.a_pr_b_r1-Psi.a_pr_b_r))


%%
tic
%2D-bsxfun(@times, 3D(\alpha*\beta*s), 1*1*s*a)-->4D[\alpha*\beta*s*a]
Fun2_a.V=       bsxfun(@times,Psi.a_pr_b_r1,reshape(1./(sqrt(-1)*Wave.k*a),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a.W=    (  -bsxfun(@times,Psi.pt_ab,  reshape(1./(sqrt(-1)*Wave.k*a),1,1,1,length(a)))); %(James-3.27,Jiaqi-77)
Fun2_a.N=        bsxfun(@times,Psi.ab_r,   reshape((sqrt(-1)*Wave.k*a),1,1,1,length(a)));
%...
%    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*a),1,1,length(GPsi.a_pr_b_r1eo_b.kappa),length(a)));%(James-3.35b)
%Fun2_a.G= specialFun(zeros(size(Geo_b.h)),Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h`h')   -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
%    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
%Fun2_a.H=       -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
%    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_a.M_2_1=    -bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a)),1,1,1,length(a))) +bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));%...
              % +specialFun(0*ones(size(Geo_b.h)),Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h(1-\kappa*h*cos\psi)');
% figure
% image( real(specialFun(zeros(size(Geo_b.h)),Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h(1-\kappa*h*cos\psi)')))
% figure
% image( imag(specialFun(zeros(size(Geo_b.h)),Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h(1-\kappa*h*cos\psi)')))

Fun2_a.M_3_1=   bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))-bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));

%Fun2_a.M  =    -Fun2_a.N+multiprod(1*Fun2_a.M_2_1,Fun2_a.V,[1,2])-multiprod(1*Fun2_a.M_3_1,Fun2_a.W,[1,2]);%(James-3.35a)
Fun2_a.M  =    -Fun2_a.N+Fun2_a.M_2_1*Fun2_a.V-Fun2_a.M_3_1*Fun2_a.W;%(James-3.35a)



Fun2_a.NM= Fun2_a.N*          Fun2_a.M       %multiprod(Fun2_a.N,Fun2_a.M,[1,2]);
%Fun2_a.L=               [-Fun2_a.G  -Fun2_a.M;...
%                          Fun2_a.N   Fun2_a.H];
n_matrix=length(Geo_b.m)*Geo_b.n;
           


% figure
% image(real(N(:,:,1,1)))
% figure
% image(imag(N(:,:,1,1)))
% figure
% image(real(M(:,:,1,1)))
% figure
% image(imag(M(:,:,1,1)))
% 
% figure
% image(real(Fun2_a.NM(:,:,1,1)))
% figure
% image(imag(Fun2_a.NM(:,:,1,1)))
% 
% 
Or=[Geo_b.n*fliplr(1:length(Geo_b.m))];
for k=1:Geo_b.n-1
Or=[Geo_b.n*fliplr(1:length(Geo_b.m))-k;Or;]
end




for kh=1:length(Geo_b.h)
    for ka=1:length(a)  
        %Fun2_a.Y(:,:,kh,ka) =                                 sqrt(-1)*Fun2_a.N_inv(:,:,kh,ka)*sqrtm(Fun2_a.NM(:,:,kh,ka));
        [Fun2_a.Vb(:,:,kh,ka),Fun2_a.Lambda(:,:,kh,ka)]=      eigs(sqrtm((Fun2_a.NM(:,:,kh,ka)*1)),n_matrix);
        %Fun2_a.Lambda(:,:,kh,ka)=sqrtm(Fun2_a.Lambda(:,:,kh,ka));
        %Fun2_a.YN(:,:,kh,ka)=Fun2_a.Y(:,:,kh,ka)*Fun2_a.N(:,:,kh,ka);
        %[Fun2_a.Wb(:,:,kh,ka),Fun2_a.Lambda1(:,:,kh,ka)]=     eigs(Fun2_a.YN(:,:,kh,ka),n_matrix);        
       % Fun2_a.Vb_inv(:,:,kh,ka)=                             inv(Fun2_a.Vb(:,:,kh,ka));
        %Fun2_a.Wb_inv(:,:,kh,ka)=                             inv(Fun2_a.Wb(:,:,kh,ka));
    end
end                      
                      
% for kh=1:length(Geo_b.h)
%     for ka=1:length(a)
%         [Fun2_a.Vb_H(:,:,kh,ka),Fun2_a.Lambda_H(:,:,kh,ka)]=    eig(( Fun2_a.L(:,:,kh,ka))); %sqrt(-1)*
% %         eigenvalues=diag(Fun2_a.Lambda_H(:,:,kh,ka));
% %         order(:,kh,ka)=find(real(eigenvalues)<0);
% %         Fun2_a.Vb_HH(:,:,kh,ka)=Fun2_a.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order(:,kh,ka),kh,ka);
% %         Fun2_a.Vb_inv_HH(:,:,kh,ka)=inv(Fun2_a.Vb_HH(:,:,kh,ka));
% %         Fun2_a.Lambda_HH(:,:,kh,ka)=diag(eigenvalues(order(:,kh,ka)));
%     end
% end
% toc






% Fun2_a.Vb_H(22:42,:,kh,ka)
% clear eigvalue order
%for hk=1:length(Geo_b.h)
 eigvalue=diag(Fun2_a.Lambda(:,:,1,1));
% eigvalue_h=diag(Fun2_a.Lambda_H(:,:,50,1));
% 
% order=(find(real(eigvalue)<0));
% temp1=eigvalue(order);
%end



[Base]=BaseJ1(Geo_b.m,Geo_b.n,1);
Eig=sqrt(Wave.k^2-Base.jmn_pm.^2)
figure;
plot(imag(i*Eig),real(i*Eig),'square');hold on
plot(-imag(i*Eig),-real(i*Eig),'square');

%% Vertification

%figure;
plot(real(eigvalue),imag(eigvalue),'x');
% hold on
% plot(imag(i*temp1),real(i*temp1),'x')
% plot(imag(eigvalue_h),real(eigvalue_h),'o');


% Duct Modes

%%
% mode=[-3:3];n_len=3;
%  N =131;Ratio=0.01; [D,r] = cheb(N,Ratio,1); 
%  Mx=0.0*ones(N+1,1);%Mx=0.9-r.^2*0.9;%Figure11 (b)
% frequency=[1000];
% rT=1
% c=343;
% [initialEigValue,mode_enlarge,cutOffLine,len]=wm2initialEigValue(N,D,r,Ratio,Mx,Wave.k,Geo_b.m,Geo_b.n-1);


%frequency*2*pi*rT/c










%%
% 
% 
% close all
% fieldnames = {'LL2.gif'};
% k=1;
% figure(3);    set(gcf,'outerposition',get(0,'screensize'));%最大化
% suptitle('Spectrum of L \kappa=[0,2],h=1,k=1.75,\tau=0,m=[-5:5],n=[1:3]');
% for kh=1
%     Lambda_H1=diag(Fun2_a.Lambda_H(:,:,kh,1));
%     Lambda_H2=diag(Fun2_a.Lambda_H(:,:,kh,2));
%     Lambda_H3=diag(Fun2_a.Lambda_H(:,:,kh,3));
%     
%     subplot(4,6,[1 2 7 8]);grid on
%     plot(real(Lambda_H1),imag(Lambda_H1),'MarkerSize',1,'Marker','x','LineWidth',1,'LineStyle','none','Color',[1,0,0]);hold on
%     %ylim([-50 50]);xlim([-20 20]);
%     title('a=1');ylabel('\lambda')
%     
%     subplot(4,6,[3 4 9 10]);grid on
%     plot(real(Lambda_H2),imag(Lambda_H2),'MarkerSize',1,'Marker','x','LineWidth',1,'LineStyle','none','Color',[1,0,0]);hold on
%     %ylim([-50 50]);xlim([-20 20]);
%     
%     subplot(4,6,[5 6 11 12]);grid on
%     plot(real(Lambda_H3),imag(Lambda_H3),'MarkerSize',1,'Marker','x','LineWidth',1,'LineStyle','none','Color',[1,0,0]);hold on
%     %ylim([-50 50]); xlim([-20 20]);
%     title('a=3')
%     
%     
%     subplot(4,6,[13 14]);grid on;
%     plot(real(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
%         'Color',[1 0 0]);hold on;ylabel('real')
%     %ylim([-10 10]);
%     subplot(4,6,[19 20]);grid on
%     plot(imag(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
%         'Color',[1 0 0]);hold on;ylabel('imag')
%     %ylim([-10 10]);
%     subplot(4,6,[15 16]);grid on
%     plot(real(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
%         'Color',[1 0 0]);hold on;ylabel('real')
%     %ylim([-10 10]);
%     subplot(4,6,[21 22]);grid on
%     plot(imag(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
%         'Color',[1 0 0]);hold on;ylabel('imag')
%     %ylim([-10 10]);
%     subplot(4,6,[17 18]);grid on
%     plot(real(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
%         'Color',[1 0 0]);hold on;ylabel('real')
%     %ylim([-10 10]);
%     subplot(4,6,[23 24]);grid on
%     plot(imag(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',1,'Marker','x','LineWidth',0.15,'LineStyle','none',...
%         'Color',[1 0 0]);hold on;ylabel('imag')
%     
%     pause(1)
% end
% 
% 
% 
% frame = getframe(gcf); % 获取整个窗口内容的图像
% im=frame2im(frame);
% [I{k},map{k}]=rgb2ind(im,256);
% imwrite(I{k},map{k},fieldnames{1},'gif','Loopcount',Inf,'DelayTime',1);
% 
% waitingTime=1-log(1:length(Geo_b.h))/log(length(Geo_b.h));
% 
% 
% for kh=2:length(Geo_b.h)
%     Lambda_H1=diag(Fun2_a.Lambda_H(:,:,kh,1));
%     Lambda_H2=diag(Fun2_a.Lambda_H(:,:,kh,2));
%     Lambda_H3=diag(Fun2_a.Lambda_H(:,:,kh,3));
%     k=k+1;
%     subplot(4,6,[1 2 7 8]);grid on
%     plot(real(Lambda_H1),imag(Lambda_H1),'MarkerSize',0.2,'Marker','x','LineWidth',0.15,'LineStyle','none','Color',[0,0,0]);hold on
%     %ylim([-50 50]);xlim([-20 20]);
%     
%     subplot(4,6,[3 4 9 10]);grid on
%     plot(real(Lambda_H2),imag(Lambda_H2),'MarkerSize',0.2,'Marker','x','LineWidth',0.15,'LineStyle','none','Color',[0,0,0]);hold on
%     %ylim([-50 50]);xlim([-20 20]);
%     
%     subplot(4,6,[5 6 11 12]);grid on
%     plot(real(Lambda_H3),imag(Lambda_H3),'MarkerSize',0.2,'Marker','x','LineWidth',0.15,'LineStyle','none','Color',[0,0,0]);hold on
%     %ylim([-50 50]); xlim([-20 20]);
%     
%     
%     subplot(4,6,[13 14]);grid on
%     plot(real(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
%         'Color',[0 0 0]);hold on;ylabel('real')
%     %ylim([-10 10]);
%     subplot(4,6,[19 20]);grid on
%     plot(imag(Lambda_H1),Geo_b.kappa(kh)*ones(size(Lambda_H1)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
%         'Color',[0 0 0]);hold on;ylabel('imag')
%     %ylim([-10 10]);
%     subplot(4,6,[15 16]);grid on
%     plot(real(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
%         'Color',[0 0 0]);hold on;ylabel('real')
%     %ylim([-10 10]);
%     subplot(4,6,[21 22]);grid on
%     plot(imag(Lambda_H2),Geo_b.kappa(kh)*ones(size(Lambda_H2)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
%         'Color',[0 0 0]);hold on;ylabel('imag')
%     %ylim([-10 10]);
%     subplot(4,6,[17 18]);grid on
%     plot(real(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
%         'Color',[0 0 0]);hold on;ylabel('real')
%     %ylim([-10 10]);
%     subplot(4,6,[23 24]);grid on
%     plot(imag(Lambda_H3),Geo_b.kappa(kh)*ones(size(Lambda_H3)),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none',...
%         'Color',[0 0 0]);hold on;ylabel('imag')
%     
%     
%     frame = getframe(gcf);% 获取整个窗口内容的图像
%     im=frame2im(frame);
%     [I{k},map{k}]=rgb2ind(im,256);
%     %追加模式
%     imwrite(I{k},map{k},fieldnames{1},'gif','WriteMode','append','DelayTime',waitingTime(k));
%     
%     %pause(0.1)
% end
% 
% % for kh=2:length(Geo_b.h)
% % subplot(1,3,1);grid on
% % Lambda_H=diag(Fun2_a.Lambda_H(:,:,kh,1));
% % plot(real(Lambda_H),imag(Lambda_H),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none');hold on
% % ylim([-50 50]);xlim([-20 20]);
% %
% % subplot(1,3,2);grid on
% % Lambda_H=diag(Fun2_a.Lambda_H(:,:,kh,2));
% % plot(real(Lambda_H),imag(Lambda_H),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none');hold on
% % ylim([-50 50]);xlim([-20 20]);
% %
% % subplot(1,3,3);grid on
% % Lambda_H=diag(Fun2_a.Lambda_H(:,:,kh,3));
% % plot(real(Lambda_H),imag(Lambda_H),'MarkerSize',0.2,'Marker','.','LineWidth',0.15,'LineStyle','none');hold on
% % ylim([-50 50]); xlim([-20 20]);
% % pause(0.1)
% % end
% 
% 
% 
