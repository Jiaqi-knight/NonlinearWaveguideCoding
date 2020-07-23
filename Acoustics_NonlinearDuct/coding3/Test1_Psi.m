%PFunction for \psi functions.
%m:circumferiential mode
%n:radial mode
%dimension:
%   2D           -  \X_{\alpha\beta}
%   3D           -  \X_{\alpha\beta\gamma}
%P Options(op):
%   ab           -  \X_{\alpha\beta}[rr]
%   pr_ab        -  \X_{[\alpha]\beta}[rr][cos\psi]
%   ps_ab        -  \X_{\{\alpha\}\beta}[rr][sin\psi]
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
%-----------------------------------------------------------------%

clc
%clear
close all
%opengl('save','software')

load('X2.mat');load('X3.mat');
subfunction_path1='.\subfunction1'
addpath(genpath(subfunction_path1));
subfunction_path2='.\chebfun-master'
addpath(genpath(subfunction_path2));
m=-5:5;
n=4;
h=0.1;tau=1;
%% X-2D-3D
op={'ab','pr_ab','ps_ab'};op_m={'1','r','r2'};
op_name2={'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}'};
op_name3={'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}'};

Op={'ab','pt_ab','ps_ab','ab_cos','ab_sin','pt_ab_cos','pt_ab_sin','ps_ab_cos','ps_ab_sin'};
Op2_name={'\Theta_{\alpha\beta}','\Theta_{(\alpha)\beta}','\Theta_{\{\alpha\}\beta}','\Theta_{\alpha\beta}[cos\phi]',...
    '\Theta_{\alpha\beta}[sin\phi]','\Theta_{(\alpha)\beta}[cos\phi]','\Theta_{(\alpha)\beta}[sin\phi]','\Theta_{\{\alpha\}\beta}[cos\phi]','\Theta_{\{\alpha\}\beta}[sin\phi]'};
Op3_name={'\Theta_{\alpha\beta\gamma}','\Theta_{(\alpha)\beta\gamma}','\Theta_{\{\alpha\}\beta\gamma}','\Theta_{\alpha\beta\gamma}[cos\phi]',...
    '\Theta_{\alpha\beta\gamma}[sin\phi]','\Theta_{(\alpha)\beta\gamma}[cos\phi]','\Theta_{(\alpha)\beta\gamma}[sin\phi]','\Theta_{\{\alpha\}\beta\gamma}[cos\phi]','\Theta_{\{\alpha\}\beta\gamma}[sin\phi]'};

%%
plotName2={'\Psi_{\alpha\beta}[r]','\Psi_{\alpha\beta}[rcos\phi]','\Psi_{\alpha\beta}[r^2cos\phi]','\Psi_{\alpha\beta}[r^2 sin\phi]','\Psi_{(\alpha)\beta}[r]','\Psi_{\{\alpha\}\beta}'}
plotName3={'\Psi_{\alpha\beta\gamma}[r]','\Psi_{\alpha\beta\gamma}[rcos\phi]','\Psi_{\alpha\beta\gamma}[r^2cos\phi]','\Psi_{\alpha\beta\gamma}[r^2 sin\phi]','\Psi_{(\alpha)\beta\gamma}[r]','\Psi_{\{\alpha\}\beta\gamma}'}

P2{1}=X2{1,2}.*Theta(m,n,2,'ab');%Fig1:\Psi_{\alpha\beta}[r]
P2{2}=X2{1,2}.*Theta(m,n,2,'ab_cos');
P2{3}=X2{1,3}.*Theta(m,n,2,'ab_cos');
P2{4}=X2{1,3}.*Theta(m,n,2,'ab_sin');
P2{5}=X2{1,2}.*Theta(m,n,2,'pt_ab');
P2{6}=X2{3,1}.*Theta(m,n,2,'ab')+X2{1,1}.*Theta(m,n,2,'ps_ab');



H=figure;colormap(jet); 
for k=1:length(plotName2)
h1=subplot(4,6,k);image(real(P2{k}),'CDataMapping','scaled');grid on;axis equal;axis off;
xlim(h1,[0 44]);ylim(h1,[0 44]);title([plotName2{k}])
set(h1,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);

h2=subplot(4,6,6+k);image(imag(P2{k}),'CDataMapping','scaled');grid on;axis equal;axis off;
xlim(h2,[0 44]);ylim(h2,[0 44]);title([plotName2{k}])
set(h2,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);

end


P3{1}=X3{1,2}.*Theta(m,n,3,'ab');%Fig1:\Psi_{\alpha\beta}[r]
P3{2}=X3{1,2}.*Theta(m,n,3,'ab_cos');
P3{3}=X3{1,3}.*Theta(m,n,3,'ab_cos');
P3{4}=X3{1,3}.*Theta(m,n,3,'ab_sin');
P3{5}=X3{1,2}.*Theta(m,n,3,'pt_ab');
P3{6}=X3{3,1}.*Theta(m,n,3,'ab')+X3{1,1}.*Theta(m,n,3,'ps_ab');
for k=1:length(plotName2)
h01=subplot(4,6,12+k);    
h1 = slice(real(P3{k}), [], [], 1:size(P3{k},3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
set(h01,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
title(plotName3{k});

h02=subplot(4,6,18+k);
h1 = slice(imag(P3{k}), [], [], 1:size(P3{k},3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
set(h02,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
title(plotName3{k});

end
colormap(jet); 
colorbar('Position',[0.05 0.35 0.02 0.35]);


