%Test for Theta function

clc;clear;close all
subfunction_path1='.\subfunction1'
addpath(genpath(subfunction_path1));
subfunction_path2='.\chebfun-master'
addpath(genpath(subfunction_path2));


%% #######Geometry########%
s =logspace(0,1,50);
h=0.1*exp(linspace(0,1.5,length(s)));
kappa=(2/3)./h;tau=0.2./h;
sw=sqrt(kappa.^2+tau.^2).*s;
x = kappa./(kappa.^2+tau.^2).*sin(sw+0);y = kappa./(kappa.^2+tau.^2).*cos(sw+0);z = tau./(kappa.^2+tau.^2).*sw;
theta_0=cumsum(tau.*[0 diff(s)]);


m=-5:5;
n=3;

op={'ab','pt_ab','a_pt_b','ps_ab','a_ps_b','ab_cos','ab_sin','pt_ab_cos','pt_ab_sin','ps_ab_cos','ps_ab_sin'};
op2_name={'\Theta_{\alpha\beta}','\Theta_{(\alpha)\beta}','\Theta_{\alpha(\beta)}','\Theta_{\{\alpha\}\beta}','\Theta_{\alpha\{\beta\}}','\Theta_{\alpha\beta}[cos\phi]',...
    '\Theta_{\alpha\beta}[sin\phi]','\Theta_{(\alpha)\beta}[cos\phi]','\Theta_{(\alpha)\beta}[sin\phi]','\Theta_{\{\alpha\}\beta}[cos\phi]','\Theta_{\{\alpha\}\beta}[sin\phi]'};
op3_name={'\Theta_{\alpha\beta\gamma}','\Theta_{(\alpha)\beta\gamma}','\Theta_{\alpha(\beta)\gamma}','\Theta_{\{\alpha\}\beta\gamma}','\Theta_{\alpha\{\beta\}\gamma}','\Theta_{\alpha\beta\gamma}[cos\phi]',...
    '\Theta_{\alpha\beta\gamma}[sin\phi]','\Theta_{(\alpha)\beta\gamma}[cos\phi]','\Theta_{(\alpha)\beta\gamma}[sin\phi]','\Theta_{\{\alpha\}\beta\gamma}[cos\phi]','\Theta_{\{\alpha\}\beta\gamma}[sin\phi]'};


figure
for k=1:length(op)
    T1=Theta(s,m,n,2,op{k},tau);
    h1=subplot(4,length(op),k);image(real(T1(:,:,1)),'CDataMapping','scaled');grid on;axis equal;axis off;title([op2_name{k}])
    xlim(h1,[0 44]);ylim(h1,[0 44]);
    set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);

    h2=subplot(4,length(op),k+length(op));image(imag(T1(:,:,1)),'CDataMapping','scaled');grid on;axis equal;axis off;title([op2_name{k}])
    xlim(h2,[0 44]);ylim(h2,[0 44]);
    set(h2,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);

end

for k=1:length(op)
    T1=Theta(s,m,n,3,op{k},tau);
    h01=subplot(4,length(op),2*length(op)+k)
    h1 = slice(real(T1(:,:,:,1)), [], [], 1:size(T1(:,:,:,1),3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
    set(h01,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
    title(op3_name{k});
    
    h02=subplot(4,length(op),3*length(op)+k)
    h2= slice(imag(T1(:,:,:,1)), [], [], 1:size(T1(:,:,:,1),3));set(h2, 'EdgeColor','none');alpha(.05);axis equal;colormap(jet);
    set(h02,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
    colorbar('Position',[0.05 0.35 0.02 0.35]);
    title(op3_name{k})
end



% if dimention==2
% figure
% h1=subplot(1,2,1);image(real(T1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')
% xlim(h1,[0.5 44.5]);ylim(h1,[0.5 44.5]);
% h2=subplot(1,2,2);image(imag(T1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')
% xlim(h2,[0.5 44.5]);ylim(h2,[0.5 44.5]);
% elseif dimention==3
% figure;
% h01=subplot(1,2,1)
% h1 = slice(real(T1), [], [], 1:size(T1,3));
% set(h1, 'EdgeColor','none')%, 'FaceColor','interp'
% alpha(.1);axis equal;
% set(h01,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
%
% h02=subplot(1,2,2)
% h2= slice(imag(T1), [], [], 1:size(T1,3));
% set(h2, 'EdgeColor','none')%, 'FaceColor','interp'
% alpha(.05);axis equal;colormap(jet);
% set(h02,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
% colorbar('Position',[0.05 0.35 0.02 0.35]);
% end