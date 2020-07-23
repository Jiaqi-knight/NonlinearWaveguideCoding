
load('X2.mat');load('X3.mat');
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
Cmn1=(sqrt(-1)).^m./(sqrt(pi)*h*BesselValue_pm.*sqrt(1-m.^2./jmn_pm.^2));%Note-23

%% X-2D
op_m={'1','r','r2'};
op={'ab','pr_ab','ps_ab'};
op_name2={'X_{\alpha\beta}','X_{[\alpha]\beta}','X_{\{\alpha\}\beta}'};
op_name3={'X_{\alpha\beta\gamma}','X_{[\alpha]\beta\gamma}','X_{\{\alpha\}\beta\gamma}'};

tic
figure
for k=1:length(op)
    for kk=1:length(op_m)
%         [X0]= deltaJ(h,Cmn1,jmn_pm,m,n,2,0,op{k},op_m{kk});
%         [X1]= deltaJ(h,Cmn1,jmn_pm,m,n,2,1,op{k},op_m{kk});
%         [X_1]= deltaJ(h,Cmn1,jmn_pm,m,n,2,-1,op{k},op_m{kk});
        %X2{k,kk}=X0+X1+X_1;
        h1=subplot(4,9,(k-1)*length(op)+kk);image(real(X2{k,kk}),'CDataMapping','scaled');grid on;axis equal;axis off;title([op_name2{k},'[',op_m{kk},']'])
        xlim(h1,[0.5 44.5]);ylim(h1,[0.5 44]);set(h1,'CLim',[-(2*pi) (2*pi)],'DataAspectRatio',[1 1 1]);
        h2=subplot(4,9,9+(k-1)*length(op)+kk);image(imag(X2{k,kk}),'CDataMapping','scaled');grid on;axis equal;axis off;title([op_name2{k},'[',op_m{kk},']'])
        xlim(h2,[0.5 44.5]);ylim(h2,[0.5 44]);set(h1,'CLim',[-(2*pi) (2*pi)],'DataAspectRatio',[1 1 1]);
        colormap(jet);
    end
end
toc
tic

for k=1:length(op)
    for kk=1:length(op_m)
%         [X0]= deltaJ(h,Cmn1,jmn_pm,m,n,3,0,op{k},op_m{kk});
%         [X1]= deltaJ(h,Cmn1,jmn_pm,m,n,3,1,op{k},op_m{kk});
%         [X_1]= deltaJ(h,Cmn1,jmn_pm,m,n,3,-1,op{k},op_m{kk});
%         X3{k,kk}=X0+X1+X_1;
        h01=subplot(4,9,18+(k-1)*length(op)+kk);
        h1 = slice(real(X3{k,kk}), [], [], 1:size(X3{k,kk},3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
        set(h01,'CLim',[-1 1],'DataAspectRatio',[1 1 1]);
        title([op_name3{k},'[',op_m{kk},']'])
        h02=subplot(4,9,27+(k-1)*length(op)+kk);
        h1 = slice(imag(X3{k,kk}), [], [], 1:size(X3{k,kk},3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
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


