%DELAT Function for \delta functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.
function [D1,O]= deltaT(m,n,dimention,k,multiply_factor) %k+ right



if dimention==2
    m1=fliplr(m);
    m2=m;
    %% 简化-等价
    %     for k1=1:length(m1)
    %         for k2=1:length(m2)
    %             D(k1,k2)=m1(k1)+m2(k2)+k;
    %         end
    %     end
    %     order_0=D==0;
    %     if nargin ==5  %multipy by m
    %         multiply_factor=m;
    %         D=order_0.*repmat(multiply_factor,length(m),1);
    %     elseif nargin ==4
    %         D=order_0;
    %     end
    %     D1=repelem(D,n,n);
    O.M1=repelem(permute(repmat(m1,[length(m1),1]),[1,2]),n,n);
    O.M2=repelem(permute(repmat(m2,[length(m1),1]),[2,1]),n,n);
    O.N1=repmat(permute(repmat([1:n],[n,1]),[1,2]),length(m1),length(m1));
    O.N2=repmat(permute(repmat([1:n],[n,1]),[2,1]),length(m1),length(m1));
    
    D1=(O.M1+O.M2+k==0);
    %     O.diag0_1=D1.*O.L1;
    %     O.diag0_2=D1.*O.L2;
    if nargin ==5  %multipy by m
        O.Multiply=repelem(permute(repmat(multiply_factor,[length(m1),1]),[1,2]),n,n);
        D1=D1.* O.Multiply;
    end
    
    
    
    
elseif dimention==3
    m1=m;
    m2=fliplr(m);
    m3=fliplr(m);
    %% 简化-等价
    %     for k1=1:length(m1)
    %         for k2=1:length(m2)
    %             for k3=1:length(m3)
    %                 D(k1,k2,k3)=m1(k1)+m2(k2)+m3(k3)+k;
    %             end
    %         end
    %     end
    %     order_0=D==0;
    %     if nargin ==5  %multipy by m
    %         multiply_factor=m;
    %         D=order_0.*repmat(multiply_factor,[length(m),1,length(m)]);    %right-hand-axis,m->i
    %     elseif nargin ==4
    %         D=order_0;
    %     end
    %D1=repelem(D,n,n,n);
    O.M1=repelem(permute(repmat(m1,[length(m1),1,length(m1)]),[1,2,3]),n,n,n);
    O.M2=repelem(permute(repmat(m2,[length(m1),1,length(m1)]),[2,1,3]),n,n,n);
    O.M3=repelem(permute(repmat(m3,[length(m1),1,length(m1)]),[3,1,2]),n,n,n);
    O.N1=repmat(permute(repmat([1:n],[n,1,n]),[1,2,3]),length(m1),length(m1),length(m1));
    O.N2=repmat(permute(repmat([1:n],[n,1,n]),[2,1,3]),length(m1),length(m1),length(m1));
    O.N3=repmat(permute(repmat([1:n],[n,1,n]),[3,1,2]),length(m1),length(m1),length(m1));
    
    D1=(O.M1+O.M2+O.M3+k==0);
    %     O.diag1=D1.*O.L1;
    %     O.diag2=D1.*O.L2;
    %     O.diag3=D1.*O.L3;
    if nargin ==5  %multipy by m
        O.Multiply = repelem(permute(repmat(multiply_factor,[length(m1),1,length(m1)]),[1,2,3]),n,n,n);
        D1=D1.* O.Multiply;
    end
    %% test for 3D
    % figure;
    % subplot(1,4,1)
    % h1 = slice(real(O.N1), [], [], 1:size(O.N1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
    % %set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
    % subplot(1,4,2)
    % h1 = slice(real(O.N2), [], [], 1:size(O.N1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
    % %set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
    % subplot(1,4,3)
    % h1 = slice(real(O.N3), [], [], 1:size(O.N1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
    % %set(h1,'CLim',[-2*pi 2*pi],'DataAspectRatio',[1 1 1]);
    % subplot(1,4,4)
    % h1 = slice(real(D1), [], [], 1:size(D1,3));set(h1, 'EdgeColor','none');alpha(.1);axis equal;
    
    
end

%% test for 2D
% figure
% subplot(2,3,3);image(real(D1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')
% subplot(2,3,4);image(imag(D1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')


%% some example of Theta function
% %% step1.2-Theta
% Theta2=2*pi*kron(speye(length(m1)),ones(n));
% Theta2_cos_phi=sparse(pi*(kron(diag(ones(length(m1)-1,1),-1),ones(n))+kron(diag(ones(length(m1)-1,1),+1),ones(n))));
% Theta2_sin_phi=sparse(-sqrt(-1)*pi*(kron(diag(ones(length(m1)-1,1),-1),ones(n))-kron(diag(ones(length(m1)-1,1),+1),ones(n))));
% Theta2_pt=2*pi*sqrt(-1)*m*kron(speye(length(m1)),ones(n));
% Theta2_ps=-tau(1)*2*pi*sqrt(-1)*kron(m1.*speye(length(m1)),ones(n));
