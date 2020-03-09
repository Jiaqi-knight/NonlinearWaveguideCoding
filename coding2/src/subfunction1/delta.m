%DELAT Function for \delta functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.
function [D1,order_0]= delta(m,n,dimention,k,multiply_factor) %k+ right



if dimention==2
    m1=m;
    m2=fliplr(m);
    for k1=1:length(m1)
        for k2=1:length(m2)
            D(k1,k2)=m1(k1)+m2(k2)+k;
        end
    end
    order_0=D==0;
    if nargin ==5  %multipy by m
        multiply_factor=m;
        D=order_0.*repmat(multiply_factor,length(m),1);
    elseif nargin ==4
        D=order_0;
    end
    D1=repelem(D,n,n);
    
    
elseif dimention==3
    m1=m;
    m2=fliplr(m);
    m3=fliplr(m);
    for k1=1:length(m1)
        for k2=1:length(m2)
            for k3=1:length(m1)
                D(k1,k2,k3)=m1(k1)+m2(k2)+m3(k3)+k;
            end
        end
    end
    order_0=D==0;
    if nargin ==5  %multipy by m
        multiply_factor=m;
        D=order_0.*repmat(multiply_factor,[length(m),1,length(m)]);    %right-hand-axis,m->i
    elseif nargin ==4
        D=order_0;
    end
    D1=repelem(D,n,n,n);
    
end
% figure
% subplot(2,3,3);image(real(D1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')
% subplot(2,3,4);image(imag(D1),'CDataMapping','scaled');grid on;axis equal;axis off;title('X_{\alpha\beta}[r]')
% figure;
% h = slice(D1, [], [], 1:size(D1,3));
% set(h, 'EdgeColor','none')%, 'FaceColor','interp'
% alpha(.1);axis equal

%
% %% step1.2-Theta
% Theta2=2*pi*kron(speye(length(m1)),ones(n));
% Theta2_cos_phi=sparse(pi*(kron(diag(ones(length(m1)-1,1),-1),ones(n))+kron(diag(ones(length(m1)-1,1),+1),ones(n))));
% Theta2_sin_phi=sparse(-sqrt(-1)*pi*(kron(diag(ones(length(m1)-1,1),-1),ones(n))-kron(diag(ones(length(m1)-1,1),+1),ones(n))));
% Theta2_pt=2*pi*sqrt(-1)*m*kron(speye(length(m1)),ones(n));
% Theta2_ps=-tau(1)*2*pi*sqrt(-1)*kron(m1.*speye(length(m1)),ones(n));
