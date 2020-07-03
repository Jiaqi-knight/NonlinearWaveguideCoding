%本模拟程序中采用REV尺度（多孔介质内部的相变传热）
  %{
  相变半径   ，格子解析度为  ，ste数 ，糊状区系数C为 
  %%@@@@@@@@修改时间2016.08.23（糊状区视为多孔介质）
  %}
clear
global tNS cxNS cyNS tT cxT cyT m n
%% %%%%%%%%%%%%%%%%%%%可以调整的一些参数%%%%%%%%%%%%%%%%%%%%%%
 %物理模型的相关参数设置
uMax     = 0.1;                                                
lx       = 250;                   %y方向的格子数
ly       = 250;                   %x方向的格子数
Nu1      = 0;                     %记录迭代次数
Thot     = 1;
Tcold    = 0;
T0       = 0;                    %初始温度
m=lx;    n=ly;                    %格子数转换
%数学模型的相关参数设置
delta_x  = 1./n;
delta_t  = uMax/n;
gr       = delta_t*delta_t/(delta_x);
buoyancy = [0,gr];
%% 
Pr       = 1.0;
Ra       = 1e6; 
St       = 5.0;                   %斯蒂芬数（在无量纲情况下显热与潜热之比，st=Cp(Th-Tc)/La）
Tms      = 0.15;                  %相变区域的中心温度-相变半径（近固体区域边界温度）
Tradius  = 0.20;
Tml      = Tms + Tradius;         %相变区域的中心温度+相变半径（近液体区域边界温度）
%%
C        = 1600;                 %范围为10^5--10^6
%%
Dam      = 1e-4;                  %小darcy数，多孔介质渗透能力的大小。达西数越大，则表示多孔介质可渗透度越大，对于流体区强度的影响越小
Dal      = 1e9;                   %大darcy数系数
K1       = 1.0;                   %相变材料固体导热系数比流体导热系数
%%
Ens      = St*Tms;               %相变区域界定开始融化的焓值
Enl      = St*Tml+1;             %相变区域界定结束融化的焓值
LaCp     = 1/St;                 %对应的是La/Cp
nu       = sqrt(Pr/Ra)*delta_t/(delta_x*delta_x);
k        = nu/Pr;                %热扩散系数
ks1      = K1*k;                 %相变固体的热扩散系数
%% 输出结果文件，图片的相关设置
maxT        = 2900000;              %总的迭代次数
tStatistics = 2000;               %输出结果文件
tPlot       = 500;                  %输出图像
maxk        = 3;                    %温度的迭代最大次数
%tao值的计算（弛豫时间）
omegaNS  = 1./(3*nu+1.5);     
omegaT   = 1./(3*k+0.5);
omegaT1  = 1./(3*ks1+0.5);

%% %%%%%%%%%%%%%%%D2Q9 LATTICE CONSTANTS速度和温度都采用的是D2Q9的模型
%速度离散
tNS   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cxNS  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cyNS  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
oppNS = [  1,   4,  5,  2,  3,    8,   9,   6,   7];
%温度离散
tT    = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cxT   = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cyT   = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
oppT  = [  1,   4,  5,  2,  3,    8,   9,   6,   7];
%% %%%%%%%%%%%%%%%%%%划分网格%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = meshgrid(1:n,1:m);

%% %%%%%%%%%%%%%%%%%%%各场的初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho  = ones(1,m,n);                  
ux   = zeros(1,m,n);
uy   = zeros(1,m,n);
T    = zeros(1,m,n);
rfl  = zeros(1,m,n);                %液相比场的初始化 

%{
%ux  = importdata('E:\LBM学习\Matlab\毕业论文计算\糊状区REV\纯相变材料\Ra1e6\Tms0.1Tml0.1\Da1e-6\36000-velocity-x.txt');
%ux  =reshape(ux,1,m,n);
%uy  = importdata('E:\LBM学习\Matlab\毕业论文计算\糊状区REV\纯相变材料\Ra1e6\Tms0.1Tml0.1\Da1e-6\36000-velocity-y.txt');
%uy  =reshape(uy,1,m,n);
%T  = importdata('E:\LBM学习\Matlab\毕业论文计算\糊状区REV\纯相变材料\Ra1e6\Tms0.1Tml0.1\Da1e-6\36000-temperature.txt');
%T  =reshape(T,1,m,n);
%}
T(:,1,:)=1;                         %高温侧温度
T(:,m,:)=0;                         %低温侧温度
rfl(:,1,:)=1;                       %高温侧温度
%平衡态速度演化方程
for i=1:9
    cuNS         = 3*(cxNS(i)*ux+cyNS(i)*uy);
    fEq(i,:,:)   = rho .* tNS(i) .* ...
                   ( 1 + cuNS + 1/2*(cuNS.*cuNS) - 3/2*(ux.^2+uy.^2) );
    fIn(i,:,:)   = fEq(i,:,:);
end

%平衡态温度演化方程
 for i=1:9
       cu          = 3*(cxT(i)*ux+cyT(i)*uy);
       tEq(i,:,:)  = T .* tT(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2)  );          %本模拟采用单温度LB，平衡态分布函数采用p50（3-49）其中 total=1
       tIn(i,:,:)  = tEq(i,:,:);
 end

 omegaTm     = ones(1,m,n);            %
 omegaTm     = omegaT.*omegaTm;        %
 Fetia       = zeros(1,m,n);
 c0          = zeros(1,m,n);           % 两个参数c1，c2，参照p47  3-26
 c1          = zeros(1,m,n);      
 Foc1        = zeros(1,m,n);
 Foc2        = zeros(1,m,n);
 
 %% %%相场位置的确定%%%%%%%%%%%
solid  = find(rfl<0.0002);                     %PCM的固体区域
mush   = find(0.0002<=rfl & rfl<1);            %PCM的糊状区区域
liquid = find(1<=rfl & rfl<2);                 %PCM的流体区域
%% %%%%%%%%%%%%%%%%%%%%%%%%主程序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 tic
 
for cycle = 1:maxT
    Nu1     =  Nu1+1;
    
    T       =  sum(tIn);
    rflt    =  rfl;                               %上一个时间步长的含液率
    rho     =  sum(fIn);   
    Tk      =  T;
    
    ux(solid)=0;   uy(solid)=0;
%% while循环语句（流体区域（）的计算过程）为单温度分布方程
%{
条件循环while end结构
%格式：while （conditions）
%                commands-1
%            end   
 %               commands-2
%①  作用
%    条件conditions成立的时候 执行commands-1 当遇到end检查条件，不满足则执行commands-2
%    注意commands-1应该有使得条件改变的句子，否则就成了死循环
%}
    kl      =  0;   
    while 1  
        kl=kl+1;
        %Tb=T;
      for i=1:9
            cu          = 3*(cxT(i)*ux+cyT(i)*uy);
            tEq(i,:,:)  = Tk .* tT(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2)  );
            Sr(i,:,:)   = -LaCp.*tT(i).*reshape((rfl-rflt),1,m,n).*(1+(1-omegaT/2.0)*cu); 
            tOut(i,:,:)  = tIn(i,:,:) - omegaT .* (tIn(i,:,:)-tEq(i,:,:))+Sr(i,:,:);
     end
    
     omegaTm(mush)  = 1./(3*(k.*rfl(mush)+ks1.*(1-rfl(mush)))+0.5);      %(k.*rfl(mush)+ks1.*(1-rfl(mush))为有效的 a
      for i=1:9
            cu (mush)       = 3*(cxT(i)*ux(mush)+cyT(i)*uy(mush));
            tEq(i,mush)  = Tk(mush) .* tT(i).* ( 1 + cu(mush) + 1/2.*(cu(mush).*cu(mush)) - 3/2*(ux(mush).^2+uy(mush).^2)  );
            Sr(i,mush)   = -LaCp.*tT(i).*((rfl(mush)-rflt(mush))).*(1+(1-omegaTm(mush)./2.0).*cu(mush)); 
            tOut(i,mush)  = tIn(i,mush) -  omegaTm(1,mush).* (tIn(i,mush)-tEq(i,mush))+Sr(i,mush);
      end
     
      for i=1:9
          %  cu        = 3*(cxT(i)*ux(solid)+cyT(i)*uy(solid));
            tEq(i,solid)  = Tk(solid) .* tT(i); %.* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux(solid).^2+uy(solid).^2)  );
            Sr(i,solid)   = -LaCp.*tT(i).*((rfl(solid)-rflt(solid)));%.*(1+(1-omegaT/2.0)*cu); 
            tOut(i,solid)  = tIn(i,solid) - omegaT1 .* (tIn(i,solid)-tEq(i,solid))+Sr(i,solid);
     end 
    
        T=sum(tOut);
        rflk=rfl;
        %（下面）确定新的液相比
        En  = St*T + reshape(rflk,1,m,n);      %焓的计算式（格子单位下）
        rfl(En<Ens) =0;                        %新的液相比为0（PCM固态）
        Earea=find(Ens<=En & En<=Enl);         %选出新的液相比为大于0，小于1的部分Earea
        rfl(Earea)=(En(Earea)-Ens)/(Enl-Ens);  %计算出糊状区的液相比
        rfl(Enl<En)=1;                         %新的液相比为1（PCM液态）
        
        solid  = find(rfl<0.0002);                     %PCM的固体区域
        mush   = find(0.0002<=rfl & rfl<1);            %PCM的糊状区区域
        liquid = find(1<=rfl & rfl<2);                 %PCM的流体区域

                %km=sum(sum(T-Tb))/sum(sum(T));
        
       if (kl>maxk)
           break;
       end
    end
    
%% %%%%%%%%%%%%%%%%%%在多孔介质部分的流动情况计算（整个流场里）Rev尺度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 全场流动均采用多孔介质内部流动方程，针对没有多孔介质时，对应的Da数取
%}
     Da=ones(1,m,n);
     Da=Dal.*Da;
     Fetia   = zeros(1,m,n);
     %mush =[mush1;mush2];
     Da(mush)=(n^2)*nu.*(rfl(mush).^3)./(C.*(1-rfl(mush)).^2);
     Kshen=Da;
     Fetia(mush)=1.75./(150.*rfl(mush).^3).^(0.5);
     c0=0.5*(1+0.5*nu.*rfl./Kshen);
     c1=0.5*Fetia.*rfl./(Kshen.^(1/2));
     vx  = reshape ( (cxNS * reshape(fIn,9,m*n)), 1,m,n) ./rho;
     vy  = reshape ( (cyNS * reshape(fIn,9,m*n)), 1,m,n) ./rho + 0.5*gr.*rfl.*T;
     ux  = vx./(c0+sqrt(c0.^2+c1.*sqrt(vx.^2+vy.^2)));
     uy  = vy./(c0+sqrt(c0.^2+c1.*sqrt(vx.^2+vy.^2)));
        ux(solid)=0;
        uy(solid)=0;
     Foc1 = -nu./Kshen.*rfl.*ux-2.*c1.*sqrt(ux.^2+uy.^2).*ux;
     Foc2 = -nu./Kshen.*rfl.*uy-2.*c1.*sqrt(ux.^2+uy.^2).*uy+gr.*rfl.*T;
     %%{
    for i=1:9
      cuNS         = 3*(cxNS(i)*ux+cyNS(i)*uy);
      fEq(i,:,:)   = rho .* tNS(i) .* ...
                       ( 1 + cuNS + 1/(2.*rfl).*(cuNS.*cuNS) - 3/(2.*rfl).*(ux.^2+uy.^2) );
      force(i,:,:) = (1-omegaNS/2.0).*tNS(i).*rho.*...
                     ((3.*cyNS(i)-(3./rfl).*uy+(9./rfl).*cxNS(i).*cyNS(i).*ux+(9./rfl).*cyNS(i).*cyNS(i).*uy).*Foc2+...
                     (3.*cxNS(i)-(3./rfl).*ux+(9./rfl).*cxNS(i).*cxNS(i).*ux+(9./rfl).*cxNS(i).*cyNS(i).*uy).*Foc1);
      fOut(i,:,:)  = fIn(i,:,:) - omegaNS .* (fIn(i,:,:)-fEq(i,:,:)) + force(i,:,:);
    end 
      %} 
     %{
      流体区域
     for i=1:9
      cuNS(liquid)         = 3*(cxNS(i)*ux(liquid)+cyNS(i)*uy(liquid));
      fEq(i,liquid)   = rho(liquid) .* tNS(i) .* ...
                       ( 1 + cuNS(liquid) + 1./(2*rfl(liquid)).*(cuNS(liquid).*cuNS(liquid)) - 3./(2*rfl(liquid)).*(ux(liquid).^2+uy(liquid).^2) );
      force(i,liquid) = (1-omegaNS/2.0).*tNS(i).*rho(liquid).*...
                     ((3.*cyNS(i)-(3./rfl(liquid)).*uy(liquid)+(9./rfl(liquid)).*cxNS(i).*cyNS(i).*ux(liquid)+(9./rfl(liquid)).*cyNS(i).*cyNS(i).*uy(liquid)).*Foc2(liquid)+...
                      (3.*cxNS(i)-(3./rfl(liquid)).*ux(liquid)+(9./rfl(liquid)).*cxNS(i).*cxNS(i).*ux(liquid)+(9./rfl(liquid)).*cxNS(i).*cyNS(i).*uy(liquid)).*Foc1(liquid));
      fOut(i,liquid)  = fIn(i,liquid) - omegaNS .* (fIn(i,liquid)-fEq(i,liquid)) + force(i,liquid);
     end  
      %糊状区域
      for i=1:9
      cuNS(mush)         = 3*(cxNS(i)*ux(mush)+cyNS(i)*uy(mush));
      fEq(i,mush)   = rho(mush) .* tNS(i) .* ...
                       ( 1 + cuNS(mush) + 1./(2*rfl(mush)).*(cuNS(mush).*cuNS(mush)) - 3./(2*rfl(mush)).*(ux(mush).^2+uy(mush).^2) );
      force(i,mush) = (1-omegaNS/2.0).*tNS(i).*rho(mush).*...
                     ((3.*cyNS(i)-(3./rfl(mush)).*uy(mush)+(9./rfl(mush)).*cxNS(i).*cyNS(i).*ux(mush)+(9./rfl(mush)).*cyNS(i).*cyNS(i).*uy(mush)).*Foc2(mush)+...
                      (3.*cxNS(i)-(3./rfl(mush)).*ux(mush)+(9./rfl(mush)).*cxNS(i).*cxNS(i).*ux(mush)+(9./rfl(mush)).*cxNS(i).*cyNS(i).*uy(mush)).*Foc1(mush));
      fOut(i,mush)  = fIn(i,mush) - omegaNS .* (fIn(i,mush)-fEq(i,mush)) + force(i,mush);
      end  
      %}
     %相变材料固相区域
     for i=1:9
         fEq(i,solid)   = rho(solid) .* tNS(i);
         fOut(i,solid)  = fIn(i,solid) - omegaNS .* (fIn(i,solid)-fEq(i,solid));
     end
    
     for i=1:9
       fOut(i,solid)  =  fIn(oppNS(i),solid);
     end
     fOut  = reshape(fOut,9,m,n);
    %% %%%%%%%%%%%%%%%%%%%%%%%%迁     移%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % STREAMING STEP FLUID
   for i=1:9
      fIn(i,:,:) = circshift(fOut(i,:,:), [0,cxNS(i),cyNS(i)]);
   end
   % STREAMING STEP FLUID
   for i=1:9
      tIn(i,:,:) = circshift(tOut(i,:,:), [0,cxT(i),cyT(i)]);
   end  
   
   tIn(4,m,:) = Tcold-tIn(1,m,:)-tIn(2,m,:)-tIn(3,m,:)-tIn(5,m,:)-tIn(6,m,:)-tIn(7,m,:)-tIn(8,m,:)-tIn(9,m,:);

   tIn(2,1,:)  = Thot-(tIn(1,1,:)+tIn(9,1,:)+tIn(3,1,:)+tIn(4,1,:)+tIn(5,1,:)+tIn(6,1,:)+tIn(7,1,:)+tIn(8,1,:));
    %% %%%%%%%%%%%%%%%%%%%%% 边  界  条  件 %%%%%%%%%%%%%%%%%%%%%%%%
  %绝热边界条件
   for i=1:9
       cuu=3*(cxT(i)*ux(:,:,1)+cyT(i)*uy(:,:,1));
       tEq(i,:,1)  = T(:,:,2) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,1).^2+uy(:,:,1).^2)  );
       tIn(i,:,1)  = tEq(i,:,1) + tIn(i,:,2)-tEq(i,:,2);
   end
   for i=1:9
       cuu=3*(cxT(i)*ux(:,:,n)+cyT(i)*uy(:,:,n));
       tEq(i,:,n)  = T(:,:,n-1) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,n).^2+uy(:,:,n).^2)  );
       tIn(i,:,n)  = tEq(i,:,n) + tIn(i,:,n-1)-tEq(i,:,n-1);
   end

    fIn  = reshape(fIn,9,m,n);
    
   for i=1:9
        fIn(i,:,1)  = fIn(oppNS(i),:,2);
        fIn(i,:,n)  = fIn(oppNS(i),:,n-1);                                                                                               
        fIn(i,1,:)  = fIn(oppNS(i),2,:);
        fIn(i,m,:)  = fIn(oppNS(i),m-1,:);
   end
  
   %% %%%%%%%%%%%%%%%%%%%数据输出部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   % VISUALIZATION结果文件的输出
   if(mod(cycle,tPlot)==0)
       Nu1                                      %迭代次数
    %% 速度、温度计算
       T =  sum(tIn);
       ux(solid)=0;   uy(solid)=0;
       ux1 = reshape(ux,m,n);     uy1 = reshape(uy,m,n); 
       u1 =sqrt(ux1.^2+uy1.^2);
       T1  = reshape(T,m,n); 
       rfl1=reshape(rfl,m,n);
       
     %% 图像输出
           subplot(2,2,1); 
           imagesc(u1(:,n:-1:1)');colorbar
           title('Fluid velocity');
           axis image; drawnow
           subplot(2,2,2); 
           imagesc(rfl1(:,n:-1:1)');colorbar
           title('Melting Ratio');
           axis image; drawnow
           subplot(2,2,3);
           imagesc(T1(:,n:-1:1)');colorbar
           title(['Temperature' num2str(Nu1)]);
           axis image; drawnow
    %% 速度温度数据的输出
        if (mod(cycle,tStatistics)==0)
            
       %速度u的文件输出
       stru = [num2str(Nu1),'-velocity.txt'];
       dlmwrite(stru,u1)
       strux = [num2str(Nu1),'-velocity-x.txt'];
       dlmwrite(strux,ux1)
       struy = [num2str(Nu1),'-velocity-y.txt'];
       dlmwrite(struy,uy1)
       %温度的文件的输出
       strT = [num2str(Nu1),'-temperature.txt'];
       dlmwrite(strT,T1)
       %含液率的输出
       strr = [num2str(Nu1),'-rfl.txt'];
       dlmwrite(strr,rfl1)
      %图片的输出
       
       end
       toc
   end
   
  end
