%{
修改说明：
时间2016.08.24（最新）
修改程序检测，试在本程序基础上修改为3d的程序（相变）
在3d相变程序上添加多孔介质（修改时间7月1日）
形式(lx:lz:ly)
                     |  /ly（Y）
                     | /
坐标显示      ------  |-------lz（H）
                     |
                     |lx（X）
将断掉的数据添加恢复计算，注意，需要改Nu1为断掉前一步的次数，
cycle为断掉的次数
%}
clear
global tNS cxNS cyNS czNS tT  cxT cyT czT m n p 
%% 物理参数设置部分
uMax         = 0.1;
%可以改变的参数（物理模型）
lx           = 50;                     %x方向格子数(必须和骨架矩阵相同)
aspect_ratio = 1;
lz           = aspect_ratio*lx;        %高度方向格子数
ly           = aspect_ratio*lx;        %垂直向里方向格子数
delta_x      = 1./(lx);
Thot         = 1;                      % Heating on bottom wall
Tcold        = 0;                      % Cooling on top wall

%%
%浮升力
delta_t      = uMax/lx;
gr           = delta_t*delta_t/delta_x;
buoyancy     = [0,gr,0];

Nu1=0;                                 %迭代计数
maxT         = 10000000;               %最大迭代次数
tPlot        = 20;                      %绘图间隔
tStatistics  = 400;                    %数据输出的间隔
%minT         = 1000;                  %最小迭代次数
m=lx;   n=lz;  p=ly;
%工况参数（按工况进行修改）
Pr       = 2;                         %Prandtl number
Ra       = 1e4;                        %Rayleigh number
St       = 2;                    %斯蒂芬数（在无量纲情况下显热与潜热之比，st=Cp(Th-Tc)/La）
%%相变温度半径
Tm       = 0.2;                  %相变中心温度
Tradius  = 0.10;                 %相变温度半径
Tms      = Tm -Tradius;         %相变区域的中心温度-相变半径（近固体区域边界温度）              
Tml      = Tm + Tradius;         %相变区域的中心温度+相变半径（近液体区域边界温度）  

K1       = 2.0;                  %相变材料固体导热系数比流体导热系数
K2       = 10;                   %固体骨架导热系数比流体导热系数
Ens      = St*Tms+0;             %相变区域界定开始融化的焓值
Enl      = St*Tml+1;             %相变区域界定结束融化的焓值
LaCp     = 1/St;                 %对应的是La/Cp
%%多孔介质参数
C        = 1600;                  %范围为10^5--10^6
Dam      = 1e-4;                  %小darcy数，多孔介质渗透能力的大小。达西数越大，则表示多孔介质可渗透度越大，对于流体区强度的影响越小
Dal      = 1e9;                   %大darcy数系数

%计算的部分参数
nu           = sqrt(Pr/Ra)*delta_t/(delta_x*delta_x);
k            = sqrt(1./(Pr*Ra))*delta_t/(delta_x*delta_x);
ks1          = K1*k;                 %相变固体的热扩散系数
ks2          = K2*k;                 %固体骨架的热扩散系数
%豫驰时间
omegaNS      = 1./(3.0*nu+0.5); 
omegaT       = 1./(3.0*k+0.5);
omegaT1      = 1./(3*ks1+0.5);
omegaT2      = 1./(3*ks2+0.5);
%温度收敛判据
%erro         = 1e-10;  
maxk         = 1e-9;

%% 速度模型
% D3Q19 LATTICE CONSTANTS 
tNS = [1/3, 1/18,1/18,1/18,1/18,1/18,1/18, 1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36]; 
cxNS = [ 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1,-1, 1,-1, 0, 0, 0, 0]; 
cyNS = [ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0 ,0, 0, 0, 1, 1,-1,-1]; 
czNS = [ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1,-1 ,1,-1]; 
oppNS = [1, 3, 2, 5, 4, 7, 6, 11, 10, 9, 8, 15, 14, 13, 12, 19, 18, 17, 16];
% D3Q19 LATTICE CONSTANTS 
tT  = [1/3, 1/18,1/18,1/18,1/18,1/18,1/18, 1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36]; 
cxT = [ 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1,-1, 1,-1, 0, 0, 0, 0]; 
cyT = [ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0 ,0, 0, 0, 1, 1,-1,-1]; 
czT = [ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1,-1 ,1,-1]; 
oppT = [1, 3, 2, 5, 4, 7, 6, 11, 10,9, 8, 15, 14, 13, 12, 19, 18, 17, 16];
%% 网格划分
[z,x,y] = meshgrid(1:p,1:m,1:n);
%% 初始参数的设置
T0           =Tml;                      %参考温度(Thot+Tms)/2
rho          = ones(1,m,n,p);
%{1
ux     = zeros(1,m,n,p);    
uz     = zeros(1,m,n,p);
uy     = zeros(1,m,n,p);
T      = zeros(1,m,n,p);
rfl    = zeros(1,m,n,p);               %液相比场的初始化 
%}

%% 以下将断掉的数据读入从新计算
%{
ux  = importdata('E:\LBMsimulation\3DSimulation\3DPhaseChange\HaveMush\100-velocity-x.txt');
ux  =reshape(ux,1,m,n,p);
uy  = importdata('E:\LBMsimulation\3DSimulation\3DPhaseChange\HaveMush\100-velocity-y.txt');
uy  =reshape(uy,1,m,n,p);
uz  = importdata('E:\LBMsimulation\3DSimulation\3DPhaseChange\HaveMush\100-velocity-y.txt');
uz  =reshape(uz,1,m,n,p);
T  = importdata('E:\LBMsimulation\3DSimulation\3DPhaseChange\HaveMush\100-temperature.txt');
T  =reshape(T,1,m,n,p);
rfl  = importdata('E:\LBMsimulation\3DSimulation\3DPhaseChange\HaveMush\100-rfl.txt');
rfl  =reshape(rfl,1,m,n,p);
%}
%% %%%%%%%%%%%%%%%固体骨架的导入（默认固体骨架的rfls==3）
gujia            = dlmread('C:\Gujia3D\0805-gj.dat',',');       %固体骨架以1,0方式导入
%rfl_1            = rot90(rfl_1,2);                             %矩阵逆时针旋转180度 （用图片导入后）
%rfl              = flipud(rfl_1);
rfls             =  find(gujia==1 | gujia==3);                  %防止未修改dat文件
rfl(rfls)        = 3;                                           %强制将文件中的固体改为3，改不改都可适应。
rfl              =  reshape(rfl,1,m,n,p);                       %真正的液相场初始化
%% 
T(:,1,:,:)=1;
T(:,m,:,:)=0;
rfl(:,1,:,:)=1;
%% 平衡态速度演化方程
for i=1:19
    cuNS         = 3*(cxNS(i)*ux+cyNS(i)*uy+czNS(i)*uz);
    fEq(i,:,:,:)   = rho .* tNS(i) .* ...
                   ( 1 + cuNS + 1/2*(cuNS.*cuNS) - 3/2*(ux.^2+uy.^2+uz.^2));
    fIn(i,:,:,:)   = fEq(i,:,:,:);
end
%平衡态温度演化方程
 for i=1:19
       cu          = 3*(cxT(i)*ux+cyT(i)*uy+czT(i)*uz);
       tEq(i,:,:,:)  = T .* tT(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2+uz.^2));
       tIn(i,:,:,:)  = tEq(i,:,:,:);
 end
 %%多孔介质内参数（结构函数，c1，c2的矩阵格式定义）
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 phi         = zeros(1,m,n,p);%%%%%%%%%%%%%%%%%%%%%%
 nu2         = nu.*ones(1,m,n,p);%%%%%%%%%%%%%%%%%%%
 omegaNS1    = omegaNS.*ones(1,m,n,p);%%%%%%%%%%%%%%
 Rfl         = zeros(1,m,n,p);%%%%%%%%%%%%%%中间rfl的临时变量
 
 omegaTm     = ones(1,m,n,p);          %
 omegaTm     = omegaT.*omegaTm;        %
 Fetia       = zeros(1,m,n,p);
 c0          = zeros(1,m,n,p);         % 两个参数c1，c2，参照p47  3-26
 c1          = zeros(1,m,n,p);      
 Foc1        = zeros(1,m,n,p);
 Foc2        = zeros(1,m,n,p);
 Foc3        = zeros(1,m,n,p);
 
 %% %%相场位置的确定%%%%%%%%%%%(有糊状区)
solid  = find(rfl<0.0002);                       %PCM的固体区域
mush1  = find(0.0002<=rfl & rfl<0.7);              %PCM的糊状区区域
mush2  = find(0.7<=rfl & rfl<1);              %PCM的糊状区区域2(悬浮液)
liquid = find(1<=rfl & rfl<2);                   %PCM的流体区域
solid2 = find(rfl==3);                         %确定固体骨架位置 
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN LOOP (TIME CYCLES)主程序%%%%%%%%%%%%%%%%%%%%%%
tic

for cycle = 1:maxT
   Nu1=Nu1+1;
   
   rho     = sum(fIn);                                                                  %(1,m,n,p)
   T       = sum(tIn);                                                                  %(1,m,n,p)
   rflt    =  rfl; 
   Tk      =  T;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%碰撞---迁移%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 while 1     %一直不停的循环
     kl     =  kl+1;
     Tb     =  T;
       %（下面）是流体区部分
       for i=1:19
            cu(liquid)     = 3*(cxT(i)*ux(liquid)+cyT(i)*uy(liquid)+czT(i)*uz(liquid));           
            tEq(i,liquid)  = Tk(liquid) .* tT(i) .* ( 1 + cu(liquid) + 1/2*(cu(liquid).*cu(liquid)) - 3/2*(ux(liquid).^2+uy(liquid).^2+uz(liquid).^2)  );
            Sr(i,liquid)   = -LaCp.*tT(i).*(rfl(liquid)-rflt(liquid)).*(1+(1-omegaT/2.0)*cu(liquid));                %源项
            tOut(i,liquid) = tIn(i,liquid) - omegaT .* (tIn(i,liquid)-tEq(i,liquid))+Sr(i,liquid);
       end
       %糊状区部分
      omegaTm(mush1)  = 1./(3*(k.*rfl(mush1)+ks1.*(1-rfl(mush1)))+0.5);
      for i=1:19
            cu(mush1)       = 3*(cxT(i)*ux(mush1)+cyT(i)*uy(mush1)+czT(i)*uz(mush1));
            tEq(i,mush1)    = Tk(mush1) .* tT(i).* ( 1 + cu(mush1) + 1/2*(cu(mush1).*cu(mush1)) - 3/2*(ux(mush1).^2+uy(mush1).^2+uz(mush1).^2)  );
            Sr(i,mush1)     = -LaCp.*tT(i).*((rfl(mush1)-rflt(mush1))).*(1+(1-omegaTm(mush1)./2.0).*cu(mush1)); 
            tOut(i,mush1)   = tIn(i,mush1) -  omegaTm(1,mush1).* (tIn(i,mush1)-tEq(i,mush1))+Sr(i,mush1);
      end
      
       omegaTm(mush2)  = 1./(3*(k.*rfl(mush2)+ks1.*(1-rfl(mush2)))+0.5);
      for i=1:19
            cu(mush2)       = 3*(cxT(i)*ux(mush2)+cyT(i)*uy(mush2)+czT(i)*uz(mush2));
            tEq(i,mush2)    = Tk(mush2) .* tT(i).* ( 1 + cu(mush2) + 1/2*(cu(mush2).*cu(mush2)) - 3/2*(ux(mush2).^2+uy(mush2).^2+uz(mush2).^2)  );
            Sr(i,mush2)     = -LaCp.*tT(i).*((rfl(mush2)-rflt(mush2))).*(1+(1-omegaTm(mush2)./2.0).*cu(mush2)); 
            tOut(i,mush2)   = tIn(i,mush2) -  omegaTm(1,mush2).* (tIn(i,mush2)-tEq(i,mush2))+Sr(i,mush2);
      end
      %相变材料固体部分计算
     for i=1:19    
            %cu(solid)        = 3*(cxT(i)*ux(solid)+cyT(i)*uy(solid)+czT(i)*uz(solid));
            tEq(i,solid)  = Tk(solid) .* tT(i); %.* ( 1 + cu(solid) + 1/2*(cu(solid).*cu(solid)) - 3/2*(ux(solid).^2+uy(solid).^2)  );
            Sr(i,solid)   = -LaCp.*tT(i).*((rfl(solid)-rflt(solid))); %.*(1+(1-omegaT/2.0)*cu(solid)); 
            tOut(i,solid)  = tIn(i,solid) - omegaT1 .* (tIn(i,solid)-tEq(i,solid))+Sr(i,solid);
     end 
     %固体骨架部分solid2
     for i=1:19    
            %cu(solid)        = 3*(cxT(i)*ux(solid)+cyT(i)*uy(solid)+czT(i)*uz(solid));
            tEq(i,solid2)  = Tk(solid2) .* tT(i); %.* ( 1 + cu(solid2) + 1/2*(cu(solid2).*cu(solid2)) - 3/2*(ux(solid2).^2+uy(solid2).^2)  );
            Sr(i,solid2)   = -LaCp.*tT(i).*((rfl(solid2)-rflt(solid2))); %.*(1+(1-omegaT/2.0)*cu(solid2)); 
            tOut(i,solid2)  = tIn(i,solid2) - omegaT2 .* (tIn(i,solid2)-tEq(i,solid2))+Sr(i,solid2);
     end 
     
         tOut  = reshape(tOut,19,m,n,p);
         T=sum(tOut);
         T(:,1,:,:)=1;                         %高温侧温度
         T(:,m,:,:)=0;                         %低温侧温度
         rflk   =  rfl;
       %（下面）确定新的液相比
        En  = St*T + reshape(rflk,1,m,n,p);    %焓的计算式（格子单位下）
        rfl(En<Ens) = 0;                       %新的液相比为0（PCM固态）
        Earea=find(Ens<=En & En<=Enl);         %选出新的液相比为大于0，小于1的部分Earea
        rfl(Earea)=(En(Earea)-Ens)/(Enl-Ens);  %计算出糊状区的液相比
        rfl(Enl<En)=1;                         %新的液相比为1（PCM液态）
        
        solid  = find(rfl<0.0002);             %PCM的固体区域
        mush1   = find(0.0002<=rfl & rfl<0.7); %PCM的糊状区区域
        mush2   = find(0.7<=rfl & rfl<1);      %PCM的糊状区区域2(悬浮液)
        liquid = find(1<=rfl & rfl<2);         %PCM的流体区域
         rfl(solid2)  =3;                      %固体骨架
         Rfl=rfl;                              %%%%%设置中间参数Rfl
        km  = abs(sum(sum(sum(T-Tb)))/sum(sum(sum(T))));
        
       if (kl>1 && km<=maxk)
           break;
       end
    end
 %% %%%%%%%%%%%%%%%%%%在多孔介质部分的流动情况计算（整个流场里）Rev尺度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     phi(mush2)=1-rfl(mush2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%固相颗粒体积分数
    nu2(mush2)=nu.*(1+0.25.*phi(mush2)+10.05.*phi(mush2)+0.00273.*exp(16.6.*phi(mush2)));%%%%%%%%%%多相流的粘度
    omegaNS1(mush2)  = 1./(3.*nu2(mush2)+0.5); %%%%%%%%%%%%%%% 多相流的弛豫时间
     rfl(mush2)=1;
     
     Da=ones(1,m,n,p);  
     Da=Dal.*Da;
     Fetia   = zeros(1,m,n,p);
     Da(mush1)=(n^2)*nu.*(rfl(mush1).^3)./(C.*(1-rfl(mush1)).^2);    %渗透率
     Kshen=Da;
     Fetia(mush1)=1.75./(150.*rfl(mush1).^3).^(0.5);
     c0=0.5*(1+0.5*nu.*rfl./Kshen);
     c1=0.5*Fetia.*rfl./(Kshen.^(1/2));
    vx  = reshape ( (cxNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho;
    vz  = reshape ( (czNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho + 0.5*gr.*rfl.*(T-T0);
    vy  = reshape ( (cyNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho;
    ux  = vx./(c0+sqrt(c0.^2+c1.*sqrt(vx.^2+vy.^2+vz.^2)));
    uz  = vz./(c0+sqrt(c0.^2+c1.*sqrt(vx.^2+vy.^2+vz.^2)));
    uy  = vy./(c0+sqrt(c0.^2+c1.*sqrt(vx.^2+vy.^2+vz.^2)));
    ux(solid)=0;       uy(solid)=0;             uz(solid)=0;
    ux(solid2) =0;     uy(solid2) =0;           uz(solid2) =0;
     Foc1 = -nu./Kshen.*rfl.*ux-2.*c1.*sqrt(ux.^2+uy.^2+uz.^2).*ux;
     Foc2 = -nu./Kshen.*rfl.*uy-2.*c1.*sqrt(ux.^2+uy.^2+uz.^2).*uy;
     Foc3 = -nu./Kshen.*rfl.*uz-2.*c1.*sqrt(ux.^2+uy.^2+uz.^2).*uz+gr.*rfl.*(T-T0);
 %%{    
  %速度
   for i=1:19
      cuNS           = 3*(cxNS(i)*ux+cyNS(i)*uy+czNS(i)*uz);
      fEq(i,:,:,:)   = rho .* tNS(i) .* ...
                       ( 1 + cuNS + 1/(2.*rfl).*(cuNS.*cuNS) - 3/(2.*rfl).*(ux.^2+uy.^2+uz.^2) );
     force(i,:,:,:) = (1-omegaNS/2.0).*tNS(i).*rho.*...
                      ((3.*cxNS(i)-(3./rfl).*ux+(3./rfl).*cuNS.*cxNS(i)).*Foc1+...
                       (3.*cyNS(i)-(3./rfl).*uy+(3./rfl).*cuNS.*cyNS(i)).*Foc2+...
                       (3.*czNS(i)-(3./rfl).*uz+(3./rfl).*cuNS.*czNS(i)).*Foc3);
      fOut(i,:,:,:)  = fIn(i,:,:,:) - omegaNS .* (fIn(i,:,:,:)-fEq(i,:,:,:)) + force(i,:,:,:);
   end
   
    for i=1:19
      cuNS(mush2)           = 3*(cxNS(i)*ux(mush2)+cyNS(i)*uy(mush2)+czNS(i)*uz(mush2));
      fEq(i,mush2)   = rho(mush2) .* tNS(i).* ...
                       ( 1 + cuNS(mush2) + 1./(2.*rfl(mush2)).*(cuNS(mush2).*cuNS(mush2)) - 3./(2.*rfl(mush2)).*(ux(mush2).^2+uy(mush2).^2+uz(mush2).^2) );
      force(i,mush2) = (1-omegaNS1(mush2)./2.0).*tNS(i).*rho(mush2).*...
                      ((3.*cxNS(i)-(3./rfl(mush2)).*ux(mush2)+(3./rfl(mush2)).*cuNS(mush2).*cxNS(i)).*Foc1(mush2)+...
                       (3.*cyNS(i)-(3./rfl(mush2)).*uy(mush2)+(3./rfl(mush2)).*cuNS(mush2).*cyNS(i)).*Foc2(mush2)+...
                       (3.*czNS(i)-(3./rfl(mush2)).*uz(mush2)+(3./rfl(mush2)).*cuNS(mush2).*czNS(i)).*Foc3(mush2));
      fOut(i,mush2)  = fIn(i,mush2) - omegaNS1(1,mush2) .* (fIn(i,mush2)-fEq(i,mush2)) + force(i,mush2);
    end
   
   for i=1:19
         fEq(i,solid)   = rho(solid) .* tNS(i);
         fOut(i,solid)  = fIn(i,solid) - omegaNS .* (fIn(i,solid)-fEq(i,solid));
   end   
     for i=1:19
         fEq(i,solid2)   = rho(solid2) .* tNS(i);
         fOut(i,solid2)  = fIn(i,solid2) - omegaNS .* (fIn(i,solid2)-fEq(i,solid2));
     end   
   
       fOut  = reshape(fOut,19,m,n,p);
     rfl  =  Rfl;
   %% 迁移部分
   for i=1:19
      fIn(i,:,:,:) = circshift(fOut(i,:,:,:), [0,cxNS(i),czNS(i),cyNS(i)]);
   end
   
   for i=1:19
      tIn(i,:,:,:) = circshift(tOut(i,:,:,:), [0,cxT(i),czT(i),cyT(i)]); 
   end
   
   
%% %%边界条件部分(绝热)非平衡外推格式
   %温度边界
   %左边界
   tIn(3,m,:,:) = Tcold -tIn(1,m,:,:)-tIn(2,m,:,:)-tIn(4,m,:,:)-tIn(5,m,:,:)-tIn(6,m,:,:)-tIn(7,m,:,:)-tIn(8,m,:,:)-tIn(9,m,:,:)...
         -tIn(10,m,:,:) -tIn(11,m,:,:)-tIn(12,m,:,:) -tIn(13,m,:,:)-tIn(14,m,:,:)-tIn(15,m,:,:)-tIn(16,m,:,:)-tIn(17,m,:,:)-tIn(18,m,:,:)-tIn(19,m,:,:);
  %右边界
   tIn(2,1,:,:)  = Thot -tIn(1,1,:,:)-tIn(3,1,:,:)-tIn(4,1,:,:)-tIn(5,1,:,:)-tIn(6,1,:,:)-tIn(7,1,:,:)-tIn(8,1,:,:)-tIn(9,1,:,:)...
         -tIn(10,1,:,:) -tIn(11,1,:,:)-tIn(12,1,:,:) -tIn(13,1,:,:)-tIn(14,1,:,:)-tIn(15,1,:,:)-tIn(16,1,:,:)-tIn(17,1,:,:)-tIn(18,1,:,:)-tIn(19,1,:,:);
   
   %下边界（:,:,1,:）lz=1
   for i=1:19
       cuu=3*(cxT(i)*ux(:,:,1,:)+cyT(i)*uy(:,:,1,:)+czT(i)*uz(:,:,1,:));
       tEq(i,:,1,:)  = T(:,:,2,:) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,1,:).^2+uy(:,:,1,:).^2+uz(:,:,1,:).^2));
       tIn(i,:,1,:)  = tEq(i,:,1,:) + tIn(i,:,2,:)-tEq(i,:,2,:);
   end
   %上边界（:,:,n,:）lz=n
   for i=1:19
       cuu=3*(cxT(i)*ux(:,:,n,:)+cyT(i)*uy(:,:,n,:)+czT(i)*uz(:,:,n,:));
       tEq(i,:,n,:)  = T(:,:,n-1,:) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,n,:).^2+uy(:,:,n,:).^2+uz(:,:,1,:).^2));
       tIn(i,:,n,:)  = tEq(i,:,n,:) + tIn(i,:,n-1,:)-tEq(i,:,n-1,:);
   end
   %前边界（:,:,:,1）ly=1
   for i=1:19
       cuu=3*(cxT(i)*ux(:,:,:,1)+cyT(i)*uy(:,:,:,1)+czT(i)*uz(:,:,:,1));
       tEq(i,:,:,1)  = T(:,:,:,2) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,:,1).^2+uy(:,:,:,1).^2+uz(:,:,:,1).^2));
       tIn(i,:,:,1)  = tEq(i,:,:,1) + tIn(i,:,:,2)-tEq(i,:,:,2);
    end
   %后边界（:,:,:,p） ly=p
    for i=1:19
       cuu=3*(cxT(i)*ux(:,:,:,p)+cyT(i)*uy(:,:,:,p)+czT(i)*uz(:,:,:,p));
       tEq(i,:,:,p)  = T(:,:,:,p-1) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,:,p).^2+uy(:,:,:,p).^2+uz(:,:,:,p).^2));
       tIn(i,:,:,p)  = tEq(i,:,:,p) + tIn(i,:,:,p-1)-tEq(i,:,:,p-1);
    end
       
  %固体处理（标准反弹格式）
     for i=1:19
       fIn(i,solid)  =  fIn(oppNS(i),solid);
     end
     
      for i=1:19
       fIn(i,solid2)  =  fIn(oppNS(i),solid2);
     end
    fIn           =  reshape(fIn,19,m,n,p);
     %速度边界     (修正反弹格式)
  for i=1:19
        fIn(i,:,1,:)  = fIn(oppNS(i),:,2,:);
        fIn(i,:,n,:)  = fIn(oppNS(i),:,n-1,:);                                                                                               
        fIn(i,1,:,:)  = fIn(oppNS(i),2,:,:);
        fIn(i,m,:,:)  = fIn(oppNS(i),m-1,:,:);
        fIn(i,:,:,1)  = fIn(oppNS(i),:,:,2);
        fIn(i,:,:,p)  = fIn(oppNS(i),:,:,p-1);
   end
     
   %% 图片及数据输出
    % VISUALIZATION
   if(mod(cycle,tPlot)==0)   
            Nu1                                                                                        %实时输出次数
    T    = sum(tIn); 
          T(:,1,:,:)=1;    T(:,m,:,:)=0;
    T1   = reshape(T,m,n,p);      
        ux(solid) =0;           uy(solid) =0;            uz(solid) =0;
        ux(solid2) =0;          uy(solid2) =0;           uz(solid2) =0;
     ux1   = reshape (ux,m,n,p);    uy1    = reshape (uy,m,n,p);    uz1   = reshape (uz,m,n,p);     
     ux1([1 m],:,:)=0;      uy1([1 m],:,:)=0;       uz1([1 m],:,:)=0; 
     ux1(:,[1 n],:)=0;      uy1(:,[1 n],:)=0 ;      uz1(:,[1 n],:)=0;
     ux1(:,:,[1 p])=0;      uy1(:,:,[1 p])=0 ;      uz1(:,:,[1 p])=0; 
     u1    = sqrt(ux1.^2+uy1.^2+uz1.^2); 
     rfl1  = reshape(rfl,m,n,p);     
     %MATLAB切片图片的输出        
          clf; 
           subplot(2,2,1); 
           imagesc(u1(:,n:-1:1,p/2).');    colorbar                                                      %切片速度的输出
           title('Fluid velocity');
           axis image; drawnow
           subplot(2,2,2);
           imagesc(rfl1(:,n:-1:1,p/2).');    colorbar                                                      %切片速度的输出
           title('Liquid content');
           axis image; drawnow
           subplot(2,2,3); 
           imagesc(T1(:,n:-1:1,p/2).');       colorbar                                                   %切片温度的输出
           title(['Temperature (The number of iterations ' num2str(Nu1) '-' ')']);
           axis image; drawnow
           %%{
     %MATLAB数据导出    
      if (mod(cycle,tStatistics)==0) 
     
       rfl1(find(rfl1==3))  =-0.1 ;      
       [Lx,Lz,Ly]=size(u1); 
       deltx=1./(Lx-1);
       deltz=1./(Lz-1);
       delty=1./(Ly-1);
       m1=0:deltx:1;
       n1=0:deltz:1;
       p1=0:delty:1;
       
       str=   [num2str(Nu1),'-','data.dat'];
       fid=fopen(str,'w'); 
       fprintf(fid,'title="3D-quad" Variables = "X" "Z" "Y" "T" "U" "Ux" "Uz" "Uy" "f" Zone I=%g,J=%g,K=%g,F=Point\n',Lx,Lz,Ly);
        for i1=1:Lx  
          for j1=1:Lz
              for v1=1:Ly
            fprintf(fid,'%g %g %g %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',m1(i1),n1(j1),p1(v1),T1(i1,j1,v1),u1(i1,j1,v1),ux1(i1,j1,v1),uz1(i1,j1,v1),uy1(i1,j1,v1),rfl1(i1,j1,v1));
             end
          end  
       end 
       fclose(fid);
      
      %%txt文件的导出，（以防外因造成的计算中断）
      rfl2  = reshape(rfl,m,n,p);
      rfl2  = reshape(rfl2,m,n*p);                                                                 %含液率
      ux2   = reshape (ux1,m,n*p);    uy2    = reshape (uy1,m,n*p);    uz2   = reshape (uz1,m,n*p); 
      u2    = reshape (u1,m,n*p);                                                                   %速度
      T2    = reshape (T1,m,n*p);                                                                   %温度
       %速度u的文件输出
      stru = [num2str(Nu1),'-velocity.txt'];
       dlmwrite(stru,u2)
      strux = [num2str(Nu1),'-velocity-x.txt'];
       dlmwrite(strux,ux2)
      struy = [num2str(Nu1),'-velocity-y.txt'];
       dlmwrite(struy,uy2)
      struy = [num2str(Nu1),'-velocity-z.txt'];
       dlmwrite(struy,uz2)
       %温度的文件的输出
       strT = [num2str(Nu1),'-temperature.txt'];
       dlmwrite(strT,T2)
       %含液率的输出
       strr = [num2str(Nu1),'-rfl.txt'];
       dlmwrite(strr,rfl2)
     
     end
  %}
 end

  toc
end 
   