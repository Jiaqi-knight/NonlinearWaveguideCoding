%{
�޸�˵����
ʱ��2016.08.24�����£�
�޸ĳ����⣬���ڱ�����������޸�Ϊ3d�ĳ�����䣩
��3d����������Ӷ�׽��ʣ��޸�ʱ��7��1�գ�
��ʽ(lx:lz:ly)
                     |  /ly��Y��
                     | /
������ʾ      ------  |-------lz��H��
                     |
                     |lx��X��
���ϵ���������ӻָ����㣬ע�⣬��Ҫ��Nu1Ϊ�ϵ�ǰһ���Ĵ�����
cycleΪ�ϵ��Ĵ���
%}
clear
global tNS cxNS cyNS czNS tT  cxT cyT czT m n p 
%% ����������ò���
uMax         = 0.1;
%���Ըı�Ĳ���������ģ�ͣ�
lx           = 50;                     %x���������(����͹Ǽܾ�����ͬ)
aspect_ratio = 1;
lz           = aspect_ratio*lx;        %�߶ȷ��������
ly           = aspect_ratio*lx;        %��ֱ���﷽�������
delta_x      = 1./(lx);
Thot         = 1;                      % Heating on bottom wall
Tcold        = 0;                      % Cooling on top wall

%%
%������
delta_t      = uMax/lx;
gr           = delta_t*delta_t/delta_x;
buoyancy     = [0,gr,0];

Nu1=0;                                 %��������
maxT         = 10000000;               %����������
tPlot        = 20;                      %��ͼ���
tStatistics  = 400;                    %��������ļ��
%minT         = 1000;                  %��С��������
m=lx;   n=lz;  p=ly;
%���������������������޸ģ�
Pr       = 2;                         %Prandtl number
Ra       = 1e4;                        %Rayleigh number
St       = 2;                    %˹�ٷ������������������������Ǳ��֮�ȣ�st=Cp(Th-Tc)/La��
%%����¶Ȱ뾶
Tm       = 0.2;                  %��������¶�
Tradius  = 0.10;                 %����¶Ȱ뾶
Tms      = Tm -Tradius;         %�������������¶�-���뾶������������߽��¶ȣ�              
Tml      = Tm + Tradius;         %�������������¶�+���뾶����Һ������߽��¶ȣ�  

K1       = 2.0;                  %�����Ϲ��嵼��ϵ�������嵼��ϵ��
K2       = 10;                   %����Ǽܵ���ϵ�������嵼��ϵ��
Ens      = St*Tms+0;             %�������綨��ʼ�ڻ�����ֵ
Enl      = St*Tml+1;             %�������綨�����ڻ�����ֵ
LaCp     = 1/St;                 %��Ӧ����La/Cp
%%��׽��ʲ���
C        = 1600;                  %��ΧΪ10^5--10^6
Dam      = 1e-4;                  %Сdarcy������׽�����͸�����Ĵ�С��������Խ�����ʾ��׽��ʿ���͸��Խ�󣬶���������ǿ�ȵ�Ӱ��ԽС
Dal      = 1e9;                   %��darcy��ϵ��

%����Ĳ��ֲ���
nu           = sqrt(Pr/Ra)*delta_t/(delta_x*delta_x);
k            = sqrt(1./(Pr*Ra))*delta_t/(delta_x*delta_x);
ks1          = K1*k;                 %�����������ɢϵ��
ks2          = K2*k;                 %����Ǽܵ�����ɢϵ��
%ԥ��ʱ��
omegaNS      = 1./(3.0*nu+0.5); 
omegaT       = 1./(3.0*k+0.5);
omegaT1      = 1./(3*ks1+0.5);
omegaT2      = 1./(3*ks2+0.5);
%�¶������о�
%erro         = 1e-10;  
maxk         = 1e-9;

%% �ٶ�ģ��
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
%% ���񻮷�
[z,x,y] = meshgrid(1:p,1:m,1:n);
%% ��ʼ����������
T0           =Tml;                      %�ο��¶�(Thot+Tms)/2
rho          = ones(1,m,n,p);
%{1
ux     = zeros(1,m,n,p);    
uz     = zeros(1,m,n,p);
uy     = zeros(1,m,n,p);
T      = zeros(1,m,n,p);
rfl    = zeros(1,m,n,p);               %Һ��ȳ��ĳ�ʼ�� 
%}

%% ���½��ϵ������ݶ�����¼���
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
%% %%%%%%%%%%%%%%%����Ǽܵĵ��루Ĭ�Ϲ���Ǽܵ�rfls==3��
gujia            = dlmread('C:\Gujia3D\0805-gj.dat',',');       %����Ǽ���1,0��ʽ����
%rfl_1            = rot90(rfl_1,2);                             %������ʱ����ת180�� ����ͼƬ�����
%rfl              = flipud(rfl_1);
rfls             =  find(gujia==1 | gujia==3);                  %��ֹδ�޸�dat�ļ�
rfl(rfls)        = 3;                                           %ǿ�ƽ��ļ��еĹ����Ϊ3���Ĳ��Ķ�����Ӧ��
rfl              =  reshape(rfl,1,m,n,p);                       %������Һ�ೡ��ʼ��
%% 
T(:,1,:,:)=1;
T(:,m,:,:)=0;
rfl(:,1,:,:)=1;
%% ƽ��̬�ٶ��ݻ�����
for i=1:19
    cuNS         = 3*(cxNS(i)*ux+cyNS(i)*uy+czNS(i)*uz);
    fEq(i,:,:,:)   = rho .* tNS(i) .* ...
                   ( 1 + cuNS + 1/2*(cuNS.*cuNS) - 3/2*(ux.^2+uy.^2+uz.^2));
    fIn(i,:,:,:)   = fEq(i,:,:,:);
end
%ƽ��̬�¶��ݻ�����
 for i=1:19
       cu          = 3*(cxT(i)*ux+cyT(i)*uy+czT(i)*uz);
       tEq(i,:,:,:)  = T .* tT(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2+uz.^2));
       tIn(i,:,:,:)  = tEq(i,:,:,:);
 end
 %%��׽����ڲ������ṹ������c1��c2�ľ����ʽ���壩
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 phi         = zeros(1,m,n,p);%%%%%%%%%%%%%%%%%%%%%%
 nu2         = nu.*ones(1,m,n,p);%%%%%%%%%%%%%%%%%%%
 omegaNS1    = omegaNS.*ones(1,m,n,p);%%%%%%%%%%%%%%
 Rfl         = zeros(1,m,n,p);%%%%%%%%%%%%%%�м�rfl����ʱ����
 
 omegaTm     = ones(1,m,n,p);          %
 omegaTm     = omegaT.*omegaTm;        %
 Fetia       = zeros(1,m,n,p);
 c0          = zeros(1,m,n,p);         % ��������c1��c2������p47  3-26
 c1          = zeros(1,m,n,p);      
 Foc1        = zeros(1,m,n,p);
 Foc2        = zeros(1,m,n,p);
 Foc3        = zeros(1,m,n,p);
 
 %% %%�ೡλ�õ�ȷ��%%%%%%%%%%%(�к�״��)
solid  = find(rfl<0.0002);                       %PCM�Ĺ�������
mush1  = find(0.0002<=rfl & rfl<0.7);              %PCM�ĺ�״������
mush2  = find(0.7<=rfl & rfl<1);              %PCM�ĺ�״������2(����Һ)
liquid = find(1<=rfl & rfl<2);                   %PCM����������
solid2 = find(rfl==3);                         %ȷ������Ǽ�λ�� 
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN LOOP (TIME CYCLES)������%%%%%%%%%%%%%%%%%%%%%%
tic

for cycle = 1:maxT
   Nu1=Nu1+1;
   
   rho     = sum(fIn);                                                                  %(1,m,n,p)
   T       = sum(tIn);                                                                  %(1,m,n,p)
   rflt    =  rfl; 
   Tk      =  T;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%��ײ---Ǩ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% whileѭ����䣨�������򣨣��ļ�����̣�Ϊ���¶ȷֲ�����
%{
����ѭ��while end�ṹ
%��ʽ��while ��conditions��
%                commands-1
%            end   
 %               commands-2
%��  ����
%    ����conditions������ʱ�� ִ��commands-1 ������end�����������������ִ��commands-2
%    ע��commands-1Ӧ����ʹ�������ı�ľ��ӣ�����ͳ�����ѭ��
%}
     kl      =  0;  
 while 1     %һֱ��ͣ��ѭ��
     kl     =  kl+1;
     Tb     =  T;
       %�����棩������������
       for i=1:19
            cu(liquid)     = 3*(cxT(i)*ux(liquid)+cyT(i)*uy(liquid)+czT(i)*uz(liquid));           
            tEq(i,liquid)  = Tk(liquid) .* tT(i) .* ( 1 + cu(liquid) + 1/2*(cu(liquid).*cu(liquid)) - 3/2*(ux(liquid).^2+uy(liquid).^2+uz(liquid).^2)  );
            Sr(i,liquid)   = -LaCp.*tT(i).*(rfl(liquid)-rflt(liquid)).*(1+(1-omegaT/2.0)*cu(liquid));                %Դ��
            tOut(i,liquid) = tIn(i,liquid) - omegaT .* (tIn(i,liquid)-tEq(i,liquid))+Sr(i,liquid);
       end
       %��״������
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
      %�����Ϲ��岿�ּ���
     for i=1:19    
            %cu(solid)        = 3*(cxT(i)*ux(solid)+cyT(i)*uy(solid)+czT(i)*uz(solid));
            tEq(i,solid)  = Tk(solid) .* tT(i); %.* ( 1 + cu(solid) + 1/2*(cu(solid).*cu(solid)) - 3/2*(ux(solid).^2+uy(solid).^2)  );
            Sr(i,solid)   = -LaCp.*tT(i).*((rfl(solid)-rflt(solid))); %.*(1+(1-omegaT/2.0)*cu(solid)); 
            tOut(i,solid)  = tIn(i,solid) - omegaT1 .* (tIn(i,solid)-tEq(i,solid))+Sr(i,solid);
     end 
     %����Ǽܲ���solid2
     for i=1:19    
            %cu(solid)        = 3*(cxT(i)*ux(solid)+cyT(i)*uy(solid)+czT(i)*uz(solid));
            tEq(i,solid2)  = Tk(solid2) .* tT(i); %.* ( 1 + cu(solid2) + 1/2*(cu(solid2).*cu(solid2)) - 3/2*(ux(solid2).^2+uy(solid2).^2)  );
            Sr(i,solid2)   = -LaCp.*tT(i).*((rfl(solid2)-rflt(solid2))); %.*(1+(1-omegaT/2.0)*cu(solid2)); 
            tOut(i,solid2)  = tIn(i,solid2) - omegaT2 .* (tIn(i,solid2)-tEq(i,solid2))+Sr(i,solid2);
     end 
     
         tOut  = reshape(tOut,19,m,n,p);
         T=sum(tOut);
         T(:,1,:,:)=1;                         %���²��¶�
         T(:,m,:,:)=0;                         %���²��¶�
         rflk   =  rfl;
       %�����棩ȷ���µ�Һ���
        En  = St*T + reshape(rflk,1,m,n,p);    %�ʵļ���ʽ�����ӵ�λ�£�
        rfl(En<Ens) = 0;                       %�µ�Һ���Ϊ0��PCM��̬��
        Earea=find(Ens<=En & En<=Enl);         %ѡ���µ�Һ���Ϊ����0��С��1�Ĳ���Earea
        rfl(Earea)=(En(Earea)-Ens)/(Enl-Ens);  %�������״����Һ���
        rfl(Enl<En)=1;                         %�µ�Һ���Ϊ1��PCMҺ̬��
        
        solid  = find(rfl<0.0002);             %PCM�Ĺ�������
        mush1   = find(0.0002<=rfl & rfl<0.7); %PCM�ĺ�״������
        mush2   = find(0.7<=rfl & rfl<1);      %PCM�ĺ�״������2(����Һ)
        liquid = find(1<=rfl & rfl<2);         %PCM����������
         rfl(solid2)  =3;                      %����Ǽ�
         Rfl=rfl;                              %%%%%�����м����Rfl
        km  = abs(sum(sum(sum(T-Tb)))/sum(sum(sum(T))));
        
       if (kl>1 && km<=maxk)
           break;
       end
    end
 %% %%%%%%%%%%%%%%%%%%�ڶ�׽��ʲ��ֵ�����������㣨���������Rev�߶�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     phi(mush2)=1-rfl(mush2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��������������
    nu2(mush2)=nu.*(1+0.25.*phi(mush2)+10.05.*phi(mush2)+0.00273.*exp(16.6.*phi(mush2)));%%%%%%%%%%��������ճ��
    omegaNS1(mush2)  = 1./(3.*nu2(mush2)+0.5); %%%%%%%%%%%%%%% �������ĳ�ԥʱ��
     rfl(mush2)=1;
     
     Da=ones(1,m,n,p);  
     Da=Dal.*Da;
     Fetia   = zeros(1,m,n,p);
     Da(mush1)=(n^2)*nu.*(rfl(mush1).^3)./(C.*(1-rfl(mush1)).^2);    %��͸��
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
  %�ٶ�
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
   %% Ǩ�Ʋ���
   for i=1:19
      fIn(i,:,:,:) = circshift(fOut(i,:,:,:), [0,cxNS(i),czNS(i),cyNS(i)]);
   end
   
   for i=1:19
      tIn(i,:,:,:) = circshift(tOut(i,:,:,:), [0,cxT(i),czT(i),cyT(i)]); 
   end
   
   
%% %%�߽���������(����)��ƽ�����Ƹ�ʽ
   %�¶ȱ߽�
   %��߽�
   tIn(3,m,:,:) = Tcold -tIn(1,m,:,:)-tIn(2,m,:,:)-tIn(4,m,:,:)-tIn(5,m,:,:)-tIn(6,m,:,:)-tIn(7,m,:,:)-tIn(8,m,:,:)-tIn(9,m,:,:)...
         -tIn(10,m,:,:) -tIn(11,m,:,:)-tIn(12,m,:,:) -tIn(13,m,:,:)-tIn(14,m,:,:)-tIn(15,m,:,:)-tIn(16,m,:,:)-tIn(17,m,:,:)-tIn(18,m,:,:)-tIn(19,m,:,:);
  %�ұ߽�
   tIn(2,1,:,:)  = Thot -tIn(1,1,:,:)-tIn(3,1,:,:)-tIn(4,1,:,:)-tIn(5,1,:,:)-tIn(6,1,:,:)-tIn(7,1,:,:)-tIn(8,1,:,:)-tIn(9,1,:,:)...
         -tIn(10,1,:,:) -tIn(11,1,:,:)-tIn(12,1,:,:) -tIn(13,1,:,:)-tIn(14,1,:,:)-tIn(15,1,:,:)-tIn(16,1,:,:)-tIn(17,1,:,:)-tIn(18,1,:,:)-tIn(19,1,:,:);
   
   %�±߽磨:,:,1,:��lz=1
   for i=1:19
       cuu=3*(cxT(i)*ux(:,:,1,:)+cyT(i)*uy(:,:,1,:)+czT(i)*uz(:,:,1,:));
       tEq(i,:,1,:)  = T(:,:,2,:) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,1,:).^2+uy(:,:,1,:).^2+uz(:,:,1,:).^2));
       tIn(i,:,1,:)  = tEq(i,:,1,:) + tIn(i,:,2,:)-tEq(i,:,2,:);
   end
   %�ϱ߽磨:,:,n,:��lz=n
   for i=1:19
       cuu=3*(cxT(i)*ux(:,:,n,:)+cyT(i)*uy(:,:,n,:)+czT(i)*uz(:,:,n,:));
       tEq(i,:,n,:)  = T(:,:,n-1,:) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,n,:).^2+uy(:,:,n,:).^2+uz(:,:,1,:).^2));
       tIn(i,:,n,:)  = tEq(i,:,n,:) + tIn(i,:,n-1,:)-tEq(i,:,n-1,:);
   end
   %ǰ�߽磨:,:,:,1��ly=1
   for i=1:19
       cuu=3*(cxT(i)*ux(:,:,:,1)+cyT(i)*uy(:,:,:,1)+czT(i)*uz(:,:,:,1));
       tEq(i,:,:,1)  = T(:,:,:,2) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,:,1).^2+uy(:,:,:,1).^2+uz(:,:,:,1).^2));
       tIn(i,:,:,1)  = tEq(i,:,:,1) + tIn(i,:,:,2)-tEq(i,:,:,2);
    end
   %��߽磨:,:,:,p�� ly=p
    for i=1:19
       cuu=3*(cxT(i)*ux(:,:,:,p)+cyT(i)*uy(:,:,:,p)+czT(i)*uz(:,:,:,p));
       tEq(i,:,:,p)  = T(:,:,:,p-1) .* tT(i) .* ( 1 + cuu + 1/2*(cuu.*cuu) - 3/2*(ux(:,:,:,p).^2+uy(:,:,:,p).^2+uz(:,:,:,p).^2));
       tIn(i,:,:,p)  = tEq(i,:,:,p) + tIn(i,:,:,p-1)-tEq(i,:,:,p-1);
    end
       
  %���崦����׼������ʽ��
     for i=1:19
       fIn(i,solid)  =  fIn(oppNS(i),solid);
     end
     
      for i=1:19
       fIn(i,solid2)  =  fIn(oppNS(i),solid2);
     end
    fIn           =  reshape(fIn,19,m,n,p);
     %�ٶȱ߽�     (����������ʽ)
  for i=1:19
        fIn(i,:,1,:)  = fIn(oppNS(i),:,2,:);
        fIn(i,:,n,:)  = fIn(oppNS(i),:,n-1,:);                                                                                               
        fIn(i,1,:,:)  = fIn(oppNS(i),2,:,:);
        fIn(i,m,:,:)  = fIn(oppNS(i),m-1,:,:);
        fIn(i,:,:,1)  = fIn(oppNS(i),:,:,2);
        fIn(i,:,:,p)  = fIn(oppNS(i),:,:,p-1);
   end
     
   %% ͼƬ���������
    % VISUALIZATION
   if(mod(cycle,tPlot)==0)   
            Nu1                                                                                        %ʵʱ�������
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
     %MATLAB��ƬͼƬ�����        
          clf; 
           subplot(2,2,1); 
           imagesc(u1(:,n:-1:1,p/2).');    colorbar                                                      %��Ƭ�ٶȵ����
           title('Fluid velocity');
           axis image; drawnow
           subplot(2,2,2);
           imagesc(rfl1(:,n:-1:1,p/2).');    colorbar                                                      %��Ƭ�ٶȵ����
           title('Liquid content');
           axis image; drawnow
           subplot(2,2,3); 
           imagesc(T1(:,n:-1:1,p/2).');       colorbar                                                   %��Ƭ�¶ȵ����
           title(['Temperature (The number of iterations ' num2str(Nu1) '-' ')']);
           axis image; drawnow
           %%{
     %MATLAB���ݵ���    
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
      
      %%txt�ļ��ĵ��������Է�������ɵļ����жϣ�
      rfl2  = reshape(rfl,m,n,p);
      rfl2  = reshape(rfl2,m,n*p);                                                                 %��Һ��
      ux2   = reshape (ux1,m,n*p);    uy2    = reshape (uy1,m,n*p);    uz2   = reshape (uz1,m,n*p); 
      u2    = reshape (u1,m,n*p);                                                                   %�ٶ�
      T2    = reshape (T1,m,n*p);                                                                   %�¶�
       %�ٶ�u���ļ����
      stru = [num2str(Nu1),'-velocity.txt'];
       dlmwrite(stru,u2)
      strux = [num2str(Nu1),'-velocity-x.txt'];
       dlmwrite(strux,ux2)
      struy = [num2str(Nu1),'-velocity-y.txt'];
       dlmwrite(struy,uy2)
      struy = [num2str(Nu1),'-velocity-z.txt'];
       dlmwrite(struy,uz2)
       %�¶ȵ��ļ������
       strT = [num2str(Nu1),'-temperature.txt'];
       dlmwrite(strT,T2)
       %��Һ�ʵ����
       strr = [num2str(Nu1),'-rfl.txt'];
       dlmwrite(strr,rfl2)
     
     end
  %}
 end

  toc
end 
   