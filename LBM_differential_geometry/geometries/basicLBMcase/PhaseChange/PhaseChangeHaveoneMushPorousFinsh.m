%��ģ������в���REV�߶ȣ���׽����ڲ�����䴫�ȣ�
  %{
  ���뾶   �����ӽ�����Ϊ  ��ste�� ����״��ϵ��CΪ 
  %%@@@@@@@@�޸�ʱ��2016.08.23����״����Ϊ��׽��ʣ�
  %}
clear
global tNS cxNS cyNS tT cxT cyT m n
%% %%%%%%%%%%%%%%%%%%%���Ե�����һЩ����%%%%%%%%%%%%%%%%%%%%%%
 %����ģ�͵���ز�������
uMax     = 0.1;                                                
lx       = 250;                   %y����ĸ�����
ly       = 250;                   %x����ĸ�����
Nu1      = 0;                     %��¼��������
Thot     = 1;
Tcold    = 0;
T0       = 0;                    %��ʼ�¶�
m=lx;    n=ly;                    %������ת��
%��ѧģ�͵���ز�������
delta_x  = 1./n;
delta_t  = uMax/n;
gr       = delta_t*delta_t/(delta_x);
buoyancy = [0,gr];
%% 
Pr       = 1.0;
Ra       = 1e6; 
St       = 5.0;                   %˹�ٷ������������������������Ǳ��֮�ȣ�st=Cp(Th-Tc)/La��
Tms      = 0.15;                  %�������������¶�-���뾶������������߽��¶ȣ�
Tradius  = 0.20;
Tml      = Tms + Tradius;         %�������������¶�+���뾶����Һ������߽��¶ȣ�
%%
C        = 1600;                 %��ΧΪ10^5--10^6
%%
Dam      = 1e-4;                  %Сdarcy������׽�����͸�����Ĵ�С��������Խ�����ʾ��׽��ʿ���͸��Խ�󣬶���������ǿ�ȵ�Ӱ��ԽС
Dal      = 1e9;                   %��darcy��ϵ��
K1       = 1.0;                   %�����Ϲ��嵼��ϵ�������嵼��ϵ��
%%
Ens      = St*Tms;               %�������綨��ʼ�ڻ�����ֵ
Enl      = St*Tml+1;             %�������綨�����ڻ�����ֵ
LaCp     = 1/St;                 %��Ӧ����La/Cp
nu       = sqrt(Pr/Ra)*delta_t/(delta_x*delta_x);
k        = nu/Pr;                %����ɢϵ��
ks1      = K1*k;                 %�����������ɢϵ��
%% �������ļ���ͼƬ���������
maxT        = 2900000;              %�ܵĵ�������
tStatistics = 2000;               %�������ļ�
tPlot       = 500;                  %���ͼ��
maxk        = 3;                    %�¶ȵĵ���������
%taoֵ�ļ��㣨��ԥʱ�䣩
omegaNS  = 1./(3*nu+1.5);     
omegaT   = 1./(3*k+0.5);
omegaT1  = 1./(3*ks1+0.5);

%% %%%%%%%%%%%%%%%D2Q9 LATTICE CONSTANTS�ٶȺ��¶ȶ����õ���D2Q9��ģ��
%�ٶ���ɢ
tNS   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cxNS  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cyNS  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
oppNS = [  1,   4,  5,  2,  3,    8,   9,   6,   7];
%�¶���ɢ
tT    = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
cxT   = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cyT   = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
oppT  = [  1,   4,  5,  2,  3,    8,   9,   6,   7];
%% %%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = meshgrid(1:n,1:m);

%% %%%%%%%%%%%%%%%%%%%�����ĳ�ʼ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho  = ones(1,m,n);                  
ux   = zeros(1,m,n);
uy   = zeros(1,m,n);
T    = zeros(1,m,n);
rfl  = zeros(1,m,n);                %Һ��ȳ��ĳ�ʼ�� 

%{
%ux  = importdata('E:\LBMѧϰ\Matlab\��ҵ���ļ���\��״��REV\��������\Ra1e6\Tms0.1Tml0.1\Da1e-6\36000-velocity-x.txt');
%ux  =reshape(ux,1,m,n);
%uy  = importdata('E:\LBMѧϰ\Matlab\��ҵ���ļ���\��״��REV\��������\Ra1e6\Tms0.1Tml0.1\Da1e-6\36000-velocity-y.txt');
%uy  =reshape(uy,1,m,n);
%T  = importdata('E:\LBMѧϰ\Matlab\��ҵ���ļ���\��״��REV\��������\Ra1e6\Tms0.1Tml0.1\Da1e-6\36000-temperature.txt');
%T  =reshape(T,1,m,n);
%}
T(:,1,:)=1;                         %���²��¶�
T(:,m,:)=0;                         %���²��¶�
rfl(:,1,:)=1;                       %���²��¶�
%ƽ��̬�ٶ��ݻ�����
for i=1:9
    cuNS         = 3*(cxNS(i)*ux+cyNS(i)*uy);
    fEq(i,:,:)   = rho .* tNS(i) .* ...
                   ( 1 + cuNS + 1/2*(cuNS.*cuNS) - 3/2*(ux.^2+uy.^2) );
    fIn(i,:,:)   = fEq(i,:,:);
end

%ƽ��̬�¶��ݻ�����
 for i=1:9
       cu          = 3*(cxT(i)*ux+cyT(i)*uy);
       tEq(i,:,:)  = T .* tT(i) .* ...
                   ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2)  );          %��ģ����õ��¶�LB��ƽ��̬�ֲ���������p50��3-49������ total=1
       tIn(i,:,:)  = tEq(i,:,:);
 end

 omegaTm     = ones(1,m,n);            %
 omegaTm     = omegaT.*omegaTm;        %
 Fetia       = zeros(1,m,n);
 c0          = zeros(1,m,n);           % ��������c1��c2������p47  3-26
 c1          = zeros(1,m,n);      
 Foc1        = zeros(1,m,n);
 Foc2        = zeros(1,m,n);
 
 %% %%�ೡλ�õ�ȷ��%%%%%%%%%%%
solid  = find(rfl<0.0002);                     %PCM�Ĺ�������
mush   = find(0.0002<=rfl & rfl<1);            %PCM�ĺ�״������
liquid = find(1<=rfl & rfl<2);                 %PCM����������
%% %%%%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 tic
 
for cycle = 1:maxT
    Nu1     =  Nu1+1;
    
    T       =  sum(tIn);
    rflt    =  rfl;                               %��һ��ʱ�䲽���ĺ�Һ��
    rho     =  sum(fIn);   
    Tk      =  T;
    
    ux(solid)=0;   uy(solid)=0;
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
    while 1  
        kl=kl+1;
        %Tb=T;
      for i=1:9
            cu          = 3*(cxT(i)*ux+cyT(i)*uy);
            tEq(i,:,:)  = Tk .* tT(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2)  );
            Sr(i,:,:)   = -LaCp.*tT(i).*reshape((rfl-rflt),1,m,n).*(1+(1-omegaT/2.0)*cu); 
            tOut(i,:,:)  = tIn(i,:,:) - omegaT .* (tIn(i,:,:)-tEq(i,:,:))+Sr(i,:,:);
     end
    
     omegaTm(mush)  = 1./(3*(k.*rfl(mush)+ks1.*(1-rfl(mush)))+0.5);      %(k.*rfl(mush)+ks1.*(1-rfl(mush))Ϊ��Ч�� a
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
        %�����棩ȷ���µ�Һ���
        En  = St*T + reshape(rflk,1,m,n);      %�ʵļ���ʽ�����ӵ�λ�£�
        rfl(En<Ens) =0;                        %�µ�Һ���Ϊ0��PCM��̬��
        Earea=find(Ens<=En & En<=Enl);         %ѡ���µ�Һ���Ϊ����0��С��1�Ĳ���Earea
        rfl(Earea)=(En(Earea)-Ens)/(Enl-Ens);  %�������״����Һ���
        rfl(Enl<En)=1;                         %�µ�Һ���Ϊ1��PCMҺ̬��
        
        solid  = find(rfl<0.0002);                     %PCM�Ĺ�������
        mush   = find(0.0002<=rfl & rfl<1);            %PCM�ĺ�״������
        liquid = find(1<=rfl & rfl<2);                 %PCM����������

                %km=sum(sum(T-Tb))/sum(sum(T));
        
       if (kl>maxk)
           break;
       end
    end
    
%% %%%%%%%%%%%%%%%%%%�ڶ�׽��ʲ��ֵ�����������㣨���������Rev�߶�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
 ȫ�����������ö�׽����ڲ��������̣����û�ж�׽���ʱ����Ӧ��Da��ȡ
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
      ��������
     for i=1:9
      cuNS(liquid)         = 3*(cxNS(i)*ux(liquid)+cyNS(i)*uy(liquid));
      fEq(i,liquid)   = rho(liquid) .* tNS(i) .* ...
                       ( 1 + cuNS(liquid) + 1./(2*rfl(liquid)).*(cuNS(liquid).*cuNS(liquid)) - 3./(2*rfl(liquid)).*(ux(liquid).^2+uy(liquid).^2) );
      force(i,liquid) = (1-omegaNS/2.0).*tNS(i).*rho(liquid).*...
                     ((3.*cyNS(i)-(3./rfl(liquid)).*uy(liquid)+(9./rfl(liquid)).*cxNS(i).*cyNS(i).*ux(liquid)+(9./rfl(liquid)).*cyNS(i).*cyNS(i).*uy(liquid)).*Foc2(liquid)+...
                      (3.*cxNS(i)-(3./rfl(liquid)).*ux(liquid)+(9./rfl(liquid)).*cxNS(i).*cxNS(i).*ux(liquid)+(9./rfl(liquid)).*cxNS(i).*cyNS(i).*uy(liquid)).*Foc1(liquid));
      fOut(i,liquid)  = fIn(i,liquid) - omegaNS .* (fIn(i,liquid)-fEq(i,liquid)) + force(i,liquid);
     end  
      %��״����
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
     %�����Ϲ�������
     for i=1:9
         fEq(i,solid)   = rho(solid) .* tNS(i);
         fOut(i,solid)  = fIn(i,solid) - omegaNS .* (fIn(i,solid)-fEq(i,solid));
     end
    
     for i=1:9
       fOut(i,solid)  =  fIn(oppNS(i),solid);
     end
     fOut  = reshape(fOut,9,m,n);
    %% %%%%%%%%%%%%%%%%%%%%%%%%Ǩ     ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %% %%%%%%%%%%%%%%%%%%%%% ��  ��  ��  �� %%%%%%%%%%%%%%%%%%%%%%%%
  %���ȱ߽�����
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
  
   %% %%%%%%%%%%%%%%%%%%%�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   % VISUALIZATION����ļ������
   if(mod(cycle,tPlot)==0)
       Nu1                                      %��������
    %% �ٶȡ��¶ȼ���
       T =  sum(tIn);
       ux(solid)=0;   uy(solid)=0;
       ux1 = reshape(ux,m,n);     uy1 = reshape(uy,m,n); 
       u1 =sqrt(ux1.^2+uy1.^2);
       T1  = reshape(T,m,n); 
       rfl1=reshape(rfl,m,n);
       
     %% ͼ�����
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
    %% �ٶ��¶����ݵ����
        if (mod(cycle,tStatistics)==0)
            
       %�ٶ�u���ļ����
       stru = [num2str(Nu1),'-velocity.txt'];
       dlmwrite(stru,u1)
       strux = [num2str(Nu1),'-velocity-x.txt'];
       dlmwrite(strux,ux1)
       struy = [num2str(Nu1),'-velocity-y.txt'];
       dlmwrite(struy,uy1)
       %�¶ȵ��ļ������
       strT = [num2str(Nu1),'-temperature.txt'];
       dlmwrite(strT,T1)
       %��Һ�ʵ����
       strr = [num2str(Nu1),'-rfl.txt'];
       dlmwrite(strr,rfl1)
      %ͼƬ�����
       
       end
       toc
   end
   
  end
