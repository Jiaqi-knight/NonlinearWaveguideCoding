%{
ʱ��2016.06.20
�޸ĳ����⣬���ڱ�����������޸�Ϊ3d�ĳ���
��ʽ(lx:lz:ly)
                     |  /ly��Y��
                     | /
������ʾ      ------  |-------lz��H��
                     |
                     |lx��X��
%}
clear
global tNS cxNS cyNS czNS tT  cxT cyT czT m n p
%% ����������ò���
uMax         = 0.1;
%���Ըı�Ĳ���������ģ�ͣ�
lx           = 50;                     %x���������
aspect_ratio = 1;
lz           = aspect_ratio*lx;        %�߶ȷ��������
ly           = aspect_ratio*lx;        %��ֱ���﷽�������
delta_x      = 1./(lx);
Thot         = 1;                      % Heating on bottom wall
Tcold        = 0;                      % Cooling on top wall
T0           = (Thot+Tcold)/2;         %�ο��¶�
%������
delta_t      = uMax/lx;
gr           =delta_t*delta_t/delta_x;
buoyancy     = [0,gr,0];

Nu1=0;                                 %��������
maxT         = 10000000;               %����������
tPlot        = 2;                     %��ͼ���
tStatistics  = 100;                   %��������ļ��
minT         = 1000;                   %��С��������
m=lx;   n=lz;  p=ly;
%���������������������޸ģ�
Pr           = 2;                     %Prandtl number
Ra           = 1e4;                     %Rayleigh number
%����Ĳ��ֲ���
nu           = sqrt(Pr/Ra)*delta_t/(delta_x*delta_x);
k            = sqrt(1./(Pr*Ra))*delta_t/(delta_x*delta_x);
%ԥ��ʱ��
omegaNS      = 1./(3.0*nu+0.5); 
omegaT       = 1./(3.0*k+0.5);
%�����о�
erro         = 1e-10;  
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
[z,x,y] = meshgrid(1:ly,1:lx,1:lz);
%% ��ʼ����������
rho    = ones(1,m,n,p);
%{1
ux     = zeros(1,m,n,p);    
uz     = zeros(1,m,n,p);
uy     = zeros(1,m,n,p);

deltax = -1.0/(m-1); 
T1chu  = 1:deltax:0; 
T2chu  = zeros(m,n,p); 
for i=1:m 
    T2chu(i,:,:) = T1chu(i); 
 end;
T      = reshape(T2chu,1,m,n,p);                                            %��ʽx����*z����*y����
%}
%% ���½��ϵ������ݶ�����¼���
%{
ux  = importdata('E:\LBMsimulation\3DSimulation\TheNumbericalSimulation\3DNaturalConvection\200-velocity-x.txt');
ux  =reshape(ux,1,m,n,p);
uy  = importdata('E:\LBMsimulation\3DSimulation\TheNumbericalSimulation\3DNaturalConvection\200-velocity-y.txt');
uy  =reshape(uy,1,m,n,p);
uz  = importdata('E:\LBMsimulation\3DSimulation\TheNumbericalSimulation\3DNaturalConvection\200-velocity-z.txt');
uz  =reshape(uz,1,m,n,p);
T  = importdata('E:\LBMsimulation\3DSimulation\TheNumbericalSimulation\3DNaturalConvection\200-temperature.txt');
T  =reshape(T,1,m,n,p);
%}
%% 
T(:,1,:,:)=1;
T(:,m,:,:)=0;

for i=1:19
    cuNS         = 3*(cxNS(i)*ux+cyNS(i)*uy+czNS(i)*uz);
    fEq(i,:,:,:)   = rho .* tNS(i) .* ...
                   ( 1 + cuNS + 1/2*(cuNS.*cuNS) - 3/2*(ux.^2+uy.^2+uz.^2));
    fIn(i,:,:,:)   = fEq(i,:,:,:);
end
 for i=1:19
       cu          = 3*(cxT(i)*ux+cyT(i)*uy+czT(i)*uz);
       tEq(i,:,:,:)  = T .* tT(i) .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2+uz.^2));
       tIn(i,:,:,:)  = tEq(i,:,:,:);
 end

tic
%% MAIN LOOP (TIME CYCLES)
for cycle = 1:maxT
   Nu1=Nu1+1;
   
   rho = sum(fIn);                                                                  %(1,m,n,p)
   T   = sum(tIn);                                                                  %(1,m,n,p)
  
   ux1    = reshape (ux,m,n,p);    uy1    = reshape (uy,m,n,p);                     %����һ��ʱ�䲽����ux��uy,uz��¼��ux1��uy1,uz1��
   uz1    = reshape (uz,m,n,p); 
   
   ux     = reshape ( (cxNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho;
   uy     = reshape ( (cyNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho;
   uz     = reshape ( (czNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho+ 0.5*gr.*(T-T0);   %
 
   ux(:,[1 m],:,:)=0;      uy(:,[1 m],:,:)=0;       uz(:,[1 m],:,:)=0; 
   ux(:,:,[1 n],:)=0;      uy(:,:,[1 n],:)=0 ;      uz(:,:,[1 n],:)=0;
   ux(:,:,:,[1 p])=0;      uy(:,:,:,[1 p])=0 ;      uz(:,:,:,[1 p])=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%��ײ---Ǩ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %�ٶ�
   for i=1:19
      cuNS           = 3*(cxNS(i)*ux+cyNS(i)*uy+czNS(i)*uz);
      fEq(i,:,:,:)   = rho .* tNS(i) .* ...
                       ( 1 + cuNS + 1/2*(cuNS.*cuNS) - 3/2*(ux.^2+uy.^2+uz.^2) );
    % force(i,:,:) = gr*(1-omegaNS/2.0).*tNS(i).*rho.*(3.*(cyNS(i)-uy)+9.*cxNS(i).*cyNS(i).*ux+9.*cyNS(i).*cyNS(i).*uy).*(T-T0);
      force(i,:,:,:) = (3*(1-omegaNS/2.0)).*tNS(i).*rho.*(((czNS(i)-uz)+3.*cuNS*czNS(i)).*gr.*(T-T0)); 
      fOut(i,:,:,:)  = fIn(i,:,:,:) - omegaNS .* (fIn(i,:,:,:)-fEq(i,:,:,:)) + force(i,:,:,:);
   end
   
  %�¶�
   for i=1:19
      cu            = 3*(cxT(i)*ux+cyT(i)*uy+czT(i)*uz);
      tEq(i,:,:,:)  = T .* tT(i) .*...
                        ( 1 + cu + 1/2*(cu.*cu) - 3/2*(ux.^2+uy.^2+uz.^2));
      tOut(i,:,:,:) = tIn(i,:,:,:) - omegaT .* (tIn(i,:,:,:)-tEq(i,:,:,:));
   end
   
   %%Ǩ�Ʋ���
   for i=1:19
      fIn(i,:,:,:) = circshift(fOut(i,:,:,:), [0,cxNS(i),czNS(i),cyNS(i)]);
   end
   
   for i=1:19
      tIn(i,:,:,:) = circshift(tOut(i,:,:,:), [0,cxT(i),czT(i),cyT(i)]); 
   end
   
%% %%�߽���������(����)��ƽ�����Ƹ�ʽ
   %�¶ȱ߽�
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
    
  %�ٶȱ߽�     (����������ʽ)
   for i=1:19
        fIn(i,:,1,:)  = fIn(oppNS(i),:,1,:);
        fIn(i,:,n,:)  = fIn(oppNS(i),:,n,:);                                                                                               
        fIn(i,1,:,:)  = fIn(oppNS(i),1,:,:);
        fIn(i,m,:,:)  = fIn(oppNS(i),m,:,:);
        fIn(i,:,:,1)  = fIn(oppNS(i),:,:,1);
        fIn(i,:,:,p)  = fIn(oppNS(i),:,:,p);
   end
   
  %��߽�
   tIn(3,m,:,:) = Tcold -tIn(1,m,:,:)-tIn(2,m,:,:)-tIn(4,m,:,:)-tIn(5,m,:,:)-tIn(6,m,:,:)-tIn(7,m,:,:)-tIn(8,m,:,:)-tIn(9,m,:,:)...
         -tIn(10,m,:,:) -tIn(11,m,:,:)-tIn(12,m,:,:) -tIn(13,m,:,:)-tIn(14,m,:,:)-tIn(15,m,:,:)-tIn(16,m,:,:)-tIn(17,m,:,:)-tIn(18,m,:,:)-tIn(19,m,:,:);
  %�ұ߽�
   tIn(2,1,:,:)  = Thot -tIn(1,1,:,:)-tIn(3,1,:,:)-tIn(4,1,:,:)-tIn(5,1,:,:)-tIn(6,1,:,:)-tIn(7,1,:,:)-tIn(8,1,:,:)-tIn(9,1,:,:)...
         -tIn(10,1,:,:) -tIn(11,1,:,:)-tIn(12,1,:,:) -tIn(13,1,:,:)-tIn(14,1,:,:)-tIn(15,1,:,:)-tIn(16,1,:,:)-tIn(17,1,:,:)-tIn(18,1,:,:)-tIn(19,1,:,:);
  
  %% ͼƬ���������
    % VISUALIZATION
   if(mod(cycle,tPlot)==0)   
          Nu1                                                                                        %ʵʱ�������
   %%%
   ux     = reshape ( (cxNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho;
   uy     = reshape ( (cyNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho;
   uz     = reshape ( (czNS * reshape(fIn,19,m*n*p)), 1,m,n,p) ./rho+ 0.5*gr.*(T-T0);   %
 
   ux(:,[1 m],:,:)=0;      uy(:,[1 m],:,:)=0;       uz(:,[1 m],:,:)=0; 
   ux(:,:,[1 n],:)=0;      uy(:,:,[1 n],:)=0 ;      uz(:,:,[1 n],:)=0;
   ux(:,:,:,[1 p])=0;      uy(:,:,:,[1 p])=0 ;      uz(:,:,:,[1 p])=0;
   
   T    = sum(tIn); 
   T(:,1,:,:)=1;    T(:,m,:,:)=0;
   T1    = reshape(T,m,n,p); 
   ux2   = reshape (ux,m,n,p);    uy2    = reshape (uy,m,n,p);    uz2   = reshape (uz,m,n,p);
   u1    = reshape(sqrt(ux2.^2+uy2.^2+uz2.^2),m,n,p); 
  %������
   I1    =sum(sum(sum((ux2(:,:,:)-ux1(:,:,:)).^2+(uy2(:,:,:)-uy1(:,:,:)).^2+(uz2(:,:,:)-uz1(:,:,:)).^2)));
   I2    =sum(sum(sum(ux2(:,:,:).^2+uy2(:,:,:).^2+uz2(:,:,:).^2)));                                      
         error=sqrt(I1)/(sqrt(I2)+1e-30); 
         error                                                                                        %���
     %MATLAB��ƬͼƬ�����        
         clf; 
           subplot(2,1,1); 
           imagesc(u1(:,n:-1:1,p/2)');    colorbar                                                      %��Ƭ�ٶȵ����
           title('Fluid velocity');
           axis image; drawnow
           subplot(2,1,2);
           imagesc(T1(:,n:-1:1,p/2)');       colorbar                                                   %��Ƭ�¶ȵ����
           title(['Temperature (The number of iterations ' num2str(Nu1) '-' ')']);
           axis image; drawnow
     %MATLAB���ݵ���    
       if (mod(cycle,tStatistics)==0)  
      [Lx,Lz,Ly]=size(T1);  
       deltx=1./(Lx-1);
       deltz=1./(Lz-1);
       delty=1./(Ly-1);
       m1=0:deltx:1;
       n1=0:deltz:1;
       p1=0:delty:1;
      str=   [num2str(Nu1),'-','data.dat'];
       fid=fopen(str,'w'); 
       fprintf(fid,'title="3D-quad" Variables = "X" "Z" "Y" "T" "U" "Ux" "Uz" "Uy"  Zone I=%g,J=%g,K=%g,F=Point\n',Lx,Lz,Ly);
       %fprintf(fid,'title="2D-quad" Variables = "X" "Y" "T" "f" Zone I=%g,J=%g,F=Point\n',ly-3,lx);
       for i1=1:Lx  
          for j1=1:Lz
              for v1=1:Ly
            fprintf(fid,'%g %g %g %8.4f %8.4f %8.4f %8.4f %8.4f\n',m1(i1),n1(j1),p1(v1),T1(i1,j1,v1),u1(i1,j1,v1),ux2(i1,j1,v1),uz2(i1,j1,v1),uy2(i1,j1,v1));
             end
          end  
       end 
       fclose(fid);
        
      %txt�ļ��ĵ��������Է�������ɵļ����жϣ�
      ux3   = reshape (ux2,m,n*p);      uz3   = reshape (uz2,m,n*p);  uy3    = reshape (uy2,m,n*p); 
      u3    = reshape (u1,m,n*p);                                                                   %�ٶ�
      T3    = reshape (T1,m,n*p);                                                                   %�¶�
       %�ٶ�u���ļ����
      stru = [num2str(Nu1),'-velocity.txt'];
       dlmwrite(stru,u3)
      strux = [num2str(Nu1),'-velocity-x.txt'];
       dlmwrite(strux,ux3)
      struy = [num2str(Nu1),'-velocity-y.txt'];
       dlmwrite(struy,uy3)
      struz = [num2str(Nu1),'-velocity-z.txt'];
       dlmwrite(struz,uz3)
       %�¶ȵ��ļ������
       strT = [num2str(Nu1),'-temperature.txt'];
       dlmwrite(strT,T3)
       
   if(error<erro && Nu1>=minT && Nu1<=maxT)                                               %�ﵽ�����о�Ҫ�����ٴ���������������ǵ�����
      Nu1                                                                                %ʵʱ�������
     error
      disp('now It is over')
     str=   [num2str(Nu1),'-','data.dat'];
       fid=fopen(str,'w'); 
       fprintf(fid,'title="3D-quad" Variables = "X" "Z" "Y" "T" "U" "Ux" "Uz" "Uy"  Zone I=%g,J=%g,K=%g,F=Point\n',Lx,Lz,Ly);
       %fprintf(fid,'title="2D-quad" Variables = "X" "Y" "T" "f" Zone I=%g,J=%g,F=Point\n',ly-3,lx);
       for i1=1:Lx  
          for j1=1:Lz
              for v1=1:Ly
            fprintf(fid,'%g %g %g %8.4f %8.4f %8.4f %8.4f %8.4f\n',m1(i1),n1(j1),p1(v1),T1(i1,j1,v1),u1(i1,j1,v1),ux2(i1,j1,v1),uz2(i1,j1,v1),uy2(i1,j1,v1));
             end
          end  
       end 
       fclose(fid);
      
       break
           return;
   end
   
  end
  toc
 end
   
end 
   
%fclose(fid);
