function [Surf_para]=MATLAB4geomTurbo
h=figure;
hold on; axis equal;view([180 -30]);

% plot3(xr*1000,yr*1000,zr*1000-170,'.');hold on
% plot3(xh*1000,yh*1000,zh*1000-170,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
%     'MarkerSize',4,...
%     'Marker','diamond',...
%     'LineStyle','none');axis equal
% plot3(xr1*1000,yr1*1000,zr*1000-170,'.');axis equal
IGV_stagger=0;IGV_lean=0;IGV_sweep=0;
rotor_stagger=45/180*pi;rotor_lean=-15/180*pi;rotor_sweep=30/180*pi;%测试估计
stator_stagger=-30/180*pi;stator_lean=0/180*pi;stator_sweep=0/180*pi;
%[fname,location]=uigetfile({'*.geomTurbo';'*.txt';'*.dat';'*.*'},'r');%选择文件
fname='DTS-02.geomTurbo';
location='E:\Jiaqi-SJTU-DOIT\Maincode\GITHUB-wjq\GreenFunction4Beforming_With_InitialEigenValue'

Data=ReadgeomTurbo(fname)   %读取文件（读取重要参数信息，读取参数data，将data转成cell格式，并存储输出）
Data.IGV_suction_rotate=rotateIGV(Data.IGV_suction,[0,1,0],IGV_stagger,[0,0,-104.015]); % rotation the IGV
Data.IGV_pressure_rotate=rotateIGV(Data.IGV_pressure,[0,1,0],IGV_stagger,[0,0,-104.015]);%%右手定则，顺时针旋转，基准为0度



grid_index=10;
[Y_IGV,Z_IGV] = meshgrid(linspace(Data.IGV_suction_rotate{1, 1}(1,2),Data.IGV_suction_rotate{1, end}(1,2),grid_index),...
    linspace(Data.IGV_suction_rotate{1,1}(1,3),Data.IGV_suction_rotate{1,1}(end,3),grid_index));
X_IGV = zeros(size(Y_IGV));
plane_IGV=surf(X_IGV,Y_IGV,Z_IGV);
[Y_rotor,Z_rotor] = meshgrid(linspace(Data.rotor_suction{1, 1}(1,2),Data.rotor_suction{1, end}(1,2),grid_index),...
    linspace(Data.rotor_suction{1,10}(1,3)-3,Data.rotor_suction{1,10}(1,3)+60,grid_index));
X_rotor = zeros(size(Y_rotor));
plane_rotor=surf(X_rotor,Y_rotor,Z_rotor);
[Y_stator,Z_stator] = meshgrid(linspace(Data.stator_suction{1, 1}(1,2),Data.stator_suction{1, end}(1,2),grid_index),...
    linspace(Data.stator_suction{1,6}(1,3),Data.stator_suction{1,6}(end,3)+2,grid_index));
X_stator = zeros(size(Y_stator))+Data.stator_suction{1,6}(1,1);
plane_stator=surf(X_stator,Y_stator,Z_stator);
set(plane_IGV,'Visible','on');set(plane_rotor,'Visible','on');set(plane_stator,'Visible','on');

Arc_s1=[828.893086090 452.272727273 0]';Arc_e1=[817.5 420.5 0]';O1=[867.5 420.5 0]';
Arc_s2=[1037.370421204 550.7 0]';Arc_e2=[828.893086090 452.272727273 0]';O2=[1037.370421204 280.7 0]';
L1=[1048.5 550.7 0]';L2=[1037.37 550.7 0]';
X1=arcPlot(Arc_s1,Arc_e1,O1);X2=arcPlot(Arc_s2,Arc_e2,O2);
X=[L1 L2 X2(:,1:end-3) X1(:,2:end-1)].';

%plot3(X(1,:),X(2,:),X(3,:));

Surf_para.bolt={[X(:,3) X(:,2)-min(X(:,2)) X(:,1)-max(X(:,1))+min(Data.hub(:,1))],[0 0 1],[0 0 0],linspace(0,2*pi,37),@surf}
Surf_bolt=rotsurf(Surf_para.bolt{1},Surf_para.bolt{2},Surf_para.bolt{3},Surf_para.bolt{4},Surf_para.bolt{5});
set(Surf_bolt, 'FaceColor',[0.68 0.92 1]);

Wall1=[185 0 -170-1400;185 0 -170];
Surf_para.duct={Wall1,[0 0 1],[0 0 0],linspace(-0.8*pi,0.8*pi,37)-1/2*pi,@surf};
Surf_wall1=rotsurf(Surf_para.duct{1},Surf_para.duct{2},Surf_para.duct{3},Surf_para.duct{4},Surf_para.duct{5});

Wall2=[sin(linspace(0+1*pi,2/4*pi+1*pi,20))'*120+305 zeros(20,1) cos(linspace(0+1*pi,2/4*pi+1*pi,20))'*120-1400-170];
Surf_para.horn={Wall2,[0 0 1],[0 0 0],linspace(-0.8*pi,0.8*pi,37)-1/2*pi,@surf};
Surf_wall2=rotsurf(Surf_para.horn{1},Surf_para.horn{2},Surf_para.horn{3},Surf_para.horn{4},Surf_para.horn{5});
set(Surf_wall1, 'FaceColor',[0.7 0.92 1]);
set(Surf_wall2, 'FaceColor',[0.7 0.92 1]);




%set(Surf_wall1, 'FaceColor',[0.68 0.92 1]);
%set(Surf_wall2, 'FaceColor',[0.68 0.92 1]);
Surf_para.shround={[Data.shround(:,2) zeros(size(Data.shround,1),1) Data.shround(:,1)],[0 0 1],[0 0 0],linspace(-0.8*pi,0.8*pi,37)-1/2*pi,@surf};
Surf_shround=rotsurf(Surf_para.shround{1},Surf_para.shround{2},Surf_para.shround{3},Surf_para.shround{4},Surf_para.shround{5});
Surf_para.hub={[Data.hub(:,2) zeros(size(Data.hub,1),1) Data.hub(:,1)],[0 0 1],[0 0 0],linspace(-1*pi,1*pi,37)-1/2*pi,@surf};
Surf_hub=rotsurf(Surf_para.hub{1},Surf_para.hub{2},Surf_para.hub{3},Surf_para.hub{4},Surf_para.hub{5});
set(Surf_shround, 'FaceColor',[0.93 0.70 0.13]);
set(Surf_hub,'LineStyle','-.',...
    'FaceColor',[0.16 0.38 0.27],...
    'EdgeColor',[0.86 0.86 0.86]);

[X11]=yexing2surf(Data.IGV_suction_rotate);
[X12]=yexing2surf(Data.IGV_pressure_rotate);
IGV(1)=surf(X11{1},X11{2},X11{3},'EdgeColor',[1 0 0]);hold on
IGV(2)=surf(X12{1},X12{2},X12{3},'EdgeColor',[1 0 0]);
Rotor{1}=plotCellData(Data.rotor_suction,'r');hold on;
Rotor{2}=plotCellData(Data.rotor_pressure,'r');hold on;
Stator{1}=plotCellData(Data.stator_suction,'r');hold on;
Stator{2}=plotCellData(Data.stator_pressure,'r');hold on;
set(IGV,'Visible','on');
set(Rotor{1},'Visible','on');set(Rotor{2},'Visible','on')
set(Stator{1},'Visible','on');set(Stator{2},'Visible','on');

%rotate_blade(IGV,[0 0 1],32,'surf'); 
%rotate_blade(Rotor{1},[0 0 1],29,'plot'); rotate_blade(Rotor{2},[0 0 1],29,'plot'); 
%rotate_blade(Stator{1},[0 0 1],37,'plot'); rotate_blade(Stator{2},[0 0 1],37,'plot'); 

%绝对坐标系
% [para.XAXIS(1) para.XAXIS(2) para.XAXIS(3)]=arrow3d( [0 0 0.5]',3,90,8,'x'); 
% [para.YAXIS(1) para.YAXIS(2) para.YAXIS(3)]=arrow3d( [0 0 0.5]',3,90,8,'y');
% [para.ZAXIS(1) para.ZAXIS(2) para.ZAXIS(3)]=arrow3d( [0 0 0.5]',3,90,8,'z');

% para.XAXIS(4)=text(110,0,0,'X','Color','r','FontSize',12);
% para.YAXIS(4)=text(0,110,0,'Y','Color','r','FontSize',12);
para.ZAXIS(4)=text(0,10,120,'Z','Color','r','FontSize',12);
%绝对相对坐标系
origin_IGV=[0 Data.IGV_suction_rotate{1, 1}(1,2:3)];
origin_rotor=[0 Data.rotor_suction{1, 1}(1,2) Data.rotor_suction{1,10}(1,3)-3];
origin_stator=[Data.stator_suction{1, 1}(1,1) Data.stator_suction{1, 1}(1,2) Data.stator_suction{1,6}(1,3)];
% 
% [para.XAXIS_IGV]=arrow3d_relative(origin_IGV',30,8,1.5); 
% [para.XAXIS_rotor]=arrow3d_relative(origin_rotor',30,8,1.5);
% [para.XAXIS_stator]=arrow3d_relative(origin_stator',30,8,1.5);

IGV_new=surf_newxyz(X_IGV,Y_IGV,Z_IGV,origin_IGV,IGV_stagger,IGV_lean,IGV_sweep,grid_index);
rotor_new1=surf_newxyz(X_rotor,Y_rotor,Z_rotor,origin_rotor,rotor_stagger,0,0,grid_index);
rotor_new2=surf_newxyz(X_rotor,Y_rotor,Z_rotor,origin_rotor,rotor_stagger,rotor_lean,0,grid_index);
rotor_new2=surf_newxyz(X_rotor,Y_rotor,Z_rotor,origin_rotor,rotor_stagger,rotor_lean,rotor_sweep,grid_index);
stator_new=surf_newxyz(X_stator,Y_stator,Z_stator,origin_stator,stator_stagger,stator_lean,stator_sweep,grid_index);


end
