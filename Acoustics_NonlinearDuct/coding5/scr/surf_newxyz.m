function New=surf_newxyz(X,Y,Z,origin,stagger,lean,sweep,grid_index)
%注意：旋转为绝对相对坐标系，最后继续转化为绝对坐标系
newxyz =([X(:), Y(:), Z(:)]-repmat(origin,length(X(:)),1))*Q_sweep(sweep)*Q_lean(lean)*Q_stagger(stagger)+repmat(origin,length(X(:)),1);%注意矩阵相乘的先后顺序
New=surf(reshape(newxyz(:,1),grid_index,grid_index),reshape(newxyz(:,2),grid_index,grid_index),reshape(newxyz(:,3),grid_index,grid_index),'EdgeColor',[1 1 1]);
end
