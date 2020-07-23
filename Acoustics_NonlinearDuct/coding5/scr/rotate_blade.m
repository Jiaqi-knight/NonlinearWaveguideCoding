function h=rotate_blade(h,azel,N,type)
%h=rotate(h,azel,alpha,origin)
%在函数rotate的基础上改进，专门适用于叶片旋转

% Determine the default origin (center of plot box).
  if ~ishghandle(h)
    error(message('MATLAB:rotate:InvalidHandle'));
  end
  ax = ancestor(h(1),'axes');
  if isempty(ax) || ax==0,
    error(message('MATLAB:rotate:InvalidHandle'));
  end
  matlab.ui.internal.UnsupportedInUifigure(ancestor(ax,'figure'));
  origin = sum([get(ax,'xlim')' get(ax,'ylim')' get(ax,'zlim')'])/2;

% find unit vector for axis of rotation
if numel(azel) == 2 % theta, phi
    theta = pi*azel(1)/180;
    phi = pi*azel(2)/180;
    u = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
elseif numel(azel) == 3 % direction vector
    u = azel(:)/norm(azel);
end

for k=1:N-1
alph = 2*pi/N*k;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

for i=1:numel(h),
  t = get(h(i),'type');
  skip = 0;
  if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'patch')
    
    % If patch, rotate vertices  
    if strcmp(t,'patch')
       verts = get(h(i),'Vertices');
       x = verts(:,1); y = verts(:,2); 
       if size(verts,2)>2
          z = verts(:,3);
       else
          z = [];
       end
       
    % If surface or line, rotate {x,y,z}data   
    else
       x = get(h(i),'xdata');
       y = get(h(i),'ydata');
       z = get(h(i),'zdata');
    end
    
    if isempty(z)
       z = -origin(3)*ones(size(y));
    end
    [m,n] = size(z);
    if numel(x) < m*n
      [x,y] = meshgrid(x,y);
    end
  elseif strcmp(t,'text')
    p = get(h(i),'position');
    x = p(1); y = p(2); z = p(3);
  elseif strcmp(t,'image')
    x = get(h(i),'xdata');
    y = get(h(i),'ydata');
    z = zeros(size(x));
  else
    skip = 1;
  end
  
  if ~skip,
    [m,n] = size(x);
    newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
    newxyz = newxyz*rot;
    newx = origin(1) + reshape(newxyz(:,1),m,n);
    newy = origin(2) + reshape(newxyz(:,2),m,n);
    newz = origin(3) + reshape(newxyz(:,3),m,n);
    
    if strcmp(type,'surf') 
        surf(newx,newy,newz,'EdgeColor',[1 1 1]);hold on
        elseif strcmp(type,'plot')
        plot3(newx,newy,newz,'Color',[1 1 1]);hold on
    end
        
        

%     if strcmp(t,'surface') || strcmp(t,'line')
%       surf(newx,newy,newz,'EdgeColor',[1 0 0]);hold on
% %     set(h(i),'xdata',newx,'ydata',newy,'zdata',newz);
%     elseif strcmp(t,'patch')
%       set(h(i),'Vertices',[newx,newy,newz]); 
%     elseif strcmp(t,'text')
%       set(h(i),'position',[newx newy newz])
%     elseif strcmp(t,'image')
%       set(h(i),'xdata',newx,'ydata',newy)
%     end
  end
end
end
end
