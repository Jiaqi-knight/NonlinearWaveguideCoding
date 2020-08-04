function plotfunction(h, rho, u, t)
% This function is designed to plot the field 
% including: 
% 1. Contour plot of velocity field
% 2. Streamlines plot of velocity field
% 3. Contour plot of pressure field
% inputs rho(Nx+1,Ny+1), u(Nx+1,Ny+1,2),
global Dimension Nx Ny delta_x delta_y delta_t c Re
%% Compute filed parameters
% position
i = 1:Nx+1;
j = 1:Ny+1;
x = i*delta_x;
y = j*delta_y;
[X,Y] = meshgrid(x,y);
startx = 1:10:Nx+1;
starty = Ny+1*ones(1,length(startx));
% time
T = t*delta_t;
% velocity field
velocity_magnitude = zeros(Nx+1,Ny+1);
for i = 1:Nx+1
    for j = 1:Ny+1
        for m = 1:Dimension
            velocity_magnitude(i,j) = velocity_magnitude(i,j) + u(i,j,m)^2;
        end
    end
end 
% vorticity magnitude field
% vorticity_magnitude = curl(X,Y,u(:,:,1),u(:,:,2));
% pressure field
p = 1./3. *c^2 .*rho ;

%% plot 
% velocity field
figure(h(1));
% subplot(2,1,1);
imagesc(rot90(velocity_magnitude));
title(['Re = ',num2str(Re),',T = ',num2str(T)]);
colormap jet;
colorbar('location','eastoutside')
drawnow;
% subplot(2,1,2);
% streamline(X,Y,u(:,:,1),u(:,:,2),startx,starty);
% box on;
% title(['Re = ',num2str(Re),',T = ',num2str(T)]);
% drawnow;
% hold on;
% hp = pcolor(X,Y,velocity_magnitude);
% set(hp,'LineStyle','none');
% hold off;
% axis equal;
% colormap jet;
% colorbar('location','southoutside')
% drawnow;
% pressure field
figure(h(2));
% hp = pcolor(X,Y,p);
imagesc(rot90(p));
colormap jet;
colorbar('location','eastoutside')
title(['Re = ',num2str(Re),',T = ',num2str(T)]);
drawnow;
% set(hp,'LineStyle','none');
% hold off;
% axis equal;

end

