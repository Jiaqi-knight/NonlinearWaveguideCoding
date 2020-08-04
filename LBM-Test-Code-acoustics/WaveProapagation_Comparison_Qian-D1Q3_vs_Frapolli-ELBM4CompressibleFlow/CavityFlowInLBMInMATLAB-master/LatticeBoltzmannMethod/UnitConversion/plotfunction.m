function plotfunction(h, rho, ux, uy, t)
% This function is designed to plot the field 
% including: 
% 1. Contour plot of velocity field
% 2. Streamlines plot of velocity field
% 3. Contour plot of pressure field
% inputs rho(Nx+1,Ny+1), u(Nx+1,Ny+1,2),
global  Nx Ny delta_x delta_y delta_t cs_Physics Re_Physics u_r rho_r L_r t_r
%% Compute filed parameters
% position
i = 1:Nx+1;
j = 1:Ny+1;
x = i*delta_x*L_r;
y = j*delta_y*L_r;
[X,Y] = meshgrid(x,y);
% startx = 1:10:Nx+1;
% starty = floor((Ny+1)/2)*ones(1,length(startx));
% time
T = t*delta_t*t_r;
% velocity magnitude
velocity_magnitude = (u_r*ux).^2 + (u_r*uy).^2;
% voticity magnitude
vorticity_magnitude = curl(X,Y,u_r*ux,u_r*uy);
% pressure field
p = cs_Physics^2 .*(rho_r*rho) ;

%% plot 
% velocity magnitude field
figure(h(1));
% subplot(2,1,1);
imagesc(rot90(velocity_magnitude));
title(['Re = ',num2str(Re_Physics),',T = ',num2str(T)]);
colormap jet;
colorbar('location','eastoutside')
drawnow;

% Voticity Magnitude field
figure(h(2));
imagesc(rot90(vorticity_magnitude));
title(['Re = ',num2str(Re_Physics),',T = ',num2str(T)]);
colormap jet;
colorbar('location','eastoutside')
drawnow;

% hlines = streamline(X,Y,ux,uy,startx,starty);
% set(hlines,'LineWidth',2,'Color','r');
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
figure(h(3));
% hp = pcolor(X,Y,p);
imagesc(rot90(p));
colormap jet;
colorbar('location','eastoutside')
title(['Re = ',num2str(Re_Physics),',T = ',num2str(T)]);
drawnow;
% set(hp,'LineStyle','none');
% hold off;
% axis equal;

end

