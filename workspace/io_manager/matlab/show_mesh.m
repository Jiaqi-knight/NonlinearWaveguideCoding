function[f]=show_mesh(D,mesh,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%   Creation date: March 31st, 2014
%     Last update: June 12th, 2014
%
%    Description: 
%          Input: 
%         Output: 
%          Notes: get(0,'ScreenSize') to get screen size
%%

switch D
    case 2
        switch mesh.D
          case 1
              f = figure('Position', [1, 1, xfigsize, yfigsize]);
              if shownodes
                  for i=1:mesh.totN
                      plot(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),strcat([pointcolor,'.']),'LineWidth',pointdim)
                      hold on
                  end
              end
              if showlines
                  for l=1:mesh.Nlinesdim1
                      plot(mesh.coordlines.dim1free(l,1),mesh.coordlines.dim1free(l,2),strcat([linecolor,'-']),'LineWidth',linedim)
                      hold on
                  end
              end
              if shownodelabels
                  for i=1:mesh.totN
                      text(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),num2str(mesh.nodes.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',nodelabelcolor)
                      hold on
                  end
              end
              if showedgelabels
                  for i=1:mesh.totE
                      text(mesh.edges.centroid(1,i),mesh.edges.centroid(2,i),num2str(mesh.edges.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',edgelabelcolor)
                      hold on
                  end
              end
              grid on
              axis([xaxismin xaxismax yaxismin yaxismax]);
              axis equal
              xlabel(xstring)
              ylabel(ystring)
              title(titlestring)
          case 2
              f = figure('Position', [1, 1, xfigsize, yfigsize]);
              if shownodes
                  for i=1:mesh.totN
                      plot(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),strcat([pointcolor,'.']),'LineWidth',pointdim)
                      hold on
                  end
              end
              if showlines
                  for j=1:mesh.Ndim2
                      for l=1:mesh.Nlinesdim1
                          plot(mesh.coordlines.dim1free(l,3*(j-1)+1),mesh.coordlines.dim1free(l,3*(j-1)+2),strcat([linecolor,'-']),'LineWidth',linedim)
                          hold on
                      end
                  end
                  for k=1:mesh.Ndim1
                      for l=1:mesh.Nlinesdim2
                          plot(mesh.coordlines.dim2free(l,3*(k-1)+1),mesh.coordlines.dim2free(l,3*(k-1)+2),strcat([linecolor,'-']),'LineWidth',linedim)
                          hold on
                      end
                  end
              end
              if shownodelabels
                  for i=1:mesh.totN
                      text(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),num2str(mesh.nodes.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',nodelabelcolor)
                      hold on
                  end
              end
              if showedgelabels
                  for i=1:mesh.totE
                      text(mesh.edges.centroid(1,i),mesh.edges.centroid(2,i),num2str(mesh.edges.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',edgelabelcolor)
                      hold on
                  end
              end
              if showfacelabels
                  for i=1:mesh.totF
                      text(mesh.faces.centroid(1,i),mesh.faces.centroid(2,i),num2str(mesh.faces.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',facelabelcolor)
                      hold on
                  end
              end
              grid on
              axis([xaxismin xaxismax yaxismin yaxismax]);
              axis equal
              xlabel(xstring)
              ylabel(ystring)
              title(titlestring)
          case 3
              disp('Dimensions do not match')
              f = 0;
        end 
    case 3
        switch mesh.D
          case 1
              f = figure('Position', [1, 1, xfigsize, yfigsize]);
              if shownodes
                  for i=1:mesh.totN
                      plot3(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),mesh.nodes.meshcoordinates(3,i),strcat([pointcolor,'.']),'LineWidth',pointdim)
                      hold on
                  end
              end
              if showlines
                  for l=1:mesh.Nlinesdim1
                      plot3(mesh.coordlines.dim1free(l,1),mesh.coordlines.dim1free(l,2),mesh.coordlines.dim1free(l,3),strcat([linecolor,'-']),'LineWidth',linedim)
                      hold on
                  end
              end
              if shownodelabels
                  for i=1:mesh.totN
                      text(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),mesh.nodes.meshcoordinates(3,i),num2str(mesh.nodes.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',nodelabelcolor)
                      hold on
                  end
              end
              if showedgelabels
                  for i=1:mesh.totE
                      text(mesh.edges.centroid(1,i),mesh.edges.centroid(2,i),mesh.edges.centroid(3,i),num2str(mesh.edges.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',edgelabelcolor)
                      hold on
                  end
              end
              grid on
              axis([xaxismin xaxismax yaxismin yaxismax zaxismin zaxismax]);
              axis equal
              xlabel(xstring)
              ylabel(ystring)
              zlabel(zstring)
              title(titlestring)
          case 2
              f = figure('Position', [1, 1, xfigsize, yfigsize]);
              if shownodes
                  for i=1:mesh.totN
                      plot3(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),mesh.nodes.meshcoordinates(3,i),strcat([pointcolor,'.']),'LineWidth',pointdim)
                      hold on
                  end
              end
              if showlines
                  for j=1:mesh.Ndim2
                      for l=1:mesh.Nlinesdim1
                          plot3(mesh.coordlines.dim1free(l,3*(j-1)+1),mesh.coordlines.dim1free(l,3*(j-1)+2),mesh.coordlines.dim1free(l,3*(j-1)+3),strcat([linecolor,'-']),'LineWidth',linedim)
                          hold on
                      end
                  end
                  for k=1:mesh.Ndim1
                      for l=1:mesh.Nlinesdim2
                          plot3(mesh.coordlines.dim2free(l,3*(k-1)+1),mesh.coordlines.dim2free(l,3*(k-1)+2),mesh.coordlines.dim2free(l,3*(k-1)+3),strcat([linecolor,'-']),'LineWidth',linedim)
                          hold on
                      end
                  end
              end
              if shownodelabels
                  for i=1:mesh.totN
                      text(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),mesh.nodes.meshcoordinates(3,i),num2str(mesh.nodes.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',nodelabelcolor)
                      hold on
                  end
              end
              if showedgelabels
                  for i=1:mesh.totE
                      text(mesh.edges.centroid(1,i),mesh.edges.centroid(2,i),mesh.edges.centroid(3,i),num2str(mesh.edges.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',edgelabelcolor)
                      hold on
                  end
              end
              if showfacelabels
                  for i=1:mesh.totF
                      text(mesh.faces.centroid(1,i),mesh.faces.centroid(2,i),mesh.faces.centroid(3,i),num2str(mesh.faces.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',facelabelcolor)
                      hold on
                  end
              end
              grid on
              axis([xaxismin xaxismax yaxismin yaxismax zaxismin zaxismax]);
              axis equal
              xlabel(xstring)
              ylabel(ystring)
              zlabel(zstring)
              title(titlestring)
          case 3
              f = figure('Position', [1, 1, xfigsize, yfigsize]);
              if shownodes
                  for i=1:mesh.totN
                      plot3(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),mesh.nodes.meshcoordinates(3,i),strcat([pointcolor,'.']),'LineWidth',pointdim)
                      hold on
                  end
              end
              if showlines
                  for i=1:mesh.Ndim3
                      for j=1:mesh.Ndim2
                          for l=1:mesh.Nlinesdim1
                              plot3(mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1)+1),mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1)+2),mesh.coordlines.dim1free(l,3*(j-1)+3*mesh.Ndim2*(i-1)+3),strcat([linecolor,'-']),'LineWidth',linedim)
                              hold on
                          end
                      end
                  end
                  for i=1:mesh.Ndim3
                      for k=1:mesh.Ndim1
                          for l=1:mesh.Nlinesdim2
                              plot3(mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1)+1),mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1)+2),mesh.coordlines.dim2free(l,3*(k-1)+3*mesh.Ndim1*(i-1)+3),strcat([linecolor,'-']),'LineWidth',linedim)
                              hold on
                          end
                      end
                  end
                  for j=1:mesh.Ndim2
                      for k=1:mesh.Ndim1
                          for l=1:mesh.Nlinesdim3
                              plot3(mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1)+1),mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1)+2),mesh.coordlines.dim3free(l,3*(k-1)+3*mesh.Ndim1*(j-1)+3),strcat([linecolor,'-']),'LineWidth',linedim)
                              hold on
                          end
                      end
                  end
              end
              if shownodelabels
                  for i=1:mesh.totN
                      text(mesh.nodes.meshcoordinates(1,i),mesh.nodes.meshcoordinates(2,i),mesh.nodes.meshcoordinates(3,i),num2str(mesh.nodes.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',nodelabelcolor)
                      hold on
                  end
              end
              if showedgelabels
                  for i=1:mesh.totE
                      text(mesh.edges.centroid(1,i),mesh.edges.centroid(2,i),mesh.edges.centroid(3,i),num2str(mesh.edges.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',edgelabelcolor)
                      hold on
                  end
              end
              if showfacelabels
                  for i=1:mesh.totF
                      text(mesh.faces.centroid(1,i),mesh.faces.centroid(2,i),mesh.faces.centroid(3,i),num2str(mesh.faces.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',facelabelcolor)
                      hold on
                  end
              end
              if showcelllabels
                  for i=1:mesh.totC
                      text(mesh.cells.centroid(1,i),mesh.cells.centroid(2,i),mesh.cells.centroid(3,i),num2str(mesh.cells.id(i)),...
	                  'VerticalAlignment','middle',...
	                  'HorizontalAlignment','left',...
	                  'FontSize',labelsize,...
                      'Color',celllabelcolor)
                      hold on
                  end
              end
              grid on
              axis([xaxismin xaxismax yaxismin yaxismax zaxismin zaxismax]);
              axis equal
              xlabel(xstring)
              ylabel(ystring)
              zlabel(zstring)
              title(titlestring)
        end
end

return