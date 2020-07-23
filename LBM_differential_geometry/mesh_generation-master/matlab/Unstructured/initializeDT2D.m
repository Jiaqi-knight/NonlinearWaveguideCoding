function[edges,triangles]=initializeDT2D(p,printflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 23rd, 2014
%    Last update: April 24th, 2014
%
%          Input: N x 4 vector points in space (indices + coordinates)
%         Output: M x 1 vector of indeces of points belonging to the hull

%% Initialize triangulation - Radial sweep

%-------------------- Step 1

N = size(p,1);

sumx = 0;
sumy = 0;

for i=1:N
    sumx = sumx + p(i,2);
    sumy = sumy + p(i,3);
end

xm = sumx/N;
ym = sumy/N;

distances = zeros(N,1);

for i=1:N
    distances(i,1) = (p(i,2)-xm)^2+(p(i,3)-ym)^2;
end

[nearests,In] = mintournamenttree(distances,false,1);

nearest = p(In(1,1),:);

clear distances nearests In

edges = [];
it = 1;
for i=1:N
    if i~=nearest(1,1)
        edges = [edges; it nearest(1,1) p(i,1) (nearest(1,2)+p(i,2))/2 (nearest(1,3)+p(i,3))/2 0 0];
        it = it + 1;
    end
end

if printflag
   f1 = figure;
   plot(nearest(1,2),nearest(1,3),'rd','LineWidth',2)
   hold on
   for i=1:N
       if i~=nearest(1,1)
           plot(p(i,2),p(i,3),'b*','LineWidth',2)
           hold on
       end
   end
   for i=1:size(edges,1)
       plot([p(edges(i,2),2);p(edges(i,3),2)],[p(edges(i,2),3);p(edges(i,3),3)],'k--')
       hold on
   end
   grid on
   xlabel('x')
   ylabel('y')
   title('Step 1: radiating edges from centroid nearest neighbour')
   pause
end

%-------------------- Step 2

prem = [p(1:nearest(1,1)-1,:);p(nearest(1,1)+1:end,:)];
prem = [prem zeros(N-1,1) zeros(N-1,1)];
for i=1:N-1
    prem(i,4) = sqrt((prem(i,2)-nearest(1,2))^2+(prem(i,3)-nearest(1,3))^2);
    if (prem(i,2)-nearest(1,2))>10^-15
        prem(i,5) = (180/pi)*atan((prem(i,3)-nearest(1,3))/(prem(i,2)-nearest(1,2)));
    elseif (prem(i,2)-nearest(1,2))<-10^-15
        prem(i,5) = 180 + (180/pi)*atan((prem(i,3)-nearest(1,3))/(prem(i,2)-nearest(1,2)));
    elseif abs((prem(i,2)-nearest(1,2)))<10^-15 && (prem(i,3)-nearest(1,3))>10^-15
        prem(i,5) = 90;
    else
        prem(i,5) = -90;
    end
end

[prem,Ipr] = mintournamenttree(prem,false,5);

tempprem = prem;
prem = [];
it = 1;
while it<=N-1
    count = 0;
    checknext = 1;
    while checknext && count<=N-1-it && it<N-1
        if tempprem(it+(count+1),5)==tempprem(it,5)
            count = count + 1;
        else
            checknext = 0;
        end
    end
    if count~=0
        [part,Ipar] = mintournamenttree(tempprem(it:it+count,:),false,4);
        prem = [prem;part];
    else
        prem = [prem;tempprem(it,:)];
    end
    it = it + (count+1);
end

clear tempprem

M = size(edges,1);
for i=1:N-1
    if i==N-1
        edges = [edges; M+i prem(i,1) prem(1,1) (prem(i,2)+prem(1,2))/2 (prem(i,3)+prem(1,3))/2 0 0];
    else
        edges = [edges; M+i prem(i,1) prem(i+1,1) (prem(i,2)+prem(i+1,2))/2 (prem(i,3)+prem(i+1,3))/2 0 0];
    end
end

triangles = [];

for i=1:N-1
    if i==N-1
        triangles = [triangles; i 1 1 i-1 (nearest(1,2)+prem(i,2)+prem(1,2))/3 (nearest(1,3)+prem(i,3)+prem(1,3))/3 nearest(1,1) prem(i,1) prem(1,1) getedge_fromnodes2D(nearest(1,1),prem(i,1),edges) getedge_fromnodes2D(prem(i,1),prem(1,1),edges) getedge_fromnodes2D(nearest(1,1),prem(1,1),edges)];
    elseif i==1
        triangles = [triangles; i i+1 i+1 N-1 (nearest(1,2)+prem(i,2)+prem(i+1,2))/3 (nearest(1,3)+prem(i,3)+prem(i+1,3))/3 nearest(1,1) prem(i,1) prem(i+1,1) getedge_fromnodes2D(nearest(1,1),prem(i,1),edges) getedge_fromnodes2D(prem(i,1),prem(i+1,1),edges) getedge_fromnodes2D(nearest(1,1),prem(i+1,1),edges)];
    else
        triangles = [triangles; i i+1 i+1 i-1 (nearest(1,2)+prem(i,2)+prem(i+1,2))/3 (nearest(1,3)+prem(i,3)+prem(i+1,3))/3 nearest(1,1) prem(i,1) prem(i+1,1) getedge_fromnodes2D(nearest(1,1),prem(i,1),edges) getedge_fromnodes2D(prem(i,1),prem(i+1,1),edges) getedge_fromnodes2D(nearest(1,1),prem(i+1,1),edges)];
    end   
end

if printflag
   f2 = figure;
   plot(nearest(1,2),nearest(1,3),'rd','LineWidth',2)
   hold on
   text(nearest(1,2),nearest(1,3),num2str(nearest(1,1)),...
	    'VerticalAlignment','middle',...
	    'HorizontalAlignment','left',...
	    'FontSize',12,...
        'Color',[0 0 1])
   hold on
   for i=1:N
       if i~=nearest(1,1)
           plot(p(i,2),p(i,3),'b*','LineWidth',2)
           hold on
           text(p(i,2),p(i,3),num2str(p(i,1)),...
	            'VerticalAlignment','middle',...
	            'HorizontalAlignment','left',...
	            'FontSize',12,...
                'Color',[0 0 1])
           hold on
       end
   end
   for i=1:size(edges,1)
       plot([p(edges(i,2),2);p(edges(i,3),2)],[p(edges(i,2),3);p(edges(i,3),3)],'k--')
       hold on
       text(edges(i,4),edges(i,5),num2str(edges(i,1)),...
	        'VerticalAlignment','middle',...
	        'HorizontalAlignment','left',...
	        'FontSize',12,...
            'Color',[0 0 0])
       hold on
   end
   for i=1:size(triangles,1)
       text(triangles(i,5),triangles(i,6),num2str(triangles(i,1)),...
	        'VerticalAlignment','middle',...
	        'HorizontalAlignment','left',...
	        'FontSize',12,...
            'Color',[1 0 0])
       hold on
   end
   grid on
   xlabel('x')
   ylabel('y')
   title('Step 2: star-shaped domain')
   pause
end

%-------------------- Step 3

[hull,A,B,C,D]=quickhull2D(p,false);

hullsize = size(hull,1);

clear hull A B C D

i = 1;
while ~isempty(prem)
    Aindex = rem(i,size(prem,1));
    Bindex = rem(i+1,size(prem,1));
    Cindex = rem(i+2,size(prem,1));
    if ~Aindex
        Aindex = size(prem,1);
    end
    if ~Bindex
        Bindex = size(prem,1);
    end
    if ~Cindex
        Cindex = size(prem,1);
    end
    A = prem(Aindex,:);
    B = prem(Bindex,:);
    C = prem(Cindex,:);
    areaABC = orient2D(A(1,2:3)',B(1,2:3)',C(1,2:3)');
    if areaABC<0
        edges = [edges; edges(end,1)+1 A(1,1) B(1,1) (A(1,2)+B(1,2))/2 (A(1,3)+B(1,3))/2 0 0; edges(end,1)+2 A(1,1) C(1,1) (A(1,2)+C(1,2))/2 (A(1,3)+C(1,3))/2 0 0; edges(end,1)+3 B(1,1) C(1,1) (B(1,2)+C(1,2))/2 (B(1,3)+C(1,3))/2 0 0];
        triangles = [triangles; triangles(end,1)+1 0 0 0 (A(1,2)+B(1,2)+C(1,2))/3 (A(1,3)+B(1,3)+C(1,3))/3 A(1,1) C(1,1) B(1,1) getedge_fromnodes2D(A(1,1),C(1,1),edges) getedge_fromnodes2D(C(1,1),B(1,1),edges) getedge_fromnodes2D(B(1,1),A(1,1),edges)];
        prem = [prem(1:Bindex-1,:);prem(Bindex+1:end,:)];
    else
        i = mod(i + 1,size(prem,1));
    end
    if size(prem,1)==hullsize
        prem = [];
    end 
end

for i=1:size(edges,1)
    [tri1,tri2] = gettriangleadjacenttoedge2D(edges(i,:),triangles);
    edges(i,6) = tri1;
    edges(i,7) = tri2;
end

for i=1:size(triangles,1)
    [tri1,tri2,tri3]=gettriangleadjacenttotriangle2D(triangles(i,:),edges,triangles);
    triangles(i,2) = tri1;
    triangles(i,3) = tri2;
    triangles(i,4) = tri3;
end

if printflag
   f3 = figure;
   plot(nearest(1,2),nearest(1,3),'rd','LineWidth',2)
   hold on
   text(nearest(1,2),nearest(1,3),num2str(nearest(1,1)),...
	    'VerticalAlignment','middle',...
	    'HorizontalAlignment','left',...
	    'FontSize',12,...
        'Color',[0 0 1])
   hold on
   for i=1:N
       if i~=nearest(1,1)
           plot(p(i,2),p(i,3),'b*','LineWidth',2)
           hold on
           text(p(i,2),p(i,3),num2str(p(i,1)),...
	            'VerticalAlignment','middle',...
	            'HorizontalAlignment','left',...
	            'FontSize',12,...
                'Color',[0 0 1])
           hold on
       end
   end
   for i=1:size(edges,1)
       plot([p(edges(i,2),2);p(edges(i,3),2)],[p(edges(i,2),3);p(edges(i,3),3)],'k--')
       hold on
       text(edges(i,4),edges(i,5),num2str(edges(i,1)),...
	        'VerticalAlignment','middle',...
	        'HorizontalAlignment','left',...
	        'FontSize',12,...
            'Color',[0 0 0])
       hold on
   end
   for i=1:size(triangles,1)
       text(triangles(i,5),triangles(i,6),num2str(triangles(i,1)),...
	        'VerticalAlignment','middle',...
	        'HorizontalAlignment','left',...
	        'FontSize',12,...
            'Color',[1 0 0])
       hold on
   end
   grid on
   xlabel('x')
   ylabel('y')
   title('Step 3: fill boundary')
   pause
end

%-------------------- Step 4

return