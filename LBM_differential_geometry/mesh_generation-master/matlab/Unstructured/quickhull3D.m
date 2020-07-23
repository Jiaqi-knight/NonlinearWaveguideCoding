function[hull,edges,A,B,C,D]=quickhull3D(p,printflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 17th, 2014
%    Last update: April 22nd, 2014
%
%          Input: N x 4 vector points in space (indices + coordinates)
%         Output: M x 1 vector of indeces of points belonging to the hull

%%

N = size(p,1);

[Xmin,Ixmin] = mintournamenttree(p,true,2);
[Xmax,Ixmax] = maxtournamenttree(p,true,2);
[Ymin,Iymin] = mintournamenttree(p,true,3);
[Ymax,Iymax] = maxtournamenttree(p,true,3);
[Zmin,Izmin] = mintournamenttree(p,true,4);
[Zmax,Izmax] = maxtournamenttree(p,true,4);

distances = [Xmin(1,:) Xmax(1,:) 0;...  % 1
             Xmin(1,:) Ymin(1,:) 0;...  % 2
             Xmin(1,:) Ymax(1,:) 0;...  % 3
             Xmin(1,:) Zmin(1,:) 0;...  % 4
             Xmin(1,:) Zmax(1,:) 0;...  % 5
             Xmax(1,:) Ymin(1,:) 0;...  % 6
             Xmax(1,:) Ymax(1,:) 0;...  % 7
             Xmax(1,:) Zmin(1,:) 0;...  % 8
             Xmax(1,:) Zmax(1,:) 0;...  % 9
             Ymin(1,:) Ymax(1,:) 0;...  % 10
             Ymin(1,:) Zmin(1,:) 0;...  % 11
             Ymin(1,:) Zmax(1,:) 0;...  % 12
             Ymax(1,:) Zmin(1,:) 0;...  % 13
             Ymax(1,:) Zmax(1,:) 0;...  % 14
             Zmin(1,:) Zmax(1,:) 0;...  % 15
             ];

for i=1:size(distances,1)
    distances(i,9) = sqrt((distances(i,6)-distances(i,2))^2+(distances(i,7)-distances(i,3))^2+(distances(i,8)-distances(i,4))^2);
end

[distmax,Idist] = maxtournamenttree(distances,true,2);

A = distmax(1,1:4);
B = distmax(1,5:8);

switch Idist
    case 1
        tocheck = [Ymin(1,:);...
                   Ymax(1,:);...
                   Zmin(1,:);...
                   Zmax(1,:) ...
                 ];
    case 2
        tocheck = [Xmax(1,:);...
                   Ymax(1,:);...
                   Zmin(1,:);...
                   Zmax(1,:) ...
                 ];
    case 3
        tocheck = [Xmax(1,:);...
                   Ymin(1,:);...
                   Zmin(1,:);...
                   Zmax(1,:) ...
                   ];
    case 4
        tocheck = [Xmax(1,:);...
                   Ymin(1,:);...
                   Ymax(1,:);...
                   Zmax(1,:) ...
                 ];
    case 5
        tocheck = [Xmax(1,:);...
                   Ymin(1,:);...
                   Ymax(1,:);...
                   Zmin(1,:) ...
                 ];
    case 6
        tocheck = [Xmin(1,:);...
                   Ymax(1,:);...
                   Zmin(1,:);...
                   Zmax(1,:) ...
                 ];
    case 7
        tocheck = [Xmin(1,:);...
                   Ymin(1,:);...
                   Zmin(1,:);...
                   Zmax(1,:) ...
                 ];
    case 8
        tocheck = [Xmin(1,:);...
                   Ymin(1,:);...
                   Ymax(1,:);...
                   Zmax(1,:) ...
                 ];
    case 9
        tocheck = [Xmin(1,:);...
                   Ymin(1,:);...
                   Ymax(1,:);...
                   Zmin(1,:) ...
                 ];
    case 10
        tocheck = [Xmin(1,:);...
                   Xmax(1,:);...
                   Zmin(1,:);...
                   Zmax(1,:) ...
                 ];
    case 11
        tocheck = [Xmin(1,:);...
                   Xmax(1,:);...
                   Ymax(1,:);...
                   Zmax(1,:) ...
                 ];
    case 12
        tocheck = [Xmin(1,:);...
                   Xmax(1,:);...
                   Ymax(1,:);...
                   Zmin(1,:)...
                 ];
    case 13
        tocheck = [Xmin(1,:);...
                   Xmax(1,:);...
                   Ymin(1,:);...
                   Zmax(1,:) ...
                 ];
    case 14
        tocheck = [Xmin(1,:);...
                   Xmax(1,:);...
                   Ymin(1,:);...
                   Zmin(1,:) ...
                 ];
    case 15
        tocheck = [Xmin(1,:);...
                   Xmax(1,:);...
                   Ymin(1,:);...
                   Ymax(1,:) ...
                 ];
end

Cindex = 0;
ABCarea = 0;
for i=1:length(tocheck)
    area = orient2D(A(1,2:end)',B(1,2:end)',tocheck(i,2:end)');
    if abs(area)>ABCarea
        ABCarea = abs(area);
        Cindex = i;
    end
end

C = tocheck(Cindex,:);

Dindex = 0;
ABCDvolume = 0;
for i=1:N
    volume = orient3D(A(1,2:end)',B(1,2:end)',C(1,2:end)',p(i,2:end)');
    if abs(volume)>ABCDvolume
        ABCDvolume = abs(volume);
        Dindex = i;
    end
end

D = p(Dindex,:);

signumABC = sign(orient3D(A(1,2:end)',B(1,2:end)',C(1,2:end)',D(1,2:end)'));
signumABD = sign(orient3D(A(1,2:end)',B(1,2:end)',D(1,2:end)',C(1,2:end)'));
signumCAD = sign(orient3D(C(1,2:end)',A(1,2:end)',D(1,2:end)',B(1,2:end)'));
signumBCD = sign(orient3D(B(1,2:end)',C(1,2:end)',D(1,2:end)',A(1,2:end)'));
tocheck = [];
for i=1:N
    if ~(signumABC*orient3D(A(1,2:end)',B(1,2:end)',C(1,2:end)',p(i,2:end)')>=0 && signumABD*orient3D(A(1,2:end)',B(1,2:end)',D(1,2:end)',p(i,2:end)')>=0 && signumCAD*orient3D(C(1,2:end)',A(1,2:end)',D(1,2:end)',p(i,2:end)')>=0 && signumBCD*orient3D(B(1,2:end)',C(1,2:end)',D(1,2:end)',p(i,2:end)')>=0)
        tocheck = [tocheck; i];
    end
end

tocheckABC = [];
for i=1:size(tocheck,1)
    if ~(signumABC*orient3D(A(1,2:end)',B(1,2:end)',C(1,2:end)',p(tocheck(i,1),2:end)')>=0)
        tocheckABC = [tocheckABC; tocheck(i,1)];
    end
end

tocheckABD = [];
for i=1:size(tocheck,1)
    if ~(signumABD*orient3D(A(1,2:end)',B(1,2:end)',D(1,2:end)',p(tocheck(i,1),2:end)')>=0)
        tocheckABD = [tocheckABD; tocheck(i,1)];
    end
end

tocheckCAD = [];
for i=1:size(tocheck,1)
    if ~(signumCAD*orient3D(C(1,2:end)',A(1,2:end)',D(1,2:end)',p(tocheck(i,1),2:end)')>=0)
        tocheckCAD = [tocheckCAD; tocheck(i,1)];
    end
end

tocheckBCD = [];
for i=1:size(tocheck,1)
    if ~(signumBCD*orient3D(B(1,2:end)',C(1,2:end)',D(1,2:end)',p(tocheck(i,1),2:end)')>=0)
        tocheckBCD = [tocheckBCD; tocheck(i,1)];
    end
end

if printflag
    figure();
    plot3(p(:,2),p(:,3),p(:,4),'*')
    hold on
    for i=1:size(tocheck,1)
        plot3(p(tocheck(i,1),2),p(tocheck(i,1),3),p(tocheck(i,1),4),'gd','LineWidth',2)
        hold on
    end
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'kd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([C(:,2);D(:,2)],[C(:,3);D(:,3)],[C(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Building the convex hull - 1')
    pause
    figure();
    for i=1:size(tocheck,1)
        plot3(p(tocheck(i,1),2),p(tocheck(i,1),3),p(tocheck(i,1),4),'gd','LineWidth',2)
        hold on
    end
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'kd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([C(:,2);D(:,2)],[C(:,3);D(:,3)],[C(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Building the convex hull - 2')
    pause
    figure();
    for i=1:size(tocheckABC,1)
        plot3(p(tocheckABC(i,1),2),p(tocheckABC(i,1),3),p(tocheckABC(i,1),4),'gd','LineWidth',2)
        hold on
    end
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'kd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([C(:,2);D(:,2)],[C(:,3);D(:,3)],[C(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Building the convex hull - 3')
    pause
    figure();
    for i=1:size(tocheckABD,1)
        plot3(p(tocheckABD(i,1),2),p(tocheckABD(i,1),3),p(tocheckABD(i,1),4),'gd','LineWidth',2)
        hold on
    end
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'kd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([C(:,2);D(:,2)],[C(:,3);D(:,3)],[C(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Building the convex hull - 4')
    pause
    figure();
    for i=1:size(tocheckCAD,1)
        plot3(p(tocheckCAD(i,1),2),p(tocheckCAD(i,1),3),p(tocheckCAD(i,1),4),'gd','LineWidth',2)
        hold on
    end
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'kd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([C(:,2);D(:,2)],[C(:,3);D(:,3)],[C(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Building the convex hull - 5')
    pause
    figure();
    for i=1:size(tocheckBCD,1)
        plot3(p(tocheckBCD(i,1),2),p(tocheckBCD(i,1),3),p(tocheckBCD(i,1),4),'gd','LineWidth',2)
        hold on
    end
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'kd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([C(:,2);D(:,2)],[C(:,3);D(:,3)],[C(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Building the convex hull - 6')
    pause
end

%------------------------------- FACE ABC ---------------------------------

if printflag
    fabc = figure;
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Face ABC')
end

tocheckABC = [tocheckABC A(1,1)*ones(length(tocheckABC),1) A(1,2)*ones(length(tocheckABC),1) A(1,3)*ones(length(tocheckABC),1) A(1,4)*ones(length(tocheckABC),1) B(1,1)*ones(length(tocheckABC),1) B(1,2)*ones(length(tocheckABC),1) B(1,3)*ones(length(tocheckABC),1) B(1,4)*ones(length(tocheckABC),1) C(1,1)*ones(length(tocheckABC),1) C(1,2)*ones(length(tocheckABC),1) C(1,3)*ones(length(tocheckABC),1) C(1,4)*ones(length(tocheckABC),1)];
indexhullABC = [];
while ~isempty(tocheckABC)
    Xindex = 0;
    ABCXvolume = 0;
    localA = tocheckABC(1,2:5);
    localB = tocheckABC(1,6:9);
    localC = tocheckABC(1,10:13);
    currentindex = 0;
    for i=1:size(tocheckABC,1)
        volume = orient3D(localA(1,2:end)',localB(1,2:end)',localC(1,2:end)',p(tocheckABC(i),2:end)');
        if abs(volume)>ABCXvolume
            currentindex = i;
            Xindex = tocheckABC(i,1);
            ABCXvolume = volume;
        end
    end
    tempcheckABC = [tocheckABC(1:currentindex-1,:); tocheckABC(currentindex+1:end,:)];
    tocheckABC = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullABC = [indexhullABC;Xindex];
        if printflag
            figure(fabc);
            plot3(X(:,2),X(:,3),X(:,4),'kd','LineWidth',2)
            hold on
            plot3([localA(:,2);localB(:,2)],[localA(:,3);localB(:,3)],[localA(:,4);localB(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localA(:,2);localC(:,2)],[localA(:,3);localC(:,3)],[localA(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localB(:,2);localC(:,2)],[localB(:,3);localC(:,3)],[localB(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localA(:,2);X(:,2)],[localA(:,3);X(:,3)],[localA(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            plot3([X(:,2);localC(:,2)],[X(:,3);localC(:,3)],[X(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localB(:,2);X(:,2)],[localB(:,3);X(:,3)],[localB(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            pause
        end
        signumABC = sign(orient3D(localA(1,2:end)',localB(1,2:end)',localC(1,2:end)',X(1,2:end)'));
        signumABX = sign(orient3D(localA(1,2:end)',localB(1,2:end)',X(1,2:end)',localC(1,2:end)'));
        signumCAX = sign(orient3D(localC(1,2:end)',localA(1,2:end)',X(1,2:end)',localB(1,2:end)'));
        signumBCX = sign(orient3D(localB(1,2:end)',localC(1,2:end)',X(1,2:end)',localA(1,2:end)'));
        for i=1:size(tempcheckABC,1)
            if tempcheckABC(i,2)==localA(1,1) && tempcheckABC(i,6)==localB(1,1) && tempcheckABC(i,10)==localC(1,1)
                if ~(signumABC*orient3D(localA(1,2:end)',localB(1,2:end)',localC(1,2:end)',p(tempcheckABC(i,1),2:end)')>=0 && signumABX*orient3D(localA(1,2:end)',localB(1,2:end)',X(1,2:end)',p(tempcheckABC(i,1),2:end)')>=0 && signumCAX*orient3D(localC(1,2:end)',localA(1,2:end)',X(1,2:end)',p(tempcheckABC(i,1),2:end)')>=0 && signumBCX*orient3D(localB(1,2:end)',localC(1,2:end)',X(1,2:end)',p(tempcheckABC(i,1),2:end)')>=0)
                    if ~(signumABX*orient3D(localA(1,2:end)',localB(1,2:end)',X(1,2:end)',p(tempcheckABC(i,1),2:end)')>=0)
                        tocheckABC = [tocheckABC; tempcheckABC(i,1) localA(1,:) localB(1,:) X(1,:)];
                    elseif ~(signumCAX*orient3D(localC(1,2:end)',localA(1,2:end)',X(1,2:end)',p(tempcheckABC(i,1),2:end)')>=0)
                        tocheckABC = [tocheckABC; tempcheckABC(i,1) localC(1,:) localA(1,:) X(1,:)];
                    else
                        tocheckABC = [tocheckABC; tempcheckABC(i,1) localB(1,:) localC(1,:) X(1,:)];
                    end
                end
            else
                tocheckABC = [tocheckABC; tempcheckABC(i,:)];
            end
        end
    end
    if printflag
        figure(fabc);
        for i=1:size(tocheckABC,1)
            plot3(p(tocheckABC(i,1),2),p(tocheckABC(i,1),3),p(tocheckABC(i,1),4),'gd','LineWidth',2)
            hold on
        end
        fprintf('There are %d points to check\n',size(tocheckABC,1));
        pause
    end
end

%------------------------------- FACE ABD ---------------------------------

if printflag
    fabd = figure;
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'rd','LineWidth',2)
    hold on
    plot3([A(:,2);B(:,2)],[A(:,3);B(:,3)],[A(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);D(:,2)],[B(:,3);D(:,3)],[B(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Face ABD')
end

tocheckABD = [tocheckABD A(1,1)*ones(length(tocheckABD),1) A(1,2)*ones(length(tocheckABD),1) A(1,3)*ones(length(tocheckABD),1) A(1,4)*ones(length(tocheckABD),1) B(1,1)*ones(length(tocheckABD),1) B(1,2)*ones(length(tocheckABD),1) B(1,3)*ones(length(tocheckABD),1) B(1,4)*ones(length(tocheckABD),1) D(1,1)*ones(length(tocheckABD),1) D(1,2)*ones(length(tocheckABD),1) D(1,3)*ones(length(tocheckABD),1) D(1,4)*ones(length(tocheckABD),1)];
indexhullABD = [];
while ~isempty(tocheckABD)
    Xindex = 0;
    ABDXvolume = 0;
    localA = tocheckABD(1,2:5);
    localB = tocheckABD(1,6:9);
    localD = tocheckABD(1,10:13);
    currentindex = 0;
    for i=1:size(tocheckABD,1)
        volume = orient3D(localA(1,2:end)',localB(1,2:end)',localD(1,2:end)',p(tocheckABD(i),2:end)');
        if abs(volume)>ABDXvolume
            currentindex = i;
            Xindex = tocheckABD(i,1);
            ABDXvolume = volume;
        end
    end
    tempcheckABD = [tocheckABD(1:currentindex-1,:); tocheckABD(currentindex+1:end,:)];
    tocheckABD = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullABD = [indexhullABD;Xindex];
        if printflag
            figure(fabd);
            plot3(X(:,2),X(:,3),X(:,4),'kd','LineWidth',2)
            hold on
            plot3([localA(:,2);localB(:,2)],[localA(:,3);localB(:,3)],[localA(:,4);localB(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localA(:,2);localD(:,2)],[localA(:,3);localD(:,3)],[localA(:,4);localD(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localB(:,2);localD(:,2)],[localB(:,3);localD(:,3)],[localB(:,4);localD(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localA(:,2);X(:,2)],[localA(:,3);X(:,3)],[localA(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            plot3([X(:,2);localD(:,2)],[X(:,3);localD(:,3)],[X(:,4);localD(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localB(:,2);X(:,2)],[localB(:,3);X(:,3)],[localB(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            pause
        end
        signumABD = sign(orient3D(localA(1,2:end)',localB(1,2:end)',localD(1,2:end)',X(1,2:end)'));
        signumABX = sign(orient3D(localA(1,2:end)',localB(1,2:end)',X(1,2:end)',localD(1,2:end)'));
        signumDAX = sign(orient3D(localD(1,2:end)',localA(1,2:end)',X(1,2:end)',localB(1,2:end)'));
        signumBDX = sign(orient3D(localB(1,2:end)',localD(1,2:end)',X(1,2:end)',localA(1,2:end)'));
        for i=1:size(tempcheckABD,1)
            if tempcheckABD(i,2)==localA(1,1) && tempcheckABD(i,6)==localB(1,1) && tempcheckABD(i,10)==localD(1,1)
                if ~(signumABD*orient3D(localA(1,2:end)',localB(1,2:end)',localD(1,2:end)',p(tempcheckABD(i,1),2:end)')>=0 && signumABX*orient3D(localA(1,2:end)',localB(1,2:end)',X(1,2:end)',p(tempcheckABD(i,1),2:end)')>=0 && signumDAX*orient3D(localD(1,2:end)',localA(1,2:end)',X(1,2:end)',p(tempcheckABD(i,1),2:end)')>=0 && signumBDX*orient3D(localB(1,2:end)',localD(1,2:end)',X(1,2:end)',p(tempcheckABD(i,1),2:end)')>=0)
                    if ~(signumABX*orient3D(localA(1,2:end)',localB(1,2:end)',X(1,2:end)',p(tempcheckABD(i,1),2:end)')>=0)
                        tocheckABD = [tocheckABD; tempcheckABD(i,1) localA(1,:) localB(1,:) X(1,:)];
                    elseif ~(signumDAX*orient3D(localD(1,2:end)',localA(1,2:end)',X(1,2:end)',p(tempcheckABD(i,1),2:end)')>=0)
                        tocheckABD = [tocheckABD; tempcheckABD(i,1) localD(1,:) localA(1,:) X(1,:)];
                    else
                        tocheckABD = [tocheckABD; tempcheckABD(i,1) localB(1,:) localD(1,:) X(1,:)];
                    end
                end
            else
                tocheckABD = [tocheckABD; tempcheckABD(i,:)];
            end
        end
    end
    if printflag
        figure(fabd);
        for i=1:size(tocheckABD,1)
            plot3(p(tocheckABD(i,1),2),p(tocheckABD(i,1),3),p(tocheckABD(i,1),4),'gd','LineWidth',2)
            hold on
        end
        fprintf('There are %d points to check\n',size(tocheckABD,1));
        pause
    end
end

%------------------------------- FACE CAD ---------------------------------

if printflag
    fcad = figure;
    plot3(A(:,2),A(:,3),A(:,4),'rd','LineWidth',2)
    hold on
    plot3(D(:,2),D(:,3),D(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3([A(:,2);D(:,2)],[A(:,3);D(:,3)],[A(:,4);D(:,4)],'--k','LineWidth',2)
    hold on
    plot3([A(:,2);C(:,2)],[A(:,3);C(:,3)],[A(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([D(:,2);C(:,2)],[D(:,3);C(:,3)],[D(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Face CAD')
end

tocheckCAD = [tocheckCAD C(1,1)*ones(length(tocheckCAD),1) C(1,2)*ones(length(tocheckCAD),1) C(1,3)*ones(length(tocheckCAD),1) C(1,4)*ones(length(tocheckCAD),1) A(1,1)*ones(length(tocheckCAD),1) A(1,2)*ones(length(tocheckCAD),1) A(1,3)*ones(length(tocheckCAD),1) A(1,4)*ones(length(tocheckCAD),1) D(1,1)*ones(length(tocheckCAD),1) D(1,2)*ones(length(tocheckCAD),1) D(1,3)*ones(length(tocheckCAD),1) D(1,4)*ones(length(tocheckCAD),1)];
indexhullCAD = [];
while ~isempty(tocheckCAD)
    Xindex = 0;
    CADXvolume = 0;
    localC = tocheckCAD(1,2:5);
    localA = tocheckCAD(1,6:9);
    localD = tocheckCAD(1,10:13);
    currentindex = 0;
    for i=1:size(tocheckCAD,1)
        volume = orient3D(localC(1,2:end)',localA(1,2:end)',localD(1,2:end)',p(tocheckCAD(i),2:end)');
        if abs(volume)>CADXvolume
            currentindex = i;
            Xindex = tocheckCAD(i,1);
            CADXvolume = volume;
        end
    end
    tempcheckCAD = [tocheckCAD(1:currentindex-1,:); tocheckCAD(currentindex+1:end,:)];
    tocheckCAD = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullCAD = [indexhullCAD;Xindex];
        if printflag
            figure(fcad);
            plot3(X(:,2),X(:,3),X(:,4),'kd','LineWidth',2)
            hold on
            plot3([localA(:,2);localD(:,2)],[localA(:,3);localD(:,3)],[localA(:,4);localD(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localA(:,2);localC(:,2)],[localA(:,3);localC(:,3)],[localA(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localD(:,2);localC(:,2)],[localD(:,3);localC(:,3)],[localD(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localA(:,2);X(:,2)],[localA(:,3);X(:,3)],[localA(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            plot3([X(:,2);localC(:,2)],[X(:,3);localC(:,3)],[X(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localD(:,2);X(:,2)],[localD(:,3);X(:,3)],[localD(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            pause
        end
        signumCAD = sign(orient3D(localC(1,2:end)',localA(1,2:end)',localD(1,2:end)',X(1,2:end)'));
        signumCAX = sign(orient3D(localC(1,2:end)',localA(1,2:end)',X(1,2:end)',localD(1,2:end)'));
        signumDCX = sign(orient3D(localD(1,2:end)',localC(1,2:end)',X(1,2:end)',localA(1,2:end)'));
        signumADX = sign(orient3D(localA(1,2:end)',localD(1,2:end)',X(1,2:end)',localC(1,2:end)'));
        for i=1:size(tempcheckCAD,1)
            if tempcheckCAD(i,2)==localC(1,1) && tempcheckCAD(i,6)==localA(1,1) && tempcheckCAD(i,10)==localD(1,1)
                if ~(signumCAD*orient3D(localC(1,2:end)',localA(1,2:end)',localD(1,2:end)',p(tempcheckCAD(i,1),2:end)')>=0 && signumCAX*orient3D(localC(1,2:end)',localA(1,2:end)',X(1,2:end)',p(tempcheckCAD(i,1),2:end)')>=0 && signumDCX*orient3D(localD(1,2:end)',localC(1,2:end)',X(1,2:end)',p(tempcheckCAD(i,1),2:end)')>=0 && signumADX*orient3D(localA(1,2:end)',localD(1,2:end)',X(1,2:end)',p(tempcheckCAD(i,1),2:end)')>=0)
                    if ~(signumCAX*orient3D(localC(1,2:end)',localA(1,2:end)',X(1,2:end)',p(tempcheckCAD(i,1),2:end)')>=0)
                        tocheckCAD = [tocheckCAD; tempcheckCAD(i,1) localC(1,:) localA(1,:) X(1,:)];
                    elseif ~(signumDCX*orient3D(localD(1,2:end)',localC(1,2:end)',X(1,2:end)',p(tempcheckCAD(i,1),2:end)')>=0)
                        tocheckCAD = [tocheckCAD; tempcheckCAD(i,1) localD(1,:) localC(1,:) X(1,:)];
                    else
                        tocheckCAD = [tocheckCAD; tempcheckCAD(i,1) localA(1,:) localD(1,:) X(1,:)];
                    end
                end
            else
                tocheckCAD = [tocheckCAD; tempcheckCAD(i,:)];
            end
        end
    end
    if printflag
        figure(fcad);
        for i=1:size(tocheckCAD,1)
            plot3(p(tocheckCAD(i,1),2),p(tocheckCAD(i,1),3),p(tocheckCAD(i,1),4),'gd','LineWidth',2)
            hold on
        end
        fprintf('There are %d points to check\n',size(tocheckCAD,1));
        pause
    end
end

%------------------------------- FACE BCD ---------------------------------

if printflag
    fbcd = figure;
    plot3(D(:,2),D(:,3),D(:,4),'rd','LineWidth',2)
    hold on
    plot3(B(:,2),B(:,3),B(:,4),'rd','LineWidth',2)
    hold on
    plot3(C(:,2),C(:,3),C(:,4),'rd','LineWidth',2)
    hold on
    plot3([D(:,2);B(:,2)],[D(:,3);B(:,3)],[D(:,4);B(:,4)],'--k','LineWidth',2)
    hold on
    plot3([D(:,2);C(:,2)],[D(:,3);C(:,3)],[D(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    plot3([B(:,2);C(:,2)],[B(:,3);C(:,3)],[B(:,4);C(:,4)],'--k','LineWidth',2)
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Face BCD')
end

tocheckBCD = [tocheckBCD B(1,1)*ones(length(tocheckBCD),1) B(1,2)*ones(length(tocheckBCD),1) B(1,3)*ones(length(tocheckBCD),1) B(1,4)*ones(length(tocheckBCD),1) C(1,1)*ones(length(tocheckBCD),1) C(1,2)*ones(length(tocheckBCD),1) C(1,3)*ones(length(tocheckBCD),1) C(1,4)*ones(length(tocheckBCD),1) D(1,1)*ones(length(tocheckBCD),1) D(1,2)*ones(length(tocheckBCD),1) D(1,3)*ones(length(tocheckBCD),1) D(1,4)*ones(length(tocheckBCD),1)];
indexhullBCD = [];
while ~isempty(tocheckBCD)
    Xindex = 0;
    BCDXvolume = 0;
    localB = tocheckBCD(1,2:5);
    localC = tocheckBCD(1,6:9);
    localD = tocheckBCD(1,10:13);
    currentindex = 0;
    for i=1:size(tocheckBCD,1)
        volume = orient3D(localB(1,2:end)',localC(1,2:end)',localD(1,2:end)',p(tocheckBCD(i),2:end)');
        if abs(volume)>BCDXvolume
            currentindex = i;
            Xindex = tocheckBCD(i,1);
            BCDXvolume = volume;
        end
    end
    tempcheckBCD = [tocheckBCD(1:currentindex-1,:); tocheckBCD(currentindex+1:end,:)];
    tocheckBCD = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullBCD = [indexhullBCD;Xindex];
        if printflag
            figure(fbcd);
            plot3(X(:,2),X(:,3),X(:,4),'kd','LineWidth',2)
            hold on
            plot3([localD(:,2);localB(:,2)],[localD(:,3);localB(:,3)],[localD(:,4);localB(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localD(:,2);localC(:,2)],[localD(:,3);localC(:,3)],[localD(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localB(:,2);localC(:,2)],[localB(:,3);localC(:,3)],[localB(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localD(:,2);X(:,2)],[localD(:,3);X(:,3)],[localD(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            plot3([X(:,2);localC(:,2)],[X(:,3);localC(:,3)],[X(:,4);localC(:,4)],'--k','LineWidth',2)
            hold on
            plot3([localB(:,2);X(:,2)],[localB(:,3);X(:,3)],[localB(:,4);X(:,4)],'--k','LineWidth',2)
            hold on
            pause
        end
        signumBCD = sign(orient3D(localB(1,2:end)',localC(1,2:end)',localD(1,2:end)',X(1,2:end)'));
        signumBCX = sign(orient3D(localB(1,2:end)',localC(1,2:end)',X(1,2:end)',localD(1,2:end)'));
        signumDBX = sign(orient3D(localD(1,2:end)',localB(1,2:end)',X(1,2:end)',localC(1,2:end)'));
        signumCDX = sign(orient3D(localC(1,2:end)',localD(1,2:end)',X(1,2:end)',localB(1,2:end)'));
        for i=1:size(tempcheckBCD,1)
            if tempcheckBCD(i,2)==localB(1,1) && tempcheckBCD(i,6)==localC(1,1) && tempcheckBCD(i,10)==localD(1,1)
                if ~(signumBCD*orient3D(localB(1,2:end)',localC(1,2:end)',localD(1,2:end)',p(tempcheckBCD(i,1),2:end)')>=0 && signumBCX*orient3D(localB(1,2:end)',localC(1,2:end)',X(1,2:end)',p(tempcheckBCD(i,1),2:end)')>=0 && signumDBX*orient3D(localD(1,2:end)',localB(1,2:end)',X(1,2:end)',p(tempcheckBCD(i,1),2:end)')>=0 && signumCDX*orient3D(localC(1,2:end)',localD(1,2:end)',X(1,2:end)',p(tempcheckBCD(i,1),2:end)')>=0)
                    if ~(signumBCX*orient3D(localB(1,2:end)',localC(1,2:end)',X(1,2:end)',p(tempcheckBCD(i,1),2:end)')>=0)
                        tocheckBCD = [tocheckBCD; tempcheckBCD(i,1) localB(1,:) localC(1,:) X(1,:)];
                    elseif ~(signumDBX*orient3D(localD(1,2:end)',localB(1,2:end)',X(1,2:end)',p(tempcheckBCD(i,1),2:end)')>=0)
                        tocheckBCD = [tocheckBCD; tempcheckBCD(i,1) localD(1,:) localB(1,:) X(1,:)];
                    else
                        tocheckBCD = [tocheckBCD; tempcheckBCD(i,1) localC(1,:) localD(1,:) X(1,:)];
                    end
                end
            else
                tocheckBCD = [tocheckBCD; tempcheckBCD(i,:)];
            end
        end
    end
    if printflag
        figure(fbcd);
        for i=1:size(tocheckBCD,1)
            plot3(p(tocheckBCD(i,1),2),p(tocheckBCD(i,1),3),p(tocheckBCD(i,1),4),'gd','LineWidth',2)
            hold on
        end
        fprintf('There are %d points to check\n',size(tocheckBCD,1));
        pause
    end
end

%--------------------------------- HULL -----------------------------------

hull = [];

for i=1:size(indexhullABC,1)
    hull = [hull;p(indexhullABC(i,1),:)];
end
for i=1:size(indexhullABD,1)
    hull = [hull;p(indexhullABD(i,1),:)];
end
for i=1:size(indexhullCAD,1)
    hull = [hull;p(indexhullCAD(i,1),:)];
end
for i=1:size(indexhullBCD,1)
    hull = [hull;p(indexhullBCD(i,1),:)];
end

[hull,Ihull] = mintournamenttree(hull,false,2);

M = size(hull,1);

edges = zeros(M,4);
for i=1:M
    ind1 = 0;
    dist1 = 0;
    ind2 = 0;
    dist2 = 0;
    ind3 = 0;
    dist3 = 0;
    for j=4:M
        if j~=i
            distance = (hull(j,2)-hull(i,2))^2+(hull(j,3)-hull(i,3))^2+(hull(j,4)-hull(i,4))^2;
            if ~ind1
                ind1 = j;
                dist1 = distance;
            elseif ~ind2
                ind2 = j;
                dist2 = distance;
            elseif ~ind3
                ind3 = j;
                dist3 = distance;
            elseif distance <= dist1
                ind3 = ind2;
                ind2 = ind1;
                ind1 = j;
                dist3 = dist2;
                dist2 = dist1;
                dist1 = distance;
            elseif distance>dist1 && distance <= dist2
                ind3 = ind2;
                ind2 = j;
                dist3 = dist2;
                dist2 = distance;
            elseif distance>dist2 && distance <= dist3
                ind3 = j;
                dist3 = distance;
            end
        end
    end
    edges(i,1) = i;
    edges(i,2) = ind1;
    edges(i,3) = ind2;
    edges(i,4) = ind3;
end

return