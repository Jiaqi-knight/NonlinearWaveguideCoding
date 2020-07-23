function[lowerhull,A,B,D]=quicklowerhull2D(p,printflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 16th, 2014
%    Last update: April 16th, 2014
%
%          Input: N x 3 vector points in plane (indices + coordinates)
%         Output: M x 1 vector of indeces of points belonging to the hull

%%

N = size(p,1);

[A,Ia] = mintournamenttree(p,true,2);
[B,Ib] = maxtournamenttree(p,true,2);

pold = p;
p = [];
for i=1:N
    if i~=Ia && i~=Ib
        p = [p;pold(i,:)];
    end
end

[p,Ip] = mintournamenttree(p,false,2);

clear Ip

%------------------------------ LOWER HULL --------------------------------

Dindex = 0;
ABDarea = 0;
tochecklower = [];
for i=1:N-2
    area = orient2D(B(1,2:end)',A(1,2:end)',p(i,2:end)');
    if area>0 && area>ABDarea
        if Dindex~=0
            tochecklower = [tochecklower;Dindex];
        end
        Dindex = i;
        ABDarea = area;
    elseif area>=0
        tochecklower = [tochecklower;i];
    end
end

D = p(Dindex,:);

if printflag
    figure();
    plot(p(:,2),p(:,3),'*')
    hold on
    for i=1:size(tochecklower,1)
        plot(p(tochecklower(i,1),2),p(tochecklower(i,1),3),'gd','LineWidth',2)
        hold on
    end
    plot(A(:,2),A(:,3),'rd','LineWidth',2)
    hold on
    plot(B(:,2),B(:,3),'rd','LineWidth',2)
    hold on
    plot(D(:,2),D(:,3),'kd','LineWidth',2)
    hold on
    plot([A(:,2);B(:,2)],[A(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    grid on
    pause
end

tochecklowerAD = [];
tochecklowerDB = [];
for i=1:size(tochecklower,1)
    area = orient2D(B(1,2:end)',A(1,2:end)',p(tochecklower(i),2:end)');
    areaAD = orient2D(D(1,2:end)',A(1,2:end)',p(tochecklower(i),2:end)');
    if (areaAD>0 && area>0) || (orient2D(B(1,2:end)',D(1,2:end)',p(tochecklower(i),2:end)')>0 && area>0)
        if areaAD>=0
            tochecklowerAD = [tochecklowerAD;tochecklower(i)];
        else
            tochecklowerDB = [tochecklowerDB;tochecklower(i)];
        end
    end
end

if printflag
    figure();
    plot(p(:,2),p(:,3),'*')
    hold on
    for i=1:size(tochecklowerAD,1)
        plot(p(tochecklowerAD(i,1),2),p(tochecklowerAD(i,1),3),'gd','LineWidth',2)
        hold on
    end
    plot(A(:,2),A(:,3),'rd','LineWidth',2)
    hold on
    plot(B(:,2),B(:,3),'rd','LineWidth',2)
    hold on
    plot(C(:,2),C(:,3),'kd','LineWidth',2)
    hold on
    plot(D(:,2),D(:,3),'kd','LineWidth',2)
    hold on
    plot([A(:,2);B(:,2)],[A(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    plot([A(:,2);D(:,2)],[A(:,3);D(:,3)],'--k','LineWidth',2)
    hold on
    grid on
    fprintf('There are %d points to check\n',size(tochecklowerAD,1));
    pause
end

tochecklowerAD = [tochecklowerAD D(1,1)*ones(length(tochecklowerAD),1) D(1,2)*ones(length(tochecklowerAD),1) D(1,3)*ones(length(tochecklowerAD),1) A(1,1)*ones(length(tochecklowerAD),1) A(1,2)*ones(length(tochecklowerAD),1) A(1,3)*ones(length(tochecklowerAD),1)];
indexhullAD = [];
while ~isempty(tochecklowerAD)
    Xindex = 0;
    ADXarea = 0;
    localD = tochecklowerAD(1,2:4);
    localA = tochecklowerAD(1,5:7);
    currentindex = 0;
    for i=1:size(tochecklowerAD,1)
        area = orient2D(localD(1,2:end)',localA(1,2:end)',p(tochecklowerAD(i),2:end)');
        if area>0 && area>ADXarea
            currentindex = i;
            Xindex = tochecklowerAD(i,1);
            ADXarea = area;
        end
    end
    tochecklowerAD = [tochecklowerAD(1:currentindex-1,:); tochecklowerAD(currentindex+1:end,:)];
    tempchecklowerAD = tochecklowerAD;
    tochecklowerAD = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullAD = [indexhullAD;Xindex];
        if printflag
            plot([X(:,2);localD(:,2)],[X(:,3);localD(:,3)],'--k','LineWidth',2)
            hold on
            plot([localA(:,2);X(:,2)],[localA(:,3);X(:,3)],'--k','LineWidth',2)
            hold on
        end
        for i=1:size(tempchecklowerAD,1)
            if tempchecklowerAD(i,2)==localD(1,1) && tempchecklowerAD(i,5)==localA(1,1)
                if orient2D(localD(1,2:end)',X(1,2:end)',p(tempchecklowerAD(i,1),2:end)')>0 || orient2D(X(1,2:end)',localA(1,2:end)',p(tempchecklowerAD(i,1),2:end)')>0
                    if orient2D(localD(1,2:end)',X(1,2:end)',p(tempchecklowerAD(i,1),2:end)')>0
                        tochecklowerAD = [tochecklowerAD; tempchecklowerAD(i,1) tempchecklowerAD(i,2:4) p(Xindex,:)];
                    else
                        tochecklowerAD = [tochecklowerAD; tempchecklowerAD(i,1) p(Xindex,:) tempchecklowerAD(i,5:7)];
                    end
                end
            else
                tochecklowerAD = [tochecklowerAD; tempchecklowerAD(i,:)];
            end
        end
    end
    if printflag
        plot([localD(:,2);localA(:,2)],[localD(:,3);localA(:,3)],'--k','LineWidth',2)
        hold on
        for i=1:size(tochecklowerAD,1)
            plot(p(tochecklowerAD(i,1),2),p(tochecklowerAD(i,1),3),'gd','LineWidth',2)
            hold on
        end
        grid on
        xlabel('x')
        ylabel('y')
        title('Convex hull')
        fprintf('There are %d points to check\n',size(tochecklowerAD,1));
        pause
    end
end

if printflag
    figure();
    plot(p(:,2),p(:,3),'*')
    hold on
    for i=1:size(tochecklowerDB,1)
        plot(p(tochecklowerDB(i,1),2),p(tochecklowerDB(i,1),3),'gd','LineWidth',2)
        hold on
    end
    plot(A(:,2),A(:,3),'rd','LineWidth',2)
    hold on
    plot(B(:,2),B(:,3),'rd','LineWidth',2)
    hold on
    plot(C(:,2),C(:,3),'kd','LineWidth',2)
    hold on
    plot(D(:,2),D(:,3),'kd','LineWidth',2)
    hold on
    plot([A(:,2);B(:,2)],[A(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    plot([D(:,2);B(:,2)],[D(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    grid on
    fprintf('There are %d points to check\n',size(tochecklowerDB,1));
    pause
end

tochecklowerDB = [tochecklowerDB B(1,1)*ones(length(tochecklowerDB),1) B(1,2)*ones(length(tochecklowerDB),1) B(1,3)*ones(length(tochecklowerDB),1) D(1,1)*ones(length(tochecklowerDB),1) D(1,2)*ones(length(tochecklowerDB),1) D(1,3)*ones(length(tochecklowerDB),1)];
indexhullDB = [];
while ~isempty(tochecklowerDB)
    Xindex = 0;
    DBXarea = 0;
    localB = tochecklowerDB(1,2:4);
    localD = tochecklowerDB(1,5:7);
    currentindex = 0;
    for i=1:size(tochecklowerDB,1)
        area = orient2D(localB(1,2:end)',localD(1,2:end)',p(tochecklowerDB(i),2:end)');
        if area>0 && area>DBXarea
            currentindex = i;
            Xindex = tochecklowerDB(i,1);
            DBXarea = area;
        end
    end
    tochecklowerDB = [tochecklowerDB(1:currentindex-1,:); tochecklowerDB(currentindex+1:end,:)];
    tempchecklowerDB = tochecklowerDB;
    tochecklowerDB = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullDB = [indexhullDB;Xindex];
        if printflag
            plot([X(:,2);localB(:,2)],[X(:,3);localB(:,3)],'--k','LineWidth',2)
            hold on
            plot([localD(:,2);X(:,2)],[localD(:,3);X(:,3)],'--k','LineWidth',2)
            hold on
        end
        for i=1:size(tempchecklowerDB,1)
            if tempchecklowerDB(i,2)==localB(1,1) && tempchecklowerDB(i,5)==localD(1,1)
                if orient2D(localB(1,2:end)',X(1,2:end)',p(tempchecklowerDB(i,1),2:end)')>0 || orient2D(X(1,2:end)',localD(1,2:end)',p(tempchecklowerDB(i,1),2:end)')>0
                    if orient2D(localB(1,2:end)',X(1,2:end)',p(tempchecklowerDB(i,1),2:end)')>0
                        tochecklowerDB = [tochecklowerDB; tempchecklowerDB(i,1) tempchecklowerDB(i,2:4) p(Xindex,:)];
                    else
                        tochecklowerDB = [tochecklowerDB; tempchecklowerDB(i,1) p(Xindex,:) tempchecklowerDB(i,5:7)];
                    end
                end
            else
                tochecklowerDB = [tochecklowerDB; tempchecklowerDB(i,:)];
            end
        end
    end
    if printflag
        plot([localB(:,2);localD(:,2)],[localB(:,3);localD(:,3)],'--k','LineWidth',2)
        hold on
        for i=1:size(tochecklowerDB,1)
            plot(p(tochecklowerDB(i,1),2),p(tochecklowerDB(i,1),3),'gd','LineWidth',2)
            hold on
        end
        grid on
        xlabel('x')
        ylabel('y')
        title('Convex hull')
        fprintf('There are %d points to check\n',size(tochecklowerDB,1));
        pause
    end
end

lowerhullAD = [];
for i=1:size(indexhullAD,1)
    lowerhullAD = [lowerhullAD;p(indexhullAD(i,1),:)];
end
lowerhullDB = [];
for i=1:size(indexhullDB,1)
    lowerhullDB = [lowerhullDB;p(indexhullDB(i,1),:)];
end

if ~isempty(lowerhullAD)
    [lowerhullAD,Iad] = mintournamenttree(lowerhullAD,false,2);
end
if ~isempty(lowerhullDB)
    [lowerhullDB,Idb] = mintournamenttree(lowerhullDB,false,2);
end

%-------------------------- LOWER CONVEX HULL -----------------------------

lowerhull = [A;lowerhullAD;D;lowerhullDB;B];

return