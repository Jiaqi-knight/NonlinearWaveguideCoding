function[upperhull,A,B,C]=quickupperhull2D(p,printflag)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: April 23rd, 2014
%    Last update: April 23rd, 2014
%
%          Input: N x 3 vector points in plane (indices + coordinates)
%         Output: M x 1 vector of indeces of points belonging to the upper hull

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

%------------------------------ UPPER HULL --------------------------------

Cindex = 0;
ABCarea = 0;
tocheckupper = [];
for i=1:N-2
    area = orient2D(A(1,2:end)',B(1,2:end)',p(i,2:end)');
    if area>0 && area>ABCarea
        if Cindex~=0
            tocheckupper = [tocheckupper;Cindex];
        end
        Cindex = i;
        ABCarea = area;
    elseif area>0
        tocheckupper = [tocheckupper;i];
    end
end

C = p(Cindex,:);

if printflag
    figure();
    plot(p(:,2),p(:,3),'*')
    hold on
    for i=1:size(tocheckupper,1)
        plot(p(tocheckupper(i,1),2),p(tocheckupper(i,1),3),'gd','LineWidth',2)
        hold on
    end
    plot(A(:,2),A(:,3),'rd','LineWidth',2)
    hold on
    plot(B(:,2),B(:,3),'rd','LineWidth',2)
    hold on
    plot(C(:,2),C(:,3),'kd','LineWidth',2)
    hold on
    plot([A(:,2);B(:,2)],[A(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    grid on
    pause
end

tocheckupperAC = [];
tocheckupperCB = [];
for i=1:size(tocheckupper,1)
    area = orient2D(A(1,2:end)',B(1,2:end)',p(tocheckupper(i),2:end)');
    areaAC = orient2D(A(1,2:end)',C(1,2:end)',p(tocheckupper(i),2:end)');
    if (areaAC>0 && area>0) || (orient2D(C(1,2:end)',B(1,2:end)',p(tocheckupper(i),2:end)')>0 && area>0)
        if areaAC>0
            tocheckupperAC = [tocheckupperAC;tocheckupper(i)];
        else
            tocheckupperCB = [tocheckupperCB;tocheckupper(i)];
        end
    end
end

if printflag
    figure();
    plot(p(:,2),p(:,3),'*')
    hold on
    for i=1:size(tocheckupperAC,1)
        plot(p(tocheckupperAC(i,1),2),p(tocheckupperAC(i,1),3),'gd','LineWidth',2)
        hold on
    end
    plot(A(:,2),A(:,3),'rd','LineWidth',2)
    hold on
    plot(B(:,2),B(:,3),'rd','LineWidth',2)
    hold on
    plot(C(:,2),C(:,3),'kd','LineWidth',2)
    hold on
    plot([A(:,2);B(:,2)],[A(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    plot([A(:,2);C(:,2)],[A(:,3);C(:,3)],'--k','LineWidth',2)
    hold on
    grid on
    fprintf('There are %d points to check\n',size(tocheckupperAC,1));
    pause
end

tocheckupperAC = [tocheckupperAC A(1,1)*ones(length(tocheckupperAC),1) A(1,2)*ones(length(tocheckupperAC),1) A(1,3)*ones(length(tocheckupperAC),1) C(1,1)*ones(length(tocheckupperAC),1) C(1,2)*ones(length(tocheckupperAC),1) C(1,3)*ones(length(tocheckupperAC),1)];
indexhullAC = [];
while ~isempty(tocheckupperAC)
    Xindex = 0;
    ACXarea = 0;
    localA = tocheckupperAC(1,2:4);
    localC = tocheckupperAC(1,5:7);
    currentindex = 0;
    for i=1:size(tocheckupperAC,1)
        area = orient2D(localA(1,2:end)',localC(1,2:end)',p(tocheckupperAC(i),2:end)');
        if area>0 && area>ACXarea
            currentindex = i;
            Xindex = tocheckupperAC(i,1);
            ACXarea = area;
        end
    end
    tocheckupperAC = [tocheckupperAC(1:currentindex-1,:); tocheckupperAC(currentindex+1:end,:)];
    tempcheckupperAC = tocheckupperAC;
    tocheckupperAC = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullAC = [indexhullAC;Xindex];
        if printflag
            plot([X(:,2);localA(:,2)],[X(:,3);localA(:,3)],'--k','LineWidth',2)
            hold on
            plot([localC(:,2);X(:,2)],[localC(:,3);X(:,3)],'--k','LineWidth',2)
            hold on
        end
        for i=1:size(tempcheckupperAC,1)
            if tempcheckupperAC(i,2)==localA(1,1) && tempcheckupperAC(i,5)==localC(1,1)
                if orient2D(localA(1,2:end)',X(1,2:end)',p(tempcheckupperAC(i,1),2:end)')>0 || orient2D(X(1,2:end)',localC(1,2:end)',p(tempcheckupperAC(i,1),2:end)')>0
                    if orient2D(localA(1,2:end)',X(1,2:end)',p(tempcheckupperAC(i,1),2:end)')>0
                        tocheckupperAC = [tocheckupperAC; tempcheckupperAC(i,1) tempcheckupperAC(i,2:4) p(Xindex,:)];
                    else
                        tocheckupperAC = [tocheckupperAC; tempcheckupperAC(i,1) p(Xindex,:) tempcheckupperAC(i,5:7)];
                    end
                end
            else
                tocheckupperAC = [tocheckupperAC; tempcheckupperAC(i,:)];
            end
        end
    end
    if printflag
        plot([localA(:,2);localC(:,2)],[localA(:,3);localC(:,3)],'--k','LineWidth',2)
        hold on
        for i=1:size(tocheckupperAC,1)
            plot(p(tocheckupperAC(i,1),2),p(tocheckupperAC(i,1),3),'gd','LineWidth',2)
            hold on
        end
        grid on
        xlabel('x')
        ylabel('y')
        title('Convex hull')
        fprintf('There are %d points to check\n',size(tocheckupperAC,1));
        pause
    end
end

if printflag
    figure();
    plot(p(:,2),p(:,3),'*')
    hold on
    for i=1:size(tocheckupperCB,1)
        plot(p(tocheckupperCB(i,1),2),p(tocheckupperCB(i,1),3),'gd','LineWidth',2)
        hold on
    end
    plot(A(:,2),A(:,3),'rd','LineWidth',2)
    hold on
    plot(B(:,2),B(:,3),'rd','LineWidth',2)
    hold on
    plot(C(:,2),C(:,3),'kd','LineWidth',2)
    hold on
    plot([A(:,2);B(:,2)],[A(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    plot([C(:,2);B(:,2)],[C(:,3);B(:,3)],'--k','LineWidth',2)
    hold on
    grid on
    fprintf('There are %d points to check\n',size(tocheckupperCB,1));
    pause
end

tocheckupperCB = [tocheckupperCB C(1,1)*ones(length(tocheckupperCB),1) C(1,2)*ones(length(tocheckupperCB),1) C(1,3)*ones(length(tocheckupperCB),1) B(1,1)*ones(length(tocheckupperCB),1) B(1,2)*ones(length(tocheckupperCB),1) B(1,3)*ones(length(tocheckupperCB),1)];
indexhullCB = [];
while ~isempty(tocheckupperCB)
    Xindex = 0;
    CBXarea = 0;
    localC = tocheckupperCB(1,2:4);
    localB = tocheckupperCB(1,5:7);
    currentindex = 0;
    for i=1:size(tocheckupperCB,1)
        area = orient2D(localC(1,2:end)',localB(1,2:end)',p(tocheckupperCB(i),2:end)');
        if area>0 && area>CBXarea
            currentindex = i;
            Xindex = tocheckupperCB(i,1);
            CBXarea = area;
        end
    end
    tocheckupperCB = [tocheckupperCB(1:currentindex-1,:); tocheckupperCB(currentindex+1:end,:)];
    tempcheckupperCB = tocheckupperCB;
    tocheckupperCB = [];
    if Xindex~=0
        X = p(Xindex,:);
        indexhullCB = [indexhullCB;Xindex];
        if printflag
            plot([X(:,2);localC(:,2)],[X(:,3);localC(:,3)],'--k','LineWidth',2)
            hold on
            plot([localB(:,2);X(:,2)],[localB(:,3);X(:,3)],'--k','LineWidth',2)
            hold on
        end
        for i=1:size(tempcheckupperCB,1)
            if tempcheckupperCB(i,2)==localC(1,1) && tempcheckupperCB(i,5)==localB(1,1)
                if orient2D(localC(1,2:end)',X(1,2:end)',p(tempcheckupperCB(i,1),2:end)')>0 || orient2D(X(1,2:end)',localB(1,2:end)',p(tempcheckupperCB(i,1),2:end)')>0
                    if orient2D(localC(1,2:end)',X(1,2:end)',p(tempcheckupperCB(i),2:end)')>0
                        tocheckupperCB = [tocheckupperCB; tempcheckupperCB(i,1) tempcheckupperCB(i,2:4) p(Xindex,:)];
                    else
                        tocheckupperCB = [tocheckupperCB; tempcheckupperCB(i,1) p(Xindex,:) tempcheckupperCB(i,5:7)];
                    end
                end
            else
                tocheckupperCB = [tocheckupperCB; tempcheckupperCB(i,:)];
            end
        end
    end
    if printflag
        plot([localC(:,2);localB(:,2)],[localC(:,3);localB(:,3)],'--k','LineWidth',2)
        hold on
        for i=1:size(tocheckupperCB,1)
            plot(p(tocheckupperCB(i,1),2),p(tocheckupperCB(i,1),3),'gd','LineWidth',2)
            hold on
        end
        grid on
        xlabel('x')
        ylabel('y')
        title('Convex hull')
        fprintf('There are %d points to check\n',size(tocheckupperCB,1));
        pause
    end
end

upperhullAC = [];
for i=1:size(indexhullAC,1)
    upperhullAC = [upperhullAC;p(indexhullAC(i,1),:)];
end
upperhullCB = [];
for i=1:size(indexhullCB,1)
    upperhullCB = [upperhullCB;p(indexhullCB(i,1),:)];
end

if ~isempty(upperhullAC)
    [upperhullAC,Iac] = maxtournamenttree(upperhullAC,false,2);
end
if ~isempty(upperhullCB)
    [upperhullCB,Icb] = maxtournamenttree(upperhullCB,false,2);
end

%------------------------------ CONVEX HULL --------------------------------

upperhull = [B;upperhullCB;C;upperhullAC;A];

return