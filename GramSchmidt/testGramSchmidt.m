function testGramSchmidt
%% A simple test for the Stabilized Gram Schmidt method for orthonormalization
% Input is a matrix including k column vector k = 2,3,...
% Output is the same matrix but the vectors are replaced with orthonormal
% vectors
%
%
% -------------------------------------------------
% code by: Reza Ahmadzadeh (reza.ahmadzadeh@iit.it
% -------------------------------------------------
%% test GramSchmidt in 2D
clc,clear,close all
v = [3 2;1 2];
h=figure;hold on
for ii=1:2
    plotVector(v(:,ii),h,'r--o');
end
W = GramSchmidt(v);
for ii=1:2
    plotVector(W(:,ii),h,'b');
end
axis equal
grid on;box on;

    function plotVector(V,h,C)
        % a function to plot a vector V on a figure given by handle h
        % C specifies the line properties
        figure(h);
        V = V(:);
        V = V / norm(V);
        plot([0,V(1,1)],[0,V(2,1)],C,'linewidth',1.5);
    end

%% test GramSchmidt in 2D
v = [3 2;1 2;-1 0];
h=figure;hold on
for ii=1:2
    plotVector3(v(:,ii),h,'r--o');
end
W = GramSchmidt(v);
for ii=1:2
    plotVector3(W(:,ii),h,'b');
end
view([123,18]);
axis equal
grid on;box on;

    function plotVector3(V,h,C)
        % a function to plot a vector V on a figure given by handle h
        % C specifies the line properties
        figure(h);
        V = V(:);
        V = V / norm(V);
        plot3([0,V(1,1)],[0,V(2,1)],[0,V(3,1)],C,'linewidth',1.5);
    end
end