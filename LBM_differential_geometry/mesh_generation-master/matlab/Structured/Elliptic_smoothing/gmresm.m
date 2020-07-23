function[x,restart]=gmresm(A,b,x0,m,tol,maxrestart)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 21st, 2014
%    Last update: May 21st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

N = size(A,1);

err = 1;
restart = 0;
while err>=tol && restart<=maxrestart
    v = zeros(N,m+1);
    h = zeros(m+1,m);
    r0 = b - A*x0;
    beta = sqrt(sum(r0.^2));
    v(:,1) = r0/beta;
    j = 1;
    innerloop = 1;
    while innerloop && j<=m
        wj = A*v(:,j);
        for i=1:j
            h(i,j) = wj'*v(:,i);
            wj = wj -h(i,j)*v(:,i);
        end
        h(j+1,j) = sqrt(sum(wj.^2));
        if h(j+1,j)==0
            innerloop = 0;
        else
            v(:,j+1) = wj/h(j+1,j);
            if j==m
                innerloop = 0;
            else
                j = j + 1;
            end
        end
    end
    if j~=m
        v = v(:,1:j);
        h = h(1:j+1,1:j);
    end
    r = zeros(j+1,1);
    r(1,1) = beta;
    for k=1:j
        P = eye(j+1);
        s = h(k+1,k)/(sqrt(h(k+1,k)^2+h(k,k)^2));
        c = h(k,k)/(sqrt(h(k+1,k)^2+h(k,k)^2));
        P(k,k) = c;
        P(k+1,k+1) = c;
        P(k,k+1) = s;
        P(k+1,k) = -s;
        h = P*h;
        r = P*r;
    end
    y = backsolve(h(1:j,:),r(1:j,:));
    x = x0 + v(:,1:j)*y;
    err = abs(r(j+1,1));
    x0 = x;
    restart = restart + 1;
end

return