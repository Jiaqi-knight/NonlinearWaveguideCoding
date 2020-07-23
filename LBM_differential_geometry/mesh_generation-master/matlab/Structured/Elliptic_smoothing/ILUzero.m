function[A,invA]=ILUzero(A)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 19th, 2014
%    Last update: May 19th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

[m n]= size(A);

if m~=n
    disp('Matrix A is rectangular. Only square matrices can be handled by this function.');
else
    for k=1:n-1
        for i=k+1:n
            if A(i,k)~=0
                if A(k,k)==0
                    disp('Error: null pivot')
                else
                    A(i,k) = A(i,k)/A(k,k);
                    for j=k+1:n
                        if A(i,j)~=0
                            A(i,j) = A(i,j)-A(i,k)*A(k,j);
                        end
                    end
                end
            end
        end
    end
end

invA = luinv(A);

return