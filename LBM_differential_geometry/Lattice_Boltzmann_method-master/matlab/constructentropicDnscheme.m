function[scheme,QD]=constructentropicDnscheme(D,scheme1D,pruneflag,prunelevels)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 27th, 2014
%    Last update: May 27th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

if nargin==2
    pruneflag = 0;
    prunelevels = -1;
end

v1D = scheme1D(:,1:end-2);
W1D = scheme1D(:,end-1);

Q1D = size(v1D,1);

scheme = [];

if pruneflag
    switch D
        case 2
            for i=1:Q1D
                for j=1:Q1D
                    prune = 0;
                    shell = round(v1D(i,1)^2+v1D(j,1)^2);
                    level = 1;
                    while ~prune && level <=size(prunelevels,1)
                        if prunelevels(level,1)==shell
                            prune = 1;
                            level = size(prunelevels,1) + 1;
                        else
                            level = level + 1;
                        end
                    end
                    if ~prune
                        scheme = [scheme; v1D(i,1) v1D(j,1) W1D(i,1)*W1D(j,1) sqrt(v1D(i,1)^2+v1D(j,1)^2)];
                    end
                end
            end
        case 3
            for i=1:Q1D
                for j=1:Q1D
                    for k=1:Q1D
                        prune = 0;
                        shell = round(v1D(i,1)^2+v1D(j,1)^2+v1D(k,1)^2);
                        level = 1;
                        while ~prune && level <=size(prunelevels,1)
                            if prunelevels(level,1)==shell
                                prune = 1;
                                level = size(prunelevels,1) + 1;
                            else
                                level = level + 1;
                            end
                        end
                        if ~prune
                            scheme = [scheme; v1D(i,1) v1D(j,1) v1D(k,1) W1D(i,1)*W1D(j,1)*W1D(k,1) sqrt(v1D(i,1)^2+v1D(j,1)^2+v1D(k,1)^2)];
                        end
                    end
                end
            end
    end
else
    switch D
        case 2
            for i=1:Q1D
                for j=1:Q1D
                    scheme = [scheme; v1D(i,1) v1D(j,1) W1D(i,1)*W1D(j,1) sqrt(v1D(i,1)^2+v1D(j,1)^2)];
                end
            end
        case 3
            for i=1:Q1D
                for j=1:Q1D
                    for k=1:Q1D
                        scheme = [scheme; v1D(i,1) v1D(j,1) v1D(k,1) W1D(i,1)*W1D(j,1)*W1D(k,1) sqrt(v1D(i,1)^2+v1D(j,1)^2+v1D(k,1)^2)];
                    end
                end
            end
    end
end

scheme = sortrows(scheme,size(scheme,2));
QD = size(scheme,1);

return