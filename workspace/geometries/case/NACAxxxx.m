function[airfoil,upperairfoil,lowerairfoil,yc]=NACAxxxx(Nx,x1,x2,x3_4,c)

%%
%        Project: Fluid-Structure Interaction on Deformable Surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 14th, 2014
%    Last update: July 16th, 2014
%
%          Input: 
%         Output: 

%%

x = (-0.5*c:c/(Nx-1):0.5*c)';

if x1==0 && x2==0
    t = x3_4/100;
    ytplus = (t./0.2).*c.*(0.2969.*(x./c + 0.5).^0.5 -0.1260.*(x./c + 0.5) -0.3516.*(x./c + 0.5).^2 +0.2843.*(x./c + 0.5).^3 -0.1036.*(x./c + 0.5).^4);
    ytminus = -(t./0.2).*c.*(0.2969.*(x./c + 0.5).^0.5 -0.1260.*(x./c + 0.5) -0.3516.*(x./c + 0.5).^2 +0.2843.*(x./c + 0.5).^3 -0.1036.*(x./c + 0.5).^4);
    upperairfoil = [x ytplus];
    lowerairfoil = [x(2:end-1) ytminus(2:end-1)];
    airfoil = [x ytplus;x ytminus];
    yc = zeros(size(x,1),1);
else
    m = x1/100;
    p = x2/10;
    t = x3_4/100;
    ytplus = (t./0.2).*c.*(0.2969.*(x./c + 0.5).^0.5 -0.1260.*(x./c + 0.5) -0.3516.*(x./c + 0.5).^2 +0.2843.*(x./c + 0.5).^3 -0.1036.*(x./c + 0.5).^4);
    ytminus = -(t./0.2).*c.*(0.2969.*(x./c + 0.5).^0.5 -0.1260.*(x./c + 0.5) -0.3516.*(x./c + 0.5).^2 +0.2843.*(x./c + 0.5).^3 -0.1036.*(x./c + 0.5).^4);
    pc = p*c;
    yc = zeros(Nx,1);
    for i=1:Nx
        if x(i)<=pc
            yc(i) = m.*x(i).*(2.*p-(x(i)./c + 0.5))./(p.^2);
        else
            yc(i) = m.*c.*(1-(x(i)./c + 0.5)).*(1-2.*p+(x(i)./c + 0.5))./((1-p).^2);
        end
    end
    upperairfoil = [x yc+ytplus];
    lowerairfoil = [x(2:end-1) yc(2:end-1)+ytplus(2:end-1)];
    airfoil = [x yc+ytplus;x yc+ytminus];
end

return