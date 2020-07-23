%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was created to be used with a 2DQ9 lattice Boltzmann scheme 
% for static boundaries aligned with the lattice grid. Given the size of 
% the density matrix in terms of number of rows 'Nr'and number of columns 
% 'Mc', and the coordinate of the points in the lattice which defines the 
% line segment associated with the solid boundaries 'crossing3'finds the 
% vectors containing the coordinates of the cells that have crossed the 
% solid boundary in each direction of propagation. This information allows
% for the swapping of cells close to the boundary according to the
% bounce-back rule (no-slip).
%
% [vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing(nl,nc,xl,yl)
%
% Outputs:
% Vextors vec_i, where i corresponds to each propagation direction
% according to the scheme below.
%
%
%                       6   2   5
%                        \  |  /
%                       3 - . - 1
%                        /  |  \
%                       7   4   8
% 
% Inputs:
% Nr and Mc are the number of lines and number of colums of the main matrix.
% 
% xl and yl are the vectors containg the coordinates of the initial
% and end points defining the boundary line segments.
%
% By Andrey R. da Silva                             McGill 06/2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s

function conditions_wall = crossing3(Nr,Mc, wall_points)
warning off

xl = wall_points{1};
yl = wall_points{2};

% begin the main for-loop
for direc=1:8

    
% directions of propagation 1 2 3 4 5 6 7 8

propx = [1 0 -1 0 1 -1 -1 1];
propy = [0 1 0 -1 1 1 -1 -1];
desto = [];
ce = [];
le = [];

% % definindo a matrix principal
% mat=ones(Nr,Mc);

for p =1:length(xl)-1
%coeficientes da reta
A = (yl(p+1)- yl(p))/(xl(p+1)-xl(p)+eps);
B = yl(p+1)-A*xl(p+1);

% definindo a linha paralela para o poligono
xv=[xl(p:p+1) fliplr(xl(p:p+1)+propx(direc))];
yv=[yl(p:p+1) fliplr(yl(p:p+1)+propy(direc))];


% gera a caixa em torno do poligono
x = floor(min(xv)):1:ceil(max(xv));
y = floor(min(yv)):1:ceil(max(yv));

% descreve os elementos dentro da caixa de forma vetorial
tamanho=length(x)*length(y);
c=zeros(1,tamanho);
l=zeros(1,tamanho);
fo=0;
for u = 1 : length(x)
for t = 1 : length(y)
     c(t+fo*length(y)) = x(1+fo);
     l(t+fo*length(y)) = y(t);
end
fo=fo+1;
end

% acha os pontos dentro do poligono
[in on] = inpolygon(c,l,xv,yv);
ji=find(in==1);
c=c(ji);l=l(ji);
norma=sqrt(2);
% acha as distancias dos pontos que cruzaram a barreira da barreira em si,
% para cada direcao

if direc ==1 
disto = c-(l-B)/(A+eps);
elseif direc == 3
    disto = c-(l-B)/(A+eps);
elseif direc == 2
    disto = l-(A*c+B);
elseif direc == 4
    disto = l-(A*c+B);
elseif direc == 5
    dx = c-(B-l+c)/(1-A);
    dy = l-(A*(l-c)-B)/(A-1);
    disto = sqrt(dx.^2+dy.^2)/norma;
elseif direc == 7
    dx = c-(B-l+c)/(1-A);
    dy = l-(A*(l-c)-B)/(A-1);
    disto = sqrt(dx.^2+dy.^2)/norma;
elseif direc == 8
    dx = c-(B-l-c)/(-A-1);
    dy = l-(A*l+A*c+B)/(1+A);
    disto = sqrt(dx.^2+dy.^2)/norma;
elseif direc == 6
    dx = c-(B-l-c)/(-A-1);
    dy = l-(A*l+A*c+B)/(1+A);
    disto = sqrt(dx.^2+dy.^2)/norma;
else
end   
desto=abs([desto disto]);
ce = [ce c];
le = [le l];
end

indice = (ce*Nr+le-Nr)+(direc-1)*(Nr*Mc);
indices(direc,1:length(indice))=indice;

end
% end of the main for-loop

vec1 = indices(1,:); vec1 = vec1(vec1>0);
vec2 = indices(2,:); vec2 = vec2(vec2>0);
vec3 = indices(3,:); vec3 = vec3(vec3>0);
vec4 = indices(4,:); vec4 = vec4(vec4>0);
vec5 = indices(5,:); vec5 = vec5(vec5>0);
vec6 = indices(6,:); vec6 = vec6(vec6>0);
vec7 = indices(7,:); vec7 = vec7(vec7>0);
vec8 = indices(8,:); vec8 = vec8(vec8>0);

conditions_wall{1} = vec1;
conditions_wall{2} = vec2;
conditions_wall{3} = vec3;
conditions_wall{4} = vec4;
conditions_wall{5} = vec5;
conditions_wall{6} = vec6;
conditions_wall{7} = vec7;
conditions_wall{8} = vec8;
