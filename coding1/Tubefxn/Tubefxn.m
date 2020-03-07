function [X,Y,Z] = tubeN(Rfxn,Nr,Paramtricfxn,tbound,Nt)
% This function is the numerial analogy to the tube function in maple. It
% generates the X,Y and Z coordinates to be used by the surf function to
% generate the intended tube figure

% Rfxn is the parametric function of Radius as a function of t
% Paramtricfxn is a cell of path parametric functions
% tbound range is the range of values of t
% Nr is the number of steps in r
% Nt is the number of steps in t
T = linspace(tbound(1),tbound(2),Nt + 1)';
R = Rfxn(T);
[Zo,Yo,Xo] = cylinder(R,Nr);
Xo = 0*Xo;
Trans = zeros(numel(T),3);
for m = 1:3
    Trans(:,m) = Paramtricfxn{m}(T);
end

Angle = DIFFhandle(Paramtricfxn,T); 
Myz  = @(thetayz) Mxyz([cos(thetayz(1)-pi/2),sin(thetayz(1)-pi/2),0],thetayz(2))*Mxyz([0,0,1],thetayz(1));
                   
                  
for n = 1:Nt + 1
    mat = Myz(Angle(n,:))*[Xo(n,:);Yo(n,:);Zo(n,:);ones(1,Nr + 1)];
    Xo(n,:) = mat(1,:);
    Yo(n,:) = mat(2,:);
    Zo(n,:) = mat(3,:);
end

X = Xo + Trans(:,1)*ones(1,Nr + 1);
Y = Yo + Trans(:,2)*ones(1,Nr + 1);
Z = Zo + Trans(:,3)*ones(1,Nr + 1);


function Theta = DIFFhandle(Paramtricfxn,T)
Tpdt    = T + 1e-10;
POS     = zeros(numel(T),3,2);

for m = 1:3
    POS(:,m,1) = Paramtricfxn{m}(T);
    POS(:,m,2) = Paramtricfxn{m}(Tpdt);
end
DEL = diff(POS,1,3);
Theta1 = atan2(DEL(:,2),DEL(:,1));
Theta2 = atan2(DEL(:,3),sqrt(DEL(:,1).^2 + DEL(:,2).^2));
findx  = find(POS(:,1,1), 1);
findy  = find(POS(:,2,1), 1);
findz  = find(POS(:,3,1), 1);
if isempty(findy) || isempty(findx)  
    Theta1  = Unwrap(Theta1);
end
if isempty(findz) || (isempty(findx) && isempty(findy))
    Theta2  = Unwrap(Theta2);
end
Theta  = [Theta1,Theta2];

function y = Unwrap(V)
y = V;
for n = 2:numel(V)
    if y(n) - y(n-1) > 3.1
        y(n:end) = y(n:end) - pi;
    elseif y(n) - y(n-1) < -3.1
        y(n:end) = y(n:end) + pi;
    end
end


function M = Mxyz(V,theta)
c          = cos(theta);
s          = sin(theta);
l          = V(1);
m          = V(2);
n          = V(3);
M          = eye(4);
M(1:3,1:3) = [l*l*(1 - c) +   c  m*l*(1 - c) - n*s  n*l*(1 - c) + m*s
              m*l*(1 - c) + n*s  m*m*(1 - c) +   c  n*m*(1 - c) - l*s
              n*l*(1 - c) - m*s  m*n*(1 - c) + l*s  n*n*(1 - c) +   c]; 
            
        
    