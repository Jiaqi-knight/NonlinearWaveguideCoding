clear all

%----
% FSELIB
%
% confirm the eigenvalues of the projection matrix
%----

N=3

al = 0.3;

P = [1-2*al al 0;
    al 1-2*al al;
    0 al 1-2*al];

Q = [4 1 0;
    1 4 1;
    0 1 4]/6.0;

PP = inv(Q)*(P+Q-eye(3));
eig(PP)

for i=1:N
 sns = sin( i/(N+1) * pi/2)^2;
 lam = ( 1- (4*al+2/3)*sns )/(1-2/3*sns)
end
