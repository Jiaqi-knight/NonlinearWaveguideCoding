%MgivenL
% Function which generates a cell array of possible magnetic quantum
% numbers given an orbital quantum number letter X.
%
% NOTE: Energy levels in Hydrogenic type atoms are characterized by quantum
% numbers n = 1, 2, 3, 4....
% Orbital quantum numbers are:
% L = 0,1,2....(n-1)
% M = 0, +/-1, +/-2,...+/- L
%
% Orbitals are typically labelled in the form nXM e.g. 3P1
% X is characterized by L.
% L = 0, X = S   min{n} = 1   max(abs(M)) = 0
% L = 1, X = P   min{n} = 2   max(abs(M)) = 1
% L = 2, X = D   min{n} = 3   max(abs(M)) = 2
% L = 3, X = F   min{n} = 4   max(abs(M)) = 3
% L = 4, X = G   min{n} = 5   max(abs(M)) = 4
% L = 5, X = H   min{n} = 6   max(abs(M)) = 5
% L = 6, X = I   min{n} = 7   max(abs(M)) = 6
%
% By convention XM orbitals are the linear combination of +abs(M) and
% -abs(M) orbitals. For example 3P-1 will combine the (angular) orbitals
% Y(L=1,M=1) - Y(L=1,M=-1). Similarly 3P1 will correspond to Y(L=1,M=1) +
% Y(L=1,M=-1).

function Ms = MgivenL(X)

%Work out orbital quantum number L from X
L = orbital2L(X);

%Work out possible values of magnetic quantum numbers M
M = -L:L;

%Create cell array
Ms = cell(1,length(M));
for n=1:length(M)
    Ms{n} = num2str(M(n));
end

%End of code