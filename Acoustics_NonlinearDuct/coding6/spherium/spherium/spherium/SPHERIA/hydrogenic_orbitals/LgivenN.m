%LgivenN
% Function which generates a cell array of possible orbital letters
% corresponding to a hydrogenic orbital with quantum number N.
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

function Ls = LgivenN(N)
L = 0:N-1;
Ls = cell(1,N-1);
for n=1:N
    if L(n)==0
        Ls{n} = 'S';
    end
    if L(n)==1
        Ls{n} = 'P';
    end
    if L(n)==2
        Ls{n} = 'D';
    end
    if L(n)==3
        Ls{n} = 'F';
    end
    if L(n)==4
        Ls{n} = 'G';
    end
    if L(n)==5
        Ls{n} = 'H';
    end
    if L(n)==6
        Ls{n} = 'I';
    end
    if L(n)==7
        Ls{n} = 'J';
    end
    if L(n)==8
        Ls{n} = 'K';
    end
    if L(n)==9
        Ls{n} = 'L';
    end
end


%End of code