function L = orbital2L(X)

%Compute orbital quantum number L
if strcmp( X, 'S' )
    L=0;
elseif strcmp( X, 'P' )
    L=1;
elseif strcmp( X, 'D' )
    L=2;
elseif strcmp( X, 'F' )
    L=3;
elseif strcmp( X, 'G' )
    L=4;
elseif strcmp( X, 'H' )
    L=5;
elseif strcmp( X, 'I' )
    L=6;
elseif strcmp( X, 'J' )
    L=7;
elseif strcmp( X, 'K' )
    L=8;
elseif strcmp( X, 'L' )
    L=9;
elseif strcmp( X, 'M' )
    L=10;
elseif strcmp( X, 'N' )
    L=11;
else
    L = NaN;
end

%End of code