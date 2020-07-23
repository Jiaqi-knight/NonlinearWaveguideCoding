function outstring = leadingzero(X,N)

outstring = [strrep(blanks(N-1 - fix(log10(X))),' ','0'),num2str(X,3)];