%=============================
% Driver for Gauss Elimination
%=============================

file1 = fopen('mat_vec.dat');

 N   = fscanf(file1,'%f',[1,1]);
 A   = fscanf(file1,'%f',[N,N]);
 rhs = fscanf(file1,'%f',[1,N]);

fclose(file1);

A = A'; % because A was read columnwise,
        % replace with the transpose

Iwlpvt=1;  % pivoting enabled (0 to disable)
Isym=0;    % system is not symmetric

[x, ...
 l,u,det, ...
 Istop ] = gel (N,A,rhs,Iwlpvt,Isym);

disp ('Solution:'); x

%-----
% done
%-----
