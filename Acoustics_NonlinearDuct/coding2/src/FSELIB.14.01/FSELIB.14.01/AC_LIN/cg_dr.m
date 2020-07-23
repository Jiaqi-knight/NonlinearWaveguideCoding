
file1 = fopen('mat_s_vec.dat');
N   = fscanf(file1,'%f',[1,1]);
A   = fscanf(file1,'%f',[N,N]);
rhs = fscanf(file1,'%f',[1,N]);
fclose(file1);

[sln] = cg (N, A, rhs);

disp ('Solution:'); sln
