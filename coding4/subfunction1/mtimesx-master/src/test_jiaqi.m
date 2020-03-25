
% A = rand(1,3,3,5)
% B = rand(3,1,2,3,5)
% 
% %A1=size(permute(A,[4,3,2,1]))
% Cx = mtimesx(A,B)
% size(Cx)
% 
% % A*B(:,:,1)
% % A*B(:,:,2)
% 
% 
% A = tensor(ones(4,2))
% B = tensor(rand(3,4,2))
% R = ttt(A,B,[1], [2])
% 
% %A.data
% 
% size(R)
% 
% size(permute(A,[1,3,2]))
% size(permute(B,[2,1,3]))
% 
mtimesx(Fun.M_2_1,Fun.V)

A = tensor(rand(4,3,2));
U = rand(3,4);
B = ttm(A,U,1);