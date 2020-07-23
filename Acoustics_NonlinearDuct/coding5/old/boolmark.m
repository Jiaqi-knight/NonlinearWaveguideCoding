clear
clc
%Here the size of 1.png is 300*900 
A=imread('1.png');  
thresh=graythresh(A)
B=im2bw(A,thresh);
B=flipud(B);%´¹Ö±·­ÞD
B=rot90(B,3);
C=~B;
imshow(C);
dlmwrite('1.dat',C,'delimiter',' ');