clear all

%---
% One-element beam buckling
%---

A = [4 2; 2 4];
B = [4 -1; -1 4]/30;
eig(A,B)
