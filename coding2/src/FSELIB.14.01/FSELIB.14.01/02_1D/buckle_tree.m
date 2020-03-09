clear all
close all

%==========================================
% FSELIB
%
% Code buckle_tree
%
% Column buckling equation
% with linear elements
%
% ne: number of elements
%==========================================

%-----------
% input data
%-----------

L = 1.0;

ne = 16; ratio = 1.0;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

%-----------------
% compact assembly
%-----------------

[at,bt,ct,ar,br,cr] = buckle_tree_sys (ne,xe,L);

A = zeros(ne,ne);
B = zeros(ne,ne);

for i=1:ne
 A(i,i) = at(i);
 B(i,i) = ar(i);
end

for i=1:ne-1
 A(i,i+1) = bt(i);
 A(i+1,i) = ct(i+1);
 B(i,i+1) = br(i);
 B(i+1,i) = cr(i+1);
end

eee = eig(A,B);
eee(1:2)

9/4*1.866^2

%-----
% done
%-----
