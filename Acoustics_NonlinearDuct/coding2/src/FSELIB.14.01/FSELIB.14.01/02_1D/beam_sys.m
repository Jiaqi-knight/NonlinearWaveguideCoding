function [stiff] = beam_sys (ne,xe)

%=============================================
% FSELIB
%
% Assembly of the global stiffness matrix
% for beam bending with cubic Hermitian elements 
%
% notation:
%
% ne:    number of elements
% xe:    position of element nodes
% stiff: global stiffness matrix
%=============================================

ng = ne+1; % number of unique global nodes

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%---------------------------------------
% initialize the global stiffness matrix
%---------------------------------------

stiff = zeros(2*ng,2*ng); 

%-----------------------
% loop over the elements
% to compute the element stiffness matrix
%-----------------------

for l=1:ne

   esm(1,1) =  12.0; esm(1,2) = 6.0*h(l);
   esm(1,3) = -12.0; esm(1,4) = 6.0*h(l);
   esm(2,1) = esm(1,2); esm(2,2) = 4.0*h(l)^2;
   esm(2,3) = -6.0*h(l); esm(2,4) = 2.0*h(l)^2;
   esm(3,1) = esm(1,3); esm(3,2) = esm(2,3); esm(3,3) = 12.0;
   esm(3,4) = -6.0*h(l);
   esm(4,1) = esm(1,4); esm(4,2) = esm(2,4); esm(4,3) = esm(3,4);
   esm(4,4) = 4.0*h(l)^2;

   esm = esm/h(l)^3;

   i1 = 2*(l-1);

   for i=1:4
     for j=1:4
       stiff(i1+i,i1+j) = stiff(i1+i,i1+j) + esm(i,j);
      end
   end

end

%-----
% done
%-----

return;
