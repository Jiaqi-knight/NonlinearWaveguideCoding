
clear all
close all

%====================================================
% proriol_ortho
%
% Verify the orthogonality of the Proriol polynomials
% by integrating over the triangle using a quadrature
%====================================================

i=2; j=1;   % examples
k=2; l=1;

%--------------------
% triangle quadrature
%--------------------

NQ=13

[xiq, etq, wq] = gauss_trgl(NQ);

%----------------------------
% integrate over the triangle
%----------------------------

olok = 0;

for q=1:NQ

  xi=xiq(q); eta=etq(q); xip = 2*xi/(1-eta)-1; etap = 2*eta-1;

  proriol_ij = jacobi(0,0,i,xip) *(1-eta)^i * jacobi(2*i+1,0,j,etap) ;
  proriol_kl = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;

  integrand =  proriol_ij * proriol_kl;
  olok = olok + integrand * wq(q);

end

olok = 0.5*olok;  % quadrature weights add up to one

olok

%--- exact for i=k and j=l:
%
% Gij = 0.5/((2*i+1)*(i+j+1))
%
%---

%-----
% done
%-----
