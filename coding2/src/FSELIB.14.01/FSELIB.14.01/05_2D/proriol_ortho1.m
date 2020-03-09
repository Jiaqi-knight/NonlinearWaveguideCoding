
clear all
close all

%======================================================
% proriol_ortho1
%
% Verify the orthogonality of the Proriol polynomials
% by integrating over the triangle
% using an integration rule
%======================================================

ipr=4; jpr=1;   % examples
kpr=0; lpr=0;

%--------------------
% integration weights
%--------------------

mint = 128;

Dxi = 1.0/mint; fc = Dxi^2; wsum = 0.0;

Ic=0;
Ic=Ic+1; w(Ic)=fc/4.0; wsum = wsum+w(Ic);

for j=2:mint-1
  Ic=Ic+1;w(Ic)=0.5*fc;wsum=wsum+w(Ic);
end

Ic=Ic+1;w(Ic)=5.0/12.0*fc;wsum=wsum+w(Ic);
Ic=Ic+1;w(Ic)=1.0/6.0*fc;wsum=wsum+w(Ic);

for i=2:mint-1
 Ic=Ic+1; w(Ic) = 0.5*fc;wsum=wsum+w(Ic);
 for j=2:mint-i
    Ic=Ic+1; w(Ic) = fc; wsum=wsum+w(Ic);
 end
 Ic=Ic+1; w(Ic)=11.0/12.0*fc; wsum=wsum+w(Ic);
 Ic=Ic+1; w(Ic)=7.0/12.0*fc; wsum=wsum+w(Ic);
end

Ic=Ic+1;w(Ic)=5.0/12.0*fc; wsum=wsum+w(Ic);
Ic=Ic+1;w(Ic)=7.0/12.0*fc; wsum=wsum+w(Ic);

Ic=Ic+1;w(Ic)=1.0/6.0*fc; wsum=wsum+w(Ic);

% wsum  % should be 0.5

%----------------------------
% integrate over the triangle
%----------------------------

olok = 0;

Ic = 0;

for i=1:mint+1   % loop over integration points
 for j=1:mint+2-i

  xi  = (i-1.0)*Dxi; eta = (j-1.0)*Dxi;
  if(eta>0.9999999) eta=0.9999999; end

  xip = 2*xi/(1-eta)-1; etap = 2*eta-1;

  proriol_ij = jacobi(0,0,ipr,xip) *(1-eta)^ipr * jacobi(2*ipr+1,0,jpr,etap) ;
  proriol_kl = jacobi(0,0,kpr,xip) *(1-eta)^kpr * jacobi(2*kpr+1,0,lpr,etap) ;

  integrand =  proriol_ij * proriol_kl;
  Ic=Ic+1;
  olok = olok + integrand * w(Ic);

end
end

olok

%--- exact:
%
% Gij = 0.5/((2*ipr+1)*(ipr+jpr+1))
%
%---

%-----
% done
%-----
