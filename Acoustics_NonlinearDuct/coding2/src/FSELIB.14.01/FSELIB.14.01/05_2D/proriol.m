function proriol = proriol (k,l,xi,eta)
                                                                                
%========================================
% evaluation of the kl Proriol polynomial
%========================================

xip = 2*xi/(1-eta)-1; etap = 2*eta-1;
proriol = jacobi(0,0,k,xip) *(1-eta)^k * jacobi(2*k+1,0,l,etap) ;

%-----
% done
%-----

return
