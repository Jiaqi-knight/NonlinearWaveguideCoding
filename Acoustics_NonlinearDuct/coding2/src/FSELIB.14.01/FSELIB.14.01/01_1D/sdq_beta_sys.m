function [ap,bp,cp,dp,ep,b] = sdq_beta_sys (ne,xe,beta,q0,fL,k,s)

%====================================================
% Compact assembly of the pentadiagonal linear system 
% for one-dimensional steady diffusion
% with quadratic elements (sdq)
% and arbitrary position of the interior node
% determined by the parameter "beta"
%====================================================

%-------------
% element size
%-------------

for l=1:ne
  h(l) = xe(l+1)-xe(l);
end

%------------------------------
% number of unique global nodes
%------------------------------

ng = 2*ne+1;

%------------------------------------
% initialize the pentadiagonal matrix
%------------------------------------

ap = zeros(ng,1); bp = zeros(ng,1); cp = zeros(ng,1);
dp = zeros(ng,1); ep = zeros(ng,1);

%-------------------------------
% initialize the right-hand side
%-------------------------------

b = zeros(ng,1); b(1) = q0/k;

%-----------------------
% loop over all elements
%-----------------------

for l=1:ne

%----
% diffusion matrix analytically:
%----

  cf = 1.0/(3.0*h(l));
  A(1,1) =  1/(1+beta)^2 * (4+3*(1+beta)^2 ) * cf;
        A(1,2) = -8/( (1+beta)*(1-beta^2) ) * cf;
            A(1,3) =  1/(1-beta^2) * (4-3*(1-beta^2))*cf;
  A(2,1) = A(1,2); 
           A(2,2) = 16/(1-beta^2)^2 * cf;
           A(2,3) = -8/((1-beta)*(1-beta^2)) * cf;
  A(3,1) = A(1,3);
  A(3,2) = A(2,3);
  A(3,3) = 1/(1-beta)^2 * (4 + 3*(1-beta)^2 )*cf;

%----
% alternative:
%
% diffusion matrix by the Vandermond matrix
% transformation
%----

% for beta=0:

 cf = 1.0/(3.0*h(l));
 A(1,1) = 7.0*cf; A(1,2) = -8.0*cf; A(1,3) = cf;
 A(2,1) = A(1,2);    A(2,2) = 16.0*cf; A(2,3) = A(1,2);
 A(3,1) = A(1,3);    A(3,2) = A(2,3);     A(3,3) = A(1,1);

%----
% vandermond matrix
%----

 vd(1,1) = 1; vd(1,2) = 0; vd(1,3) = 0;
 vd(2,1) = beta*(beta-1)/2; vd(2,2) = 1-beta^2;
 vd(2,3) = beta*(beta+1)/2;
 vd(3,1) = 0; vd(3,2)=0; vd(3,3)=1;

 invvd = inv(vd);

%----
% transform
%----

 A = invvd'*A*invvd;

%----------------------------------------------
% mass matrix by the 4-point lobatto quadrature:
%----------------------------------------------

% lobatto base points and weights

  xi(1)=-1.0; xi(2) =-1/sqrt(5); xi(3)=1/sqrt(5); xi(4)=1.0;
  w(1)=1/6; w(2)=5/6; w(3)=5/6; w(4)=1/6;

% initialize:

  for i=1:3
    for j=1:3
      B(i,j) = 0;
    end
  end

% lobatto quadrature:

  for q=1:4

    psi(1) = (xi(q)-beta)*(xi(q)-1)/(2*(1+beta));
    psi(2) = (1-xi(q)^2)/(1-beta^2);
    psi(3) = (xi(q)-beta)*(xi(q)+1)/(2*(1-beta));

    B(1,1) = B(1,1) + psi(1)*psi(1)*w(q);
    B(1,2) = B(1,2) + psi(1)*psi(2)*w(q);
    B(1,3) = B(1,3) + psi(1)*psi(3)*w(q);
    B(2,2) = B(2,2) + psi(2)*psi(2)*w(q);
    B(2,3) = B(2,3) + psi(2)*psi(3)*w(q);
    B(3,3) = B(3,3) + psi(3)*psi(3)*w(q);

  end

 B(2,1)=B(1,2); B(3,1)=B(1,3); B(3,2)=B(2,3); 

 B = h(l)/2.0*B;

%---------------
% mass matrix by the Vandermond matrix
% transformation
%---------------

% for beta=0:

  cf = h(l)/30.0;
  B(1,1) = 4.0*cf; B(1,2) =  2.0*cf;  B(1,3) = -cf;
  B(2,1) = B(1,2); B(2,2) = 16.0*cf;  B(2,3) = B(1,2);
  B(3,1) = B(1,3); B(3,2) = B(2,3);   B(3,3) = B(1,1);

%----
% transform
%----

  B = invvd'*B*invvd;

%---------
% assemble:
%---------

  cl1 = 2*l-1; cl2 = 2*l; cl3 = 2*l+1;

  ap(cl1) = ap(cl1) + A(1,1);
  bp(cl1) = bp(cl1) + A(1,2);
  cp(cl1) = cp(cl1) + A(1,3);

  dp(cl2) = dp(cl2) + A(2,1);
  ap(cl2) = ap(cl2) + A(2,2);
  bp(cl2) = bp(cl2) + A(2,3);

  ep(cl3) = ep(cl3) + A(3,1);
  dp(cl3) = dp(cl3) + A(3,2);
  ap(cl3) = ap(cl3) + A(3,3);

  b(cl1) = b(cl1) + (B(1,1)*s(cl1) + B(1,2)*s(cl2) + B(1,3)*s(cl3))/k;
  b(cl2) = b(cl2) + (B(2,1)*s(cl1) + B(2,2)*s(cl2) + B(2,3)*s(cl3))/k;
  b(cl3) = b(cl3) + (B(3,1)*s(cl1) + B(3,2)*s(cl2) + B(3,3)*s(cl3))/k;

end

%----------------------------------
% implement the Dirichlet condition
%----------------------------------

b(ng-2) = b(ng-2) - cp(ng-2)*fL;
b(ng-1) = b(ng-1) - bp(ng-1)*fL;

%-----
% done
%-----

return;
