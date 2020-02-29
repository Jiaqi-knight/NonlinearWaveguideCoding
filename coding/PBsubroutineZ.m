function PBsubroutine();
  opengl('save','software')

  format long e
  format long
  format compact
  fzoptions = optimset('display','off');
  warning off MATLAB:MKDIR:DirectoryExists
  clear all
  omega = 20;
  k0    = -17+70i;
  Z1    = 0.4+1.5000i;
 %Z1    = 1e9 ;
  Z0    = 0.4+1.5000i;
 %Z0    = 1e9 ;
  Np    = 3;
  yp    = [0:Np-1]/(Np-1);       % note: calculation interval is [yp(1),yp(Np)], Np = length(yp)
  nn    =  60;
  NN    = 160;
 %tic
  [k,fp,dfp,count] = PB2D(k0,yp,omega,@U0,@C0,nn,NN,Z0,Z1);
 %toc
  display([sprintf('k = %+1.12G%+1.12Gi, ',real(k),imag(k)),sprintf(' count = %d',count)])
end

function u = U0(y);
         tau   = 0.5;
         sigma = 0.3;
         u = tau + sigma*y  ;
end
function u = C0(y);
         u = ones(1,length(y));% 0.5+0.5*4*y.*(1-y);
end

function [k,fp,dfp,cnt] = PB2D(k0,yp,om,uu,cc,n,N,z0,z1)
                          % k0 initial estimate of eigen value k
                          % fp, dfp the found function and its derivative along the plot points
                          % cnt number of iterations
                          % interval [a,b], a = yp(1), b = yp(Np)
                          % plot points yp, Np = length(yp)
                          % om frequency
                          % uu and cc are function handlers to U0 en C0
                          % nn number of basis functions
                          % N  number of Gauss Legendre points
                          % z0, z1 impedancy at y=a, resp y=b
 %
 %global
  Np       = length(yp);   % number of plot points, >1
  a        = yp(1);
  b        = yp(Np);
 [ygl,Wgl] = legzo(N,a,b);                   % GL points in [a,b]
  ytot     = [ygl,yp];                       % N+Np = total number y-points
  N0       = N+1;                            % index boundary point a
  N1       = N+Np;                           % index boundary point b
  errM     = 1e-10;

  %%%%%%%%%%%
  %%%%%%%%%%% start Make Chebyshev polynomials
  %%%%%%%%%%%
     t      = (2*ytot -b-a)/(b-a);     % we use the GL and plot-points as if on [-1,1]
     T      = zeros(n+1,N+Np);
    DT      = zeros(n+1,N+Np);
     Tgl    = zeros(n+1,N);
    DTgl    = zeros(n+1,N);
     Tp     = zeros(n+1,Np);
    DTp     = zeros(n+1,Np);
     T0     = zeros(n+1,1);
     T1     = zeros(n+1,1);

     U      = zeros(n+1,N+Np);
     T(1,:) = ones(   1,N+Np);         %  T0 = 1
     U(1,:) = ones(   1,N+Np);         %  U0 = 1
    DT(1,:) = zeros(  1,N+Np);         % DT0 = 0
     T(2,:) = t;
     U(2,:) = 2*t;
    DT(2,:) = ones(   1,N+Np);
   for jj = 3:n+1
        T(jj,:) = 2*t.*T(jj-1,:) - T(jj-2,:);  % Chebyshev T  of order jj-1 ( so T(1,:) is T_0(x) )
        U(jj,:) = 2*t.*U(jj-1,:) - U(jj-2,:);  % Chebyshev U
       DT(jj,:) = (jj-1)*U(jj-1,:)          ;  % Chebyshev T'
   end
   DT  =  DT*2/(b-a);          % correct derivative
   T0  =  T(:,N0);
   T1  =  T(:,N1);
   for jj = 1:n+1
          Tp(jj,:) =  T(jj,N+1:N+Np);      % phi_n = T_{n}
         DTp(jj,:) = DT(jj,N+1:N+Np);      % phi_n = T_{n}
         Tgl(jj,:) =  T(jj,1:N);
        DTgl(jj,:) = DT(jj,1:N);
   end
   TTgl   = zeros(n+1,n+1,N);
  DTTgl   = zeros(n+1,n+1,N);
  for jj = 1:n+1
     for mm = jj:n+1
           TTgl(jj,mm,:) =  Tgl(jj,:).*Tgl(mm,:);    % in view of intensive re-use while making M
          DTTgl(jj,mm,:) = DTgl(jj,:).*DTgl(mm,:);   % careful: maybe apply reshape
     end
  end
  %%%%%%%%%%%
  %%%%%%%%%%% end   Make Chebyshev polynomials
  %%%%%%%%%%%
  mu0  = k0;
  rest = 1;
  cnt  = 0;
  ugl  = uu(ygl)       ;
  cgl  = cc(ygl)             ;
 [M0,DM0] = MakeM(mu0);
  while (rest> errM & cnt<50 &abs(mu0)<1000)
    cnt     = cnt+1;
   [V,E]    = eig(M0,DM0);
    E       = diag(E);
   [q,p]    = min(abs(E));
    mu0     = mu0-E(p);
   [M0,DM0] = MakeM(mu0);
    rest    = norm(M0*V(:,p));
  end
           k = mu0;
          fp = V(:,p).'*Tp;
         dfp = V(:,p).'*DTp;
   %
   % normalisering van fp. Let op: dit zou tot axiaal niet-gladde eigenfuncties kunnen leiden
   %
 [fmax,ymax] = max(abs(fp));                       %
         Afp = exp(-i*angle(fp(ymax)))/fmax;
          fp = Afp*fp                         ;    % fp is reel in geval van harde wanden
         dfp = Afp*dfp                        ;    %
  if real(fp(1)) <0, fp=-fp; dfp=-dfp; end         % dit gaat automatisch goed als fp is complex


  function [M,DM] = MakeM(mu)    % heeft legzo nodig
     Nyg = [1,N];                % dit is natuurlijk N, maar vanwege reshape moet het [1,N] zijn
     Om1 = (om - mu*ugl)./cgl;
     Om2 = Om1.*Om1;
     Om3 = Om2.*Om1;
     M   = zeros(n+1);
     DM  = zeros(n+1);
     for jj = 1:n+1
       for mm = jj:n+1
           int1gl   =         -reshape(DTTgl(jj,mm,:),Nyg)./Om2 + (1-mu^2./Om2).*reshape(TTgl(jj,mm,:),Nyg);
           int2gl   = -2*(ugl.*reshape(DTTgl(jj,mm,:),Nyg)      +        om*mu.*reshape(TTgl(jj,mm,:),Nyg))./cgl./Om3;
           M(jj,mm) = sum(Wgl.*int1gl);
          DM(jj,mm) = sum(Wgl.*int2gl);
           if abs(z0)<1E9
                  M(jj,mm) = M(jj,mm) +  T0(jj)*T0(mm)/(i*om*z0);
           end
           if abs(z1)<1E9
                  M(jj,mm) = M(jj,mm) +  T1(jj)*T1(mm)/(i*om*z1);
           end
           if jj ~= mm, M(mm,jj) = M(jj,mm); DM(mm,jj) = DM(jj,mm); end
       end
     end
  end  % eind MakeM

end  % eind PB2D


function [x,w]=legzo(n, a, b)
%       =========================================================
%       Purpose : Compute the zeros of Legendre polynomial Pn(x)
%                 in the interval [a,b], and the corresponding
%                 weighting coefficients for Gauss-Legendre
%                 integration
%       Input :   n    --- Order of the Legendre polynomial
%                 a    --- Lower boundary (optional)
%                 b    --- Upper boundary (optional)
%       Output:   x(n) --- Zeros of the Legendre polynomial
%                 w(n) --- Corresponding weighting coefficients
%       =========================================================
  if nargin == 1
      a = -1;
      b =  1;
  end;
  x = zeros(1, n);
  w = zeros(1, n);
  m = (n+1)/2;
  h = b-a;

  for ii=1:m
      z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
      z1 = z+1;
      while abs(z-z1)>eps
          p1 = 1;
          p2 = 0;
          for jj = 1:n
              p3 = p2;
              p2 = p1;
              p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
          end
          pp = n*(z*p1-p2)/(z^2-1); % The L.P. derivative.
          z1 = z;
          z = z1-p1/pp;
      end
      x(ii) = z; % Build up the abscissas.
      x(n+1-ii) = -z;
      w(ii) = h/((1-z^2)*(pp^2)); % Build up the weights.
      w(n+1-ii) = w(ii);
  end

  if a ~= -1 || b ~= 1
      x = (x+1)*(h/2) + a;
  end
end  % einde legzo
