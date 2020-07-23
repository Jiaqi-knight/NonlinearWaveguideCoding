function PBGalerkin2DL_Ztrace();
  clear all
  opengl('save','software')

  format long e
  format long
  format compact
  global omega   tau sigma n y ygl ugl cgl N N0 N1 TTgl DTTgl T0 T1 Wgl dispth
  fzoptions = optimset('display','off');
  warning off MATLAB:MKDIR:DirectoryExists

  EigInitVal  = 0;
  plotvanfile = 0;
  volgsigma   = 0;
  onlysigma   = 0;
  smallfrqapp = 0;
  makemovie   = 1;         % only if plotoffile=1
  profiel     = 'par';     % 'flat', 'linear' or 'parabolic' or 'bou' .....
  dispth      = 0.010;
  dpsbug      = 0;         % dpsimplify bug
  softside    = 2;         % =1 if y=1, =0  if y=0 , =2 if both the same impedance

  RZ     = 0.5;      % Re(Z);
  IZ     = 60 ;
  dZ     = 0.020;
  Zsign  = -1   ;    % + = opward, - = down
  n      =  20;      % highest order basis functions. Matrix is (n+1)x(n+1)
  N      =  40;      % nr of nodes for Gauss-Legendre
  N0     = 1;        % because of BC at 0, then N GL-points, and then the N+2-th point N1 at 1
  N1     = N+2;      % because of BC at 1
  omega  =  5   ;

  tau    = 0.5 ;
  sigma1 = 0.5 ;
  Nsg    = 10  ;
  Nsg1   = Nsg+1;
  Merr   = 1e-9;
  mu0brk = 1000;
  dperr  = 1e-6;   % 1e-6
  Rs     = 10 ;
  Nmu    = 10  ;             % even, to get the same number of left and right modes
  Nmu5   = floor(Nmu/2);
  if (sigma1 ~= 0) & (profiel == 'fla')
        error 'sigma<>0 and flat profile'
  end
  if sigma1 == 0
      Nsg     = 0;
      ssgma   = 0.0;
      profile = 'par';
  else
      ssgma  = sigma1*[0:Nsg]./Nsg;   % Nsg intervals, Nsg+1 end- and inter-points for sigma = 0 .. sigma_eind
  end
 %Zarray   = RZ + Zsign*[-1E10,[-IZ:dZ:-IZ+5*dZ]]*i;  % if Zsign=+1
 %Zarray   = RZ + Zsign*[   [IZ-5*dZ:dZ:IZ],1E10]*i;  % if Zsign=-1
  Zarray   = RZ + Zsign*[-1E10,[-IZ:dZ:IZ],1E10]*i;  % the last value is problematic if the mode disappears to infinity. Therefore we don't use it
  NZa      = length(Zarray);
  xmin     = -25 ;
  xmax     =  20;
  ymin     = -25   ;
  ymax     =  25   ;

  Nmovie   =   1;
  IZstart  =  15;             % absolute value is enough
  IZstart  = -Zsign*IZstart;

  %%%%%%%%%%%%%
  % Boundingbox correction
  dx1 = -5;
  dx2 =  5;
  dy1 = -5;
  dy2 =  5;
  %%%%%%%%%%%%%%

  %%%%%%%%%%
  if ~plotvanfile
   if EigInitVal
       % om=5, tau=0.5, sig =0.5:
       mu00 = [+3.703,+2.583,-2.119-3.322,-2.222-8.194,-2.191-12.08,-2.191+12.08,-2.222+8.194,-2.119+3.322,-5.756,-7.846];
      %mu00  = sort(mu00,'ascend');
       Nmu   = length(mu00);
       display([sprintf('Nmu = %d',Nmu)]);
       Nsg   = 0;
       sigma = sigma1;
   else
       Nmu  = 2*Nmu5;
       mu00 = zeros(1,Nmu);
       for j=1:Nmu5
           mu00(      j) = (-omega*tau-i*sqrt(-omega^2+(1-tau^2)*((j-1)*pi)^2))/(1-tau^2); % harde wand op y=0 en y=1
           mu00(Nmu+1-j) = (-omega*tau+i*sqrt(-omega^2+(1-tau^2)*((j-1)*pi)^2))/(1-tau^2); %
       end
   end
   mu     = zeros(NZa,Nmu);
  end

  %%%%%%%%%%
  % Naamgeving
  if omega <1, casenm=strcat('w',num2str(floor(10*mod(omega,1)),'.%d')); else, casenm=strcat('w',num2str(omega,'%02.f')); end

  casenm = strcat(casenm,'-t',num2str(100*tau,'%02.f'),'s',num2str(100*sigma1,'%02.f'));
  switch softside
      case 2,    casenm   = strcat(casenm,'-B',num2str(RZ,'%0.2f'));
      case 1,    casenm   = strcat(casenm,'-R',num2str(RZ,'%0.2f'));
      case 0,    casenm   = strcat(casenm,'-L',num2str(RZ,'%0.2f'));
      otherwise, display('verkeerde softside'), return
  end
  casedir = strcat('.\data\',casenm,'\');
  casenm  = strcat(casenm,'-n',int2str(n),'N',int2str(N));
  switch profiel(1:3)
    case 'fla', casenm = strcat(casenm,'-fla');
    case 'lin', casenm = strcat(casenm,'-lin');
    case 'par', casenm = strcat(casenm,'-par');
    case 'bou', casenm = strcat(casenm,'-bou');
                casenm = strcat(casenm,'-d',num2str(100*dispth,'%02.f'));
    otherwise , error 'Invalid option for profiel', return
  end
  savefile = strcat(casedir,casenm,'.dat');
  if exist(savefile, 'file') & ~plotvanfile
     warningMessage = sprintf('Warning: file exists:\n%s', savefile);
     uiwait(msgbox(warningMessage));
     return
  end
  if ~exist(savefile, 'file') & plotvanfile
     warningMessage = sprintf('Warning: file does not exists:\n%s', savefile);
     uiwait(msgbox(warningMessage));
     return
  end
  parmfile = strcat(casedir,casenm,'.par');
  epsfile  = strcat(casedir,casenm,'.eps');
  jpgfile  = strcat(casedir,casenm,'.jpg');
  pngfile  = strcat(casedir,casenm,'.png');
  sgmfile  = strcat(casedir,casenm,'.sgm');
  giffile  = strcat(casedir,casenm,'.gif');
  movfile  = strings(1,Nmovie);
  for j=1:Nmovie
    movfile(j)  = strcat(casedir,casenm,num2str(j,'mv%02d'),'.png');
  end
  mkdir(casedir);
  %%%%%%%%%%%
  if ~plotvanfile
     fid = fopen(parmfile,'w');
     cl  = onCleanup(@()fclose(fid));
     fprintf(fid,'case name  = %s\n', casenm);
     fprintf(fid,'m-file     = %s.m\n',mfilename);  % mfilename geeft geen extensie

     fprintf(fid,'omega      = %+1.12G\n', omega);
     fprintf(fid,'tau        = %+1.12G\n', tau);
     fprintf(fid,'sigma1     = %+1.12G\n', sigma1);
     fprintf(fid,'softside   = %1d\n'    , softside);
     fprintf(fid,'RZ         = %+1.12G\n', RZ);
     fprintf(fid,'n          = %1d\n'    , n      );
     fprintf(fid,'N          = %1d\n'    , N      );
     fprintf(fid,'profiel    = %s\n'     , profiel);
     fprintf(fid,'IZ         = %+1.12G\n', IZ);
     fprintf(fid,'dZ         = %+1.12G\n', dZ);
     fprintf(fid,'Zsign      = %+1d\n'   , Zsign  );
     fprintf(fid,'Nsg        = %1d\n'    , Nsg    );
     fprintf(fid,'Rs         = %+1.12G\n', Rs   );
     fprintf(fid,'Nmu        = %1d\n'    , Nmu    );
     fprintf(fid,'EigInitVal = %1d\n'    , EigInitVal);
     if EigInitVal
      if Nmu==1
        fprintf(fid,'mu00       = [%+1.12G%+1.12Gi]\n', real(mu00(1)),imag(mu00(1)));
      else
        fprintf(fid,'mu00       = [%+1.12G%+1.12Gi, ', real(mu00(1)),imag(mu00(1)));
        for j=2:Nmu-1
             fprintf(fid,'%+1.12G%+1.12Gi, ', real(mu00(j)),imag(mu00(j)));
        end;
        fprintf(fid,'%+1.12G%+1.12Gi]\n', real(mu00(Nmu)),imag(mu00(Nmu)));
      end
     end
     if profiel(1:3) == 'bou'
        fprintf(fid,'dispth     = %+1.12G\n', dispth);
     end;
     fprintf(fid,'xmin       = %+1.12G\n', xmin );
     fprintf(fid,'xmax       = %+1.12G\n', xmax );
     fprintf(fid,'ymin       = %+1.12G\n', ymin );
     fprintf(fid,'ymax       = %+1.12G\n', ymax );
     %fclose(fid);
  end
  %%%%%%%%%%%
  if plotvanfile % to make a plot from earlier saved points in .dat-file
     %fid = fopen(savefile,'r');
     %cl  = onCleanup(@()fclose(fid));
      A   = dlmread(savefile);
      [NZa,Nmu] = size(A);
      Nmu = (Nmu-2)/2;
      Z   = zeros(NZa,1);
      Z   = A(:,1)+A(:,2);
      mu  = zeros(NZa,Nmu);
      for j=1:Nmu
        mu(:,j) = A(:,2*j+1)+A(:,2*j+2);
      end
      %fclose(fid);
      sigma = sigma1;
      if smallfrqapp
         kappr(1,:) = SmallFreqApprox(RZ+Zsign*1E10*i,omega,tau,sigma);
         for nz = 2:NZa-1
              kappr(nz,:) = SmallFreqApprox(Z(nz),omega,tau,sigma);
         end
         kappr(NZa,:) = kappr(1,:);
      end
      figure(100)
      clf(100,'reset')
      hold on
      grid on
      axis([xmin xmax ymin ymax]);
      plot([xmin,xmax],[0,0],'--','Color',[0.5,0.5,0.5],'Linewidth',0.6)
      plot([0,0],[ymin,ymax],'--','Color',[0.5,0.5,0.5],'Linewidth',0.6)
      if volgsigma
       plot(real(mu00),imag(mu00),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',1)
      end
      plot(real(mu(1,:)),imag(mu(1,:)),'o','MarkerEdgeColor','[0,0,0]','MarkerFaceColor','[0,0,0]','MarkerSize',1.2)
      for j=1:Nmu
        switch (j<=Nmu5)*(mod(j-1,8)+1)+(j>Nmu5)*(mod(Nmu-j,8)+1)
          case 1,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-k','LineWidth',                          0.5)
          case 2,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-b','LineWidth',                          0.5)
          case 3,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-r','LineWidth',                          0.5)
          case 4,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-','Color',[0,0.7,0.7],'LineWidth',       0.5)          %donkercyan
          case 5,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-m','LineWidth',                          0.5)
          case 6,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-','Color',[200/255,173/255,      0],'LineWidth',0.5)   %okergeel
          case 7,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-','Color',[ 59/255,116/255, 75/255],'LineWidth',0.5)   %grijsgroen
          case 8,ps = [real(mu(:,j)), imag(mu(:,j))];plot(ps(:,1),ps(:,2),'-','Color',[156/255,130/255,189/255],'LineWidth',0.5)   %grijspaars
        end
      end
      if smallfrqapp
       plot(real(kappr(:,1)),imag(kappr(:,1)),'--','Color',[200/255,173/255,0],'LineWidth',0.5)
       plot(real(kappr(:,2)),imag(kappr(:,2)),'--','Color',[200/255,173/255,0],'LineWidth',0.5)
      end
      if profiel(1:3) == 'bou'
         title(['\omega = ',sprintf('%0.1f',omega),  ...
                                                     ...
                ', \tau = ',sprintf('%0.2f',tau),    ...
                ', \sigma = ',sprintf('%0.2f',sigma),...
                ', \delta = ',sprintf('%0.3G',dispth),...
                ', Re(Z) = ',sprintf('%0.2f',RZ)     ...
                ],'fontweight','normal','fontsize',10)
      else
         title(['\omega = ',sprintf('%0.1f',omega),  ...
                                                     ...
                ', \tau = ',sprintf('%0.2f',tau),    ...
                ', \sigma = ',sprintf('%0.2f',sigma),...
                ', Re(Z) = ',sprintf('%0.2f',RZ)     ...
                ],'fontweight','normal','fontsize',10)
      end
      %%%%%%%%%
      box on
      set(gcf,'color',[1 1 1])
      set(gcf,'inverthardcopy','off')
      set(gcf,'papertype','a4letter')
      set(gcf,'paperunits','centimeters')
      set(gcf,'paperposition',[0  0  10  10 ])
      set(gcf,'clipping','off')
      %eval(['print -depsc -cmyk ' epsfile]);
      print(epsfile,'-depsc','-cmyk','-r600')
      fixepsbbox(epsfile,dx1,dy1,dx2,dy2)
     %print2eps(epsfile)
      print(jpgfile,'-djpeg','-r600')
      print(pngfile,'-dpng' ,'-r600')
      if makemovie
         set(gcf,'paperposition',[0  0  10  10 ])
         set(gca,'FontSize',5 );
         ftot = gcf;
         savefig(ftot,'tempfig','compact');
         [U,I1] = min(abs(imag(Z(:))-IZstart));
         [U,I2] = min(abs(imag(Z(:))+IZstart));
         Istep  = round((I2-I1)/(Nmovie-1));
         Imovie = [I1:Istep:I1+(Nmovie-1)*Istep];
         Zmovie = Z(Imovie);
         mumvie = mu(Imovie,:);
         for j=1:Nmovie
             openfig('tempfig');
             title(['\omega = ',sprintf('%0.1f',omega),  ...
                                                         ...
                    ', \tau = ',sprintf('%0.2f',tau),    ...
                    ', \sigma = ',sprintf('%0.2f',sigma),...
                    ', Z = ',sprintf('%0.2f',RZ),        ...
                      sprintf('%+06.2f',imag(Zmovie(j))),'i',     ...
                    ],'fontweight','bold','fontsize',6,'FontName','FixedWidth')
             hold on
             grid on
             box on
             g0=plot(real(mumvie(j,:)),imag(mumvie(j,:)),'ro');
             set(g0,'LineWidth',0.5,'Markersize',2, 'MarkerEdgeColor','r','MarkerFaceColor','r');
             g1=plot(xmax,imag(Zmovie(j)),'ko');
             set(g1,'LineWidth',0.5,'Markersize',2, 'MarkerEdgeColor','k','MarkerFaceColor','k');
             print(movfile(j),'-dpng' ,'-r300')
             hold off
             close(gcf);
             pause(0.10);
         end
         if ~isfile(giffile)
            fopen(giffile,'w');
         end
      end
      return
  end % plotvanfile

  %%%%%%   make Chebyshev and derivatives
  MakeT2D(0,1);     % interval is [0,1]
  %%%%%%%%%%%  sigma iteration
  tic
  if Nsg == 0
     for j=1:Nmu
       mu0 = mu00(j);
       sigma = sigma1;
       Z  = Zarray(1);
       rest = 1;
       cnt  = 0;
       [M0,DM0]   = MakeM(mu0);
       while (rest>Merr & cnt<50)
         cnt      = cnt+1;
         [V,E]    = eig(M0,DM0);
         E        = diag(E);
         [q,p]    = min(abs(E));
         mu0      = mu0-E(p);
         [M0,DM0] = MakeM(mu0);
         rest     = norm(M0*V(:,p));     % note: eigenvector V is normalised
       end
       mu(1,j) = mu0;
     end;  % over j  (Nmu)
  else
    tic
    mus = zeros(Nsg1,Nmu);
    mus(1,:) = mu00;
    Zs  = 1./(1e-10+1/Rs*4*[0:Nsg].*[Nsg:-1:0]./Nsg^2);
    for j=1:Nmu
      fprintf('j=%2i, ',j);
      if volgsigma,
          fprintf('\n');
      else
         if mod(j,10)==0, fprintf('\n'), end
      end
      Z     = Zs(1);
      sigma = ssgma(1);
      mu0  = mu00(j);
      cnt  = 0;
      if volgsigma
         switch true
             case abs(Z)<9.9999,                   display([sprintf('s=%0.3f, Z = %0.2f%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
             case abs(Z)>=9.9999 && abs(Z)<99.999, display([sprintf('s=%0.3f, Z = %0.1f%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
             case abs(Z)>=99.999 && abs(Z)<999.99, display([sprintf('s=%0.3f, Z = %0.f.%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
             case abs(Z)>=999.99 && abs(Z)<9999.9, display([sprintf('s=%0.3f, Z = %4.f%+4.fi,   r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
             case abs(Z)>=9999.9 && abs(Z)<1E10,   display([sprintf('s=%0.3f, Z = %2.0fe%+d%+2.0fe%+di, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z)/10^(floor(log10(abs(real(Z))))),floor(log10(abs(real(Z)))),imag(Z)/10^(floor(log10(abs(imag(Z))))),floor(log10(abs(imag(Z)))), cnt,real(mu0),imag(mu0))])
             otherwise,                            display([sprintf('s=%0.3f, Z = inf +inf.i, r=%1i, mu0 = %4.12f %+4.12fi',sigma,cnt,real(mu0),imag(mu0))])
         end
      end
      for ns = 2:Nsg1      % eerst Nsg om sigma te varieren
          Z     = Zs(ns);
          sigma = ssgma(ns);
          rest = 1;
          cnt  = 0;
          [M0,DM0]   = MakeM(mu0);
          while (rest>1e-10 & cnt<50)
            cnt      = cnt+1;
            %display([rest,cnt,norm(M0),norm(DM0)])
            [V,E]    = eig(M0,DM0);
            E        = diag(E);
            [q,p]    = min(abs(E));
            mu0      = mu0-E(p);
            [M0,DM0] = MakeM(mu0)             ;
            rest     = norm(M0*V(:,p));
          end
          mus(ns,j) = mu0;

          if volgsigma
             switch true
                 case abs(Z)<9.9999,                   display([sprintf('s=%0.3f, Z = %0.2f%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
                 case abs(Z)>=9.9999 && abs(Z)<99.999, display([sprintf('s=%0.3f, Z = %0.1f%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
                 case abs(Z)>=99.999 && abs(Z)<999.99, display([sprintf('s=%0.3f, Z = %0.f.%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
                 case abs(Z)>=999.99 && abs(Z)<9999.9, display([sprintf('s=%0.3f, Z = %4.f%+4.fi,   r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
                 case abs(Z)>=9999.9 && abs(Z)<1E10,   display([sprintf('s=%0.3f, Z = %2.0fe%+d%+2.0fe%+di, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z)/10^(floor(log10(abs(real(Z))))),floor(log10(abs(real(Z)))),imag(Z)/10^(floor(log10(abs(imag(Z))))),floor(log10(abs(imag(Z)))), cnt,real(mu0),imag(mu0))])
                 otherwise,                            display([sprintf('s=%0.3f, Z = inf +inf.i, r=%1i, mu0 = %4.12f %+4.12fi',sigma,cnt,real(mu0),imag(mu0))])
             end
          end
          mu0 = 2*mus(ns,j)-mus(ns-1,j); % ns=1 komt niet voor want die kennen we exact
      end % over ns
    end % over Nmu
    mu(1,:) = mus(Nsg1,:);
    telapseds = toc;
  end % if Nsg > 0
  toc
  %%%%%%%%%%
  %%%%%%%%%%%%%
  if volgsigma
     figure(101)
     clf(101,'reset')
     hold on;
     grid on;
     axis([xmin xmax ymin ymax]);
     plot([xmin,xmax],[0,0],'--','Color',[0.5,0.5,0.5],'Linewidth',0.5)
     plot([0,0],[ymin,ymax],'--','Color',[0.5,0.5,0.5],'Linewidth',0.5)
     if volgsigma
      plot(real(mu00),imag(mu00),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',1)
     end
     plot(real(mu(1,:)),imag(mu(1,:)),'o','MarkerEdgeColor','[0.5,0.5,0.5]','MarkerFaceColor','[0.5,0.5,0.5]','MarkerSize',2)
     pause(0.1)
     %% test voor sigma iteratie:
     for j = 1:Nmu
       switch (j<=Nmu5)*(mod(j-1,8)+1)+(j>Nmu5)*(mod(Nmu-j,8)+1)
         case 1, plot(real(mus(:,j)),imag(mus(:,j)),'-k','LineWidth',0.5)
         case 2, plot(real(mus(:,j)),imag(mus(:,j)),'-b','LineWidth',0.5)
         case 3, plot(real(mus(:,j)),imag(mus(:,j)),'-r','LineWidth',0.5)
         case 4, plot(real(mus(:,j)),imag(mus(:,j)),'-c','LineWidth',0.5)
         case 5, plot(real(mus(:,j)),imag(mus(:,j)),'-m','LineWidth',0.5)
         case 6, plot(real(mus(:,j)),imag(mus(:,j)),'Color',[200/255,173/255,0],'LineWidth',0.5)
         case 7, plot(real(mus(:,j)),imag(mus(:,j)),'Color',[0.2314,0.4549,0.2941],'LineWidth',0.5)
         case 8, plot(real(mus(:,j)),imag(mus(:,j)),'Color',[0.61,0.51,0.74],'LineWidth',0.5)
       end
     end
     %% test voor sigma iteratie:
     if onlysigma
        fid = fopen(sgmfile,'w');
        cl  = onCleanup(@()fclose(fid));
        for ns=1:Nsg1
          for j=1:Nmu
                switch true
                  case abs(real(mus(ns,j)))<9.9999,                                 fprintf(fid,[blanks(3),'%+4.12f'], real(mus(ns,j)));
                  case abs(real(mus(ns,j)))>=9.9999 && abs(real(mus(ns,j)))<99.999, fprintf(fid,[blanks(2),'%+4.12f'], real(mus(ns,j)));
                  otherwise,                                                        fprintf(fid,[blanks(1),'%+4.12f'], real(mus(ns,j)));
                end
                switch true
                  case abs(imag(mus(ns,j)))<9.9999,                                 fprintf(fid,[blanks(3),'%+4.12fi'], imag(mus(ns,j)));
                  case abs(imag(mus(ns,j)))>=9.9999 && abs(imag(mus(ns,j)))<99.999, fprintf(fid,[blanks(2),'%+4.12fi'], imag(mus(ns,j)));
                  otherwise,                                                        fprintf(fid,[blanks(1),'%+4.12fi'], imag(mus(ns,j)));
                end
          end;
          fprintf(fid,['\n']);
        end;
        quit
     end
  end
  figure(100)
  clf(100,'reset')
  hold on;
  grid on;
  axis([xmin xmax ymin ymax]);
  plot([xmin,xmax],[0,0],'--','Color',[0.5,0.5,0.5],'Linewidth',0.6)
  plot([0,0],[ymin,ymax],'--','Color',[0.5,0.5,0.5],'Linewidth',0.6)
  plot(real(mu00),imag(mu00),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',1)
  plot(real(mu(1,:)),imag(mu(1,:)),'o','MarkerEdgeColor','[0,0,0]','MarkerFaceColor','[0,0,0]','MarkerSize',1.6)
  if profiel(1:3) == 'bou'
     title(['\omega = ',sprintf('%0.1f',omega),  ...
                                                 ...
            ', \tau = ',sprintf('%0.2f',tau),    ...
            ', \sigma = ',sprintf('%0.2f',sigma),...
            ', \delta = ',sprintf('%0.3G',dispth),...
            ', Re(Z) = ',sprintf('%0.2f',RZ)     ...
            ],'fontweight','normal','fontsize',10)
  else
     title(['\omega = ',sprintf('%0.1f',omega),  ...
                                                 ...
            ', \tau = ',sprintf('%0.2f',tau),    ...
            ', \sigma = ',sprintf('%0.2f',sigma),...
            ', Re(Z) = ',sprintf('%0.2f',RZ)     ...
            ],'fontweight','normal','fontsize',10)
  end
  pause(0.1)
  tic
  for j=1:Nmu
    fprintf('j=%2i, ',j)
    if mod(j,10)==0, fprintf('\n'), end
    mu0  = mu(1,j);
    sigma = sigma1;
    if (j==1 & smallfrqapp & profiel=='lin')
            kappr(1,:) = SmallFreqApprox(RZ+Zsign*1E10*i,omega,tau,sigma);
    end
    for nz = 2:NZa-1   % NZa om Z te varieren
        Z  = Zarray(nz);
        rest = 1;
        cnt  = 0;
        [M0,DM0]   = MakeM(mu0);
        while (rest>1e-10 & cnt<50)
          cnt      = cnt+1;
          [V,E]    = eig(M0,DM0);
          E        = diag(E);
          [q,p]    = min(abs(E));
          mu0      = mu0-E(p);
          [M0,DM0] = MakeM(mu0);
          rest     = norm(M0*V(:,p));
        end
        mu(nz,j) = mu0;
        if (j==1 & smallfrqapp)
            kappr(nz,:) = SmallFreqApprox(Z,omega,tau,sigma);
        end

        if volgsigma
        if abs(imag(mu0))<1e-10
           %title(['\mu = ',sprintf('%0.10f',real(mu0))],'fontweight','normal')
           switch true
               case abs(imag(Z))<10,                           display([sprintf('s=%0.2f, Z = %0.2f%+0.2fi, r=%1i, mu0 = %4.12f',sigma,real(Z), imag(Z), cnt,real(mu0))])
               case abs(imag(Z))>=10 && abs(imag(Z))<100,      display([sprintf('s=%0.2f, Z = %0.2f%+0.1fi, r=%1i, mu0 = %4.12f',sigma,real(Z), imag(Z), cnt,real(mu0))])
               case abs(imag(Z))>=100 && abs(imag(Z))<1000,    display([sprintf('s=%0.2f, Z = %0.2f%+4.f.i, r=%1i, mu0 = %4.12f',sigma,real(Z), imag(Z), cnt,real(mu0))])
               case abs(imag(Z))>=1000 && abs(imag(Z))<10000,  display([sprintf('s=%0.2f, Z = %0.2f%+4.fi,  r=%1i, mu0 = %4.12f',sigma,real(Z), imag(Z), cnt,real(mu0))])
               case abs(imag(Z))>=10000 && abs(imag(Z))<1E10,  display([sprintf('s=%0.2f, Z = %0.2f%+2.0fe%+di, r=%1i, mu0 = %4.12f',sigma,real(Z), imag(Z)/10^(floor(log10(abs(imag(Z))))),floor(log10(abs(imag(Z)))), cnt,real(mu0))])
               otherwise,                                      display([sprintf('s=%0.2f, Z = %0.2f+inf.i,  r=%1i, mu0 = %4.12f',sigma, real(Z),cnt,real(mu0))])
           end
        else
           %title(['\mu = ',sprintf('%0.10f%+0.10fi',real(mu0),imag(mu0))],'fontweight','normal')
           switch true
               case abs(imag(Z))<10,                           display([sprintf('s=%0.2f, Z = %0.2f%+0.2fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
               case abs(imag(Z))>=10 && abs(imag(Z))<100,      display([sprintf('s=%0.2f, Z = %0.2f%+0.1fi, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
               case abs(imag(Z))>=100 && abs(imag(Z))<1000,    display([sprintf('s=%0.2f, Z = %0.2f%+4.f.i, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
               case abs(imag(Z))>=1000 && abs(imag(Z))<10000,  display([sprintf('s=%0.2f, Z = %0.2f%+4.fi,  r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z), cnt,real(mu0),imag(mu0))])
               case abs(imag(Z))>=10000 && abs(imag(Z))<1E10,  display([sprintf('s=%0.2f, Z = %0.2f%+2.0fe%+di, r=%1i, mu0 = %4.12f %+4.12fi',sigma,real(Z), imag(Z)/10^(floor(log10(abs(imag(Z))))),floor(log10(abs(imag(Z)))), cnt,real(mu0),imag(mu0))])
               otherwise,                                      display([sprintf('s=%0.2f, Z = %0.2f+inf.i,  r=%1i, mu0 = %4.12f %+4.12fi',sigma, real(Z),cnt,real(mu0),imag(mu0))])
           end
        end
        end
        mu0 = 2*mu(nz,j)-mu(nz-1,j);  % we hoeven niet nz=1 apart te houden omdat we starten met nz=2
    end % over nz
    [U,I] = min(abs(mu(NZa-1,j)-mu(1,:)));                            % we sluiten aan op de dichtstbijzijnde HW mode
    if U<1, mu(NZa,j) = mu(1,I); else, mu(NZa,j) = mu(NZa-1,j); end
    if (j==1 & smallfrqapp), kappr(NZa,:) = kappr(1,:); end
    switch (j<=Nmu5)*(mod(j-1,8)+1)+(j>Nmu5)*(mod(Nmu-j,8)+1)
      case 1, plot(real(mu(:,j)),imag(mu(:,j)),'-k','LineWidth',                              0.5)
      case 2, plot(real(mu(:,j)),imag(mu(:,j)),'-b','LineWidth',                              0.5)
      case 3, plot(real(mu(:,j)),imag(mu(:,j)),'-r','LineWidth',                              0.5)
      case 4, plot(real(mu(:,j)),imag(mu(:,j)),'-','Color',[0,0.7,0.7],'LineWidth',           0.5)
      case 5, plot(real(mu(:,j)),imag(mu(:,j)),'-m','LineWidth',                              0.5)
      case 6, plot(real(mu(:,j)),imag(mu(:,j)),'-','Color',[200/255,173/255,0],'LineWidth',   0.5)
      case 7, plot(real(mu(:,j)),imag(mu(:,j)),'-','Color',[0.2314,0.4549,0.2941],'LineWidth',0.5)
      case 8, plot(real(mu(:,j)),imag(mu(:,j)),'-','Color',[0.61,0.51,0.74],'LineWidth',      0.5)
    end
    if smallfrqapp
     plot(real(kappr(:,1)),imag(kappr(:,1)),'--r','LineWidth',0.5)
     plot(real(kappr(:,2)),imag(kappr(:,2)),'--r','LineWidth',0.5)
    end
    pause(0.01)
  end;  % over j  (Nmu)
  telapsedz = toc;
  fid = fopen(parmfile,'a');
  cl  = onCleanup(@()fclose(fid));
  if Nsg == 0
     fprintf(fid,'time       = %1.12G seconds\n', telapsedz);
     display([sprintf('\ntime = %1.12G seconds', telapsedz)]);;
  else
     fprintf(fid,'time       = %1.12G + %1.12G seconds\n', telapseds, telapsedz);
     display([sprintf('\ntime = %1.12G + %1.12G seconds', telapseds, telapsedz)]);
  end;
  %%%%%%%%%
  box on
  set(gcf,'color',[1 1 1])
  set(gcf,'inverthardcopy','off')
  set(gcf,'papertype','a4letter')
  set(gcf,'paperunits','centimeters')
  set(gcf,'paperposition',[0  0  10  10 ])
  set(gcf,'clipping','off')
  %eval(['print -depsc -cmyk ' epsfile]);
  print(epsfile,'-depsc','-cmyk','-r600')
  fixepsbbox(epsfile,dx1,dy1,dx2,dy2)
  print(jpgfile,'-djpeg','-r600')
  print(pngfile,'-dpng' ,'-r600')
  %%%%%%%%%%
  fid = fopen(savefile,'w');
  cl  = onCleanup(@()fclose(fid));
  for nz = 1:NZa
    Z = Zarray(nz);
    switch true
      case abs(real(Z))<9.9999,                         fprintf(fid,[blanks(3),'%+4.4f'], real(Z));
      case abs(real(Z))>=9.9999 && abs(real(Z))<99.999, fprintf(fid,[blanks(2),'%+4.4f'], real(Z));
      case abs(real(Z))>=99.999 && abs(real(Z))<999.99, fprintf(fid,[blanks(1),'%+4.4f'], real(Z));
      otherwise,                                        fprintf(fid,[blanks(1),'%+9.1e'], real(Z));
    end
    switch true
      case abs(imag(Z))<9.9999,                         fprintf(fid,[blanks(3),'%+4.4fi'], imag(Z));
      case abs(imag(Z))>=9.9999 && abs(imag(Z))<99.999, fprintf(fid,[blanks(2),'%+4.4fi'], imag(Z));
      case abs(imag(Z))>=99.999 && abs(imag(Z))<999.99, fprintf(fid,[blanks(1),'%+4.4fi'], imag(Z));
      otherwise,                                        fprintf(fid,[blanks(1),'%+9.1ei'], imag(Z));
    end
    for j=1:Nmu
          switch true
            case abs(real(mu(nz,j)))<9.9999,                                fprintf(fid,[blanks(3),'%+4.12f'], real(mu(nz,j)));
            case abs(real(mu(nz,j)))>=9.9999 && abs(real(mu(nz,j)))<99.999, fprintf(fid,[blanks(2),'%+4.12f'], real(mu(nz,j)));
            otherwise,                                                      fprintf(fid,[blanks(1),'%+4.12f'], real(mu(nz,j)));
          end
          switch true
            case abs(imag(mu(nz,j)))<9.9999,                                fprintf(fid,[blanks(3),'%+4.12fi'], imag(mu(nz,j)));
            case abs(imag(mu(nz,j)))>=9.9999 && abs(imag(mu(nz,j)))<99.999, fprintf(fid,[blanks(2),'%+4.12fi'], imag(mu(nz,j)));
            otherwise,                                                      fprintf(fid,[blanks(1),'%+4.12fi'], imag(mu(nz,j)));
          end
    end;
    fprintf(fid,['\n']);
  end;
  %fclose(fid);

  function u = U0(yy);
           switch profiel(1:3)
             case 'fla',  u = tau;
             case 'lin',  u = tau + sigma*yy;
             case 'par',  u = tau - sigma*(2*yy-1).^2;
             case 'bou'
                 a = dispth;
                 u = tau + sigma*( -1 + tanh((1-yy)/a) + (1-tanh(1/a))*( 1 + (1 + (1+tanh(1/a))/a)*yy).*(1-yy) );
             otherwise,   error 'Invalid option for profiel'; return
           end
  end
  function c = C0(yy);
           c = ones(size(yy));
  end

  function kp = SmallFreqApprox(Z,omg,ta,sg)
           if (profiel == 'lin')
            switch softside
                case 0,
                    zeta0 = Z*omg;
                    kp(1) =   omg/(sqrt(1/(1-i/zeta0)+sg^2/4) + (ta+sg/2));
                    kp(2) = - omg/(sqrt(1/(1-i/zeta0)+sg^2/4) - (ta+sg/2));
                case 1,
                    zeta1 = Z*omg;
                    kp(1) =   omg/(sqrt(1/(1-i/zeta1)+sg^2/4) + (ta+sg/2));
                    kp(2) = - omg/(sqrt(1/(1-i/zeta1)+sg^2/4) - (ta+sg/2));
                case 2,
                   zeta0 = Z*omg;
                   zeta1 = Z*omg;
                   kp(1) =   omg/(sqrt(1/(1-i*(1/zeta0+1/zeta1))+sg^2/4) + (ta+sg/2));
                   kp(2) = - omg/(sqrt(1/(1-i*(1/zeta0+1/zeta1))+sg^2/4) - (ta+sg/2));
%                     zeta0 = Z/omg;
%                     zeta1 = Z/omg;
%                     kp(1) = omg*(1-i*(zeta0+zeta1))/(ta+sg/2+sqrt(i*(zeta0+zeta1)*((ta+sg/2)^2+sg^2/12)-sg^2/12));
%                     kp(2) = omg*(1-i*(zeta0+zeta1))/(ta+sg/2-sqrt(i*(zeta0+zeta1)*((ta+sg/2)^2+sg^2/12)-sg^2/12));
                otherwise
                    error 'Invalid option for SmallFreqApprox', return
            end
           else
             kp=[0,0];
           end
  end

  function [M,DM] = MakeM(mu)
     ugl = U0(ygl)             ;
     cgl = C0(ygl)             ;
     Om1 = (omega - mu*ugl)./cgl;
     Om2 = Om1.*Om1;
     Om3 = Om2.*Om1;
     M   = zeros(n+1);
     DM  = zeros(n+1);
     for nn = 1:n+1
       for mm = nn:n+1
           int1gl     = -reshape(DTTgl(nn,mm,:),size(ygl))./Om2 + (1-mu^2./Om2).*reshape(TTgl(nn,mm,:),size(ygl));
           int2gl     = -2*(ugl.*reshape(DTTgl(nn,mm,:),size(ygl)) + omega*mu.*reshape(TTgl(nn,mm,:),size(ygl)))./cgl./Om3;
           M(nn,mm) = sum(Wgl.*int1gl);
          DM(nn,mm) = sum(Wgl.*int2gl);
           if abs(Z)<1E9
              switch softside
                  case  2,   M(nn,mm) = M(nn,mm) + (T0(nn)*T0(mm)+T1(nn)*T1(mm))/(i*omega*Z);
                  case  1,   M(nn,mm) = M(nn,mm) +  T1(nn)*T1(mm)/(i*omega*Z);
                  case  0,   M(nn,mm) = M(nn,mm) +  T0(nn)*T0(mm)/(i*omega*Z);
                  otherwise, display('verkeerde softside'); return
              end
           end
           if nn ~= mm, M(mm,nn) = M(nn,mm); DM(mm,nn) = DM(nn,mm); end
       end
     end
  end

  function MakeT2D(a,b)
   [tgl,Wgl] = legzo(N)           ;  % N punten
     t       = [-1, tgl, 1]       ;  % 1 + N + 1 = N+2=N1 punten; eindpunten vanwege rvw
     y       = (b-a)/2*t  +(b+a)/2;  % vertaal interval [-1,1] naar [a,b], y wordt nog niet gebruikt
     ygl     = y(2:N+1)           ;
     Wgl     = Wgl*(b-a)/2        ;

     T       = zeros(n+1,N1);
     U       = zeros(n+1,N1);
     DT      = zeros(n+1,N1);
     T(1,:)  = ones(1,N1);
     U(1,:)  = ones(1,N1);
     DT(1,:) = zeros(1,N1);
     T(2,:)  = t;
     U(2,:)  = 2*t;
     DT(2,:) = ones(1,N1);
     for nn = 3:n+1
        T(nn,:)  = 2*t.*T(nn-1,:) - T(nn-2,:);  % Chebyshev T  van orde nn-1 ( dus T(1,:) is T_0(x) )
        U(nn,:)  = 2*t.*U(nn-1,:) - U(nn-2,:);  % Chebyshev U
        DT(nn,:) = (nn-1)*U(nn-1,:)          ;  % Chebyshev T'
     end
     DT     = DT*2/(b-a);                       % corrigeer
      Tgl   =  T(:,2:N+1);
     DTgl   = DT(:,2:N+1);
      T0    =  T(:,N0);
      T1    =  T(:,N1);
    %DT0    = DT(n+1,1);                       % misschien nodig bij andere randvoorwaarden
    %DT1    = DT(n+1,N1);
    TTgl      = zeros(n+1,n+1,N);
    DTTgl     = zeros(n+1,n+1,N);
    for nn = 1:n+1
       for mm = nn:n+1
           TTgl(nn,mm,:)   =   Tgl(nn,:).*Tgl(mm,:);   % vanwege intensief hergebruik bij maken van M
           DTTgl(nn,mm,:)  =  DTgl(nn,:).*DTgl(mm,:);   % voorzichtig: misschien reshape toepassen
       end
    end
    %if nn ~= mm, TTgl(mm,nn,:) = TTgl(nn,mm,:); DTTgl(mm,nn,:) = DTTgl(nn,mm,:); end % niet echt nodig
  end
end


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
end

function [ps,ix] = dpsimplify(p,tol)

% Recursive Douglas-Peucker Polyline Simplification, Simplify
%
% [ps,ix] = dpsimplify(p,tol)
%
% dpsimplify uses the recursive Douglas-Peucker line simplification
% algorithm to reduce the number of vertices in a piecewise linear curve
% according to a specified tolerance. The algorithm is also know as
% Iterative Endpoint Fit. It works also for polylines and polygons
% in higher dimensions.
%
% In case of nans (missing vertex coordinates) dpsimplify assumes that
% nans separate polylines. As such, dpsimplify treats each line
% separately.
%
% For additional information on the algorithm follow this link
% http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
%
% Input arguments
%
%     p     polyline n*d matrix with n vertices in d
%           dimensions.
%     tol   tolerance (maximal euclidean distance allowed
%           between the new line and a vertex)
%
% Output arguments
%
%     ps    simplified line
%     ix    linear index of the vertices retained in p (ps = p(ix))
%
% Examples
%
% 1. Simplify line
%
%     tol    = 1;
%     x      = 1:0.1:8*pi;
%     y      = sin(x) + randn(size(x))*0.1;
%     p      = [x' y'];
%     ps     = dpsimplify(p,tol);
%
%     plot(p(:,1),p(:,2),'k')
%     hold on
%     plot(ps(:,1),ps(:,2),'r','LineWidth',2);
%     legend('original polyline','simplified')
%
% 2. Reduce polyline so that only knickpoints remain by
%    choosing a very low tolerance
%
%     p = [(1:10)' [1 2 3 2 4 6 7 8 5 2]'];
%     p2 = dpsimplify(p,eps);
%     plot(p(:,1),p(:,2),'k+--')
%     hold on
%     plot(p2(:,1),p2(:,2),'ro','MarkerSize',10);
%     legend('original line','knickpoints')
%
% 3. Simplify a 3d-curve
%
%     x = sin(1:0.01:20)';
%     y = cos(1:0.01:20)';
%     z = x.*y.*(1:0.01:20)';
%     ps = dpsimplify([x y z],0.1);
%     plot3(x,y,z);
%     hold on
%     plot3(ps(:,1),ps(:,2),ps(:,3),'k*-');
%
%
%
% Author: Wolfgang Schwanghart, 13. July, 2010.
% w.schwanghart[at]unibas.ch


if nargin == 0
    help dpsimplify
    return
end

error(nargchk(2, 2, nargin))

% error checking
if ~isscalar(tol) || tol<0;
    error('tol must be a positive scalar')
end


% nr of dimensions
nrvertices    = size(p,1);
dims    = size(p,2);

% anonymous function for starting point and end point comparision
% using a relative tolerance test
compare = @(a,b) abs(a-b)/max(abs(a),abs(b)) <= eps;

% what happens, when there are NaNs?
% NaNs divide polylines.
Inan      = any(isnan(p),2);
% any NaN at all?
Inanp     = any(Inan);

% if there is only one vertex
if nrvertices == 1 || isempty(p);
    ps = p;
    ix = 1;

% if there are two
elseif nrvertices == 2 && ~Inanp;
    % when the line has no vertices (except end and start point of the
    % line) check if the distance between both is less than the tolerance.
    % If so, return the center.
    if dims == 2;
        d    = hypot(p(1,1)-p(2,1),p(1,2)-p(2,2));
    else
        d    = sqrt(sum((p(1,:)-p(2,:)).^2));
    end

    if d <= tol;
        ps = sum(p,1)/2;
        ix = 1;
    else
        ps = p;
        ix = [1;2];
    end

elseif Inanp;

    % case: there are nans in the p array
    % --> find start and end indices of contiguous non-nan data
    Inan = ~Inan;
    sIX = strfind(Inan',[0 1])' + 1;
    eIX = strfind(Inan',[1 0])';

    if Inan(end)==true;
        eIX = [eIX;nrvertices];
    end

    if Inan(1);
        sIX = [1;sIX];
    end

    % calculate length of non-nan components
    lIX = eIX-sIX+1;
    % put each component into a single cell
    c   = mat2cell(p(Inan,:),lIX,dims);

    % now call dpsimplify again inside cellfun.
    if nargout == 2;
        [ps,ix]   = cellfun(@(x) dpsimplify(x,tol),c,'uniformoutput',false);
        ix        = cellfun(@(x,six) x+six-1,ix,num2cell(sIX),'uniformoutput',false);
    else
        ps   = cellfun(@(x) dpsimplify(x,tol),c,'uniformoutput',false);
    end

    % write the data from a cell array back to a matrix
    ps = cellfun(@(x) [x;nan(1,dims)],ps,'uniformoutput',false);
    ps = cell2mat(ps);
    ps(end,:) = [];

    % ix wanted? write ix to a matrix, too.
    if nargout == 2;
        ix = cell2mat(ix);
    end


else


% if there are no nans than start the recursive algorithm
ixe     = size(p,1);
ixs     = 1;

% logical vector for the vertices to be retained
I   = true(ixe,1);

% call recursive function
p   = simplifyrec(p,tol,ixs,ixe);
ps  = p(I,:);

% if desired return the index of retained vertices
if nargout == 2;
    ix  = find(I);
end

end

% _________________________________________________________
function p  = simplifyrec(p,tol,ixs,ixe)

    % check if startpoint and endpoint are the same
    % better comparison needed which included a tolerance eps

    c1 = num2cell(p(ixs,:));
    c2 = num2cell(p(ixe,:));

    % same start and endpoint with tolerance
    sameSE = all(cell2mat(cellfun(compare,c1(:),c2(:),'UniformOutput',false)));


    if sameSE;
        % calculate the shortest distance of all vertices between ixs and
        % ixe to ixs only
        if dims == 2;
            d    = hypot(p(ixs,1)-p(ixs+1:ixe-1,1),p(ixs,2)-p(ixs+1:ixe-1,2));
        else
            d    = sqrt(sum(bsxfun(@minus,p(ixs,:),p(ixs+1:ixe-1,:)).^2,2));
        end
    else
        % calculate shortest distance of all points to the line from ixs to ixe
        % subtract starting point from other locations
        pt = bsxfun(@minus,p(ixs+1:ixe,:),p(ixs,:));

        % end point
        a = pt(end,:)';

        beta = (a' * pt')./(a'*a);
        b    = pt-bsxfun(@times,beta,a)';
        if dims == 2;
            % if line in 2D use the numerical more robust hypot function
            d    = hypot(b(:,1),b(:,2));
        else
            d    = sqrt(sum(b.^2,2));
        end
    end

    % identify maximum distance and get the linear index of its location
    [dmax,ixc] = max(d);
    ixc  = ixs + ixc;

    % if the maximum distance is smaller than the tolerance remove vertices
    % between ixs and ixe
    if dmax <= tol;
        if ixs ~= ixe-1;
            I(ixs+1:ixe-1) = false;
        end
    % if not, call simplifyrec for the segments between ixs and ixc (ixc
    % and ixe)
    else
        p   = simplifyrec(p,tol,ixs,ixc);
        p   = simplifyrec(p,tol,ixc,ixe);

    end

end
end

function fixepsbbox(filename,dx1,dy1,dx2,dy2)
% function fixepsbbox(filename)
%
% matlab seems to compute a bounding box on eps files which is too
% large in the x-direction
%
% this script fixes the bounding box
% it is 99% stolen from fixeps.m, located here:
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4818&objectType=file
% the only change is that this changes the bbox numbers
% seriously, the only change is the addition of lines 22,23,33,34
%
% boundingbox has form of:
% %%BoundingBox:    x1   y1   x2   y2
% where (x1,y1) is lower-left and (x2,y2) is upper-right
%
% matlab computes x1 too small and x2 too large
% changes lines to:
% %%BoundingBox:    x1+dx1 y1 x2+dx2 y2
% default amount to change bbox - found this fixed my plots just fine
%dx1 = -10;	% amount to move x1
%dx2 =  25;	% amount to move x2
fid = fopen(filename,'r+');
k=0;
while k <2                                                  % 2 locations to replace.
    tline = fgetl(fid);                                     % get one line text
    stridx=strfind(tline,'Box:');
	if isempty(stridx)==0
        len=length(tline);                                  % the original line length
		bb=sscanf(tline(stridx+4:end),'%i');                % read the numbers
		bb(1) = bb(1) + dx1;								% change x1
		bb(2) = bb(2) + dy1;								% change y1
		bb(3) = bb(3) + dx2;								% change x2
		bb(4) = bb(4) + dy2;								% change y2
		bbstr=sprintf('%g %g %g %g',bb);                    % write bb numbers to string
        tline=tline(1:stridx+3);                             % keep the "%%(page)boundingbox" string (with starting '%%')
		spaces(1:len-length(tline)-length(bbstr)-1)=' ';    % add trailing spaces as to overwrite old line completely
		tline=[tline ' ' bbstr spaces];                     % concate numbers and blank spaces to "%%(page)boundingbox"
        fseek(fid,-len-2,'cof');                            % before using fprintf search to correct position
		count = fprintf(fid,'%s',tline);
        fseek(fid,2,'cof');                                 % seek to beginning of line (for windows text file) on
                                                            % for linux: change '2' to '1' I think
        k=k+1;
	end
end
fclose(fid);
end
