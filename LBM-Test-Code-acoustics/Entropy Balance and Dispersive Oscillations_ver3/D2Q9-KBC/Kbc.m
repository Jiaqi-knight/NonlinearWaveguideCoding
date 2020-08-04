%%    populations_mirr=populations+repmat(alpha,1,1,Q).*(popequilibriums-populations);
%    populations = LBMPropagate(populations,nx,ny,scheme); %Propagate the populations
function [populations]=Kbc(beta,nx,ny,Q,rho,scheme,populations,popequilibriums,kbc_coluring_scheme)
%Finding the moments KBC
dM_xy=zeros(nx,ny);  % the difference in moments M-M_eq
dM_yy=zeros(nx,ny);
dM_xx=zeros(nx,ny);
dM_xyy=zeros(nx,ny);
dM_xxy=zeros(nx,ny);
dM_xxyy=zeros(nx,ny);
for i=1:9
    dM_xy =dM_xy+scheme(i,1)*scheme(i,2)*(populations(:,:,i)-popequilibriums(:,:,i));
    dM_xx =dM_xx+scheme(i,1)*scheme(i,1)*(populations(:,:,i)-popequilibriums(:,:,i));
    dM_yy =dM_yy+scheme(i,2)*scheme(i,2)*(populations(:,:,i)-popequilibriums(:,:,i));
    dM_xyy =dM_xyy+scheme(i,1)*scheme(i,2)*scheme(i,2)*(populations(:,:,i)-popequilibriums(:,:,i));
    dM_xxy =dM_xxy+scheme(i,1)*scheme(i,1)*scheme(i,2)*(populations(:,:,i)-popequilibriums(:,:,i));
    dM_xxyy=dM_xxyy+scheme(i,1)*scheme(i,1)*scheme(i,2)*scheme(i,2)*(populations(:,:,i)-popequilibriums(:,:,i));
end

dM_xy=dM_xy./rho;
dM_yy=dM_yy./rho;
dM_xx=dM_xx./rho;
dM_xxy=dM_xxy./rho;
dM_xyy=dM_xyy./rho;
dM_xxyy=dM_xxyy./rho;


delS=zeros(nx,ny,Q);
%unsigned int kbc_coluring_scheme=0;
%kbc_coluring_scheme=1;
switch (kbc_coluring_scheme)
    case 0
        %KBC minimalistic grouping % lecture 7 p.8)
        delS(:,:,1) = 0.0;
        delS(:,:,2) = 0.5*rho.*0.5.*(dM_xx-dM_yy);
        delS(:,:,3) = 0.5*rho.*0.5.*-1.0 .*(dM_xx-dM_yy);
        delS(:,:,4) = 0.5*rho.*0.5.*(dM_xx-dM_yy);
        delS(:,:,5) = 0.5*rho.*0.5*-1.0 .*(dM_xx-dM_yy);
        delS(:,:,6) = 0.25*rho.*scheme(5,1)*scheme(5,2).*dM_xy;
        delS(:,:,7) = 0.25*rho.*scheme(6,1)*scheme(6,2).*dM_xy;
        delS(:,:,8) = 0.25*rho.*scheme(7,1)*scheme(7,2).*dM_xy;
        delS(:,:,9) = 0.25*rho.*scheme(8,1)*scheme(8,2).*dM_xy;
        
    case 1
        % classical
        delS(:,:,1) = -1.0*rho.*(dM_xx+dM_yy);
        delS(:,:,2) = 0.5*rho.*0.5 .*(2*dM_xx);
        delS(:,:,3) = 0.5*rho.*0.5 .*(2*dM_yy);
        delS(:,:,4) = 0.5*rho*0.5 .*(2*dM_xx);
        delS(:,:,5) = 0.5*rho*0.5 .*(2*dM_yy);
        delS(:,:,6) = 0.25*rho*scheme(5,1)*scheme(5,2).*dM_xy;
        delS(:,:,7) = 0.25*rho*scheme(6,1)*scheme(6,2).*dM_xy;
        delS(:,:,8) = 0.25*rho*scheme(7,1)*scheme(7,2).*dM_xy;
        delS(:,:,9) = 0.25*rho*scheme(8,1)*scheme(8,2).*dM_xy;
end

% fi = ki + si + hi
%calculating delH = fi - fi_eq - delS

delH = populations-popequilibriums-delS;


entScalarProd_dSdH=zeros(nx,ny);
entScalarProd_dHdH=zeros(nx,ny);

% Calculating gamma
for i=1:9
    entScalarProd_dSdH =entScalarProd_dSdH+ delS(:,:,i).*delH(:,:,i)./popequilibriums(:,:,i);
    entScalarProd_dHdH =entScalarProd_dHdH+ delH(:,:,i).*delH(:,:,i)./popequilibriums(:,:,i);
end

alpha = (1.0/beta) - (2.0 - 1.0/beta)*(entScalarProd_dSdH./entScalarProd_dHdH);
alpha(alpha  <= 1 | entScalarProd_dHdH == 0) = 2; % gamma = 1/beta; ->       regularised     % gamma = 2 -> LBGK


populations=popequilibriums+(1-2.0*beta)*delS+repmat(1.0-alpha*beta,1,1,Q).*delH;

end