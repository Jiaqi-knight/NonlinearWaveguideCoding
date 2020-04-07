%EFun Function for intial based  functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.



function [Out]= EFun(Geo_b,Psi,Wave,rk) 

a=Wave.a;
b=Wave.b;
a_b=Wave.a_b;
k=Wave.k; 
gamma=Wave.gamma;
n_matrix=length(Geo_b.m)*Geo_b.n;

tic
disp('EFun..Begin...')
%2D-bsxfun(@times, 3D(\alpha*\beta*s), 1*1*s*a)-->4D[\alpha*\beta*s*a]
Fun2_a.I2=      bsxfun(@times,Psi.ab_r,reshape(ones(size(a)),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Out.I2=Fun2_a.I2;
Fun2_a.V=       bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_b.V=       bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*b),1,1,1,length(b)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a_b.V=     bsxfun(@times,Psi.a_pr_b_r,reshape(1./(sqrt(-1)*k*a_b),1,1,1,length(b),length(a)));    %(James-3.26,Jiaqi-76)  --expand,{alpha*beta*a*s}
Fun2_a.W=       -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*a),1,1,1,length(a))); %(James-3.27,Jiaqi-77)
Fun2_b.W=       -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*b),1,1,1,length(b))); %(James-3.27,Jiaqi-77)
Fun2_a_b.W=     -bsxfun(@times,Psi.pt_ab,reshape(1./(sqrt(-1)*k*a_b),1,1,1,length(b),length(a))); %(James-3.27,Jiaqi-77)
Fun2_a.N=       bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*a),1,1,length(Geo_b.kappa),length(a)));%(James-3.35b)
Fun2_a_b.N=     bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*a_b),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*permute(a_b,[3,1,2])),1,1,length(Geo_b.kappa),length(b),length(a)));%(James-3.35b)
Fun2_b.N=       bsxfun(@times,Psi.ab_r,reshape((sqrt(-1)*k*b),1,1,1,length(b)))...
    -bsxfun(@times,Psi.ab_r2_cos,reshape((sqrt(-1)*k*Geo_b.kappa.'.*b),1,1,length(Geo_b.kappa),length(b )));%(James-3.35b)
for i=1:length(Geo_b.s);for j=1:length(a);Fun2_a.N_inv(:,:,i,j)=Fun2_a.N(:,:,i,j)^(-1);end;end;
for i=1:length(Geo_b.s);for j=1:length(b); for k=1:length(a);if a(k)-b(j)==0;Fun2_a_b.N_inv(:,:,i,j,k)=NaN+NaN*sqrt(-1);else;Fun2_a_b.N_inv(:,:,i,j,k)=Fun2_a_b.N(:,:,i,j,k)^(-1);end;end;end;end;
for i=1:length(Geo_b.s);for j=1:length(b);Fun2_b.N_inv(:,:,i,j)=Fun2_b.N(:,:,i,j)^(-1);end;end;

Fun2_a.G=       -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_b.G=       -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(b)),1,1,1,length(b)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(b)),1,1,1,length(b)));
Fun2_a_b.G=     -bsxfun(@times,Psi.ps_ab_r_1,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.ps_ab_r_2,reshape(ones(size(a_b)),1,1,1,length(b),length(a)));
Fun2_a.H=       -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a)),1,1,1,length(a)));%(James-3.35d)
Fun2_b.H=       -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(b)),1,1,1,length(b)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(b)),1,1,1,length(b)));%(James-3.35d)
Fun2_a_b.H=     -bsxfun(@times,Psi.a_ps_b_r_1,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.a_ps_b_r_2,reshape(ones(size(a_b)),1,1,1,length(b),length(a)));%(James-3.35d)
Fun2_a.M_2_1=   bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));
Fun2_a_b.M_2_1= bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*permute(ones(size(a_b)),[3,1,2]),1,1,length(Geo_b.kappa),length(b),length(a)));
Fun2_b.M_2_1=   bsxfun(@times,Psi.pr_ab_r,reshape(ones(size(b)),1,1,1,length(b)))...
    -bsxfun(@times,Psi.pr_ab_r2_cos,reshape(Geo_b.kappa.'.*ones(size(b)),1,1,length(Geo_b.kappa),length(b)));

Fun2_a.M_3_1=   bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));
Fun2_a_b.M_3_1= bsxfun(@times,Psi.pt_ab,reshape(ones(size(a_b)),1,1,1,length(b),length(a)))...
    -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*permute(ones(size(a_b)),[3,1,2]),1,1,length(Geo_b.kappa),length(b),length(a)));
Fun2_b.M_3_1=   bsxfun(@times,Psi.pt_ab,reshape(ones(size(b)),1,1,1,length(b)))...
    -bsxfun(@times,Psi.pt_ab_r_cos,reshape(Geo_b.kappa.'.*ones(size(b)),1,1,length(Geo_b.kappa),length(b)));

Fun2_a.M  =     -Fun2_a.N-multiprod(Fun2_a.M_2_1,Fun2_a.V,[1,2])-multiprod(Fun2_a.M_3_1,Fun2_a.W,[1,2]);%(James-3.35a)
Fun2_a_b.M=     -Fun2_a_b.N-multiprod(Fun2_a_b.M_2_1,Fun2_a_b.V,[1,2])-multiprod(Fun2_a_b.M_3_1,Fun2_a_b.W,[1,2]);%(James-3.35a)
Fun2_b.M=       -Fun2_b.N-multiprod(Fun2_b.M_2_1,Fun2_b.V,[1,2])-multiprod(Fun2_b.M_3_1,Fun2_b.W,[1,2]);%(James-3.35a)

Fun2_a.A_2_1_1= bsxfun(@times,Psi.ab_s1,reshape(ones(size(a)),1,1,1,length(a)));
Fun2_a.A_2_1_2= bsxfun(@times,Psi.ab_s2,reshape(ones(size(a)),1,1,1,length(a)));
Fun2_a.A_2_1_3= bsxfun(@times,Psi.ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.ab_r_cos,reshape((Geo_b.kappa.'.*ones(size(a))),1,1,length(Geo_b.kappa),length(a)));
Fun2_a.A_2_1_4= Fun2_a.M_2_1;
Fun2_a.A_3_1_2= bsxfun(@times,Psi.pt_ab,reshape(ones(size(a)),1,1,1,length(a)))...
    -bsxfun(@times,Psi.pt_ab_cos,reshape(Geo_b.kappa.'.*ones(size(a)),1,1,length(Geo_b.kappa),length(a)));


%3D-bsxfun(@times, 4D(\alpha*\beta*\gamma*s), 1*1*1*s*b*a)-->6D[\alpha*\beta*gamma*s*b*a]
Fun3_ab.I3_1=   bsxfun(@times,Psi.abc,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.I3_r=   bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));

Fun3_ab.A_1=    bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
    -bsxfun(@times,Psi.abc_r2_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(sqrt(-1)*k*b),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.A_2=    permute(bsxfun(@times,Fun2_a.M_2_1,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.N_inv=  permute(bsxfun(@times,Fun2_a.N_inv,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_2_1=  permute(bsxfun(@times,Fun2_a.A_2_1_1+Fun2_a.A_2_1_2+Fun2_a.A_2_1_3+Fun2_a.A_2_1_4,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_2_2=  bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.A_3=    permute(bsxfun(@times,Fun2_a.M_3_1,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_3_1_2=permute(bsxfun(@times,Fun2_a.A_3_1_2,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);
Fun3_ab.A_3_2=  bsxfun(@times,Psi.abc_r_sin,reshape(bsxfun(@times,Geo_b.kappa.'*ones(size(b)),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.A =     -Fun3_ab.A_1...
    -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-multiprod(Fun3_ab.I3_r,Fun3_ab.A_2_1,[2,3])-Fun3_ab.A_2_2),[1,2])...
    -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(multiprod(Fun3_ab.I3_r,Fun3_ab.A_3_1_2,[2,3])+Fun3_ab.A_3_2),[1,2]);%(James-3.35e)
Fun3_ab.B_1=    bsxfun(@times,Psi.abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*((gamma-1)/2*sqrt(-1)*k*ones(size(b))),reshape(a,1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
    -bsxfun(@times,Psi.abc_r2_cos,reshape(bsxfun(@times,Geo_b.kappa.'*((gamma-1)/2*sqrt(-1)*k*ones(size(b))),reshape(a,1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.B_2=    Fun3_ab.A_1;
Fun3_a_b.V=     permute(Fun2_a_b.V,[1,2,6,3,4,5]);
Fun3_b.V=       permute(bsxfun(@times,Fun2_b.V,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_ab.B_3=    multiprod(multiprod(Fun3_ab.B_2,Fun3_a_b.V,[1,2]),Fun3_b.V,[2,3]);
Fun3_a_b.W=     permute(Fun2_a_b.W,[1,2,6,3,4,5]);
Fun3_b.W=       permute(bsxfun(@times,Fun2_b.W,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.G=     permute(Fun2_a_b.G,[1,2,6,3,4,5]);
Fun3_b.G=       permute(bsxfun(@times,Fun2_b.G,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.H=     permute(Fun2_a_b.H,[1,2,6,3,4,5]);
Fun3_b.H=       permute(bsxfun(@times,Fun2_b.H,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_ab.B_4=    multiprod(multiprod(Fun3_ab.B_2,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
Fun3_a_b.M=     permute(Fun2_a_b.M,[1,2,6,3,4,5]);
Fun3_b.M=       permute(bsxfun(@times,Fun2_b.M,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);

Fun3_a_b.I=     permute(bsxfun(@times,Fun2_a.I2,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_b.I=       permute(bsxfun(@times,Fun2_a.I2,reshape(ones(size(b)),1,1,1,1,length(b))),[6,1,2,3,5,4]);

Fun3_ab.B_5_1=  multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_2=  multiprod(multiprod(Fun3_ab.B_1/((gamma-1)/2),Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_31= bsxfun(@times,Psi.pr_abc_r,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
    -bsxfun(@times,Psi.pr_abc_r2_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.B_5_3=  multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.V,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_41= bsxfun(@times,Psi.pt_abc,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
    -bsxfun(@times,Psi.pt_abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.B_5_4=  multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.W,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_5_51= bsxfun(@times,Psi.abc,reshape(bsxfun(@times,ones(size(Geo_b.s.'))*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)))...
    -bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.B_5_5=  multiprod(multiprod(Fun3_ab.B_5_51,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);

Fun3_ab.B_6_1=  multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_2=  multiprod(multiprod(Fun3_ab.B_1/((gamma-1)/2),Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_3=  multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.V,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_4=  multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.W,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_5=  multiprod(multiprod(Fun3_ab.B_5_51,Fun3_a_b.W,[1,2]),Fun3_b.V,[2,3]);

Fun3_ab.B =     -Fun3_ab.B_1-Fun3_ab.B_2-Fun3_ab.B_3-Fun3_ab.B_4...
    -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(Fun3_ab.B_5_1+Fun3_ab.B_5_2+Fun3_ab.B_5_3+Fun3_ab.B_5_4+Fun3_ab.B_5_5),[1,2])...
    -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(Fun3_ab.B_6_1+Fun3_ab.B_6_2+Fun3_ab.B_6_3+Fun3_ab.B_6_4+Fun3_ab.B_6_5),[1,2]);

Fun3_ab.C_1=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.I,[1,2]),Fun3_b.M,[2,3]);
Fun3_ab.C_2=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.M,[1,2]),Fun3_b.I,[2,3]);
Fun3_ab.C_3=    Fun3_ab.B_1/((gamma-1)/2);
Fun3_ab.C_4=    multiprod(multiprod(Fun3_ab.B_5_31,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.C_5=    multiprod(multiprod(Fun3_ab.B_5_41,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.B_6_1=  bsxfun(@times,Psi.abc_r_cos,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.C_6=    multiprod(multiprod(Fun3_ab.B_6_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.B_7_1=  bsxfun(@times,Psi.abc_r_sin,reshape(bsxfun(@times,Geo_b.kappa.'*(ones(size(b))),reshape(ones(size(a)),1,1,[])),1,1,1,length(Geo_b.s),length(b),length(a)));
Fun3_ab.C_7=    multiprod(multiprod(Fun3_ab.B_7_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.C=      Fun3_ab.C_1+Fun3_ab.C_2+Fun3_ab.C_3+Fun3_ab.C_4+Fun3_ab.C_5+Fun3_ab.C_6-Fun3_ab.C_7; %(James-3.35g)

Fun3_ab.D_1=    bsxfun(@times,Psi.ps_abc_ps_r, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
Fun3_ab.D_2=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.G,[1,2]),Fun3_b.I,[2,3]);
Fun3_ab.D_3=    multiprod(multiprod(Fun3_ab.I3_r,Fun3_a_b.I,[1,2]),Fun3_b.G,[2,3]);
Fun3_ab.D_4=    bsxfun(@times,Psi.ps_abc_r, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
Fun3_ab.D=      -Fun3_ab.D_1+Fun3_ab.D_2+Fun3_ab.D_3+Fun3_ab.D_4;

Fun3_ab.E_1_1=  multiprod(multiprod(Fun3_ab.D_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.E_1_2=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.G,[1,2]),Fun3_b.V,[2,3]);
Fun3_ab.E_1_3=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.I,[1,2]),multiprod(Fun3_b.G,Fun3_b.V,[2,3]),[2,3]);
Fun3_ab.E_1_4_1=bsxfun(@times,Psi.ps_abc, reshape(ones(size(a_b)),1,1,1,1,length(b),length(a)));
Fun3_ab.E_1_4=  multiprod(multiprod(Fun3_ab.E_1_4_1,Fun3_a_b.I,[1,2]),Fun3_b.V,[2,3]);

Fun3_ab.E_2_1=  multiprod(multiprod(Fun3_ab.D_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.E_2_2=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.G,[1,2]),Fun3_b.W,[2,3]);
Fun3_ab.E_2_3=  multiprod(multiprod(Fun3_ab.I3_1,Fun3_a_b.I,[1,2]),multiprod(Fun3_b.H,Fun3_b.W,[2,3]),[2,3]);
Fun3_ab.E_2_4=  multiprod(multiprod(Fun3_ab.E_1_4_1,Fun3_a_b.I,[1,2]),Fun3_b.W,[2,3]);

Fun3_ab.E =     -multiprod(permute(multiprod(Fun3_ab.A_2,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-Fun3_ab.E_1_1+Fun3_ab.E_1_2+Fun3_ab.E_1_3+Fun3_ab.E_1_4),[1,2])...
    -multiprod(permute(multiprod(Fun3_ab.A_3,Fun3_ab.N_inv,[2,3]),[2,3,1,4,5,6]),(-Fun3_ab.E_2_1+Fun3_ab.E_2_2-Fun3_ab.E_2_3+Fun3_ab.E_2_4),[1,2]);
disp('Warning:a-b!=0, be NaN, need to be deleted.')
disp('EFun..Finish!')
toc



%% #######Boudnary########%

if rk==1
tic
disp('BFun..Begin...')
Fun2_a.NM=              multiprod(Fun2_a.N,Fun2_a.M,[1,2]);
Fun2_a_b.NM=            multiprod(Fun2_a_b.N,Fun2_a_b.M,[1,2]);
Fun2_b.NM=              multiprod(Fun2_b.N,Fun2_b.M,[1,2]);
Fun2_a.L=               [-Fun2_a.G  -Fun2_a.M;...
    Fun2_a.N   Fun2_a.H];
Fun2_b.L=               [-Fun2_b.G  -Fun2_b.M;...
    Fun2_b.N   Fun2_b.H];
Fun2_a_b.L=               [-Fun2_a_b.G  -Fun2_a_b.M;...
    Fun2_a_b.N   Fun2_a_b.H];
Fun2_a.L2=              multiprod(Fun2_a.L,Fun2_a.L,[1,2]);

%Case1: Torsion Free Outlet
%a

for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        Fun2_a.Y(:,:,kh,ka) =                                 sqrt(-1)*Fun2_a.N_inv(:,:,kh,ka)*sqrtm(Fun2_a.NM(:,:,kh,ka));
        [Fun2_a.Vb(:,:,kh,ka),Fun2_a.Lambda(:,:,kh,ka)]=      eigs(sqrt(-1)*sqrtm(Fun2_a.NM(:,:,kh,ka)),n_matrix);
        Fun2_a.YN(:,:,kh,ka)=Fun2_a.Y(:,:,kh,ka)*Fun2_a.N(:,:,kh,ka);
        [Fun2_a.Wb(:,:,kh,ka),Fun2_a.Lambda1(:,:,kh,ka)]=     eigs(Fun2_a.YN(:,:,kh,ka),n_matrix);        
        Fun2_a.Vb_inv(:,:,kh,ka)=                             inv(Fun2_a.Vb(:,:,kh,ka));
        Fun2_a.Wb_inv(:,:,kh,ka)=                             inv(Fun2_a.Wb(:,:,kh,ka));
    end
end
%b
for kh=1:length(Geo_b.h)
    for kb=1:length(b)
        Fun2_b.Y(:,:,kh,kb) = sqrt(-1)*Fun2_b.N_inv(:,:,kh,kb)*sqrtm(Fun2_b.NM(:,:,kh,kb));
        [Fun2_b.Vb(:,:,kh,kb),Fun2_b.Lambda(:,:,kh,kb)]=eigs(sqrt(-1)*sqrtm(Fun2_b.NM(:,:,kh,kb)),n_matrix);
        [Fun2_b.Wb(:,:,kh,kb),temp]=eigs(Fun2_b.Y(:,:,kh,kb)*Fun2_b.N(:,:,kh,kb),n_matrix);
        Fun2_b.Vb_inv(:,:,kh,kb)=inv(Fun2_b.Vb(:,:,kh,kb));
        Fun2_b.Wb_inv(:,:,kh,kb)=inv(Fun2_b.Wb(:,:,kh,kb));
    end
end
%a_b
temp1=(NaN+NaN*sqrt(-1))*ones(size(Fun2_a_b.N(:,:,1,1,1)));
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        for kb=1:length(b)
            if a(ka)-b(kb)==0
                Fun2_a_b.Y(:,:,kh,kb,ka)=         temp1;
                Fun2_a_b.Vb(:,:,kh,kb,ka)=        temp1;
                Fun2_a_b.Lambda(:,:,kh,kb,ka)=    temp1;
                Fun2_a_b.Wb(:,:,kh,kb,ka)=        temp1;
                Fun2_a_b.Vb_inv(:,:,kh,kb,ka)=    temp1;
                Fun2_a_b.Wb_inv(:,:,kh,kb,ka)=    temp1;
            else
                Fun2_a_b.Y(:,:,kh,kb,ka) = sqrt(-1)*Fun2_a_b.N_inv(:,:,kh,kb,ka)*sqrtm(Fun2_a_b.NM(:,:,kh,kb,ka));
                [Fun2_a_b.Vb(:,:,kh,kb,ka),Fun2_a_b.Lambda(:,:,kh,kb,ka)]=eigs(sqrt(-1)*sqrtm(Fun2_a_b.NM(:,:,kh,kb,ka)),n_matrix);
                [Fun2_a_b.Wb(:,:,kh,kb,ka),temp]=eigs(Fun2_a_b.Y(:,:,kh,kb,ka)*Fun2_a_b.N(:,:,kh,kb,ka),n_matrix);
                Fun2_a_b.Vb_inv(:,:,kh,kb,ka)=inv(Fun2_a_b.Vb(:,:,kh,kb,ka));
                Fun2_a_b.Wb_inv(:,:,kh,kb,ka)=inv(Fun2_a_b.Wb(:,:,kh,kb,ka));
            end
        end
    end
end
%2D->3D
Fun2_a.Y_minus=        -Fun2_a.Y;
Fun3_a_b.Vb=           permute(Fun2_a_b.Vb,[1,2,6,3,4,5]);
Fun3_b.Vb=             permute(bsxfun(@times,Fun2_b.Vb,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.Vb_inv=       permute(Fun2_a_b.Vb_inv,[1,2,6,3,4,5]);
Fun3_b.Vb_inv=         permute(bsxfun(@times,Fun2_b.Vb_inv,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.YV=           permute(multiprod(Fun2_a_b.Y,Fun2_a_b.Vb,[1,2]),[1,2,6,3,4,5]);
Fun3_b.YV=             permute(bsxfun(@times,multiprod(Fun2_b.Y,Fun2_b.Vb,[1,2]),reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a.Wb=             permute(bsxfun(@times,Fun2_a.Wb,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Wb_inv=         permute(bsxfun(@times,Fun2_a.Wb_inv,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Y=              permute(bsxfun(@times,Fun2_a.Y,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_ab.II=            ones(size(Fun3_ab.A));
Fun3_a.Lambda=         permute(bsxfun(@times,Fun2_a.Lambda,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a_b.Lambda=       permute(Fun2_a_b.Lambda,[1,2,6,3,4,5]);
Fun3_b.Lambda=         permute(bsxfun(@times,Fun2_b.Lambda,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3.Lambda_alpha_a=   multiprod(Fun3_a.Lambda,multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.I,[2,3]),[1,2]);
Fun3.Lambda_beta_a_b=  multiprod(multiprod(Fun3_ab.II,Fun3_a_b.Lambda,[1,2]),Fun3_b.I,[2,3]);
Fun3.Lambda_gamma_b=   multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.Lambda,[2,3]);
Fun3_ab.Yp=            (multiprod(Fun3_a.Wb_inv,multiprod(multiprod(Fun3_ab.A,Fun3_a_b.YV,[1,2]),Fun3_b.YV,[2,3]),[1,2])...
                       +multiprod(Fun3_a.Wb_inv,multiprod(multiprod(Fun3_ab.B,Fun3_a_b.Vb,[1,2]),Fun3_b.Vb,[2,3]),[1,2])...
                       -multiprod(multiprod(Fun3_a.Wb_inv,Fun3_a.Y,[1,2]),multiprod(multiprod(Fun3_ab.C,Fun3_a_b.YV,[1,2]),Fun3_b.Vb,[2,3]),[1,2]))...
                        ./(Fun3.Lambda_alpha_a+Fun3.Lambda_beta_a_b+Fun3.Lambda_gamma_b);
Fun3_ab.YY=            multiprod(Fun3_a.Wb,multiprod(multiprod(Fun3_ab.Yp,Fun3_a_b.Vb_inv,[1,2]),Fun3_b.Vb_inv,[2,3]),[1,2]);
Fun3_ab.YY_minus=      -Fun3_ab.YY;

%Warning:a-b!=0, be NaN, need to be deleted.



%Case2: Torsion Helical Duct
%a
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        [Fun2_a.Vb_H(:,:,kh,ka),Fun2_a.Lambda_H(:,:,kh,ka)]= eigs(( Fun2_a.L(:,:,kh,ka)),n_matrix*2); %sqrt(-1)*
        eigenvalues_a=diag(Fun2_a.Lambda_H(:,:,kh,ka));
        order_a=find(real(eigenvalues_a)<0);
        Fun2_a.Vb_HH(:,:,kh,ka)=Fun2_a.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order_a,kh,ka);
        Fun2_a.Vb_inv_HH(:,:,kh,ka)=inv(Fun2_a.Vb_HH(:,:,kh,ka));
        Fun2_a.Lambda_HH(:,:,kh,ka)=diag(eigenvalues_a(order_a));
    end
end
%b
for kh=1:length(Geo_b.h)
    for kb=1:length(b)
        [Fun2_b.Vb_H(:,:,kh,kb),Fun2_b.Lambda_H(:,:,kh,kb)]= eigs(( Fun2_b.L(:,:,kh,kb)),n_matrix*2); %sqrt(-1)*
        eigenvalues_b=diag(Fun2_b.Lambda_H(:,:,kh,kb));
        order_b=find(real(eigenvalues_b)<0);
        Fun2_b.Vb_HH(:,:,kh,kb)=Fun2_b.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order_b,kh,kb);
        Fun2_b.Vb_inv_HH(:,:,kh,kb)=inv(Fun2_b.Vb_HH(:,:,kh,kb));
        Fun2_b.Lambda_HH(:,:,kh,kb)=diag(eigenvalues_b(order_b));
    end
end
%a_b
temp1=(NaN+NaN*sqrt(-1))*ones(size(Fun2_a_b.N(:,:,1,1,1)));
temp2=(NaN+NaN*sqrt(-1))*ones(size(Fun2_a_b.L(:,:,1,1,1)));

for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        for kb=1:length(b)
            if a(ka)-b(kb)==0
                Fun2_a_b.Lambda_H(:,:,kh,kb,ka)=  temp2;
                Fun2_a_b.Vb_H(:,:,kh,kb,ka)=      temp2;
                Fun2_a_b.Vb_HH(:,:,kh,kb,ka)=     temp1;
                Fun2_a_b.Vb_inv_HH(:,:,kh,kb,ka)= temp1;
                Fun2_a_b.Vb_inv(:,:,kh,kb,ka)=    temp1;
                Fun2_a_b.Lambda_HH(:,:,kh,kb,ka)= temp1;
            else
                [Fun2_a_b.Vb_H(:,:,kh,kb,ka),Fun2_a_b.Lambda_H(:,:,kh,kb,ka)]= eigs((Fun2_a_b.L(:,:,kh,kb,ka)),n_matrix*2); %sqrt(-1)*
                eigenvalues_a_b=diag(Fun2_a_b.Lambda_H(:,:,kh,kb,ka));
                order_a_b=find(real(eigenvalues_a_b)<0);
                Fun2_a_b.Vb_HH(:,:,kh,kb,ka)=Fun2_a_b.Vb_H(length(Geo_b.m)*Geo_b.n+1:end,order_a_b,kh,kb,ka);
                Fun2_a_b.Vb_inv_HH(:,:,kh,kb,ka)=inv(Fun2_a_b.Vb_HH(:,:,kh,kb,ka));
                Fun2_a_b.Lambda_HH(:,:,kh,kb,ka)=diag(eigenvalues_a_b(order_a_b));
            end
        end
    end
end
Fun2_a.Y_H=    multiprod(Fun2_a.N_inv,multiprod(multiprod(Fun2_a.Vb_HH,Fun2_a.Lambda_HH,[1,2]),Fun2_a.Vb_inv_HH,[1,2])-Fun2_a.H,[1,2]);
Fun2_b.Y_H=    multiprod(Fun2_b.N_inv,multiprod(multiprod(Fun2_b.Vb_HH,Fun2_b.Lambda_HH,[1,2]),Fun2_b.Vb_inv_HH,[1,2])-Fun2_b.H,[1,2]);
Fun2_a_b.Y_H=  multiprod(Fun2_a_b.N_inv,multiprod(multiprod(Fun2_a_b.Vb_HH,Fun2_a_b.Lambda_HH,[1,2]),Fun2_a_b.Vb_inv_HH,[1,2])-Fun2_a_b.H,[1,2]);


Fun2_a.YNG=multiprod(Fun2_a.Y_H,Fun2_a.N,[1,2])+Fun2_a.G;     %James-3.52
Fun2_b.YNG=multiprod(Fun2_b.Y_H,Fun2_b.N,[1,2])+Fun2_b.G;
Fun2_a_b.YNG=multiprod(Fun2_a_b.Y_H,Fun2_a_b.N,[1,2])+Fun2_a_b.G;
for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        [Fun2_a.Wb_H(:,:,kh,ka),Fun2_a.Lambda1_H(:,:,kh,ka)]=eigs(Fun2_a.YNG(:,:,kh,ka),n_matrix);
         Fun2_a.Wb_inv_H(:,:,kh,ka)=inv(Fun2_a.Wb_H(:,:,kh,ka));
    end
end
for kh=1:length(Geo_b.h)
    for kb=1:length(b)
        [Fun2_b.Wb_H(:,:,kh,kb),Fun2_b.Lambda1_H(:,:,kh,kb)]=eigs(Fun2_b.YNG(:,:,kh,kb),n_matrix);
         Fun2_b.Wb_inv_H(:,:,kh,kb)=inv(Fun2_b.Wb_H(:,:,kh,kb));
    end
end

for kh=1:length(Geo_b.h)
    for ka=1:length(a)
        for kb=1:length(b)
            if a(ka)-b(kb)==0
                Fun2_a_b.Wb_H(:,:,kh,kb,ka)=temp1;
                Fun2_a_b.Lambda1_H(:,:,kh,kb,ka)=temp1;               
            else
                [Fun2_a_b.Wb_H(:,:,kh,kb,ka),Fun2_a_b.Lambda1_H(:,:,kh,kb,ka)]=eigs(Fun2_a_b.YNG(:,:,kh,kb,ka),n_matrix);
                 Fun2_a_b.Wb_inv_H(:,:,kh,kb,ka)=inv(Fun2_a_b.Wb_H(:,:,kh,kb,ka));
            end
        end
    end
end

%2D->3D
Fun3_a_b.Vb_H=           permute(Fun2_a_b.Vb_HH,[1,2,6,3,4,5]);

Fun3_b.Vb_H=             permute(bsxfun(@times,Fun2_b.Vb_HH,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a_b.Vb_inv_H=       permute(Fun2_a_b.Vb_inv_HH,[1,2,6,3,4,5]);
Fun3_b.Vb_inv_H=         permute(bsxfun(@times,Fun2_b.Vb_inv_HH,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);

Fun3_a_b.YV_H=            permute(multiprod(Fun2_a_b.Y_H,Fun2_a_b.Vb_HH,[1,2]),[1,2,6,3,4,5]);%!!!!!!!!!¾ùÎª0£¬bug
Fun3_b.YV_H=              permute(bsxfun(@times,multiprod(Fun2_b.Y_H,Fun2_b.Vb_HH,[1,2]),reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a.Wb_H=              permute(bsxfun(@times,Fun2_a.Wb_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Wb_inv_H=          permute(bsxfun(@times,Fun2_a.Wb_inv_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Y_H=               permute(bsxfun(@times,Fun2_a.Y_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a.Lambda_HH=         permute(bsxfun(@times,Fun2_a.Lambda_HH,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);
Fun3_a_b.Lambda_HH=       permute(Fun2_a_b.Lambda_HH,[1,2,6,3,4,5]);
Fun3_b.Lambda_HH=         permute(bsxfun(@times,Fun2_b.Lambda_HH,reshape(ones(size(a)),1,1,1,1,length(a))),[6,1,2,3,4,5]);
Fun3_a.Lambda1_H=         permute(bsxfun(@times,Fun2_a.Lambda1_H,reshape(ones(size(b)),1,1,1,1,length(b))),[1,2,6,3,5,4]);

Fun3.Lambda_H_alpha_a=    multiprod(Fun3_a.Lambda1_H,multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.I,[2,3]),[1,2]);
Fun3.Lambda_H_beta_a_b=   multiprod(multiprod(Fun3_ab.II,Fun3_a_b.Lambda_HH,[1,2]),Fun3_b.I,[2,3]);
Fun3.Lambda_H_gamma_b=    multiprod(multiprod(Fun3_ab.II,Fun3_a_b.I,[1,2]),Fun3_b.Lambda_HH,[2,3]);


Fun3_ab.Yp_H=            (multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.A,Fun3_a_b.YV_H,[1,2]),Fun3_b.YV_H,[2,3]),[1,2])...
                         +multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.B,Fun3_a_b.Vb_H,[1,2]),Fun3_b.Vb_H,[2,3]),[1,2])...
                         -multiprod(multiprod(Fun3_a.Wb_inv_H,Fun3_a.Y_H,[1,2]),multiprod(multiprod(Fun3_ab.C,Fun3_a_b.YV_H,[1,2]),Fun3_b.Vb_H,[2,3]),[1,2])...
                         -multiprod(multiprod(Fun3_a.Wb_inv_H,Fun3_a.Y_H,[1,2]),multiprod(multiprod(Fun3_ab.D,Fun3_a_b.YV_H,[1,2]),Fun3_b.YV_H,[2,3]),[1,2])...
                         +multiprod(Fun3_a.Wb_inv_H,multiprod(multiprod(Fun3_ab.E,Fun3_a_b.YV_H,[1,2]),Fun3_b.Vb_H,[2,3]),[1,2]))...
                         ./(Fun3.Lambda_H_alpha_a+Fun3.Lambda_H_beta_a_b+Fun3.Lambda_H_gamma_b);
Fun3_ab.YY_H=            multiprod(Fun3_a.Wb_H,multiprod(multiprod(Fun3_ab.Yp_H,Fun3_a_b.Vb_inv_H,[1,2]),Fun3_b.Vb_inv_H,[2,3]),[1,2]);


%Warning:a-b!=0, be NaN, need to be deleted.
disp('BFun..Finish!')
toc
%% OUTPUT:Y0,YY0
Out.Y0_a=Fun2_a.Y_H;
Out.Y0_b=Fun2_b.Y_H;
Out.Y0_a_b=Fun2_a_b.Y_H;
Out.YY0_ab=Fun3_ab.YY_H;
Out.I_b=Fun3_b.I(:,:,:,1,:,:);
Out.I_a_b=Fun3_a_b.I(:,:,:,1,:,:);

end
%% OUTPUT:N_a,M_a,H_a,G_a,N_a_b,N_b,H_a_b,H_b,A,B,C,D,E
Out.M_a=Fun2_a.M;
Out.M_b=Fun2_b.M;
Out.M_a_b=Fun2_a_b.M;
Out.N_a=Fun2_a.N;
Out.N_b=Fun2_b.N;
Out.N_a_b=Fun2_a_b.N;
Out.H_a=Fun2_a.H;
Out.H_b=Fun2_b.H;
Out.H_a_b=Fun2_a_b.H;
Out.G_a=Fun2_a.G;
Out.G_b=Fun2_b.G;
Out.G_a_b=Fun2_a_b.G;
Out.A=Fun3_ab.A;
Out.B=Fun3_ab.B;
Out.C=Fun3_ab.C;
Out.D=Fun3_ab.D;
Out.E=Fun3_ab.E;

Out.C(find(isnan(Out.C)))=0;
Out.D(find(isnan(Out.D)))=0;

end

