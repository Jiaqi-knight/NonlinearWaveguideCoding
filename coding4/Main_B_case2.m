%Main-Boundary
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.
%-----------------------------------------------------------------%

% clc
% clear
% close all
% subfunction_path1='.\subfunction1';
% addpath(genpath(subfunction_path1));
% load('Data_m3_n1_a2_b6.mat');


%% #######Boudnary########%
%Infinite Uniform Duct Outlet
%Case2: Torsion Helical Duct

%which will suitable for the Case 1
%with H=G=0;
%L=(-G -M;N H);
%L^2=((G^2-MN) (GM-MH);(HN-NG) (H^2-NM))


%question1: how to reconstruct the block matrix of nD matrix
size([Fun2_a.N Fun2_a.N;Fun2_a.N Fun2_a.N])




Bdry_a.N=Fun2_a.N(:,:,end,:);
Bdry_a_b.N=Fun2_a_b.N(:,:,end,:,:);
Bdry_b.N=Fun2_b.N(:,:,end,:);

Bdry_a.N_inv=Fun2_a.N_inv(:,:,end,:);
Bdry_a_b.N_inv=Fun2_a_b.N_inv(:,:,end,:,:);
Bdry_b.N_inv=Fun2_b.N_inv(:,:,end,:);

Bdry_a.NM=multiprod(Fun2_a.N(:,:,end,:),Fun2_a.M(:,:,end,:),[1,2]);
Bdry_a_b.NM=multiprod(Fun2_a_b.N(:,:,end,:,:),Fun2_a_b.M(:,:,end,:,:),[1,2]);
Bdry_b.NM=multiprod(Fun2_b.N(:,:,end,:),Fun2_b.M(:,:,end,:),[1,2]);



G=Fun2_a.G(:,:,end,:);
M=Fun2_a.M(:,:,end,:)
N=Fun2_a.N(:,:,end,:)
H=Fun2_a.H(:,:,end,:)

Bdry_a.L=[-G -M;...
           N, H];
Bdry_a.L2=[multiprod(G,G)-multiprod(M,N) multiprod(G,M)-multiprod(M,H);...
          multiprod(H,N)-multiprod(N,G) multiprod(H,H)-multiprod(N,M)];%Bdry_a.L2=multiprod(Bdry_a.L,Bdry_a.L,[1,2]);

for ka=1:length(a)
  
[Bdry_a.V_H(:,:,1,ka),Bdry_a.Lambda_H(:,:,1,ka)]=eig(sqrt(-1)*Bdry_a.L(:,:,1,ka));
      
end

figure;
Lambda_H=diag(Bdry_a.Lambda_H(:,:,1,1));
plot(real(Lambda_H),imag(Lambda_H),'.')
