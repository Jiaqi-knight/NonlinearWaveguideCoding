function buckle_wrinkle (eigenvalue,ne,ng,c,p,gfl,gfln,ngm1,mat)

%====================================
% compute and plot the eigenfunctions
%====================================

mat1 = mat - eigenvalue*eye(ngm1);

for i=1:ngm1-1
 rhs2(i) = - mat1(i,ngm1);
end

%----
% formulate the eigenvector system
%---

for i=1:ngm1-1
 for j=1:ngm1-1
  mat2(i,j) = mat1(i,j);
 end
end

eigenvector = rhs2/mat2';
eigenvector(ngm1) = 1.0;

%---
% extract transverse displacement
% at the vertex nodes
%---

Ic = 0;

for i=1:ng

 if(gfl(i,1) == 0) % interior node
   Ic = Ic+1;
   f(i) = eigenvector(Ic);
   if(gfln(i) == 0)  % vertex node
     Ic=Ic+2;
   end
 else
   f(i) = 0; 
 end

end

%-------------------------------
% interpolate the midside values
%-------------------------------

for i=1:ne
  f(c(i,4))=0.5*( f(c(i,1))+f(c(i,2)) );
  f(c(i,5))=0.5*( f(c(i,2))+f(c(i,3)) );
  f(c(i,6))=0.5*( f(c(i,3))+f(c(i,1)) );
end

%-----
% plot
%-----
                                                                              
% plot_6 (ne,ng,p,c,f);

%-----
% generate the extended connectivity matrix
% for 3-node sub-triangles
%-----

Ic=0;

for i=1:ne
 Ic=Ic+1;
 c3(Ic,1)=c(i,1); c3(Ic,2)=c(i,4); c3(Ic,3)=c(i,6);
 Ic=Ic+1;
 c3(Ic,1)=c(i,4); c3(Ic,2)=c(i,2); c3(Ic,3)=c(i,5);
 Ic=Ic+1;
 c3(Ic,1)=c(i,5); c3(Ic,2)=c(i,3); c3(Ic,3)=c(i,6);
 Ic=Ic+1;
 c3(Ic,1)=c(i,4); c3(Ic,2)=c(i,5); c3(Ic,3)=c(i,6);
end
                                                                              
trimesh (c3,p(:,1),p(:,2),f);

%-----
% done
%-----

return
