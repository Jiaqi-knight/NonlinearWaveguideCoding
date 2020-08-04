%A function to propagate the populations
function populations = LBMPropagate(oldpopulations,nx,ny,scheme)
populations = oldpopulations;
c_x=scheme(:,1);
c_y=scheme(:,2);

for k=1:length(c_x)
    if c_x(k)>0
        %         up=populations(nx-c_x(k)+1:nx,:,k); %periodic condition
        populations(c_x(k)+1:nx,:,k)   = populations(1:nx-c_x(k),:,k);
        %         populations(1:c_x(k),:,k)  =  left;
    elseif c_x(k)<0
        %         right=populations(1:-c_x(k),:,k);
        populations(1:nx+c_x(k),:,k) = populations(1-c_x(k):nx,:,k);
        %         populations(nx+c_x(k)+1:nx,:,k) = right;
    end
end
for k=1:length(c_y)
    if c_y(k)>0
        %         left=populations(:,nx-c_y(k)+1:nx,k); %periodic condition
        populations(:,c_y(k)+1:ny,k)   = populations(:,1:ny-c_y(k),k);
        %         populations(:,1:c_y(k),k)  =  left;
    elseif c_y(k)<0
        %         right=populations(:,1:-c_y(k),k);
        populations(:,1:ny+c_y(k),k) = populations(:,1-c_y(k):ny,k);
        %         populations(:,nx+c_y(k)+1:nx,k) = right;
    end
end
% void advect(){
% //first initialize buffers
% 
% % horizontal faces
% for (int ix=0;ix< static_cast<int>(l.nx);++ix){
%     for (unsigned int m=0; m<velocity_set().size; ++m){
%         l.f[m][l.index(ix,-1)] = l.f[m][l.index(ix,0)];
%         l.f[m][l.index(ix,l.ny)] = l.f[m][l.index(ix,l.ny-1)];
%         }
%         }
%         % vertical faces faces
%         for (int iy=0;iy < static_cast<int>(l.ny);++iy){
%             for (unsigned int m=0; m < velocity_set().size; ++m){
%                 l.f[m][l.index(-1,iy)] = l.f[m][l.index(l.nx-1,iy)];
%                 l.f[m][l.index(l.nx,iy)] = l.f[m][l.index(0,iy)];
%                 }
%                 }
%                 /*      % corners
%                 for (unsigned int m=0; m < velocity_set().size; ++m){
%                     l.f[m][l.index(-1,-1)] = l.f[m][l.index(l.nx-1,l.ny-1)];
%                     l.f[m][l.index(l.nx,l.ny)] = l.f[m][l.index(0,0)];
%                     l.f[m][l.index(l.nx,-1)] = l.f[m][l.index(0,l.ny-1)];
%                     l.f[m][l.index(-1,l.ny)] = l.f[m][l.index(l.nx-1,0)];
%                     }*/
%                     //shift populations in the opposite directions to the velocities so that not the same population is copied across the whole lattice
%                     for (unsigned int m=0; m<velocity_set().size; ++m){
%                         if (shift[m] > 0){
%                             for (unsigned int n = l.index(l.nx-1,l.ny-1); n>=l.index(0,0);--n){
%                                 l.f[m][n] = l.f[m][n-shift[m]];
%                                 }
%                                 }
%                                 else if (shift[m] < 0){
%                                         for (unsigned int n = l.index(0,0);n<=l.index(l.nx-1,l.ny-1);++n){
%                                             l.f[m][n] = l.f[m][n-shift[m]];
%                                             }
%                                             }
%                                             }
%                                             
%                                             return ;
%                                             }
