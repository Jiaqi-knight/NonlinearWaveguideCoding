function[f]=stream3D(N,Q,fold,lattice)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: May 30th, 2014
%    Last update: June 24th, 2014
%
%    Description: 
%          Input: 
%         Output: 
%     f(:,:,1) = [f(:,1:2,1) f(:,2:Mc-1,1)];
%     f(:,:,2) = [f(1:2,:,2);f(2:Nr-1,:,2)];
%     f(:,:,3) = [f(:,2:Mc-1,3) f(:,Mc-1:Mc,3)];
%     f(:,:,4) = [f(2:Nr-1,:,4);f(Nr-1:Nr,:,4)];
%     f(:,:,5) = [f(:,1:2,5) f(:,2:Mc-1,5)];
%     f(:,:,5) = [f(1:2,:,5);f(2:Nr-1,:,5)];
%     f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
%     f(:,:,6) = [f(1:2,:,6);f(2:Nr-1,:,6)];
%     f(:,:,7) = [f(:,2:Mc-1,7) f(:,Mc-1:Mc,7)];
%     f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)];
%     f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
%     f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)];
%%

     
            