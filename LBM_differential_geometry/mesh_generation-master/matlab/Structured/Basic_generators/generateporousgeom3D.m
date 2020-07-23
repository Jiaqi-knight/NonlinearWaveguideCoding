function[conf]=generateporousgeom3D(Nx,Ny,Nz,D,lattice,T,J,h,initconf,tmax,dtsave,date,baseoutfilename,latticefolder,meshfolder)

%%
%         Course: Computational Statistical Physics
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 30th, 2014
%    Last update: July 30th, 2014
% 
%%

N = Nx*Ny*Nz;

beta = 1/T;

coordnum = 2*D;
neighnconf = 2^coordnum;

[indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8] = getindices3D(Nx,Ny,Nz);

periodicity = [1;2;3];
flagintbounds = 0;
indicesintbounds = 0;
typeintbounds = 0;
neighbours = build_neighbourhoods3D(N,Nx,Ny,Nz,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

clear indicesbulk indicesinternalbulk indicesF1 indicesF2 indicesF3 indicesF4 indicesF5 indicesF6 indicesinternalF1 indicesinternalF2 indicesinternalF3 indicesinternalF4 indicesinternalF5 indicesinternalF6 indicesE1 indicesE2 indicesE3 indicesE4 indicesE5 indicesE6 indicesE7 indicesE8 indicesE9 indicesE10 indicesE11 indicesE12 indicesinternalE1 indicesinternalE2 indicesinternalE3 indicesinternalE4 indicesinternalE5 indicesinternalE6 indicesinternalE7 indicesinternalE8 indicesinternalE9 indicesinternalE10 indicesinternalE11 indicesinternalE12 indicesC1 indicesC2 indicesC3 indicesC4 indicesC5 indicesC6 indicesC7 indicesC8 indicesinternalC1 indicesinternalC2 indicesinternalC3 indicesinternalC4 indicesinternalC5 indicesinternalC6 indicesinternalC7 indicesinternalC8

neighboursconf = zeros(neighnconf,coordnum);
for i=1:neighnconf
    tempconf = dectobin(coordnum,i-1);
    for j=1:coordnum
        neighboursconf(i,j) = (-1)^(1-tempconf(j));
    end
end

deltaE = [2*(-J*sum(neighboursconf,2)-h) -2*(-J*sum(neighboursconf,2)-h)];

clear neighboursconf

paccept = ones(neighnconf,2);

[rows,cols] = find(deltaE>0);
for i=1:length(rows)
    paccept(rows(i),cols(i)) = exp(-beta*deltaE(rows(i),cols(i)));
end

conf = sparse(dectobin(N,initconf));

for i=1:tmax
    for j=1:10
        randi(N,1);
    end
    site = randi(N,1);
    spin = (-1).^(1-conf(site,1));
    if spin==-1
        index = 1;
    else
        index = 2;
    end
    for j=1:10
        rand(1);
    end
    p = rand(1);
    neighbourindex = 1+bintodec(conf(neighbours(site,2:7),1));
    if p<paccept(neighbourindex,index)
        conf(site) = 1-conf(site);
    end
    if mod(i,dtsave)==0
        flags = [1;1;1];
        saveporousgeom3Dtovtk(date,i,tmax,baseoutfilename,latticefolder,meshfolder,flags,N,Nx,Ny,Nz,lattice,(-1).^(1-conf));      
        
%         figure();
%         plot3(lattice(find(conf==0),1),lattice(find(conf==0),2),lattice(find(conf==0),3),'.b')
%         hold on
%         plot3(lattice(find(conf==1),1),lattice(find(conf==1),2),lattice(find(conf==1),3),'.g')
%         hold on
%         grid on
%         xlabel('$\xi_{1}$','Interpreter','LaTex')
%         ylabel('$\xi_{2}$','Interpreter','LaTex')
%         zlabel('$\xi_{3}$','Interpreter','LaTex')
%         title(strcat('Lattice domain - computational space ',num2str(i),' iterations'))
%         axis equal
% 
%         figure();
%         plot3(lattice(find(conf==0),7),lattice(find(conf==0),8),lattice(find(conf==0),9),'.b')
%         hold on
%         plot3(lattice(find(conf==1),7),lattice(find(conf==1),8),lattice(find(conf==1),9),'.g')
%         hold on
%         grid on
%         xlabel('x  [m]')
%         ylabel('y  [m]')
%         zlabel('z  [m]')
%         title(strcat('Physical domain ',num2str(i),' iterations'))
%         axis equal
    end
end

return