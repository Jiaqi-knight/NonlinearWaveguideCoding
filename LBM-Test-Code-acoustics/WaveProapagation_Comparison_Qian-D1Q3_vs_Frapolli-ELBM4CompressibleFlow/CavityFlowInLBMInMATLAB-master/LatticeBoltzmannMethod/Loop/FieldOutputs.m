function file_name = FieldOutputs(rho, u, t, name)
% This function designed to outputs the whole field in Tecplot format
% Input field varibles rho(Nx+1,Ny+1), u(Nx+1,Ny+1,2), and time step t
% and same file name 
global Nx Ny Lx Ly
%% Output each file's name
name_temp = [name,'_',num2str(t),'.dat'] ;  % output file name
fp = fopen(name_temp,'w') ;                   % file pointer
%% print file's head
fprintf(fp,'Titile = "LBM Lid Drived Flow" \n');
fprintf(fp,'VARIABLES = "X", "Y", "Rho", "U", "V" \n');
fprintf(fp,'zone I=%d, J=%d, F = POINT \n',Nx+1,Ny+1);

%% print data
% This loop is for 2 dimension case
for i = 1:Nx+1
     for j = 1:Ny+1
           fprintf(fp,'%.15E   %.15E   %.15E   %.15E   %.15E   \n',...
                   i/Lx,j/Ly,rho(i,j),u(i,j,1),u(i,j,2));
     end
end

fclose(fp);            % close file pointer
file_name = name_temp; % output file's full name to main
end