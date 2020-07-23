function f = Stream(f,N,M,K,xi)

% Create temporary grid to prevent overwriting useful populations
fnew = zeros(size(f));

% Stream one lattice site at a time
for i = 1:N
    for j = 1:M
        for k = 1:K

            for n = 1:6
            
                % Compute destination coordinates
                dest_x = 1 + mod(i-1+xi(1,n)+N,N);
                dest_y = 1 + mod(j-1+xi(2,n)+M,M);
                dest_z = 1 + mod(k-1+xi(3,n)+K,K);

                % Stream values
                fnew(n,dest_z,dest_y,dest_x) = f(n,k,j,i);
                
            end

        end
    end
end

fnew(7,:,:,:) = f(7,:,:,:); % Rest particle distributions
f = fnew;

end