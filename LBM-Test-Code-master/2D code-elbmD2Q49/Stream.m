function f = Stream(f,N,M,c)

% Create temporary grid to prevent overwriting useful populations
fnew = zeros(size(f));

% Stream one lattice site at a time
for i = 1:N
    for j = 1:M

        for n = 2:25
            
            % Compute destination coordinates
            dest_x = 1 + mod(i-1+c(n,1)+N,N);
            dest_y = 1 + mod(j-1+c(n,2)+M,M);
            
            % Stream values
            fnew(dest_y,dest_x,n) = f(j,i,n);


        end
    end
end

fnew(:,:,1) = f(:,:,1); % Rest particle distributions
f = fnew;

end