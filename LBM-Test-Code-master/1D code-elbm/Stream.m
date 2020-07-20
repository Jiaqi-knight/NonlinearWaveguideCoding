function f = Stream(f,N,c)

% Create temporary grid to prevent overwriting useful populations
fnew = zeros(size(f));

% Stream one lattice site at a time
for i = 1:N
        for n = 2:5
            
            % Compute destination coordinates
            dest_x = 1 + mod(i-1+c(n)+N,N);            
            % Stream values
            fnew(dest_x,n) = f(i,n);

        end
end

fnew(1:10,:) = f(1:10,:);
fnew(end-9:end,:) = f(end-9:end,:);

fnew(:,1) = f(:,1); % Rest particle distributions
f = fnew;

end