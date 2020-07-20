function [feq] = Equilibrium(P,w,xi,J,c,cs,feq,N,M,K,type)
% Loop through the lattice points to compute equilibrium based on
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4
if type==0 %moment matching
    for i = 1:N
        for j = 1:M
            for k = 1:K
                for v = 7
                    feq(v,k,j,i) = P(k,j,i) + (c^2/ cs^2) * P(k,j,i) * (w(v) -1);
                end
                for v = 1:6
                    
                    % Compute c_ia * u_a which is actually the dot product of xi
                    % and J
                    A = dot( xi(:,v), J(:,k,j,i) );
                    
                    % Compute Q_iab u_a u_b =
                    %          (c_a^2 - cs^2)u_a^2 + 2c_a c_b u_a u_b +
                    %                   (c_b^2 - cs^2)u_b^2
                    %                 B = (xi(1,v)^2 - cs^2) * J(1,k,j,i)^2 + ...
                    %                     (xi(2,v)^2 - cs^2) * J(2,k,j,i)^2 + ...
                    %                     (xi(3,v)^2 - cs^2) * J(3,k,j,i)^2 + ...
                    %                     2*xi(1,v)*xi(2,v)*J(1,k,j,i)*J(2,k,j,i) + ...
                    %                     2*xi(1,v)*xi(3,v)*J(1,k,j,i)*J(3,k,j,i) + ...
                    %                     2*xi(2,v)*xi(3,v)*J(2,k,j,i)*J(3,k,j,i);
                    
                    % Compute f^eq
                    feq(v,k,j,i) = (w(v)/cs^2)*( c^2 * P(k,j,i) + A ) ;
                    
                end
                
            end
        end
    end
    
elseif type==1 %Hermite Polynomials
    for i = 1:N
        for j = 1:M
            for k = 1:K
                for v = 7
                    feq(v,k,j,i) = P(k,j,i) -(2.5*P(k,j,i)-1.5*c^2/ cs^2*P(k,j,i))*(1-w(v))+ 1.5*(c^2/ cs^2) * P(k,j,i) * (c^2-cs^2);
                end
                for v = 1:6
                    
                    % Compute c_ia * u_a which is actually the dot product of xi
                    % and J
                    A = dot( xi(:,v), J(:,k,j,i) );
                    B = dot( xi(:,v), xi(:,v) );
                    % Compute Q_iab u_a u_b =
                    %          (c_a^2 - cs^2)u_a^2 + 2c_a c_b u_a u_b +
                    %                   (c_b^2 - cs^2)u_b^2
                    %                 B = (xi(1,v)^2 - cs^2) * J(1,k,j,i)^2 + ...
                    %                     (xi(2,v)^2 - cs^2) * J(2,k,j,i)^2 + ...
                    %                     (xi(3,v)^2 - cs^2) * J(3,k,j,i)^2 + ...
                    %                     2*xi(1,v)*xi(2,v)*J(1,k,j,i)*J(2,k,j,i) + ...
                    %                     2*xi(1,v)*xi(3,v)*J(1,k,j,i)*J(3,k,j,i) + ...
                    %                     2*xi(2,v)*xi(3,v)*J(2,k,j,i)*J(3,k,j,i);
                    
                    % Compute f^eq
                    feq(v,k,j,i) = w(v)*(P(k,j,i) + A/cs^2 +P(k,j,i)/2/cs^4*(c^2-cs^2)*(B-3*cs^2));
                    
                end
                
            end
        end
    end
end

end

