function [fnew, feq] = Collide(T0,rho,w,c,u,cs,feq,N,M,f,omega)
% Loop through the lattice points to compute the new ditribution functions.
% Equilibrium based on:
%       rho * w * (1 + c_ia u_a / cs^2 + Q_iab u_a u_b / 2*cs^4

for i = 1:N
    for j = 1:M
        for k = 1:25
            % Compute c_ia * u_a which is actually the dot product of c and
            % u
            A = dot(c(k,:).',reshape(u(j,i,:),2,1));
            B=0;
            for alpha=1:2
                for beta=1:2
                    B=B+u(j,i,alpha)*u(j,i,beta)*(c(k,alpha)*c(k,beta)-T0*(alpha==beta));
                end
            end
            C=0;
            for alpha=1:2
                for beta=1:2
                    for gamma=1:2
                        C=C+u(j,i,alpha)*u(j,i,beta)*u(j,i,gamma)*c(k,gamma)*(c(k,alpha)*c(k,beta)-3*T0*(alpha==beta));
                    end
                end
            end
            
           
            
            % Compute f^eq
            feq(j,i,k) = rho(j,i) * w(k) * ...
                (1 + (A / T0) + (B / (2*T0^2))+ (C / (6*T0^3)));
            
            % Alternative form for feq I have seen published
            %             feq2(j,i,k) = rho(j,i) * w(k) * ...
            %                 ( 1 + ...
            %                 ((3/2)*dot(c(:,k),reshape(u(j,i,:),2,1))) + ...
            %                 ((9/2)*(dot(c(:,k),reshape(u(j,i,:),2,1))^2)) - ...
            %                 ((3/2)*norm(reshape(u(j,i,:),2,1))^2) );
            
        end
        
    end
end

% Recompute distribution function f
fnew = (omega * (feq -f ) ) + f ;

end
