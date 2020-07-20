function [fnew, feq] = Collide(P,w,xi,J,c,cs,feq,N,M,K,f,omega)
% Compute the new ditribution functions using LBKG SRT collision

% Update equilibrium
feq = Equilibrium(P,w,xi,J,c,cs,feq,N,M,K,0);

% Recompute distribution function f
fnew = (-omega * (f - feq) ) + f ; % + force_i;
    
end
    
