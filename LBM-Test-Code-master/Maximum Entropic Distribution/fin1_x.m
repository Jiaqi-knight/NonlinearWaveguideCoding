function fin=fin1_x(x);
% This is the external function which calculates
% the fin(x) in the special case of the Gamma distribution.
% This is to be used with ME_dens1.
M=3;
fin=zeros(length(x),M);
fin(:,1)=ones(size(x));
fin(:,2)=x;
fin(:,3)=log(x);
return
end
