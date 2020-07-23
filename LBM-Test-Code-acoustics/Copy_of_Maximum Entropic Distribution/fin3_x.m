function fin=fin3_x(x,M);
% This is the external function which calculates
% the fin(x) in the special case of the Fourier moments.
% This is to be used with ME_DENS3.
%
x=x(:); lx=length(x); % x axis
xmin=x(1); xmax=x(lx); dx=x(2)-x(1);
%
fin=zeros(lx,M); %
fin(:,1)=ones(size(x)); % fi0(x)=1
w0=2*pi/(xmax-xmin);jw0x=(sqrt(-1)*w0)*x;
for n=2:M
fin(:,n)=exp(-(n-1)*jw0x);
end
return
end