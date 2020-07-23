function [Base]=BaseJ1(m,n,h)

%% use +m to generate -m
% m0=m((length(m)+1)/2:end);
% m>=0

for k=1:length(m)
    temp1=roots(diff(chebfun(@(t) besselj(sign(m(k)).*m(k),t),[0,300]))+0.00001);
    temp2= sign(m(k)).^(sign(m(k)).*m(k))* besselj(sign(m(k)).*m(k),temp1);
    Base.jmn_pm(:,k)=temp1(1:n);
    BesselValue_pm(:,k)=temp2(1:n);
end

% repmat((-1).^m,n,1)
% Base.jmn_pm=[fliplr(jmn) jmn(:,2:end)];
% BesselValue_pm=[fliplr(repmat((-1).^m0,n,1).*BesselValue) BesselValue(:,2:end)];

Base.Cmn1=bsxfun(@times,(sqrt(-1)).^m./(sqrt(pi)*BesselValue_pm.*sqrt(1-m.^2./Base.jmn_pm.^2)),reshape(1./h,1,1,[]));%£¨Jiaqi-130£©4*11*50

end