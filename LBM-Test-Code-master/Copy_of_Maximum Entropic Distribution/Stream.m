function fnew = Stream(Q,f,N,c)
for k=1:Q
fnew(:,k)=circshift(f(:,k),c(k));
end
fnew(1:2,:)=f(1:2,:);
fnew(end-1:end,:)=f(end-1:end,:);
end