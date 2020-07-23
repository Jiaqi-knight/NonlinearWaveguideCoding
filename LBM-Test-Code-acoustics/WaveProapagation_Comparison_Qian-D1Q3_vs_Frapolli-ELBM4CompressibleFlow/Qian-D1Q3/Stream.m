function fnew = Stream(f,N,c)

fnew(:,1)=circshift(f(:,1),c(1));
fnew(:,2)=circshift(f(:,2),c(2));
fnew(:,3)=circshift(f(:,3),c(3));
% fnew(1:2,:)=f(1:2,:);
% fnew(end-1:end,:)=f(end-1:end,:);
end