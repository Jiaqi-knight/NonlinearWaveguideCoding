function H=Hf(scheme,f)
H=sum(f.*log(f./permute(scheme(:,3),[2,3,1])),3);

end