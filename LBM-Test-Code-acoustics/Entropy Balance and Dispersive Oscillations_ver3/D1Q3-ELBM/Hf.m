function H=Hf(scheme,f)
H=sum(f.*log(f./scheme(:,2)),1);

end