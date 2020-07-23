function [X]=yexing2surf(Blade)
for k=1:length(Blade)
s(k)=size(Blade{1, k},1);
end
s1=sort(s);
s_max=s1(end-1);
s_min=min(s);
xx1=[];xx2=[];xx3=[];
for k=1:length(Blade)
    xyz=[];
    xyz=resample(Blade{1,k},s_min,size(Blade{1,k},1));
    xx1=[xx1 xyz(:,1)];
    xx2=[xx2 xyz(:,2)];
    xx3=[xx3 xyz(:,3)];
end
X={xx1;xx2;xx3};
end
