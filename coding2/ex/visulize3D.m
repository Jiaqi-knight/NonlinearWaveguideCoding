clc
clear
close all
%% step1.3-Theta 
m1=-5:5;
m2=fliplr(m1);
m3=fliplr(m1);
for k1=1:length(m1)
    for k2=1:length(m2)
        for k3=1:length(m1)
            D(k1,k2,k3)=m1(k1)+m2(k2)+m3(k3);
        end
    end
end

order_0=D==0;order_N0=D~=0;
D(order_0)=1;D(order_N0)=0;

figure;h = slice(D, [], [], 1:size(D,3));
set(h, 'EdgeColor','none')%, 'FaceColor','interp'
alpha(.2);axis equal

D1=repelem(D,5,5,5);

figure;h = slice(D1, [], [], 1:size(D1,3));
set(h, 'EdgeColor','none')%, 'FaceColor','interp'
alpha(.2);axis equal

