%% erfenfa 


function [alpha]=erfenfa(scheme,populations,popequilibriums,bot,top,err)
%输入：f为要求解的方程函数表达式，bot为求解区间的下界，top为求解区间的上界，err为所要求的误差范围
%输出：i为二分法求解的次数，x为最终的根，fx为方程的根x对应的函数值（精度要求内接近于零或等于零）

Finf_=F_eval(scheme,populations,popequilibriums,bot);
Fsup_=F_eval(scheme,populations,popequilibriums,top);
val_=Finf_.*Fsup_;

alpha(val_>0)=2;
index=find(val_<=0);
i=ceil((log(top-bot)- log(err))/log(2))-1; %n为二分法运算总的次数；ceil是上取整 ,相对应的floor是下取整
for k=1:length(index)
top1=top; bot1=bot;   
while abs(top1-bot1)>err
    x=(bot1+top1)/2;
    fx=F_eval(scheme,populations(:,index(k)),popequilibriums(:,index(k)),x);
    if fx==0
        bot1=x; top1=x;
    elseif F_eval(scheme,populations(:,index(k)),popequilibriums(:,index(k)),bot1)*fx<0
        top1=x;
    else
        bot1=x;
    end
end
alpha(index(k))=x;
end
% fprintf('\nThe result:\n二分法运算次数i=%d;方程的根x=%.4f;f(x)=%.4f\n',i,x,fx);%保留4位小数
end