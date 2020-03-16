% Sn1
% function dy=myode(t,y)
% N=200;
% n=26;
% Sn01=xlsread('Sn.xlsx', 1, 'Z1:Z201');
% for i=1:N+1           
%             be_n=2*n-1;
%                     
%             dy=zeros(4,1);
%             dy(1)=y(2);
%             dy(2)=y(3);
%             dy(3)=y(4);
%             dy(4)=Sn01(i)-2*y(4)./t+(1+2*be_n^2)*y(3)./t^2-...
%                 (1+2*be_n^2)*y(2)./t^3-(be_n^4-4*be_n^2)*y(1)./t^4;
%  
% end R=rk4('f','g','k','p',0.5,1,0,0,0,0,0.005)

function R=rk4(f,g,k,p,t0,tf,xa,ya,za,ua,h)
% clc;
% clear all;

b=1;c=sqrt(2);

n=26;

Sn01=xlsread('Sn.xlsx', 2, 'Z1:Z401');

be_n=2*n-1;

N=200;

f=@(t,x,y,z,u)(y);
g=@(t,x,y,z,u)(z);
k=@(t,x,y,z,u)(u);
p=@(t,x,y,z,u)(-2./t.*u+(1+2*be_n^2)./t^2.*z-(1+2*be_n^2)./t^3.*y-(be_n^4-4*be_n^2)./t^4.*x);

T=zeros(1,N);
X=zeros(1,N);
Y=zeros(1,N);
Z=zeros(1,N);
U=zeros(1,N);

T=b:(c-b)/N:c;

X(1)=xa;
Y(1)=ya;
Z(1)=za;
U(1)=ua;

for j=1:N
    J=2*j;
    JJ=2*j+1;
    

    f1=feval(f,T(j),X(j),Y(j),Z(j),U(j));
    g1=feval(g,T(j),X(j),Y(j),Z(j),U(j));
    k1=feval(k,T(j),X(j),Y(j),Z(j),U(j));
    p1=Sn01(j)+feval(p,T(j),X(j),Y(j),Z(j),U(j));
    
    f2=feval(f,T(j)+h/2,X(j)+h/2*f1,Y(j)+h/2*g1,Z(j)+h/2*k1,U(j)+h/2*p1);
    g2=feval(g,T(j)+h/2,X(j)+h/2*f1,Y(j)+h/2*g1,Z(j)+h/2*k1,U(j)+h/2*p1);
    k2=feval(k,T(j)+h/2,X(j)+h/2*f1,Y(j)+h/2*g1,Z(j)+h/2*k1,U(j)+h/2*p1);
    p2=Sn01(J)+feval(p,T(j)+h/2,X(j)+h/2*f1,Y(j)+h/2*g1,Z(j)+h/2*k1,U(j)+h/2*p1);
    
    f3=feval(f,T(j)+h/2,X(j)+h/2*f2,Y(j)+h/2*g2,Z(j)+h/2*k2,U(j)+h/2*p2);
    g3=feval(g,T(j)+h/2,X(j)+h/2*f2,Y(j)+h/2*g2,Z(j)+h/2*k2,U(j)+h/2*p2);
    k3=feval(k,T(j)+h/2,X(j)+h/2*f2,Y(j)+h/2*g2,Z(j)+h/2*k2,U(j)+h/2*p2);
    p3=Sn01(J)+feval(p,T(j)+h/2,X(j)+h/2*f2,Y(j)+h/2*g2,Z(j)+h/2*k2,U(j)+h/2*p2);
    
    f4=feval(f,T(j)+h,X(j)+h*f3,Y(j)+h*g3,Z(j)+h*k3,U(j)+h*p3);
    g4=feval(g,T(j)+h,X(j)+h*f3,Y(j)+h*g3,Z(j)+h*k3,U(j)+h*p3);
    k4=feval(k,T(j)+h,X(j)+h*f3,Y(j)+h*g3,Z(j)+h*k3,U(j)+h*p3);
    p4=Sn01(JJ)+feval(p,T(j)+h,X(j)+h*f3,Y(j)+h*g3,Z(j)+h*k3,U(j)+h*p3);
    
    X(j+1)=X(j)+h*(f1+2*f2+2*f3+f4)/6;
    Y(j+1)=Y(j)+h*(g1+2*g2+2*g3+g4)/6;
    Z(j+1)=Z(j)+h*(k1+2*k2+2*k3+k4)/6;
    U(j+1)=U(j)+h*(p1+2*p2+2*p3+p4)/6;

end
figure(1)
plot(T,X,'r',T,Y,'c',T,Z,'m',T,U,'y');
figure(2)
plot(T,X,'r',T,Y,'c',T,Z,'m');
figure(3)
plot(T,X,'r',T,Y,'c');
figure(4)
plot(T,X,'r');
grid on;

xlswrite('Tnp-2.xlsx',X',1,'Z1');
xlswrite('Tnp-2.xlsx',Y',2,'Z1');
xlswrite('Tnp-2.xlsx',Z',3,'Z1');
xlswrite('Tnp-2.xlsx',U',4,'Z1');
