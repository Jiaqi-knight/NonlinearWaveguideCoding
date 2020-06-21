function X=arcPlot(Arc_s,Arc_e,O)
R1=Arc_s-O;
R2=Arc_e-O;
theta=subspace(R1,R2);
R=[];
R=[R1,R];
X=[];
W=cross(R1,R2);
th=linspace(0,theta,20);
dt=(th(2)-th(1))/norm(W);
for m=1:length(th)
    Rota=R(:,end)+cross(W,R(:,end))*dt;
    R=[R,Rota];
    X=[X,O+Rota];
end
X=[Arc_s,X,Arc_e];
end
