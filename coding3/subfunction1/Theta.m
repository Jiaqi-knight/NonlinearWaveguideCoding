%THETA Function for \Theta functions.
%m:circumferiential mode
%n:radial mode
%dimension:
%   2D           -  \Theta_{\alpha\beta}
%   3D           -  \Theta_{\alpha\beta\gamma}
%THETA Options(op):
%   ab           -  \Theta_{\alpha\beta}
%   pt_ab        -  \Theta_{(\alpha)\beta}
%   ps_ab        -  \Theta_{\{\alpha\}\beta}
%   ab_cos    -  \Theta_{\alpha\beta}[cos\phi]
%   ab_sin    -  \Theta_{\alpha\beta}[sin\phi]
%   pt_ab_cos    -  \Theta_{\alpha\beta}[cos\phi]
%   pt_ab_sin    -  \Theta_{\alpha\beta}[sin\phi]
%   ps_ab_cos    -  \Theta_{\alpha\beta}[cos\phi]
%   ps_ab_sin    -  \Theta_{\alpha\beta}[sin\phi]
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.

function [T1]= Theta(Geo,dimention,op);
%tau=1;
tau=Geo.tau;
s=Geo.s;
m=Geo.m;
n=Geo.n;
if dimention==2
    switch op
        case 'ab'
            T1=bsxfun(@times,2*pi*deltaT(m,n,dimention,0),reshape(ones(size(s)),1,1,[]));
        case 'pt_ab'
            multiply_factor=m;
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(ones(size(s)),1,1,[]));
        case 'a_pt_b'
            multiply_factor=fliplr(m);
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(ones(size(s)),1,1,[]));
        case   'ps_ab'
            multiply_factor=m;
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(-tau,1,1,[]));
        case   'a_ps_b'
            multiply_factor=fliplr(m);
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(-tau,1,1,[]));
        case   'ab_cos'
            T1=bsxfun(@times,pi*(deltaT(m,n,dimention,1)+deltaT(m,n,dimention,-1)),reshape(ones(size(s)),1,1,[]));
        case   'ab_sin'
            T1=bsxfun(@times,-sqrt(-1)*pi*(deltaT(m,n,dimention,1)-deltaT(m,n,dimention,-1)),reshape(ones(size(s)),1,1,[]));
        case   'pt_ab_cos'
            multiply_factor=m;
            T1=bsxfun(@times,sqrt(-1)*pi*(deltaT(m,n,dimention,1,multiply_factor)+deltaT(m,n,dimention,-1,multiply_factor)),reshape(ones(size(s)),1,1,[]));
        case   'pt_ab_sin'
            multiply_factor=m;
            T1=bsxfun(@times,pi*(deltaT(m,n,dimention,1,multiply_factor)-deltaT(m,n,dimention,-1,multiply_factor)),reshape(ones(size(s)),1,1,[]));
        case   'ps_ab_cos'
            multiply_factor=m;
            T1=bsxfun(@times,sqrt(-1)*pi*(deltaT(m,n,dimention,1,multiply_factor)+deltaT(m,n,dimention,-1,multiply_factor)),reshape(-tau,1,1,[]));
        case   'ps_ab_sin'
            multiply_factor=m;
            T1=bsxfun(@times,pi*(deltaT(m,n,dimention,1,multiply_factor)-deltaT(m,n,dimention,-1,multiply_factor)),reshape(-tau,1,1,[]));
    end
elseif dimention==3
    switch op
        case 'ab'
            T1=bsxfun(@times,2*pi*deltaT(m,n,dimention,0),reshape(ones(size(s)),1,1,1,[]));
        case 'pt_ab'
            multiply_factor=m;
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(ones(size(s)),1,1,1,[]));
        case 'a_pt_b'
            multiply_factor=fliplr(m);
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(ones(size(s)),1,1,1,[]));
        case   'ps_ab'
            multiply_factor=m;
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(-tau,1,1,1,[]));
        case   'a_ps_b'
            multiply_factor=fliplr(m);
            T1=bsxfun(@times,2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor),reshape(-tau,1,1,1,[]));
        case   'ab_cos'
            T1=bsxfun(@times,pi*(deltaT(m,n,dimention,1)+deltaT(m,n,dimention,-1)),reshape(ones(size(s)),1,1,1,[]));
        case   'ab_sin'
            T1=bsxfun(@times,-sqrt(-1)*pi*(deltaT(m,n,dimention,1)-deltaT(m,n,dimention,-1)),reshape(ones(size(s)),1,1,1,[]));
        case   'pt_ab_cos'
            multiply_factor=m;
            T1=bsxfun(@times,sqrt(-1)*pi*(deltaT(m,n,dimention,1,multiply_factor)+deltaT(m,n,dimention,-1,multiply_factor)),reshape(ones(size(s)),1,1,1,[]));
        case   'pt_ab_sin'
            multiply_factor=m;
            T1=bsxfun(@times,pi*(deltaT(m,n,dimention,1,multiply_factor)-deltaT(m,n,dimention,-1,multiply_factor)),reshape(ones(size(s)),1,1,1,[]));
        case   'ps_ab_cos'
            multiply_factor=m;
            T1=bsxfun(@times,sqrt(-1)*pi*(deltaT(m,n,dimention,1,multiply_factor)+deltaT(m,n,dimention,-1,multiply_factor)),reshape(-tau,1,1,1,[]));
        case   'ps_ab_sin'
            multiply_factor=m;
            T1=bsxfun(@times,pi*(deltaT(m,n,dimention,1,multiply_factor)-deltaT(m,n,dimention,-1,multiply_factor)),reshape(-tau,1,1,1,[]));
    end
end