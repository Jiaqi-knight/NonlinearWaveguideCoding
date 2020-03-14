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

function [T1]= Theta(m,n,dimention,op);
tau=1
switch op
    case 'ab'
        T1=2*pi*deltaT(m,n,dimention,0);
    case 'pt_ab'
        multiply_factor=m;
        T1=2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor);
    case   'ps_ab'
        %global tau
        multiply_factor=m;
        T1=-tau*2*pi*sqrt(-1)*deltaT(m,n,dimention,0,multiply_factor);
    case   'ab_cos'
        T1=pi*(deltaT(m,n,dimention,1)+deltaT(m,n,dimention,-1));
    case   'ab_sin'
        T1=-sqrt(-1)*pi*(deltaT(m,n,dimention,1)-deltaT(m,n,dimention,-1));
    case   'pt_ab_cos'
        multiply_factor=m;    
        T1=sqrt(-1)*pi*(deltaT(m,n,dimention,1,multiply_factor)+deltaT(m,n,dimention,-1,multiply_factor));
    case   'pt_ab_sin'
        multiply_factor=m;
        T1=pi*(deltaT(m,n,dimention,1,multiply_factor)-deltaT(m,n,dimention,-1,multiply_factor));
    case   'ps_ab_cos'
        %global tau
        multiply_factor=m;
        T1=-tau*sqrt(-1)*pi*(deltaT(m,n,dimention,1,multiply_factor)+deltaT(m,n,dimention,-1,multiply_factor));
    case   'ps_ab_sin'
        %global tau
        multiply_factor=m;
        T1=-tau*pi*(deltaT(m,n,dimention,1,multiply_factor)-deltaT(m,n,dimention,-1,multiply_factor));
        
end

end