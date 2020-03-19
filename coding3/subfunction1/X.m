%X Function for \Theta functions.
%m:circumferiential mode
%n:radial mode
%dimension:
%   2D           -  \X_{\alpha\beta}
%   3D           -  \X_{\alpha\beta\gamma}
%THETA Options(op):
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.


function [X2]= X(s,h,Cmn1,jmn_pm,m,n,dimention,op,op_m);
        X2=deltaJ(s,h,Cmn1,jmn_pm,m,n,dimention,0,op,op_m)...
            +deltaJ(s,h,Cmn1,jmn_pm,m,n,dimention,1,op,op_m)...
            +deltaJ(s,h,Cmn1,jmn_pm,m,n,dimention,-1,op,op_m);
end