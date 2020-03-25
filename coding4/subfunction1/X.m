%X Function for \Theta functions.
%m:circumferiential mode
%n:radial mode
%dimension:
%   2D           -  \X_{\alpha\beta}
%   3D           -  \X_{\alpha\beta\gamma}
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.


function [X2]= X(Base,Geo,dimention,op,op_m);

        X2=   deltaJ(Geo.h,Geo.h_diff,Base.Cmn1,Base.jmn_pm,Geo.m,Geo.n,dimention,0,op,op_m)...
             +deltaJ(Geo.h,Geo.h_diff,Base.Cmn1,Base.jmn_pm,Geo.m,Geo.n,dimention,1,op,op_m)...
             +deltaJ(Geo.h,Geo.h_diff,Base.Cmn1,Base.jmn_pm,Geo.m,Geo.n,dimention,-1,op,op_m);
end