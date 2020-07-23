%X Function for \Theta functions.
%m:circumferiential mode
%n:radial mode
%dimension:
%   2D           -  \X_{\alpha\beta}
%   3D           -  \X_{\alpha\beta\gamma}
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%   Copyright 2020, SJTU.


function [X2]= X(Geo_b,dimention,op,op_m);

        X2=   deltaJ_h(Geo_b.h,Geo_b.h_diff,Geo_b.m,Geo_b.n,dimention,0,op,op_m)...
             +deltaJ_h(Geo_b.h,Geo_b.h_diff,Geo_b.m,Geo_b.n,dimention,1,op,op_m)...
             +deltaJ_h(Geo_b.h,Geo_b.h_diff,Geo_b.m,Geo_b.n,dimention,-1,op,op_m);
%         X2_1=   deltaJ(Geo_b.h,Geo_b.h_diff,Geo_b.m,Geo_b.n,dimention,0,op,op_m)...
%              +deltaJ(Geo_b.h,Geo_b.h_diff,Geo_b.m,Geo_b.n,dimention,1,op,op_m)...
%              +deltaJ(Geo_b.h,Geo_b.h_diff,Geo_b.m,Geo_b.n,dimention,-1,op,op_m);
%          max(max( X2-X2_1))
end