%TFun Function for intial based  functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.



function [T2,T3]= TFun(Geo_b) %k+ right


%% different h can be similarly deal with
tic
disp('TFun..Begin...')

T2.ab=         Theta(Geo_b,2,'ab');
T2.ps_ab=      Theta(Geo_b,2,'ps_ab');
T2.pt_ab=      Theta(Geo_b,2,'pt_ab');
T2.ab_cos=     Theta(Geo_b,2,'ab_cos');
T2.a_ps_b=     Theta(Geo_b,2,'a_ps_b');
T2.pt_ab_cos=  Theta(Geo_b,2,'pt_ab_cos');

T3.abc=        Theta(Geo_b,3,'ab');
T3.pt_abc=     Theta(Geo_b,3,'pt_ab');
T3.abc_cos=    Theta(Geo_b,3,'ab_cos');
T3.abc_sin=    Theta(Geo_b,3,'ab_sin');
T3.pt_abc_cos= Theta(Geo_b,3,'pt_ab_cos');
disp('TFun..Finish!')
toc

end