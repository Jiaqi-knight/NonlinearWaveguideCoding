%Fun Function for intial based  functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.



function [X2,X3]= XFun(Geo_b) %k+ right


%% different h can be similarly deal with
tic
disp('XFun..Begin...')
X2.ab=         X(Geo_b,2,'ab','1');        %1/h
X2.ab_r=       X(Geo_b,2,'ab','r');        %1
X2.ab_r2 =     X(Geo_b,2,'ab','r2');       %h
X2.pr_ab_r=    X(Geo_b,2,'pr_ab','r');     %1/h
X2.pr_ab_r2=   X(Geo_b,2,'pr_ab','r2');    %1
X2.a_pr_b_r =  X(Geo_b,2,'a_pr_b','r');    %1
X2.a_ps_b_r =  X(Geo_b,2,'a_ps_b','r');    %h'/h
X2.ps_ab_r =   X(Geo_b,2,'ps_ab','r');     %h'/h
toc;tic
X3.abc=        X(Geo_b,3,'ab','1');        %1/h^2
X3.abc_r=      X(Geo_b,3,'ab','r');        %1/h
X3.abc_r2=     X(Geo_b,3,'ab','r2');       %1
X3.ps_abc_r=   X(Geo_b,3,'ps_ab','r');     %h'/h^2
X3.ps_abc=     X(Geo_b,3,'ps_ab','1');     %h'/h^3
X3.pr_abc_r=   X(Geo_b,3,'pr_ab','r');     %1/h^2
X3.pr_abc_r2=  X(Geo_b,3,'pr_ab','r2');    %1/h
X3.ps_abc_ps_r=X(Geo_b,3,'ps(ab)','r');    %h'/h^2
disp('XFun..Finish!')
toc

end
