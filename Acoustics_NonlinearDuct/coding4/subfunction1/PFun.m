%PsiFun Function for intial based  functions.
%   ref:https://github.com/Jiaqi-knight/NonlinearWaveguideCoding
%   Email:Jiaqi_Wang@sjtu.edu.cn
%Copyright 2020, SJTU.



function [Psi]= PFun(Geo_b) 


%% different h can be similarly deal with
[X2,X3]= XFun(Geo_b)
[T2,T3]= TFun(Geo_b)

tic
disp('PFun..Begin...')
Psi.ab_r=           X2.ab_r         .*T2.ab;
Psi.a_pr_b_r =      X2.a_pr_b_r     .*T2.ab;
Psi.ab_r2_cos =     X2.ab_r2        .*T2.ab_cos;
Psi.ps_ab_r_1 =     X2.ps_ab_r      .*T2.ab;
Psi.ps_ab_r_2 =     X2.ab_r         .*T2.ps_ab;
Psi.a_ps_b_r_1 =    X2.a_ps_b_r     .*T2.ab;
Psi.a_ps_b_r_2 =    X2.ab_r         .*T2.a_ps_b;
Psi.pr_ab_r=        X2.pr_ab_r      .*T2.ab;
Psi.pr_ab_r2_cos=   X2.pr_ab_r2     .*T2.ab_cos;
Psi.pt_ab=          X2.ab           .*T2.pt_ab;
Psi.pt_ab_r_cos=    X2.ab_r         .*T2.pt_ab_cos;
Psi.ab=             X2.ab           .*T2.ab;
Psi.ab_r_cos=       X2.ab_r         .*T2.ab_cos;
Psi.pt_ab_cos=      X2.ab           .*T2.pt_ab_cos;
Psi.ab_s1= specialFun(Geo_b.theta_0,Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'hh`^2/[1-\kappa*h*cos\psi]');
Psi.ab_s2= specialFun(Geo_b.theta_0,Geo_b.h_diff,Geo_b.h,Geo_b.kappa,Geo_b.m,Geo_b.n,'h(1-\kappa*h*cos\psi)');


Psi.abc=                X3.abc          .*T3.abc;
Psi.abc_r=              X3.abc_r        .*T3.abc;
Psi.abc_r2_cos=         X3.abc_r2       .*T3.abc_cos;
Psi.abc_r_cos=          X3.abc_r        .*T3.abc_cos;
Psi.abc_r_sin=          X3.abc_r        .*T3.abc_sin;
Psi.pr_abc_r=           X3.pr_abc_r     .*T3.abc;
Psi.pr_abc_r2_cos=      X3.pr_abc_r2    .*T3.abc_cos;
Psi.pt_abc=             X3.abc          .*T3.pt_abc;
Psi.pt_abc_r_cos=       X3.abc_r        .*T3.pt_abc_cos;
Psi.abc=                X3.abc          .*T3.abc;
Psi.abc_r_cos=          X3.abc_r        .*T3.abc_cos;
Psi.abc_r_sin=          X3.abc_r        .*T3.abc_sin;
Psi.ps_abc_ps_r=        X3.ps_abc_ps_r  .*T3.abc;  %(Jiaqi-153)
Psi.ps_abc_r=           X3.ps_abc_r     .*T3.abc;
Psi.ps_abc=             X3.ps_abc       .*T3.abc;
disp('PFun..Finish!')
toc

end