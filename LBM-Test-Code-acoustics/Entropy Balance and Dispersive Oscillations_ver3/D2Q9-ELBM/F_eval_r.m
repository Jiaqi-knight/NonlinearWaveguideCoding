%%-----------------------------------------------------------
%%    Function that evaluate the following function:
%%      F(alfa)=H(f_i)-H(f_i+alfa*(f_i-f_i^eq))
%%-----------------------------------------------------------
function [fn_loc]=F_eval_r(x,fn,df)
{
    fn_loc=0;df_loc=0;
    Delta_i=0;F_alfa=0;
    %%-------------------------------------------------------
    %%      H(f+alfa*Delta) - H(f) = 0
    %%-------------------------------------------------------
    app_log=0;
    sum_f=0;
    for k=1:9
    Delta_i=(equilib.link[k]-node_loc.link[k]);
    F_alfa=(node_loc.link[k]+x*Delta_i);
    app_log=(log(log_[k]*fabs(F_alfa)));
    %%---------------------------------------------------
    sum_f=sum_f+node_loc.link[k]*(log(log_[k]*node_loc.link[k]));
    fn_loc=fn_loc+(F_alfa*app_log);
    df_loc=df_loc+(Delta_i*(app_log+1)); %%  (2)

    end
    
    fn_loc=fn_loc-sum_f;

end


