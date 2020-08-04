function [alpha]=NewtonMod(Finf_,Fsup_,index,)

x1=1.6;x2=2.2;
%%-------------------------------------------------------
%%ATTENTION: the values Finf_ and Fsup_ (that are global variables)
%       have been evaluated in function alpha_().
if (Finf_==0.0)
x=x1;
end
if (Fsup_==0.0)
    x=x2;
end
%%-------------------------------------------------------
if(Fsup_> 0.0)
    %%Orient the search so that f(xl) < 0.
    xl=x1;
    xh=x2;
else
    xh=x1;
    xl=x2;
end
%%-------------------------------------------------------
rts=2; %%this is the starting point for the root finder
%%the optimization strategies have to be applied here
(*FD_eval)(rts,&f,&df);
if (fabs(f)<=TOL) 
    x= rts;
end
    %%
    dxold=fabs(x2-x1);   %%(1)
    %%   dxold=f/df;         %%(2)
    dx=dxold;            %%and the last step.
    %%-------------------------------------------------------
    for j=1:MAXIT
        %%---------------------------------------------------
        %%Loop over allowed iterations.
        if((((rts-xh)*df-f)*((rts-xl)*df-f)>0.0)||(fabs(2.0*f)>fabs(dxold*df)))
            %%Bisect if Newton out of range,or not decreasing fast enough.
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts)
                x= rts; %%Change in root is negligible.
            end
        else
            %%Newton step acceptable. Take it.
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts =rts- dx;
            if (temp == rts)
                x= rts;
            end
        end
        %%---------------------------------------------------
        if (fabs(dx)<TOLL)
            x= rts;%%Convergence criterion.
        end
        %%---------------------------------------------------
        (*FD_eval)(rts,&f,&df);
        %%---------------------------------------------------
        %%The one new function evaluation per iteration.
        if (f < 0.0) %%Maintain the bracket on the root.
            xl=rts;
        else
            xh=rts;
        end
        printf("Error. We exceed the maximum number of iterations!\n");
        x= 0.0;
    end
end
