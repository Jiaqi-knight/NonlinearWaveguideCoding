% Function to find the constant entropy parameter for ELBM, in the case of
% LBKG this is simply 2
function [alpha, convergence] = LBMEntropicParameter(scheme,populations, ...
    popequilibriums,methodChoice,nx,ELBMEpsilon,Norm)
%Choice of method:
%0 - LBGK polynomial equilibria,
%1 - LBGK entropic equilibria,
%2 - ELBM Parabola iterations,
%3 - ELBM Bisection method
if (methodChoice == 0)
    alpha = 2*ones(1,nx);
    convergence = zeros(1,nx);
elseif (methodChoice == 1)
    alpha = 2*ones(1,nx);
    convergence = zeros(1,nx);
elseif (methodChoice == 2)
    convergence = zeros(1,nx);
    alpha = 2*ones(1,nx);
    Sf = LBMEvaluateS(scheme,populations,popequilibriums,zeros(1,nx));%eq-4.7
    SEquil = LBMEvaluateS(scheme,popequilibriums,popequilibriums,zeros(1,nx));%S(f)
    remaining = find(abs(SEquil - Sf) > 10^-15);
    popnorm(remaining) = LBMNorms(populations(:,remaining), ...
        popequilibriums(:,remaining),length(remaining),Norm);
    iteration = 1;
    while (isempty(remaining) == 0)
        alpha(remaining) = LBMInteriorParabola(scheme,populations(:,remaining),...
            popequilibriums(:,remaining),Sf(remaining),alpha(remaining));
        Spop = LBMEvaluateS(scheme,populations(:,remaining), ...
            popequilibriums(:,remaining),alpha(remaining));
        Sdiff = LBMEvaluateDiffS(scheme,populations(:,remaining), ...
            popequilibriums(:,remaining),alpha(remaining));
        dAlpha(remaining) = (Spop - Sf(remaining))./Sdiff;
        nowdone = find( popnorm(remaining).*abs(dAlpha(remaining)) ...
            <= ELBMEpsilon); %eq-4.10
        nowremaining = find( popnorm(remaining).* ...
            abs(dAlpha(remaining)) > ELBMEpsilon);
        done = remaining(nowdone);
        remaining = setdiff(remaining,done);
        SDone = LBMEvaluateS(scheme,populations(:,done), ...
            popequilibriums(:,done),alpha(done)) - Sf(done);
        negFind = find(SDone < 0);
        negativeEntropy = done(negFind);
        alpha(negativeEntropy) = alpha(negativeEntropy) ...
            - 2*dAlpha(negativeEntropy);
        SDone = LBMEvaluateS(scheme,populations(:,done), ...
            popequilibriums(:,done),alpha(done)) - Sf(done);
        if(SDone < 0)
            keyboard
        end
        convergence(done) = iteration;
        iteration = iteration + 1;
        if(~isreal(alpha))
            keyboard
        end
        if(iteration > 50)
%             keyboard
           alpha(remaining)=2;
        end
    end
elseif (methodChoice == 3)
    convergence = zeros(1,nx);
    alpha = 2*ones(1,nx);
    Sf = LBMEvaluateS(scheme,populations,popequilibriums,zeros(1,nx));
    SEquil = LBMEvaluateS(scheme,popequilibriums,popequilibriums,zeros(1,nx));
    remaining = find(abs(SEquil - Sf) > 10^-14);
    popnorm(remaining) = LBMNorms(populations(:,remaining), ...
        popequilibriums(:,remaining),length(remaining),Norm);
    alpha(remaining) = LBMInteriorParabola(scheme,populations(:,remaining), ...
        popequilibriums(:,remaining),Sf(remaining),alpha(remaining)) ;
    Spop = LBMEvaluateS(scheme,populations(:,remaining), ...
        popequilibriums(:,remaining),alpha(remaining)) - Sf(remaining);
    leftSide = find(Spop > 0);
    rightSide = find(Spop <= 0);
    left(remaining(leftSide)) = alpha(remaining(leftSide));
    right(remaining(rightSide)) = alpha(remaining(rightSide));
    Spop = LBMEvaluateS(scheme,populations(:,remaining), ...
        popequilibriums(:,remaining),alpha(remaining));
    Sdiff = LBMEvaluateDiffS(scheme,populations(:,remaining), ...
        popequilibriums(:,remaining),alpha(remaining));
    dAlpha(remaining) = (Spop - Sf(remaining))./Sdiff;
    right(remaining(leftSide)) = alpha(remaining(leftSide)) ...
        - 2*dAlpha(remaining(leftSide));
    left(remaining(rightSide)) = alpha(remaining(rightSide)) ...
        + 2*dAlpha(remaining(rightSide));
    if(right < left)
        keyboard
    end
    iteration = 1;
    mid = left + (right-left)/2;
    while (isempty(remaining) == 0)
        mid(remaining) = left(remaining) ...
            + (right(remaining)-left(remaining))/2;
        Sleft = LBMEvaluateS(scheme,populations(:,remaining), ...
            popequilibriums(:,remaining),left(remaining));
        Sright = LBMEvaluateS(scheme,populations(:,remaining), ...
            popequilibriums(:,remaining),right(remaining));
        Smid = LBMEvaluateS(scheme,populations(:,remaining), ...
            popequilibriums(:,remaining),mid(remaining));
        nowdone = find( popnorm(remaining).*(right(remaining) ...
            - left(remaining)) <= ELBMEpsilon);
        nowremaining = find( popnorm(remaining).*(right(remaining) ...
            - left(remaining)) > ELBMEpsilon);
        done = remaining(nowdone);
        remaining = setdiff(remaining,done);
        alpha(done) = left(done);
        convergence(done) = iteration;
        for itr = 1:length(remaining)
            if (Smid(nowremaining(itr)) > Sf(remaining(itr)))
                left(remaining(itr)) = mid(remaining(itr));
            else
                right(remaining(itr)) = mid(remaining(itr));
            end
        end
        iteration = iteration + 1;
        if(~isreal(alpha))
            keyboard
        end
        if(iteration > 100)
            keyboard
        end
    end
end

