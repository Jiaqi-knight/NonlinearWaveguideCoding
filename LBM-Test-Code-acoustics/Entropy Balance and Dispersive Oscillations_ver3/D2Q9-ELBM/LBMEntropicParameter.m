% Function to find the constant entropy parameter for ELBM, in the case of
% LBKG this is simply 2
function [alpha, convergence] = LBMEntropicParameter(scheme,populations, ...
    popequilibriums,nx,ELBMEpsilon,Norm)





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

