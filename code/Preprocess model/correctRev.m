function model = correctRev(model)
% Correctes model.rev where bounds are set to irreversible but rev
% indicates reversible or vice versa.

    wrong_rev = [];
    for i=1:length(model.rxns)
        if model.lb(i)==0 && model.ub(i)>0 && model.rev(i)==1
            wrong_rev = [wrong_rev; i];
        elseif model.lb(i)<0 && model.ub(i)==0 && model.rev(i)==1
            wrong_rev = [wrong_rev; i];
        end
    end


    ShouldbeRev = [];
    for i=1:length(model.rxns)
        if model.lb(i)<0 && model.ub(i)>0 && model.rev(i) == 0
            ShouldbeRev = [ShouldbeRev; i];
        end
    end

     model.rev(wrong_rev) = 0;
     model.rev(ShouldbeRev) = 1;
end
