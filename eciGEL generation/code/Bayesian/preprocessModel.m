function updatedModel = preprocessModel(model, rxnIdxs, gam, ngam)
    
    % unconstrain all active reactions to take care of some odd constraints
    constrainedReverseReactions = find(model.lb < 0 & model.lb > -1000);
    model.lb(constrainedReverseReactions) = -1000;
    constrainedForwardReactions = find(model.ub > 0 & model.ub < 1000);
    model.lb(constrainedForwardReactions) = 1000;

    % Turn off xylose uptake
    model.lb(rxnIdxs.Xylose) = 0;  % Turns off xylose uptake
    model.ub(rxnIdxs.Xylose) = 0;  % Turns off xylose uptake


    % Turn on glucose uptake at batch rate
    model.lb(rxnIdxs.Glucose) = -1000;

    % Cap oxygen uptake to the highest chemostat value (higher than batch
    % for now)
     %model.lb(rxnIdxs.o2) = -51.61;

    % Change ngam
    model.lb(rxnIdxs.ngam) = ngam;  % ngam = 65

    % Change gam (originally 40)
    model = changeGAM(model, gam, rxnIdxs);  % gam = 45


    % Add glycogen exchange to model exopolysaccharides, and fix rate to
    % observed value
    % exopol_qp = 2.78;  % mmol/g, h
    % exopol_exchange_bound = exopol_qp*0.162;  % Normalize glycogen exchange to mols of glucose based on model reaction stoichiometry
    % addReaction(model, 'EX_Glycogen', 'reactionFormula', 'Glycogen ->', 'lowerBound', exopol_exchange_bound, 'upperBound', exopol_exchange_bound);
    % 

     updatedModel = model;
end