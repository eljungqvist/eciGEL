function rxnIdxs = findRxnIdxs(model)
    
    rxnIdxs.Glucose = find(strcmp(model.rxns, 'EX_Glucose'));
    rxnIdxs.Xylose = find(strcmp(model.rxns, 'EX_Xylose'));
    rxnIdxs.o2 = find(strcmp(model.rxns, 'EX_Oxygen'));
    rxnIdxs.co2 = find(strcmp(model.rxns, 'EX_CO2'));
    rxnIdxs.Acetate = find(strcmp(model.rxns, 'EX_Acetate'));
    rxnIdxs.biomass = find(strcmp(model.rxns, 'Biomass'));
    rxnIdxs.Glycerol = find(strcmp(model.rxns, 'EX_C00116[e]'));  % Glycerol exchange reaction index
    rxnIdxs.Lactate = find(strcmp(model.rxns, 'EX_Lactate')); 
    rxnIdxs.ComplexI = find(strcmp(model.rxns, 'Complex_I')); 
    rxnIdxs.ComplexII = find(strcmp(model.rxns, 'Complex_II'));
    rxnIdxs.ComplexIII = find(strcmp(model.rxns, 'Complex_III')); 
    rxnIdxs.ATPsyn = find(strcmp(model.rxns, 'ATP_Synthase')); 
    rxnIdxs.H_ext_metidx = find(strcmp(model.mets, 'C00080_e'));
    
    % Need to fix the ones below later
    rxnIdxs.Galactose = find(strcmp(model.rxns, 'EX_Galactose'));
    rxnIdxs.Maltose = find(strcmp(model.rxns, 'EX_C00208[e]'));
    rxnIdxs.Sucrose = find(strcmp(model.rxns, 'EX_C00089[e]'));
    rxnIdxs.Cellobiose = find(strcmp(model.rxns, 'EX_Cellobiose'));
    rxnIdxs.Glycogen =      'temp';
    rxnIdxs.Mannitol = find(strcmp(model.rxns, 'EX_C00392[e]'));

    % GAM- and NGAM-related
    rxnIdxs.ngam = find(contains(model.rxns, 'R00086'));

    atp_metidx = find(strcmp(model.mets, 'C00002_c'));
    biomass_rxnidx = find(strcmp(model.rxns, 'Biomass'));
    rxnIdxs.gamATP = [atp_metidx, biomass_rxnidx];
    
    adp_metidx = find(strcmp(model.mets, 'C00008_c'));
    rxnIdxs.gamADP = [adp_metidx, biomass_rxnidx];

    pi_metidx = find(strcmp(model.mets, 'C00009_c'));
    rxnIdxs.gamPi = [pi_metidx, biomass_rxnidx];

    h2o_metidx = find(strcmp(model.mets, 'C00001_c'));
    rxnIdxs.gamH2O = [h2o_metidx, biomass_rxnidx];

    h_metidx = find(strcmp(model.mets, 'C00080_c'));
    rxnIdxs.gamH = [h_metidx, biomass_rxnidx];
