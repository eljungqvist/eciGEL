function updatedModel = changeGAM(model, gam, rxnIdxs)

    model.S(rxnIdxs.gamATP(1), rxnIdxs.gamATP(2)) = -gam;
    model.S(rxnIdxs.gamH2O(1), rxnIdxs.gamH2O(2)) = -gam;
    
    model.S(rxnIdxs.gamH(1), rxnIdxs.gamH(2)) = gam;
    model.S(rxnIdxs.gamPi(1), rxnIdxs.gamPi(2)) = gam;
    model.S(rxnIdxs.gamADP(1), rxnIdxs.gamADP(2)) = gam;

    updatedModel = model;

end