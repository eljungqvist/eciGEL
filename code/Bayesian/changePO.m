function updatedModel = changePO(model, ComplexI_H_transported, ComplexII_H_transported, ComplexIII_H_transported, ATPsynthase_H_required, rxnIdxs)

    % Complex I. Changes how many Protons are pumped over the membrane in
    % Complex I of the electron transport chain. Default X = 4
    % NADH + Q + (X+1) H+ -> NAD+ + QH2 + X H[e]+

    model.S(rxnIdxs.H_ext_metidx, rxnIdxs.ComplexI) = ComplexI_H_transported;
    model.S(rxnIdxs.gamH(1), rxnIdxs.ComplexI) = -(ComplexI_H_transported+1);
   

    % Complex II. Changes how many Protons are pumped over the membrane in
    % Complex I of the electron transport chain. Default X = 1.
    % FADH2 + Q + X H => FAD+ + QH2 + X H+[e]
    model.S(rxnIdxs.H_ext_metidx, rxnIdxs.ComplexII) = ComplexII_H_transported;
    model.S(rxnIdxs.gamH(1), rxnIdxs.ComplexII) = -ComplexI_H_transported;
    
    % Complex III. Changes how many Protons are pumped over the membrane in
    % Complex IIII of the electron transport chain. Default X = 2
    % 0.5 O2 + QH2 + X H2 -> X H+[e] + Q + H2O

    model.S(rxnIdxs.H_ext_metidx, rxnIdxs.ComplexIII) = ComplexIII_H_transported;
    model.S(rxnIdxs.gamH(1), rxnIdxs.ComplexIII) = -ComplexIII_H_transported;
    
    % ATP Synthase. Changes how many external Protons are needed to generate one ATP. Default X = 2
    % ADP + Pi + X H+[e] -> ATP + (X-1) H+ + H2O

    model.S(rxnIdxs.H_ext_metidx, rxnIdxs.ATPsyn) = -ATPsynthase_H_required;
    model.S(rxnIdxs.gamH(1), rxnIdxs.ATPsyn) = (ATPsynthase_H_required - 1);

    updatedModel = model;

end
