function ProtConc = ExtractProtConcs(model, solution)

% Extract predicted protein concentrations from solution of ecGEM.
% Returns a table with Protein ID and protein concentration 
% (in mg protein / g CDW) predicted by FBA.
% Usage:
% ProteinConcentrations = ExtractProtConcs(model, solution)


Enzymes = model.ec.enzymes;
Concentration =[];
usage_str = 'usage_prot_';

for i=1:length(Enzymes)

    UsageRxnIdx = find(strcmp(model.rxns, strcat('usage_prot_', Enzymes{i})));
    ProteinID(i) = Enzymes(i);
    Concentration(i) = solution.x(UsageRxnIdx)*-1;
end

ProtConc = table(ProteinID', Concentration');
ProtConc.Properties.VariableNames = ["ProteinID", "Concentration"];
end
