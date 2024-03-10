
load('../../../Models/eciGEL627_posterior.mat');
sol = solveLP(ecModel);
Proteomics_DataStruct = loadProtData(4,'../../../Databases/Proteomics_fractions.tsv');
% Proteomics_DataStruct.abundances = Proteomics_DataStruct.abundances .* 510;
Proteomics_Data = struct2table(Proteomics_DataStruct);
Proteomics_Data = sortrows(Proteomics_Data);

Predicted_prots = ExtractProtConcs(ecModel, sol);
Predicted_prots = sortrows(Predicted_prots);

for i=1:height(Predicted_prots)
    if any(strcmp(Proteomics_Data.uniprotIDs, Predicted_prots.ProteinID(i)))
        idx = find(strcmp(Proteomics_Data.uniprotIDs, Predicted_prots.ProteinID(i)));
        Comparison.Protein(i) = Predicted_prots.ProteinID(i);
        Comparison.Predicted(i) = Predicted_prots.Concentration(i);
        Comparison.Measured(i) = Proteomics_Data.abundances(idx);
    else
        Comparison.Protein(i) = Predicted_prots.ProteinID(i);
        Comparison.Predicted(i) = Predicted_prots.Concentration(i);
        Comparison.Measured(i) = 0;        
    end
end

regression = fitlm(Comparison.Predicted, Comparison.Measured);
%disp(regression.Rsquared.Ordinary)
rsqrd = num2str(regression.Rsquared.Ordinary);
txt = "R2 = " + rsqrd;
figure
scatter(Comparison.Predicted, Comparison.Measured)
hold on
lsline;
text(5,5,0,txt,'FontSize',14)
hold off
xlabel('Predicted conc')
ylabel('Measured Conc')
% [~,PosIdx] = intersect(Proteomics_Data.uniprotIDs, Predicted_prots.ProteinID);
% Metabolic_Enzymes = Proteomics_Data(PosIdx,:);
