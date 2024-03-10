function [cordova_fluxes, Model_fluxes_scaled, Cordova_flux_fraction, Model_flux_fraction] = Compare_to_Cordova(model, solution)

%Relevant_Rxns = [];
Model_fluxes_scaled = [];
Flux_Sum = [];
Model_flux_fraction = [];
cordova_fluxes = readtable('../../../Databases/Cordova_fluxes.xlsx');  % Imports Cordova's glucose 13C MFA fluxes
cordova_fluxes = cordova_fluxes.Var1';

Comparisonfluxes=["R00771","R00756","R01068","R01015","R01061","R01512",...
                  "R01518","R00658","R00200","R00209","R00344","R00351",...
                  "R01325","R00267","R08549","R00405","R02164","R01082",...
                  "R00342","R02035","R01528","R01529",...
                  "R01056","R01641","R08575","R01067","R01858","R01138",...
                  "R02734"];

for i=1:length(Comparisonfluxes)
    All_rxns_with_idx = model.rxns(find(contains(model.rxns, Comparisonfluxes(i))));
    
    if length(All_rxns_with_idx)>1

        for j=1:length(All_rxns_with_idx)

            if contains(All_rxns_with_idx(j), 'REV')
                solution.x(find(ismember(model.rxns, All_rxns_with_idx(j))))=-1*solution.x(find(ismember(model.rxns, All_rxns_with_idx(j))));
            end
        end
    end

    Flux_Sum = [Flux_Sum; abs(sum(solution.x(find(contains(model.rxns, All_rxns_with_idx)))))];

    
end

R00200=find(contains(Comparisonfluxes, "R00200"));
R01858=find(contains(Comparisonfluxes, "R01858"));
R01138=find(contains(Comparisonfluxes, "R01138"));
PTS = solution.x(find(ismember(model.rxns, 'PTS')));
R00206 = solution.x(find(ismember(model.rxns, 'R00206')));
Flux_Sum(R00200) = Flux_Sum(R00200) + Flux_Sum(R01858) + Flux_Sum(R01138) + PTS - R00206;


R00405=find(contains(Comparisonfluxes, "R00405"));
R02734=find(contains(Comparisonfluxes, "R02734"));
Flux_Sum(R00405) = Flux_Sum(R00405) + Flux_Sum(R02734);
Flux_Sum(R02734) = [];
Flux_Sum(R01138) = [];
Flux_Sum(R01858) = [];

% Calculate how close Modeled fluxes are to Cordova measured fluxes. To
% normalize each flux, calculate the fraction of measured flux that the
% modeled flux fulfills.

Model_flux_fraction = abs(Flux_Sum'./cordova_fluxes);
Cordova_flux_fraction = cordova_fluxes./cordova_fluxes;

Model_fluxes_scaled = abs(Flux_Sum./solution.x(find(contains(model.rxns, 'EX_Glucose'))))*100;
Model_fluxes_scaled = Model_fluxes_scaled';
