function [c13flux, model_flux, C13_flux_fraction, Model_flux_fraction] = Compare_to_Cordova(model, solution)
% Compares fluxes of central carbon metabolism between eciGEL prediction and
% 13C-MFA analysis performed by Cordova et al. 2015.


c13fluxes = readtable('../../../Databases/Cordova_fluxes_rxns.xlsx', 'ReadVariableNames', true);  % Imports Cordova's glucose 13C MFA fluxes


% 13C-MFA measures total flux between two metabolites, i.e. forward flux
% minus reverse flux. Some glycolysis and TCA reactions KEGG defines in the 
% reverse direction, so that forward in MFA is reverse in the model.
for i=1:length(c13fluxes.RxnID)
    rxn = extractBefore(c13fluxes.RxnID(i), 7);
    All_rxns_with_idx = model.rxns(find(contains(model.rxns, rxn)));

    if length(All_rxns_with_idx)>1 
        if contains(c13fluxes.RxnID(i), 'REV') % REV direction is the 'correct' direction 
          fwdrxns = All_rxns_with_idx(contains(All_rxns_with_idx, 'REV'));
          revrxns = All_rxns_with_idx(~contains(All_rxns_with_idx, 'REV'));
          model_flux(i) = sum(solution.x(contains(model.rxns, fwdrxns))) - sum(solution.x(contains(model.rxns, revrxns)));
        else % fwd direction is correct direction
          fwdrxns = All_rxns_with_idx(~contains(All_rxns_with_idx, 'REV'));
          revrxns = All_rxns_with_idx(contains(All_rxns_with_idx, 'REV'));
          model_flux(i) = sum(solution.x(contains(model.rxns, fwdrxns))) - sum(solution.x(contains(model.rxns, revrxns)));            
        end
    else
        model_flux(i) = solution.x(contains(model.rxns, c13fluxes.RxnID(i)));
    end
end
model_flux = model_flux';

% We're only measuring the total flux between two metabolites. If there are
% isozymes we want to compare their total to MFA. Total is calculated in
% loop above. Here remove isozyme reactions so they are not compared twice.
rowsToDelete = [];
for i = 1:height(c13fluxes)
    % Check if the entry contains 'EXP_2'
    if contains(c13fluxes.RxnID{i}, 'EXP_2')
        rowsToDelete = [rowsToDelete; i]; % Mark this row for deletion
    % Check if the entry contains 'EXP_1' and modify it if true
    elseif contains(c13fluxes.RxnID{i}, 'EXP_1')
        c13fluxes.RxnID{i} = c13fluxes.RxnID{i}(1:6); % Keep only the first 6 characters
    end
end
% Delete the marked rows in one operation
c13fluxes(rowsToDelete, :) = [];
model_flux(rowsToDelete) = [];

% PEP -> Pyr conversion performed by several reactions, introducing a lot
% of uncertainty. Remove from comparison. 
PEPPYR = find(strcmp(c13fluxes.RxnID, 'R00200_REV'));
c13fluxes(PEPPYR, :) = [];
model_flux(PEPPYR) = [];

% Cordova fluxes are normalized to a substrate uptake of 100. Do the same for
% model fluxes
PTSflux = solution.x(find(strcmp(model.rxns, 'PTS')));
model_flux = model_flux.*(100/PTSflux);

c13flux = c13fluxes.Flux;
% Calculate how close Modeled fluxes are to Cordova measured fluxes. To
% normalize each flux, calculate the fraction of measured flux that the
% modeled flux fulfills.
Model_flux_fraction = model_flux./c13fluxes.Flux;
C13_flux_fraction = c13fluxes.Flux./c13fluxes.Flux;

%Model_fluxes_scaled = abs(Flux_Sum./solution.x(find(contains(model.rxns, 'EX_Glucose'))))*100;
%Model_fluxes_scaled = Model_fluxes_scaled';
