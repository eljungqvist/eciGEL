%% Setup parameters

adapterLocation = fullfile(findGECKOroot, 'eciGEL', 'iGEL604_adapter.m');
ModelAdapterManager.setDefault(adapterLocation);

%% Create uniprotConversion.tsv

DB = loadDatabases('both');




%% Load model
model_conv = loadConventionalGEM();
rxnIdxs = findRxnIdxs(model_conv);

% Manually update GAM and NGAM here based on current Pirt parameter estimates
%model_conv = changeGAM(model_conv, 55, rxnIdxs);

%model_conv.lb(rxnIdxs.ngam) = 20;
%model_conv.ub(rxnIdxs.ngam) = 20;
%%
sol = solveLP(model_conv)
printFluxes(model_conv, sol.x, false)
%% Do some stuff to format model

% Fix compartment names
%model_conv.compNames = {'Cytosol'; 'Extracellular'};

% Import ec codes from mat file
% load(fullfile(modelAdapter.params.path,'models','iGEL628.mat'));
% model_conv.eccodes = erase(model.rxnECNumbers, 'ec-code/');  % Fetch ec-codes from iGEL604, for reasons not copied over to the xml version during sbml conversion...

% Clear compartment from metabolite names
%model_conv.metNames = erase(model_conv.metNames, ' [Cytosol]');
%model_conv.metNames = erase(model_conv.metNames, ' [Extracellular]');
%%
% Remove some odd gene associations
rnaRxnID = find(contains(model_conv.rxns, 'RNA'));
model_conv.grRules(rnaRxnID) = {''}; % For some reason mapped to a gene for glycogen synthesis now
model_conv.rxnGeneMat(rnaRxnID, :) = 0;

dnaRxnID = find(contains(model_conv.rxns, 'DNA'));
model_conv.grRules(dnaRxnID) = {''};  % Seems correctly mapped to DNA synthase, but kcat is odd with the mmols nucleotides needed
model_conv.rxnGeneMat(dnaRxnID, :) = 0;

NGAMRxnID = find(strcmp(model_conv.rxns, 'R00086'));
model_conv.grRules(NGAMRxnID) = {''}; % Mapped to a proposed phosphohydrolase. Creates an unnatural demand for this protein.
model_conv.rxnGeneMat(NGAMRxnID, :) = 0;

% Annotate names for exchange reactions
exRxnIdxs = find(contains(model_conv.rxns, 'EX_')); % Finds exchange reactions
for i = 1:size(exRxnIdxs)
    model_conv.rxnNames(exRxnIdxs(i)) = replace(model_conv.rxns(exRxnIdxs(i)), 'EX_', '');
    model_conv.rxnNames(exRxnIdxs(i)) = append(model_conv.rxnNames(exRxnIdxs(i)), ' exchange');
end

%% Annotate names for pseudoreactions
% Cell wall
tmpIdx = find(strcmp(model_conv.rxns, 'CellWall'));
model_conv.rxnNames(tmpIdx) = {'Cell wall pseudoreaction'};

% Cell membrane
tmpIdx = find(strcmp(model_conv.rxns, 'CellMembrane'));
model_conv.rxnNames(tmpIdx) = {'Cell membrane pseudoreaction'};

% Teichuronic acid
tmpIdx = find(strcmp(model_conv.rxns, 'Teichuronic_acid'));
model_conv.rxnNames(tmpIdx) = {'Teichuronic acid pseudoreaction'};

% Peptidoglycan
tmpIdx = find(strcmp(model_conv.rxns, 'Peptidoglycan'));
model_conv.rxnNames(tmpIdx) = {'Peptidoglycan pseudoreaction'};

% Protein synthesis
tmpIdx = find(strcmp(model_conv.rxns, 'aa_to_protein'));
model_conv.rxnNames(tmpIdx) = {'Protein biosynthesis pseudoreaction'};

% RNA
tmpIdx = find(strcmp(model_conv.rxns, 'RNA'));
model_conv.rxnNames(tmpIdx) = {'RNA pseudoreaction'};

% DNA
tmpIdx = find(strcmp(model_conv.rxns, 'DNA'));
model_conv.rxnNames(tmpIdx) = {'DNA pseudoreaction'};

% Glycogen 
tmpIdx = find(strcmp(model_conv.rxns, 'Glycogen'));
model_conv.rxnNames(tmpIdx) = {'Glycogen pseudoreaction'};

% Cofactors
tmpIdx = find(strcmp(model_conv.rxns, 'Cofactor_Pool'));
model_conv.rxnNames(tmpIdx) = {'Cofactor pool pseudoreaction'};

% Biomass
tmpIdx = find(strcmp(model_conv.rxns, 'Biomass'));
model_conv.rxnNames(tmpIdx) = {'Biomass pseudoreaction'};

% NGAM
tmpIdx = find(strcmp(model_conv.rxns, 'R00086'));
model_conv.rxnNames(tmpIdx) = {'NGAM pseudoreaction'};

%% Make EC model

% Set all EX_rxn lbs to 0. Otherwise the reactions get flipped when EcModel
% is made. Then make sure to turn on the uptake of required substrates.
EXRxnIdx = find(contains(model_conv.rxns, 'EX_'));
UptakeRxnsIdx = contains(model_conv.rxns, 'EX') & model_conv.lb < 0;
UptakeRxns = model_conv.rxns(UptakeRxnsIdx);
model_conv.lb(UptakeRxnsIdx) = 0;
[ecModel, noUniprot] = makeEcModel(model_conv, false, iGEL604_adapter);
ecUptakeRxns = ismember(ecModel.rxns, UptakeRxns);
ecModel.lb(ecUptakeRxns) = -1000;

% Gives warning on genes not found in uniprot.tsv, with geneIDs showing up
% blank. Seems to be an artefact of the .xml conversion of the model where
% model.genes have been extended with 15 blank genes. These are not found
% in the original .mat version of the model, and are not connected to any
% reaction in the .xml model. 
%%
sol = solveLP(ecModel)
printFluxes(ecModel, sol.x, false);
%% Get complex data - not available for LC300, might do manually sometime

%complexInfo = getComplexData('', modelAdapter);                         
%[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel, complexInfo);

%% Save model for now
saveEcModel(ecModel, modelAdapter, 'yml', 'ecModel_emil')

%% Add kcats to model
ecModel = getECfromGEM(ecModel);  % Gets the EC codes from the model
noEC = cellfun(@isempty, ecModel.ec.eccodes);  % Finds reactions lacking EC 

%%
ecModel = getECfromDatabase(ecModel, noEC);  % Find EC for missing reactions from database

%% funkar inte?
% Fix N/A in eccodes
naEC = find(contains(ecModel.ec.eccodes, 'N/A'));
ecModel.ec.eccodes(naEC) = {''};

% kcatList_fuzzy = fuzzyKcatMatching(ecModel); % Use this to try to match
% kcats to Brenda instead of using only DLKcat

%% Find SMILES
[ecModel, noSMILES] = findMetSmiles(ecModel);  

%% Prepare for DLkcat
writeDLKcatInput(ecModel);  % Outputs file containing input for DLKcat prediction

%% Machine-learning DLkcat

% Requires DLKcat to be runnable from Matlab. Can be done manually as well outside Matlab.
% Current repository already contains this output, so unless the model is
% changed, or it is desired, there is no need to run DLKcat again.

runDLKcat();  

%% Import and merge kcat data
kcatList_DLKcat = readDLKcatOutput(ecModel);  % Reads predicted kcats

%%
%kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat,
%kcatList_fuzzy);  % MErges Brenda and predicted kcats if fuzzy kcat
%matching is used above
%ecModel = selectKcatValue(ecModel, kcatList_merged);

ecModel = selectKcatValue(ecModel, kcatList_DLKcat); % Directly pick kcats from only DLKcat


%% Add kcats to isozymes missing kcat value (optional)
ecModel = getKcatAcrossIsozymes(ecModel);

% The current Bayesian script uses the ecModel from this stage as input, so
% if running the Bayesian fitting is desired, right-click and save
% "ecModel" in the workspace at this stage.

% After this point, the rest of the code is pre-Bayesian fitting
% implementation
save('models/eciGEL629',"ecModel");
