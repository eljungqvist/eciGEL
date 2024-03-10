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
%% Load Proteomics data and constrain model
Proteomics_Data = loadProtData(4);
ecModel = fillProtConcs(ecModel, Proteomics_Data);
ecModel = constrainProtConcs(ecModel);

%% Add "standard" kcat and MW to reactions missing gene information to constrain them
%[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);  % May need to redo this as maintenance and biomass now get kcat values associated

%% Set some constants for ecModel

ecModel_modified = manualModifications(ecModel, modelAdapter); 



%% Constrain model with kcats
%ecModel_constrained = applyKcatConstraints(ecModel_modified);

ecModel_constrained = applyKcatConstraints(ecModel_modified);
ecModel_constrained = setProtPoolSize(ecModel_constrained, [], modelAdapter.params.f);


%ecModel_constrained = setProtPoolSize(ecModel_constrained, [], 0.07);

%%  Test model vs experimental data
%changeCobraSolver('glpk', 'lp')
modsolution = testModelGrowth(ecModel_constrained)

%% use sensitivityTuning instead - not working?
%ecModel_constrained = applyKcatConstraints(ecModel);
%ecModel_constrained = setProtPoolSize(ecModel_constrained, [], modelAdapter.params.f);
cSourceIDx = find(contains(ecModel_constrained.rxns, modelAdapter.params.c_source));
ecModel_constrained.lb(cSourceIDx) = -1000;

[ecModel_tuned, tunedKcats] = sensitivityTuning(ecModel_constrained, 2.2, modelAdapter, 2);

%% Test tuned model
changeCobraSolver('gurobi', 'lp')
modsolution = testModelGrowth(ecModel_tuned);

%% Retune for acetate as well
model_HAc = changeRxnBounds(ecModel_tuned, modelAdapter.params.c_source, 0, 'b');
model_HAc = changeRxnBounds(model_HAc, 'EX_Acetate', -1000, 'l');
[ecModel_tuned2, tunedKcats2] = sensitivityTuning(model_HAc, 0.69, modelAdapter, 2);

%%
modsolution2 = testModelGrowth(ecModel_tuned);

%% --------------------------------------------------
%          BELOW IS UNUSED CODE FOR NOW
% ---------------------------------------------------

%% test modified model

%[rxns, percentdecrease] = testKcats(ecModel_modified, modelAdapter);
[rxns, percentdecrease] = testKcats(ecModel, modelAdapter);

%%
significantKcats = find(percentdecrease > 0.1);  % Rxns with more than 10% growth rate decrease from cutting kcat by 10x
significantReactions = rxns(significantKcats);  % Reactions with significance (not considering isozyme reactions as well...)
  %{'ATP_Synthase'}
  %  {'Complex_I'   }
  % {'Complex_III' }
   % {'R00239'      }
   % {'R01224_REV'  }
   % {'R01512_REV'  }
   % {'R01626'      }
   % {'R01777'      }
   % {'R03083'      }
   
%% Find most used enzymes
%ecModel_constrained.lb(3106) = -2.65
sol = optimizeCbModel(ecModel_constrained);
sol.f
usageData = enzymeUsage(ecModel_constrained, sol.x);
usageReport = reportEnzymeUsage(ecModel_constrained, usageData, 0.01) %> 90%


%% --------------------------------
% Function definitions below here

%%
function modsolution = testModelGrowth(ecModel)
% Takes an ec-Constrained iGEL-model and compares growth predictions with actual data

hac_id = find(strcmp(ecModel.rxns, 'EX_Acetate'));
o2_id = find(strcmp(ecModel.rxns, 'EX_Oxygen'));
ecModel.lb(hac_id) = 0;

disp('  ')
disp('Carbon source  Exp_Mu  Mod_mu  Exp_qs Mod_qs Exp_qHac  Mod_qHac  Mod_qO2')


% Glucose data

glucose_reaction = 'EX_Glucose';
glc_mu_expect = 2.2;  % Emil
glc_qs_expect = 35.7;  % Emil
glc_qHAc_expect = 10.1;  % Emil
glucose_id = find(strcmp(ecModel.rxns, 'EX_Glucose'));

ecModel.lb(glucose_id) = -1000;
glc_solution = optimizeCbModel(ecModel);
glc_mu_mod = glc_solution.f;
qsGlc_mod = -glc_solution.x(glucose_id);
glc_qHAc_mod = glc_solution.x(hac_id);
glc_qO2_mod = glc_solution.x(o2_id);

disp(sprintf('Glucose %f %f %f %f %f %f %f',  ...
    glc_mu_expect, glc_mu_mod, glc_qs_expect, ...
    qsGlc_mod, glc_qHAc_expect, glc_qHAc_mod, glc_qO2_mod));

modsolution = glc_solution;

% Xylose data
xyl_mu_expect = 1.25;  % Jeanett
xyl_qs_expect = 5.24/0.150;  % Cordova
xyl_qHAc_expect = 0;  % Cordova
xylose_id = find(strcmp(ecModel.rxns, 'EX_Xylose'));

ecModel.lb(glucose_id) = 0;
ecModel.ub(xylose_id) = 1000; %13.61;  % Switch to xylose uptake
xyl_solution = optimizeCbModel(ecModel);
xyl_mu_mod = xyl_solution.f;
qsXyl_mod = xyl_solution.x(xylose_id);
xyl_qHAc_mod = xyl_solution.x(hac_id);
xyl_qO2_mod = xyl_solution.x(o2_id);

disp(sprintf('Xylose %f %f %f %f %f %f %f',  ...
    xyl_mu_expect, xyl_mu_mod, xyl_qs_expect, ...
    qsXyl_mod, xyl_qHAc_expect, xyl_qHAc_mod, xyl_qO2_mod));

% Galactose 
galactose_id = find(strcmp(ecModel.rxns, 'EX_Galactose'));
ecModel.ub(xylose_id) = 0;
ecModel.lb(galactose_id) = -1000;
gal_mu_expect = 0.83;  % Cordova

galactose_sol = optimizeCbModel(ecModel);
gal_mu_mod = galactose_sol.f;
qsGal_mod = galactose_sol.x(galactose_id);
gal_qHAc_mod = galactose_sol.x(hac_id);
gal_qO2_mod = galactose_sol.x(o2_id);

disp(sprintf('Galactose %f %f N/A %f N/A %f %f',  ...
    gal_mu_expect, gal_mu_mod, ...
    qsGal_mod, gal_qHAc_mod, gal_qO2_mod));

% Acetate 
ecModel.lb(galactose_id) = 0;
ecModel.lb(hac_id) = -1000;
hac_mu_expect = 0.69;  % Jeanett

hac_sol = optimizeCbModel(ecModel);
hac_mu_mod = hac_sol.f;
qsHAc_mod = hac_sol.x(hac_id);
hac_qHAc_mod = hac_sol.x(hac_id);
hac_qO2_mod = hac_sol.x(o2_id);

disp(sprintf('Acetate %f %f N/A %f N/A %f %f',  ...
    hac_mu_expect, hac_mu_mod, ...
    qsHAc_mod, hac_qHAc_mod, hac_qO2_mod));

% Glycerol
ecModel.lb(hac_id) = 0;
glycerol_id = find(strcmp(ecModel.rxns, 'EX_C00116[e]'));
ecModel.lb(glycerol_id) = -1000;
glycerol_mu_expect = 1.37; % Jeanett

glycerol_sol = optimizeCbModel(ecModel);
glycerol_mu_mod = glycerol_sol.f;
qsGlycerol_mod = glycerol_sol.x(hac_id);
glycerol_qHAc_mod = glycerol_sol.x(hac_id);
glycerol_qO2_mod = glycerol_sol.x(o2_id);

disp(sprintf('Glycerol %f %f N/A %f N/A %f %f',  ...
    glycerol_mu_expect, glycerol_mu_mod, ...
    qsGlycerol_mod, glycerol_qHAc_mod, glycerol_qO2_mod));
% Mannose (not working now)
% mannose_id = find(strcmp(ecModel.rxns, 'EX_Mannose'));
% ecModel.lb(galactose_id) = 0;
% ecModellb(mannose_id) = -1000;
% man_mu_expect = 0.69; % 1.15;  % Cordova
% 
% mannose_sol = optimizeCbModel(ecModel);
% man_mu_mod = mannose_sol.f;
% qsMan_mod = mannose_sol.x(mannose_id);
% man_qHAc_mod = mannose_sol.x(hac_id);
% man_qO2_mod = mannose_sol.x(o2_id);
% 
% disp(sprintf('Mannose %f %f N/A %f N/A %f %f',  ...
%     man_mu_expect, man_mu_mod, ...
%     qsMan_mod, man_qHAc_mod, man_qO2_mod));

end

%% Manual modifications
function updatedEcModel = manualModifications(ecModel, modelAdapter)
activeUb = find(ecModel.ub > 0);
ecModel.ub(activeUb) = 1000; % Unconstrain all forward reactions

etoh_id = find(strcmp(ecModel.rxns, 'EX_Ethanol'));  % Need to remove this one from model!
ecModel.ub(etoh_id) = 0;

lac_id = find(strcmp(ecModel.rxns, 'EX_Lactate'));  % Should not be needed!
ecModel.ub(lac_id) = 0;

pyr_id = find(strcmp(ecModel.rxns, 'EX_Pyruvate'));
ecModel.ub(pyr_id) = 0;

r00207_id = find(strcmp(ecModel.rxns, 'R00207'));
ecModel.ub(r00207_id) = 0;  % Pyruvate + o2+pi = AcP + H2O2 + CO2


%Update MW and kcat manually
%pdh
pdh_genes = {'IB49_08605'; 'IB49_15580'; 'IB49_15575'; 'IB49_15570'; 'IB49_03685'};
for i = 1:length(pdh_genes)
   pdh_id = find(strcmp(ecModel.ec.genes, pdh_genes(i)));
   ecModel.ec.mw(pdh_id) = 216000;  % Assumes 1 active site/monomer in E1, E2, E3, and MWs from E coli
end

% atp synthase
%atps_id = find(strcmp(ecModel.ec.rxns, 'ATP_Synthase'));
%ecModel.ec.kcat(atps_id) = 16;  % E coli WT enzyme is 16.

% Respiration
%cI_id = find(strcmp(ecModel.ec.rxns, 'Complex_I'));
%cIII_id = find(strcmp(ecModel.ec.rxns, 'Complex_III'));
%ecModel.ec.kcat(cI_id) = 1;  % Arbitrarily low for now
%ecModel.ec.kcat(cIII_id) = 1;  % Arbitrarily low for now

% PTS Glucose
%ptsG_id = find(strcmp(ecModel.ec.rxns, 'PTS'));
%ecModel.ec.kcat(ptsG_id) = 5; %semi-arbitrary from brenda



%highKcats = find(ecModel.ec.kcat > 100);  %     finds kcats greater than 100
%ecModel.ec.kcat(highKcats) = 100;  % Sets max kcat to 100

%low_MWs = find(ecModel.ec.mw < 30);
%ecModel.ec.mw(low_MWs) = 30;  % Sets smallest allowed enzyme MW to 20

updatedEcModel = defineGrowthMedium(ecModel, modelAdapter);

end

function updatedModel = defineGrowthMedium(ecModel, modelAdapter)
    
    excludedVitamins = {'EX_PyridoxinHydrochloride', 'EX_Thiamine', 'EX_Nicotinicacid', 'EX_paminobenzoate', 'EX_Folate'}
    updatedModel = changeRxnBounds(ecModel, excludedVitamins, 0, 'b');
    updatedModel = ecModel;
    updatedModel.ub(73) = 0;  % set xylose uptake to 0
    c_source_idx = find(contains(updatedModel.rxns, modelAdapter.params.c_source));
    updatedModel.lb(c_source_idx) = -1000;   %-36;  % Set Emil's glucopse uptake

end

%% Test kcats having an influence
function [rxns, percent_decrease] = testKcats(ecModel, modelAdapter)
    c_source_idx = find(contains(ecModel.rxns, modelAdapter.params.c_source));
    ecModel.lb(c_source_idx) = -1000;  % Set unrestricted glucose uptake

    constrainedModel = applyKcatConstraints(ecModel);
    constrainedModel = setProtPoolSize(constrainedModel, [], 0.5);

    origsol = optimizeCbModel(constrainedModel);

    kcatsWithEffect = [];
    percent_decrease = [];

    for i= 1:length(ecModel.ec.kcat)
        origKcat = ecModel.ec.kcat(i);
        ecModel.ec.kcat(i) = origKcat/10;
        newconsmodel = applyKcatConstraints(ecModel);
        newconsmodel = setProtPoolSize(newconsmodel, [], 0.5);
        newsol = optimizeCbModel(newconsmodel);
        if newsol.f < origsol.f
            kcatsWithEffect = [kcatsWithEffect; ecModel.ec.rxns(i)];
            percent_decrease = [percent_decrease; (1-(newsol.f/origsol.f))];
        end
        ecModel.ec.kcat(i) = origKcat;

    end
    rxns = kcatsWithEffect;
    

end


