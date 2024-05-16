% Script to run the DL-kcat Bayesian methodology
% Reads an unconstrained ec model of LC300 created by Gecko.
% Creates a prior kcat distribution based on DL-kcat values in the model
% And runs one generation comparing growth rate of kcat population with
% growth rates and fluxes from experimental data
% 2023-09-08 Created by Martin Gustavsson
% 2024-04-23 Current version by Emil Ljungqvist

%% Set up Cobra
initCobraToolbox
changeCobraSolver('gurobi', 'lp')  % Preferably use this

%%  Load EC model  - has kcats, protweights, etc
cd(fileparts(matlab.desktop.editor.getActiveFilename));  % Change folder to directory of this script

% Loads the unconstrained ecModel created by Gecko 3
load('../../../Models/eciGEL638.mat');
%%
rxnIdxs = findRxnIdxs(ecModel);

% Applies manual modifications to model prior to running the Bayesian fitting
ecModel = manualModifications(ecModel); 

 % Set GAM and NGAM
ecModel = preprocessModel(ecModel, rxnIdxs, 50, 50);

% Set PO-ratio
ecModel = changePO(ecModel,4,0,4,4,rxnIdxs);

ecModel = applyKcatConstraints(ecModel);

% Set size of Protein pool with 0.51 g/gDW protein, and 0.43 g/gprot
% metabolic enzymes
ecModel = setProtPoolSize(ecModel, 0.51, 0.43, 0.5); 

% Reads table with kapps and finds their indexes in model ec-structure
kapps = readtable('../../../Databases/kapp.tsv', 'FileType','text','Delimiter',{'\t'});
for i=1:height(kapps)
    Idxs(i) = find(ismember(ecModel.ec.rxns, kapps.Reactions(i)));
end
kapps.Idxs = Idxs';
ecModel.ec.source(kapps.Idxs) = {'kapp'};
tempmodel=ecModel;
%%  Load experimental data
% Imports all experimental data
expdata_all = readtable('../../../Databases/Cultivation_data/Cultivation_data.csv');  

% Conversion factors to normalize fluxes and residuals to Cmol in order to
% avoid weighting high-carbon compounds higher than low-carbon compounds
carbonNum.Glucose = 6;  % Cmol per mol glucose
carbonNum.Xylose = 5;  % Cmol per mol xylose
carbonNum.Glycerol = 3;  % Glycerol
carbonNum.Galactose = 6;  % Cmol per mol glucose
carbonNum.Maltose = 12;  % Cmol per mol xylose
carbonNum.Sucrose = 12;  % Glycerol
carbonNum.Cellobiose = 12;  % Cmol per mol xylose
carbonNum.Mannitol = 6;
carbonNum.Acetate = 2;  % Cmol per mol HAc
carbonNum.Lac = 3;  % Cmol per mol lactate
carbonNum.CO2 = 1; 
carbonNum.biomass = 1000/26.74; % Gives Cmmol/g, given than 1 Cmol biomass = 26.74g.

%% Scales fluxes to cmmol.
expdata_all.q_Ac = expdata_all.q_Ac*carbonNum.Acetate;  % Scale qHAc to cmmol
expdata_all.q_Lac = expdata_all.q_Lac*carbonNum.Lac;  % Scale qLac to cmmol
expdata_all.Mu = expdata_all.Mu*carbonNum.biomass;  % Scale growth to cmmol biomass


%% Concert qs values to cmmol for all substrates of interest to have comparable data for different csources
rowsToChange = contains(expdata_all.Csource, 'Glucose');
% Scale qs to account for missing carbon
expdata_all.q_s_cmmol(rowsToChange) = expdata_all.q_s(rowsToChange) * carbonNum.Glucose;
expdata_all.qs_scaled = expdata_all.q_s .* (expdata_all.C_Balance./100);
expdata_all.qs_cmmol_scaled = expdata_all.q_s_cmmol .* (expdata_all.C_Balance./100);
 
rowsToChange = contains(expdata_all.Csource, 'Xylose');
expdata_all.q_s_cmmol(rowsToChange) = expdata_all.q_s(rowsToChange) * carbonNum.Xylose;

rowsToChange = contains(expdata_all.Csource, 'Glycerol');
expdata_all.q_s_cmmol(rowsToChange) = expdata_all.q_s(rowsToChange) * carbonNum.Glycerol;

% Set scaled qs to NaN for all Csources without qs data
expdata_all.q_s_cmmol(expdata_all.q_s_cmmol==0)=NaN;

% Split dataset in two depending on cultivation mode
expdata_batch = expdata_all(find(strcmp(expdata_all.Cultmode,'Batch')),:); % Batch data
expdata_chemostat = expdata_all(find(strcmp(expdata_all.Cultmode,'Chemostat')),:);  % Chemostat data


%%  Create prior kcat distributions and run one cycle

kcats = ecModel.ec.kcat;  % Original kcats in the draft ecModel, derived from DLKcat
numKcats = 150;  % Number of random kcats in the population
sigma = 1;  % Standard deviation

% Bound max and min kcats (log10 values input)
lb = -3;
ub = 8;

% Turn all '0 kcats' to the average kcat instead
kcats(kcats==0) = mean(kcats(kcats~=0));  


% Create 150 random kcats for each reaction taken from a normal
% distribution centered on the original DLKcat estimate
random_kcats_all = arrayfun(@getRandomKcats, kcats, (ones(height(kcats), 1)*sigma), (ones(height(kcats), 1)*numKcats), 'UniformOutput', false);
random_kcats_all = cell2mat(random_kcats_all);


% Model growth for all 150 models generated, each containing a combination
% of random kcats for each reaction
[modeledGrowthRates, modeledGrowthAndFluxes, expBatchGrowthRates, expGrowthAndFlux, mod_CCF, C13_CCF] = modelGrowth(tempmodel, random_kcats_all, expdata_batch, expdata_chemostat, rxnIdxs, carbonNum, kapps);


modeledGrowthRates(isnan(modeledGrowthRates)) = 0;  % If model can't grow, set growth to 0

% Calculate the residual between model growth or flux and experimental data
% for each point
[residuals, residuals_batch, residuals_chemostat, residuals_ccf] = calculateRmse(expBatchGrowthRates, modeledGrowthRates, expGrowthAndFlux, modeledGrowthAndFluxes, mod_CCF, C13_CCF);


% Find the indices of the 100 kcat combinations with smallest residual
[~,lowestResiduals_idxs] = mink(residuals, 100);
disp(mink(residuals,1))

% Find the best 100 kcat combinations
best_100 = random_kcats_all(:,lowestResiduals_idxs);

% Store the kcats
kcatDists{1} = best_100;

i=2;
disp(i)
Storedresiduals=[];
%% Bayesian fitting
% Iteratively update posterior distribution to improve fit to model data

% Set this value to a suitable one for iterating. It is possible to first 
% test with e.g. 10, and then increase and rerun this section, as the 
% current generation number "i" is kept in the workspace.
% As a guideline, 100 generations with 150 kcats takes around 4 hours or so
% on my machine.

maxGenerations = 150; 

tic
for i = i:maxGenerations
    
    % Find mean and standard deviation for the 100 best kcats for reach
    % reaction from the previous iteration
    [mu, sigma] = updatePriorDistribution(best_100);
    
    % Generate 150 new kcats for each reaction using distributions based on
    % the best 100 from the last run
    random_kcats_all = arrayfun(@getRandomKcats, mu, sigma, repmat(150, height(mu), 1), 'UniformOutput', false);
    random_kcats_all = cell2mat(random_kcats_all);
   
    % Set kapps for reactions where they are applicable
    for j = 1:height(kapps)
        rowIndex = kapps.Idxs(j); % Get the index for the row in kappsDist
        kappValue = kapps.kapp(j); % Get the kapp value to set
        random_kcats_all(rowIndex, :) = kappValue; % Set all columns in that row to the kappValue
    end

    % Model growth using the new random kcats
    [modeledGrowthRates, modeledGrowthAndFluxes, expBatchGrowthRates, expGrowthAndFlux, mod_CCF, C13_CCF] = modelGrowth(tempmodel, random_kcats_all, expdata_batch, expdata_chemostat, rxnIdxs, carbonNum, kapps);

    modeledGrowthRates(isnan(modeledGrowthRates)) = 0;

    % Calculate residuals and find the best 100 kcat combinations
    [residuals, residuals_batch, residuals_chemostat, residuals_ccf] = calculateRmse(expBatchGrowthRates, modeledGrowthRates, expGrowthAndFlux, modeledGrowthAndFluxes, mod_CCF, C13_CCF);
    [~,lowestResiduals_idxs] = mink(residuals, 100);
    Storedresiduals = [Storedresiduals, mink(residuals,1)]; %stores the best RMSE from each generation in an array
    % Store the best kcat combinations and iterate
    best_100 = random_kcats_all(:,lowestResiduals_idxs);
    kcatDists{i} = best_100;
    disp(i)
end

% As a backup, save workspace to file here
save('Workspaces/date_workspace.mat')

disp(toc)

%% Get the mean kcat for each reaction based on the 100 best distributions after iterations above
% and use these to create a "posterior model". After this point the model
% is completed. Below are some model analyses and function definitions
posteriorMeanKcats = getPosteriorMeanKcats(best_100);
posteriorModel = manualModifications(ecModel); 
posteriorModel = preprocessModel(posteriorModel, rxnIdxs, 50, 50); 
posteriorModel = changePO(posteriorModel,4,0,4,4,rxnIdxs); 
posteriorModel = setProtPoolSize(posteriorModel, 0.51, 0.43, 0.5);  

posteriorModel.ec.kcat = posteriorMeanKcats;
posteriorModel.ec.kcat(kapps.Idxs) = kapps.kapp;

posteriorModel = applyKcatConstraints(posteriorModel);

sol = optimizeCbModel(posteriorModel)

flux2csv(posteriorModel, sol.x, 'test.txt')

%%  Model growth on different carbon sources after each generation of Bayesian fitting
for i = 1:width(kcatDists)
    kcats = kcatDists{i};
    for j = 1:height(kcats)
        meanKcats(j,1) = mean(kcats(j,:));

    end
    
    meanKcatCell{i,1} = meanKcats;
end

 posteriorModel.lb(rxnIdxs.Glucose) = 0;
 posteriorModel.lb(rxnIdxs.Acetate) = 0;
 posteriorModel.lb(rxnIdxs.Glycerol) = 0;
 posteriorModel.lb(rxnIdxs.Xylose) = 0;

for i = 1:height(meanKcatCell)
    % Apply kcat constraint from each generation
    posteriorModel.ec.kcat = meanKcatCell{i};
    posteriorModel = applyKcatConstraints(posteriorModel);
    
    posteriorModel.lb(rxnIdxs.Glucose) = 0;
    posteriorModel.lb(rxnIdxs.Acetate) = 0;
    posteriorModel.lb(rxnIdxs.Glycerol) = 0;

    posteriorModel.lb(rxnIdxs.Xylose) = 0;


    %glucose growth
    posteriorModel.lb(rxnIdxs.Glucose) = -1000;
    solglc = optimizeCbModel(posteriorModel);
    glucoseMod.Mu(i) = solglc.f;
    glucoseMod.qs(i) = solglc.x(rxnIdxs.Glucose);
    glucoseMod.qHAc(i) = solglc.x(rxnIdxs.Acetate);
    glucoseMod.qO2(i) = solglc.x(rxnIdxs.o2);
    glucoseMod.qCO2(i) = solglc.x(rxnIdxs.co2);
    glucoseMod.qLac(i) = solglc.x(rxnIdxs.Lactate);
    
  
    %acetate growth
    % posteriorModel.lb(rxnIdxs.Glucose) = 0;
    % posteriorModel.lb(rxnIdxs.Acetate) = -1000;
    % solhac = optimizeCbModel(posteriorModel);
    % hacMod.Mu(i) = solhac.f;
    % hacMod.qs(i) = solhac.x(rxnIdxs.Acetate);
    % hacMod.qHAc(i) = solhac.x(rxnIdxs.Acetate);
    % hacMod.qO2(i) = solhac.x(rxnIdxs.o2);
    % hacMod.qCO2(i) = solhac.x(rxnIdxs.co2);
    % hacMod.qLac(i) = solhac.x(rxnIdxs.Lactate);


    %glycerol growth
    posteriorModel.lb(rxnIdxs.Acetate) = 0;
    posteriorModel.lb(rxnIdxs.Glucose) = 0;
    posteriorModel.lb(rxnIdxs.Glycerol) = -1000;
    solgly = optimizeCbModel(posteriorModel);
    glycerolMod.Mu(i) = solgly.f;
    glycerolMod.qs(i) = solgly.x(rxnIdxs.Glycerol);
    glycerolMod.qHAc(i) = solgly.x(rxnIdxs.Acetate);
    glycerolMod.qO2(i) = solgly.x(rxnIdxs.o2);
    glycerolMod.qCO2(i) = solgly.x(rxnIdxs.co2);
    glycerolMod.qLac(i) = solgly.x(rxnIdxs.Lactate);

    %Xylose growth
    posteriorModel.lb(rxnIdxs.Glycerol) = 0;
    posteriorModel.lb(rxnIdxs.Xylose) = -1000;

    solxyl = optimizeCbModel(posteriorModel);
    xyloseMod.Mu(i) = solxyl.f;
    xyloseMod.qs(i) = solxyl.x(rxnIdxs.Xylose);
    xyloseMod.qHAc(i) = solxyl.x(rxnIdxs.Acetate);
    xyloseMod.qO2(i) = solxyl.x(rxnIdxs.o2);
    xyloseMod.qCO2(i) = solxyl.x(rxnIdxs.co2);
    xyloseMod.qLac(i) = solxyl.x(rxnIdxs.Lactate);

    posteriorModel.ub(rxnIdxs.Xylose) = 0;



end
%%  Plot model convergence to experimental data
% plot glucose
figure
subplot(2,3,1)
plot(glucoseMod.Mu)
hold on
plot(ones(i)*2.2)
hold off
title('Batch glucose growth')
xlabel('Number of iterations')
ylabel('Growth rate')

subplot(2,3,2)
plot(glucoseMod.qHAc)
hold on
plot(ones(i)*7.86)
hold off
title('Batch glucose qHAc')
xlabel('Number of iterations')
ylabel('qHAc (mmol/g,h)')

subplot(2,3,3)
plot(-glucoseMod.qO2)
hold on
plot(ones(i)*44.42)
hold off
title('Batch glucose qO_2')
xlabel('Number of iterations')
ylabel('qO_2 (mmol/g, h)')

subplot(2,3,4)
plot(glucoseMod.qCO2)
hold on
plot(ones(i)*45.70)
hold off
title('Batch glucose qCO_2')
xlabel('Number of iterations')
ylabel('qCO_2 (mmol/g, h')

subplot(2,3,5)
plot(glucoseMod.qs)
hold on
plot(ones(i)*-30.06)
hold off
title('Batch glucose q_s')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')

subplot(2,3,6)
plot(glucoseMod.qLac)
hold on
plot(ones(i)*1.09)
hold off
title('Batch glucose qLac')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')

% Plot xylose
figure
subplot(2,3,1)
plot(xyloseMod.Mu)
hold on
plot(ones(i)*1.52)
hold off
title('Batch xylose growth')
xlabel('Number of iterations')
ylabel('Growth rate')

subplot(2,3,2)
plot(xyloseMod.qHAc)
%hold on
%plot(ones(i)*10.1)
%hold off
title('Batch xylose qHAc')
xlabel('Number of iterations')
ylabel('qHAc (mmol/g,h)')

subplot(2,3,3)
plot(-xyloseMod.qO2)
%hold on
%plot(ones(i)*1.95)
%hold off
title('Batch xylose qO_2')
xlabel('Number of iterations')
ylabel('qO_2 (mmol/g, h)')

subplot(2,3,4)
plot(xyloseMod.qCO2)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch xylose qCO_2')
xlabel('Number of iterations')
ylabel('qCO_2 (mmol/g, h')

subplot(2,3,5)
plot(xyloseMod.qs)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch xylose q_s')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')

subplot(2,3,6)
plot(xyloseMod.qLac)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch xylose qLac')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')


% Plot glycerol
figure
subplot(2,3,1)
plot(glycerolMod.Mu)
hold on
plot(ones(i)*1.95)
hold off
title('Batch glycerol growth')
xlabel('Number of iterations')
ylabel('Growth rate')

subplot(2,3,2)
plot(glycerolMod.qHAc)
%hold on
%plot(ones(i)*10.1)
%hold off
title('Batch glycerol qHAc')
xlabel('Number of iterations')
ylabel('qHAc (mmol/g,h)')

subplot(2,3,3)
plot(-glycerolMod.qO2)
%hold on
%plot(ones(i)*1.95)
%hold off
title('Batch glycerol qO_2')
xlabel('Number of iterations')
ylabel('qO_2 (mmol/g, h)')

subplot(2,3,4)
plot(glycerolMod.qCO2)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch glycerol qCO_2')
xlabel('Number of iterations')
ylabel('qCO_2 (mmol/g, h')

subplot(2,3,5)
plot(glycerolMod.qs)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch glycerol q_s')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')

subplot(2,3,6)
plot(glycerolMod.qLac)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch glycerol qLac')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')

% Plot HAc
figure
subplot(2,3,1)
plot(hacMod.Mu)
hold on
plot(ones(i)*0.48)
hold off
title('Batch HAc growth')
xlabel('Number of iterations')
ylabel('Growth rate')

subplot(2,3,2)
plot(hacMod.qHAc)
%hold on
%plot(ones(i)*10.1)
%hold off
title('Batch HAc qHAc')
xlabel('Number of iterations')
ylabel('qHAc (mmol/g,h)')

subplot(2,3,3)
plot(-hacMod.qO2)
%hold on
%plot(ones(i)*1.95)
%hold off
title('Batch HAc qO_2')
xlabel('Number of iterations')
ylabel('qO_2 (mmol/g, h)')

subplot(2,3,4)
plot(hacMod.qCO2)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch HAc qCO_2')
xlabel('Number of iterations')
ylabel('qCO_2 (mmol/g, h')

subplot(2,3,5)
plot(hacMod.qs)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch HAc q_s')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')

subplot(2,3,6)
plot(hacMod.qLac)
%hold on
%plot(ones(i)*1.35)
%hold off
title('Batch HAc qLac')
xlabel('Number of iterations')
ylabel('q_s (mmol/g, h')


%% Plot residuals over generations
figure
scatter(1:149, Storedresiduals)
xlabel('Generations')
ylabel('R.M.S.E')

%% Compare fluxes to Cordova

sol=optimizeCbModel(posteriorModel);
[rxns, cordova_fluxes, Model_fluxes_scaled, Cordova_flux_fraction, Model_flux_fraction] = Compare_to_Cordova(posteriorModel, sol);
regression = fitlm(cordova_fluxes, Model_fluxes_scaled);
%disp(regression.Rsquared.Ordinary)
rsqrd = num2str(regression.Rsquared.Ordinary);
txt = "R2 = " + rsqrd;
figure
scatter(cordova_fluxes, Model_fluxes_scaled)
hold on
lsline;
text(10,140,0,txt,'FontSize',14)
hold off
xlabel('13C MFA fluxes')
ylabel('eciGEL fluxes')

%% ------------------------------------
% Function definitions below this point
% -------------------------------------


% Manual modifications, scripts taken from iGEL604_ec_script 2023-09-08
function updatedEcModel = manualModifications(ecModel)
activeUb = find(ecModel.ub > 0);
ecModel.ub(activeUb) = 1000; % Unconstrain all forward reactions

% Allow production of all amino acids
aaExchanges = {'EX_Alanine', 'EX_Arginine', 'EX_Asparagine', 'EX_Aspartate', 'EX_Cysteine' ,... 
    'EX_Glutamate', 'EX_Glutamine', 'EX_Glycine', 'EX_Histidine', 'EX_Isoleucine','EX_Leucine',...
    'EX_Lysine', 'EX_Methionine','EX_Phenylalanine', 'EX_Proline' , 'EX_Serine', 'EX_Threonine',...
    'EX_Tryptophan', 'EX_Tyrosine', 'EX_Valine'}; % Defines names of amino acid exchange reactions

ecModel.ub(contains(ecModel.rxns, aaExchanges)) = 1000;

updatedEcModel = defineGrowthMedium(ecModel);

end

function updatedModel = defineGrowthMedium(ecModel)
    
    excludedVitamins = {'EX_PyridoxinHydrochloride', 'EX_Thiamine', 'EX_Nicotinicacid', 'EX_paminobenzoate', 'EX_Folate'};
    %updatedModel = changeRxnBounds(ecModel, excludedVitamins, 0, 'b');
    updatedModel = setParam(ecModel, 'lb', excludedVitamins, 0);
    %updatedModel = ecModel;
    updatedModel.lb(find(ismember(updatedModel.rxns, 'EX_Xylose'))) = 0;  % set xylose uptake to 0
    c_source_idx = find(contains(updatedModel.rxns, 'EX_Glucose'));
    updatedModel.lb(c_source_idx) = -1000;   %-36;  % Set Emil's glucose uptake

end

%% Sampling kcats
function randomKcats = getRandomKcats(mu, sigma, numberOfKcats)
    mu = log10(mu);
    lb = -3; % Bound kcats within reasonable values
    ub = 8; % Bound kcats within reasonable values
    
    priorDist = makedist('normal', 'mu', mu, 'sigma', sigma);
    randomKcats = random(priorDist, 1, numberOfKcats);
    randomKcats(randomKcats < lb) = lb;
    randomKcats(randomKcats > ub) = ub;
    randomKcats = 10.^(randomKcats);
end

function[batchMod, growthFluxMod, batchExp, growthFluxExp, mod_CCF_fraction, C13_CCF_fraction] = modelGrowth(model, random_kcats_all, expdata_batch,  expdata_chemostat, rxnIdxs, carbonNum, kapps)
    batchMod = {};  % Stores growth rates for tested carbon sources from the generated models
    batchExp = {};
    batchCCF = {};
    growthFluxMod = {};  % Stores growth rates for tested carbon sources from the generated models
    growthFluxExp = {};
    models = {};

    numKcats = width(random_kcats_all);


    for i = 1:numKcats  % Generate and test 150 random kcat combinations
        model.ec.kcat = random_kcats_all(:,i);
        models{i} = applyKcatConstraints(model);
    end

    for i = 1:length(models)
        % Batch growth
        [modeledGrowth, expGrowth] = model_batch_growth(models{i}, expdata_batch, rxnIdxs, carbonNum);
        batchMod{i} = modeledGrowth;
        batchExp{i} = expGrowth;  

        % Central carbon fluxes
        [~, ~, C13_fraction, modeled_fraction] = get_central_carbon_fluxes(models{i}, rxnIdxs);
        C13_CCF_fraction{i} = C13_fraction;
        mod_CCF_fraction{i} = modeled_fraction;
    
        % Chemostat growth
        [modeledGrowth, expGrowth] = model_chemostat_growth(models{i}, expdata_chemostat, rxnIdxs, carbonNum);
        
        growthFluxMod{i} = modeledGrowth;
        growthFluxExp{i} = expGrowth;

    end

    batchMod = cell2mat(batchMod');
    growthFluxMod = cell2mat(growthFluxMod');

    mod_CCF_fraction = cell2mat(mod_CCF_fraction');

    batchExp = cell2mat(batchExp');
    growthFluxExp = cell2mat(growthFluxExp');
    C13_CCF_fraction = cell2mat(C13_CCF_fraction');

end

function [C13_fluxes, modeled_fluxes, C13_fraction, modeled_fraction] = get_central_carbon_fluxes(model, rxnIdxs)

growthRates = [];
expRates = [];
Central_carbon_fluxes = [];

csource = 'Glucose';
csource_rxn = rxnIdxs.(csource); 
model.lb(csource_rxn) = -1000;
sol = optimizeCbModel(model);

    if sol.stat > 0
        [~, C13_fluxes, modeled_fluxes, C13_fraction, modeled_fraction] = Compare_to_Cordova(model, sol);
        C13_fraction = C13_fraction';
        modeled_fraction = modeled_fraction';
    else
        modeled_fluxes(1:28) = NaN;
        C13_fluxes(1:28) = NaN;
        C13_fraction(1:28) = NaN;
        modeled_fraction(1:28) = NaN;
    end

model.lb(csource_rxn) = 0;
end

function [growthRates, expRates] = model_batch_growth(model, expdata_batch, rxnIdxs, carbonNum)
%Loop through experiental batch data and create matching models  


growthRates = [];
expRates = [];

for i = 1:height(expdata_batch)
    
    csource = string(expdata_batch.Csource(i)); % Gets name of carbon source
    csource_rxn = rxnIdxs.(csource); % Finds index of csource in model
    

    if ~strcmp(csource_rxn, 'temp')  % Limit to only the 4 c-sources that have complete pathways in iGEL604

        model.lb(csource_rxn) = -1000;
        sol = optimizeCbModel(model);
        if sol.stat > 0  % Check if model has a solution and store modeled fluxes
            modMu(i) = sol.f*carbonNum.biomass;
            modHAc(i) = sol.x(rxnIdxs.Acetate)*carbonNum.Acetate;
            modLac(i) = sol.x(rxnIdxs.Lactate)*carbonNum.Lac;
            modCO2(i) = sol.x(rxnIdxs.co2)*carbonNum.CO2;
            modO2(i) = abs(sol.x(rxnIdxs.o2));
            modqs(i) = abs(sol.x(csource_rxn)*carbonNum.(csource));
        else  % If no solution, set growth to 0
            modMu(i) = 0;
            modHAc(i) = NaN;
            modLac(i) = NaN;
            modCO2(i) = NaN;
            modO2(i) = NaN;
            modqs(i) = NaN;
        end 
        model.lb(csource_rxn) = 0;
    end
        
    growthRates = [modMu, modHAc, modLac, modCO2, modO2, modqs]; 
    expRates = [expdata_batch.Mu', expdata_batch.q_Ac', expdata_batch.q_Lac', expdata_batch.q_CO2', expdata_batch.q_O2', expdata_batch.q_s_cmmol'];
    
end
end


function [modelData, expData] = model_chemostat_growth(model, expdata_with_flux, rxnIdxs, carbonNum)
    
    expQs = expdata_with_flux.q_s;

    expData = [];
    modelData = [];
    modMu = [];
    modHAc = [];
    modLac = [];
    modCO2 = [];
    modO2 = [];
    for i = 1:height(expdata_with_flux)
        
        csource = string(expdata_with_flux.Csource(i)); % Gets name of carbon source
        csource_rxn = rxnIdxs.(csource); % Finds index of csource in model
        
        model.lb(csource_rxn) = -expdata_with_flux.q_s(i);
        sol = optimizeCbModel(model);
        
        if sol.stat > 0
            modMu(i) = sol.f*carbonNum.biomass;
            modHAc(i) = sol.x(rxnIdxs.Acetate)*carbonNum.Acetate;
            modLac(i) = sol.x(rxnIdxs.Lactate)*carbonNum.Lac;
            modCO2(i) = sol.x(rxnIdxs.co2)*carbonNum.CO2;
            modO2(i) = abs(sol.x(rxnIdxs.o2));
            modGlc(i) = abs(sol.x(rxnIdxs.Glucose)*carbonNum.Glucose);
        else
            modMu(i) = 0;
            modHAc(i) = NaN;
            modLac(i) = NaN;
            modCO2(i) = NaN;
            modO2(i) = NaN;
            modGlc(i) = NaN;
        end

    end
    
    % Pack output data
    modelData = [modMu, modHAc, modLac, modCO2, modO2, modGlc]; 
    expData =[expdata_with_flux.Mu', expdata_with_flux.q_Ac', expdata_with_flux.q_Lac', expdata_with_flux.q_CO2', expdata_with_flux.q_O2', expdata_with_flux.q_s_cmmol'];
end



function [residuals, residuals_batch, residuals_chemostat, residuals_ccf] = calculateRmse(experimentalData_batch, modeledData_batch, experimentalData_chemostat, modeledData_chemostat, modCCF, C13CCF)
    % Calculates rmse between model data and experimental data
    % Used to rank best kcat distributions
    % All data including production rates other than mu should be
    % normalized to Cmmol in the input

    residuals_batch = [];
    residuals_chemostat = [];
    residuals_ccf = [];
    

    % Find residuals to the batch data. 
    for i = 1:height(experimentalData_batch)
        residual = rmse(experimentalData_batch(i,:),modeledData_batch(i,:), 'omitnan');
        residuals_batch = [residuals_batch; residual];
    end

    %Find residuals between modelled central carbon fluxes and C13 data
    %from cordova. Only for glucose as C-source
    for i=1:height(modCCF)
        residual = rmse(C13CCF(i,:), modCCF(i,:), 'omitnan');
        residuals_ccf = [residuals_ccf; residual];
    end

    % Find residuals to the chemostat data. Takes mu, qHac, qLac, qCO2 and
    % qO2 into account. Residual is the sum of the residuals of each
    % factor.
    for i = 1:height(experimentalData_chemostat)
        residual = rmse(experimentalData_chemostat(i,:), modeledData_chemostat(i,:), 'omitnan');
        residuals_chemostat = [residuals_chemostat; residual];
    end
    
  
   residuals = sqrt((residuals_batch.^2 + residuals_chemostat.^2 + residuals_ccf.^2)/3);  % Gets the average residuals from the three data types
   
end

function randomKcats = updateRandomKcats(kcats, numberToGenerate)
    
    priorDist = fitdist(log10(kcats),'Normal');
    
    lb = -3; % Bound kcats within reasonable values
    ub = 8; % Bound kcats within reasonable values
    
    randomKcats = random(priorDist, 1, numberToGenerate);
    randomKcats(randomKcats < lb) = lb;
    randomKcats(randomKcats > ub) = ub;
    randomKcats = (10.^(randomKcats));
     

end


function [mu, sigma] = updatePriorDistribution(best_100)
    % Updates prior distribution by fitting normal distribution to the 100
    % best kcats for each reaction
    
    mu = [];
    sigma = [];
    for i = 1:height(best_100)
        logtransformed_kcats = log10(best_100(i,:));
        pd = fitdist(logtransformed_kcats', 'Normal');
        mu = [mu; 10^(pd.mu)];
        sigma = [sigma; pd.sigma];
    end
end

function kcats = getPosteriorMeanKcats(best_100)
    for i = 1:height(best_100)
        kcats(i) = mean(best_100(i,:));

    end
    kcats = kcats';
end
