
[modelData, expData] = model_chemostat_growth(posteriorModel, expdata_chemostat, rxnIdxs, carbonNum); 

modelMu = modelData(2:29);  % Mu
modelHAc = modelData(32:59); % HAc
modelLac = modelData(62:89); % Lac
modelCO2 = modelData(92:119); % CO2
modelO2 = -modelData(122:149); % O2
modelqs = -modelData(152:179); % qS

% Unscale experimental data
expData(1:8) = expData(1:8)./carbonNum.biomass;
expData(9:16) = expData(9:16)./carbonNum.Acetate;
expData(17:24) = expData(17:24)./carbonNum.Lac;
expData(25:32) = expData(25:32)./carbonNum.CO2;
expData(41:48) = expData(41:48)./carbonNum.Glucose;
%%
figure
plot(modelMu, modelHAc)
hold on
scatter(expData(1:8), expData(9:16))
hold off
xlabel('Dilution rate')
ylabel('qHAc')

figure
plot(modelMu, modelO2)
hold on
scatter(expData(1:8), expData(33:40))
hold off
xlabel('Dilution rate')
ylabel('qO2')

figure
plot(modelMu, modelCO2)
hold on
scatter(expData(1:8), expData(25:32))
hold off
xlabel('Dilution rate')
ylabel('qCO2')

figure
plot(modelMu, modelqs)
hold on
scatter(expdata_chemostat.Mu/carbonNum.biomass, expdata_chemostat.q_s.*expdata_chemostat.C_Balance/100)
hold off
xlabel('Dilution rate')
ylabel('qS')

%%
function [modelData, expData] = model_chemostat_growth(model, expdata_chemostat, rxnIdxs, carbonNum)
    
    expQs = expdata_chemostat.q_s;

    expData = [];
    modelData = [];
    modMu = [];
    modHAc = [];
    modLac = [];
    modCO2 = [];
    modGlc = [];

    for i = 1:30  % Iteratively model chemostat growth at all experimentally tested q_s values, and report fluxes
        
        model.lb(rxnIdxs.Glucose) = -i;  % Set glucose uptake to experimental value
        sol = optimizeCbModel(model);  % Solve model for maximized growth
        
        if sol.stat > 0  % Check if model has a solution and store modeled fluxes
            modMu(i) = sol.f;
            modHAc(i) = sol.x(rxnIdxs.Acetate);
            modLac(i) = sol.x(rxnIdxs.Lactate);
            modCO2(i) = sol.x(rxnIdxs.co2);
            modO2(i) = sol.x(rxnIdxs.o2);
            modGlc(i) = sol.x(rxnIdxs.Glucose);
        else  % If no solution, set growth to 0
            modMu(i) = 0;
            modHAc(i) = NaN;
            modLac(i) = NaN;
            modCO2(i) = NaN;
            modO2(i) = NaN;
        end 

    end


    % Pack output data
    modelData = [modMu, modHAc, modLac, modCO2, modO2, modGlc]; 
    expData =[expdata_chemostat.Mu', expdata_chemostat.q_Ac', expdata_chemostat.q_Lac', expdata_chemostat.q_CO2', expdata_chemostat.q_O2', expdata_chemostat.q_s'];  % Not needed to do this many times...

end