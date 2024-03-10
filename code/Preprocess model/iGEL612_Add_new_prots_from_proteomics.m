

%load('../../../Models/iGEL612.mat');
model=ravenCobraWrapper(model);

%remove weird gene
model=removeGenes(model,'IB49_1815');
%% 
% Proteins found in proteomics that were not assigned in model. Added to
% reactions according to their annotation.
genes.genes = {'IB49_03915', 'IB49_05305', 'IB49_05310', 'IB49_04735', 'IB49_03065', 'IB49_05515',...
               'IB49_09375', 'IB49_16340', 'IB49_00730', 'IB49_17385', 'IB49_12915', 'IB49_11305',...
               'IB49_12715', 'IB49_12970', 'IB49_14520', 'IB49_09935', 'IB49_05065', 'IB49_13625', 'IB49_06865'};
model = addGenesRaven(model, genes);

% Reactions that already exist

%'AKU25739.1'    'IB49_03915'    'Glucokinase'      R00299
model = changeGrRules(model, 'R00299', 'IB49_03915', false);

%'AKU25985.1'	'IB49_05305'	'electron transfer flavoprotein subunit alpha'  Complex_I
%'AKU25986.1'	'IB49_05310'	'electron transfer flavoprotein subunit beta'   Complex_I
%'AKU25881.1'	'IB49_04735'	'NADH dehydrogenase'                            Complex_I
model = changeGrRules(model, 'Complex_I', ' and IB49_05305 and IB49_05310 and IB49_04735 and IB49_06865', false);

%'AKU25606.1'	'IB49_03065'	'phosphoglycerate mutase'                       R01518
model = changeGrRules(model, 'R01518', 'IB49_03065', false);

%'AKU26624.1'	'IB49_09375'	'ATP synthase subunit I'    ATP_Synthase
%'AKU27702.1'	'IB49_16340'	'ATP synthase'              ATP_Synthase
%'AKU25281.1'	'IB49_00730'	'ATP synthase subunit B'    ATP_Synthase
model = changeGrRules(model, 'ATP_Synthase', ' and IB49_09375 and IB49_16340 and IB49_00730', false);

%'AKU27863.1'	'IB49_17385'	'malate:quinone oxidoreductase'     R00342
model = changeGrRules(model, 'R00342', 'IB49_17385', false);

%'AKU26924.1'	'IB49_11305'	'glycerol phosphate lipoteichoic acid synthase' Teichuronic_acid
model = changeGrRules(model, 'Teichuronic_acid', 'IB49_11305', false);

%'AKU27134.1'	'IB49_12715'	'L-asparaginase'        R00485
model = changeGrRules(model, 'R00485', 'IB49_12715', false);

%'AKU27174.1'	'IB49_12970'	'glutamate synthase'    R00114 R00248 R00256
model = changeGrRules(model, 'R00114', 'IB49_17540 and IB49_12970', true);

% Glutamate dehydrogenase
model = changeGrRules(model, 'R00248', 'IB49_02960', true); 

%Glutamate aminohydrolase
model = changeGrRules(model, 'R00256', 'IB49_02335 and IB49_11725', true);

%'AKU27408.1'	'IB49_14520'	'UDP-glucose epimerase'  R00291
model = changeGrRules(model, 'R00291', 'IB49_14520', false);

% New reaction
%'AKU26716.1'	'IB49_09935'	'mannose-6-phosphate isomerase'
reaction.rxns = 'R00772';
reaction.equations = 'C00275_c <=> C00085_c';
reaction.rxnNames = 'mannose-6-phosphate isomerase';
reaction.lb = -1000;
reaction.ub = 1000;
reaction.eccodes = '5.3.1.8';
reaction.grRules = 'IB49_09935';
model=addRxns(model,reaction,1,'c',false,true);

%'AKU26024.1'	'IB49_05515'	'S-adenosylmethionine decarboxylase'    R00178
reaction.rxns = 'R00178';
reaction.equations = 'C00019_c + C00080_c => C01137_c + C00011_c';
reaction.rxnNames = 'S-adenosylmethionine decarboxylase';
reaction.lb = 0;
reaction.ub = 1000;
reaction.eccodes = '4.1.1.50';
reaction.grRules = 'IB49_05515';
model=addRxns(model,reaction,1,'c',false,true);


%'AKU27165.1'	'IB49_12915'	'formate dehydrogenase' R09481 
reaction.rxns = 'R09481';
reaction.equations = 'C00058_c + C00003_c => C00011_c + C00004_c + C00080_c';
reaction.rxnNames = 'formate dehydrogenase';
reaction.lb = 0;
reaction.ub = 1000;
reaction.eccodes = '1.17.1.9';
reaction.grRules = 'IB49_12915';
model=addRxns(model,reaction,1,'c',false,true);


model.grRules(find(contains(model.rxns, 'ATP_Synthase')))={'IB49_09335 and IB49_09340 and IB49_09345 and IB49_09350 and IB49_09355 and IB49_09360 and IB49_09365 and IB49_09370 and IB49_09375 and IB49_16340 and IB49_00730'};
model.grRules(find(ismember(model.rxns, 'Complex_I')))={'IB49_09280 and IB49_09330 and IB49_09325 and IB49_09315 and IB49_09310 and IB49_09305 and IB49_09300 and IB49_09295 and IB49_09290 and IB49_09285 and IB49_09280 and IB49_05305 and IB49_05310 and IB49_04735 and IB49_06865'};


metabolite.mets = {'C06148_c', 'C04895_c'};
metabolite.metNames = {'2,5-Diamino-6-(5-triphosphoryl-3,4-trihydroxy-2-oxopentyl)-amino-4-oxopyrimidine' ,'7,8-Dihydroneopterin 3-triphosphate'};
metabolite.compartments = {'c', 'c'};
metabolite.metFormulas = {'C9H18N5O14P3', 'C9H16N5O13P3'};
model = addMets(model, metabolite);

reaction.rxns = {'R02237', 'R04620', 'R04639', 'R05048', 'R00258'};
reaction.equations = {'C00002_c + C00921_c + C00025_c => C00008_c + C00009_c + C00415_c','C04895_c + 3 C00001_c => C04874_c + 3 C00009_c', 'C04895_c + C00001_c <=> C06148_c', 'C05923_c <=> C06148_c', 'C00041_c + C00026_c <=> C00022_c + C00025_c'};
reaction.rxnNames = {'dihydrofolate synthase', 'alkaline phosphatase', '2-amino-4-hydroxy-6-(erythro-1,2,3-trihydroxypropyl) dihydropteridine triphosphate hydrolase', '2,5-diaminopyrimidine nucleoside triphosphate mutase', 'Alanine transaminase'};
reaction.lb = [0, 0, -1000, -1000, -1000];
reaction.ub = [1000, 1000, 0, 1000, 1000];
reaction.eccodes = {'6.3.2.17', '3.1.3.1', '3.5.4.16', '3.5.4.16', '2.6.1.2'};
reaction.grRules = {'IB49_05065', 'IB49_05390', 'IB49_02845', 'IB49_02845', 'IB49_13625'};
model=addRxns(model,reaction,1,'c',false,true);

model = changeRxns(model, 'Cofactor_Pool', '1 C00003_c + 0.054 C00006_c + 0.025 C00019_c + 0.013 C00016_c + 0.022 C00061_c + 0.04 C00018_c + 0.013 C00010_c + 0.40 C00120_c + 0.22 C00504_c + 0.16 C00032_c => CofactorPool_c', 1, 'c', false);
model = setParam(model, 'lb', 'R00228', 0);
model.rev(find(contains(model.rxns, 'R00228'))) = 0;


% Modifications to allow exportModel to work. Otherwise gives error related
% to formatting of subsystems
for i = 1:length(model.subSystems)
    if isnumeric(model.subSystems{i})
        model.subSystems{i} = num2str(model.subSystems{i}); % Convert numeric values to strings
    elseif islogical(model.subSystems{i})
        model.subSystems{i} = num2str(model.subSystems{i}); % Convert logical values to strings ("0" or "1")
    elseif iscell(model.subSystems{i})
        % Handle cell arrays, possibly by joining cell array strings or converting numbers to strings
        cellArray = model.subSystems{i};
        strArray = cellfun(@num2str, cellArray, 'UniformOutput', false); % Convert numeric cells to strings
        model.subSystems{i} = strjoin(strArray, ', '); % Join cell array elements into a single string
    elseif ~ischar(model.subSystems{i})
        % For any other non-character type, convert to a generic placeholder string or handle specifically
        model.subSystems{i} = 'UnknownSubsystem'; % Placeholder for unknown or unhandled types
    end
    
    % No action needed if already a character vector
end


model=setParam(model, 'ub', 'R00006',0);
model.rev(find(ismember(model.rxns, 'R00006')))=0;
model = correctRev(model);
%%
rxnIdxs = findRxnIdxs(model);
model = preprocessModel(model, rxnIdxs, 35, 30); %30, 15
model = changePO(model,3,1,2,4,rxnIdxs); %4,2,6


