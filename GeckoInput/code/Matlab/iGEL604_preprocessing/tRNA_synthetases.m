% 2023-05-19 Martin Gustavsson
% Script that adds tRNA synthetase reactions and the associated genes to
% iGEL604 in order to explicitly model the enzyme cost of tRNA activation
% when transitioning into a Gecko-format model of G. LC300. Manually adds
% 20 tRNA-charging reactions and updates the "aa_to_protein" reaction to
% use aminoacyl-tRNAs instead of aas directly. 

%% Find directory paths  - can be replaced by making this a function instead

mfilePath = mfilename('fullpath');
if(contains(mfilePath,'LiveEditorEvaluationHelper'));
    mfilePath = matlab.desktop.editor.getActiveFilename;
end

codeDir = fileparts(mfilePath);
modelDir = fullfile(codeDir, '../Models/');

%% Load parent model

modelName = 'iGEL604.mat';

load(fullfile(modelDir, modelName));
   %modelDir = ]
% Dummy code

%% Add 20 reactions including gene association

% For now, makes aa + ATP -> AMP + ppi + aa-tRNA, excluding the consumption
% and later reformation of tRNA since this will anyway always be
% steady-state and just adds extra metabolites (especially if each aa
% should have it's unique tRNAs)

% Asn
updatedModel = addReaction(model, 'R03648', 'reactionFormula', 'C00152[c] + C00002[c] -> C00020[c] + C00013[c] + C03402[c]', 'geneRule', 'IB49_02625', 'reactionName', 'Asparaginyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.22'};

% Ala
updatedModel = addReaction(updatedModel, 'R03038', 'reactionFormula', 'C00041[c] + C00002[c] -> C00020[c] + C00013[c] + C00886[c]', 'geneRule', 'IB49_04645', 'reactionName', 'Alanyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.7'};

% Asp
updatedModel = addReaction(updatedModel, 'R05577', 'reactionFormula', 'C00049[c] + C00002[c] -> C00020[c] + C00013[c] + C02984[c]', 'geneRule', 'IB49_04725', 'reactionName', 'Aspartyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.12'};

% His
updatedModel = addReaction(updatedModel, 'R03655', 'reactionFormula', 'C00135[c] + C00002[c] -> C00020[c] + C00013[c] + C02988[c]', 'geneRule', 'IB49_04730', 'reactionName', 'Histidyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.21'};

% Val
updatedModel = addReaction(updatedModel, 'R03665', 'reactionFormula', 'C00183[c] + C00002[c] -> C00020[c] + C00013[c] + C02554[c]', 'geneRule', 'IB49_05070', 'reactionName', 'Valyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.9'};

% Thr
updatedModel = addReaction(updatedModel, 'R03663', 'reactionFormula', 'C00188[c] + C00002[c] -> C00020[c] + C00013[c] + C02992[c]', 'geneRule', 'IB49_05480', 'reactionName', 'Threonyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.3'};

% Tyr
updatedModel = addReaction(updatedModel, 'R02918', 'reactionFormula', 'C00082[c] + C00002[c] -> C00020[c] + C00013[c] + C02839[c]', 'geneRule', 'IB49_05870', 'reactionName', 'Tyrosyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.1'};

% Leu
updatedModel = addReaction(updatedModel, 'R03657', 'reactionFormula', 'C00123[c] + C00002[c] -> C00020[c] + C00013[c] + C02047[c]', 'geneRule', 'IB49_06080', 'reactionName', 'Leucyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.4'};

% Arg
updatedModel = addReaction(updatedModel, 'R03646', 'reactionFormula', 'C00062[c] + C00002[c] -> C00020[c] + C00013[c] + C02163[c]', 'geneRule', 'IB49_09560', 'reactionName', 'Arginyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.19'};

% Gly
updatedModel = addReaction(updatedModel, 'R03654', 'reactionFormula', 'C00037[c] + C00002[c] -> C00020[c] + C00013[c] + C02412[c]', 'geneRule', 'IB49_09715', 'reactionName', 'Glycyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.14'};

% Ser
updatedModel = addReaction(updatedModel, 'R03662', 'reactionFormula', 'C00065[c] + C00002[c] -> C00020[c] + C00013[c] + C02553[c]', 'geneRule', 'IB49_10235', 'reactionName', 'Seryl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.11'};

% Met
updatedModel = addReaction(updatedModel, 'R03659', 'reactionFormula', 'C00073[c] + C00002[c] -> C00020[c] + C00013[c] + C02430[c]', 'geneRule', 'IB49_10355', 'reactionName', 'Methionyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.10'};

% Lys
updatedModel = addReaction(updatedModel, 'R03658', 'reactionFormula', 'C00047[c] + C00002[c] -> C00020[c] + C00013[c] + C01931[c]', 'geneRule', 'IB49_10590', 'reactionName', 'Lysyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.6'};

% Glu
updatedModel = addReaction(updatedModel, 'R05578', 'reactionFormula', 'C00025[c] + C00002[c] -> C00020[c] + C00013[c] + C02987[c]', 'geneRule', 'IB49_10700', 'reactionName', 'Nondiscriminating glutamyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.24'}; % This enzyme is annotated as 
% nondiscriminating and can make Glutamyl-tRNA with glutamine tRNA instead.

% Cys
updatedModel = addReaction(updatedModel, 'R03650', 'reactionFormula', 'C00097[c] + C00002[c] -> C00020[c] + C00013[c] + C03125[c]', 'geneRule', 'IB49_10710', 'reactionName', 'Cysteinyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.16'};

% Trp
updatedModel = addReaction(updatedModel, 'R03664', 'reactionFormula', 'C00078[c] + C00002[c] -> C00020[c] + C00013[c] + C03512[c]', 'geneRule', 'IB49_14365', 'reactionName', 'Tryptophanyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.2'};

% Ile
updatedModel = addReaction(updatedModel, 'R03656', 'reactionFormula', 'C00407[c] + C00002[c] -> C00020[c] + C00013[c] + C03127[c]', 'geneRule', 'IB49_15960', 'reactionName', 'Isoleucine-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.5'};

% Pro
updatedModel = addReaction(updatedModel, 'R03661', 'reactionFormula', 'C00148[c] + C00002[c] -> C00020[c] + C00013[c] + C02702[c]', 'geneRule', 'IB49_16515', 'reactionName', 'Prolyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.15'};

% Phe
updatedModel = addReaction(updatedModel, 'R03660', 'reactionFormula', 'C00079[c] + C00002[c] -> C00020[c] + C00013[c] + C03511[c]', 'geneRule', 'IB49_05415 and IB49_05420', 'reactionName', 'Phenylalanyl-tRNA synthetase');
updatedModel.rxnECNumbers(end) = {'6.1.1.20'};

% Gln
% Annotated as going from Glu-tRNA to Gln-Trna by transfer of amide
% from Gln, as in Bacillus subtilis
% Two annotated subunits in LC300 genome, missing gatB. Manual blast finds
% a truncated gatB sequence, need to manually add this with gene name and
% sequence later, but for now this is left out (making the protein complex
% smaller than it actually is, which gives a lower protein cost in the EC model).

updatedModel = addReaction(updatedModel, 'R03905', 'reactionFormula', 'C06112[c] + C00064[c] + C00002[c] + C00001[c] -> C02282[c] + C00025[c] + C00009[c] + C00008[c]', 'geneRule', 'IB49_11825 and 11830', 'reactionName', 'Glutaminyl-tRNA synthase (glutamine-hydrolysing)');
updatedModel.rxnECNumbers(end) = {'6.3.5.7'};

% Glu-tRNA(Gln)
% Nonspecific reaction where Glu is instead charged onto tRNA(Gln) for
% later use in the amidation reaction above (R03905)
updatedModel = addReaction(updatedModel, 'R03651', 'reactionFormula', 'C00025[c] + C00002[c] -> C00020[c] + C00013[c] + C06112[c]', 'geneRule', 'IB49_10700', 'reactionName', 'L-Glutamate:tRNA(Gln) ligase (AMP-forming)');
updatedModel.rxnECNumbers(end) = {'6.1.1.24'}; % This enzyme is annotated as 
% nondiscriminating and can make Glutamyl-tRNA with glutamine tRNA instead.


%% Update "aa_to_protein" reaction

% Should the water still be in here on the reactant side?
updatedModel = addReaction(updatedModel, 'aa_to_protein', 'reactionFormula', '8.455 C00001[c] + 0.5669 C02987[c] + 0.8787 C02412[c] + 1.075 C00886[c] + 18.455 C00044[c] + 0.5813 C01931[c] + 0.4037 C02984[c] + 0.4253 C02163[c] + 0.5674 C02282[c] + 0.4 C02553[c] + 0.2087 C02430[c] + 0.1181 C03512[c] + 0.2884 C03511[c] + 0.2455 C02839[c] + 0.1317 C03125[c] + 0.7333 C02047[c] + 0.1727 C02988[c] + 0.3865 C02702[c] + 0.4042 C03402[c] + 0.6649 C02554[c] + 0.4961 C02992[c] + 0.4789 C03127[c] 	->	18.455 C00009[c] + 18.455 C00035[c] + 1.146 Protein[c] ', 'reactionName', 'Proteins synthesis from amino acids');
%% Save model
outfile = fullfile(modelDir, 'iGEL626.mat')
save(outfile, "updatedModel")