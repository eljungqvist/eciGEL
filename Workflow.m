% Workflow for generating eciGEL from iGEL604

% 1. Load conventional GEM iGEL604
% load("/Users/emilljungqvist/Modeling/iGEL604.mat")

% 2. Manual Curation (code/Preprocessing)
%       - Update available carbon sources from Vitamin manuscript (expand_carbon_sources.m)
%       - Add new gene-associations based on proteomics. Converts from Cobra to Raven structure (Add_genes_from_proteomics.m)
%       - Correct central carbon metabolism gene associations (Correct_grRules_CCC.m)
% 
% 3. Update Pirt Parameters manually (code/Pirt Parameters) 

% 4. Export model to .xml with exportModel. Gives errors that are fake
%    news.

% 5. Create ecModel (code/GeckoInput/iiGEL604_ec_script.m)

% 6. Calculate k_apps for central carbon metabolism reactions (code/Proteomics/calculate_kapp.m)

% 6. Update kcats (code/Bayesian/Bayesian_fitting.m)

% 7. Implement proteomics


