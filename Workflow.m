% Workflow for generating eciGEL from iGEL604

% 1. Load conventional GEM iGEL604
% load("/Users/emilljungqvist/Modeling/iGEL604.mat")

% 2. Manual Curation (ecLC300/Scripts/Preprocessing)
%       - Update available carbon sources from Vitamin manuscript (expand_carbon_sources)
%       - Add new gene-associations based on proteomics. Converts from Cobra to Raven structure (Add_genes_from_proteomics)
%       - Correct central carbon metabolism gene associations (Correct_grRules_CCC)
% 
% 3. Update Pirt Parameters (ecLC300/Scripts/Pirt Parameters)

% 4. Export model to .xml with exportModel. Gives errors that are fake
%    news.

% 5. Import Proteomics data and calculate f-factor. Modify model adapter
%    with caluclated value in step 6.
% 
% 6. Create ecModel (Gecko-3.1.3/eciGEL/code/iiGEL604_ec_script.m)

% 7. Update kcats (ecLC300/Scripts/Bayesian/Bayesian_fitting.m)

% 8. Implement proteomics


