%% Script to add reactions to iGEL604, to make it capable of growing on the experimentally proven carbon sources
load("/Users/emilljungqvist/Modeling/iGEL604.mat")
changemetnames
%%
% Maltose import

metID='C00208[e]';
metName='Maltose[e]';
Formula='C12H22O11';
KEGGID='C00208';
PubChemID='6255';
InChi='InChI=1S/C12H22O11/c13-1-3-5(15)6(16)9(19)12(22-3)23-10-4(2-14)21-11(20)8(18)7(10)17/h3-20H,1-2H2/t3-,4-,5-,6+,7-,8-,9-,10-,11?,12-/m1/s1';
Charge=0;
model=addMetabolite(model, metID, 'metName', metName, 'metFormula', Formula,...
    'KEGGId', KEGGID, 'PubChemID', PubChemID, 'InChi', InChi, 'Charge', Charge)
clear metID metName Formula KEGGID PubChemID InChi Charge

model = addExchangeRxn(model, 'C00208[e]', -1000, 0)

% Maltose import through ABC transporter. Blasts to b. sub. maltodextrin
% transporter
ReactionID='Maltose_ABC';
ReactionName='Maltose abc transport';    
ReactionFormula='C00208[c] + 2 C00009[c] + 2 C00008[c] + 2 C00080[c] => C00208[e] + 2 C00002[c] + 2 C00001[c]';
Reversibility=false;
lb=-1000;
ub=0;
grRules='IB49_13790 and IB49_13795 and IB49_13800 and IB49_13425';
model=addReaction(model, ReactionID, 'reactionName',ReactionName,...
    'reactionFormula', ReactionFormula, 'reversible', Reversibility,...
    'lowerBound', lb, 'upperBound', ub, 'geneRule', grRules);


% adding reactions for gene IB49_13340 (alpha-glucosidase EC 3.2.1.20)
ReactionID='R00028';
ReactionName='maltose glucohydrolase';    
ReactionFormula='C00208[c] + C00001[c] <=> 2 C00031[c]';
Reversibility=false;
lb=0;
ub=1000;
grRules='IB49_13340';
model=addReaction(model, ReactionID, 'reactionName',ReactionName,...
    'reactionFormula', ReactionFormula, 'reversible', Reversibility,...
    'lowerBound', lb, 'upperBound', ub, 'geneRule', grRules);

%% Sucrose metabolism

%Fructokinase
ReactionID='R00760';
ReactionName='fructokinase';    
ReactionFormula='C00002[c] + C00095[c] <=> C00008[c] + C00085[c]';
Reversibility=true;
lb=-1000;
ub=1000;
grRules='IB49_08420';
model=addReaction(model, ReactionID, 'reactionName',ReactionName,...
    'reactionFormula', ReactionFormula, 'reversible', Reversibility,...
    'lowerBound', lb, 'upperBound', ub, 'geneRule', grRules);
model=changeRxnBounds(model, 'EX_Xylose', 0, 'l');

%% Glycerol metabolism

% Glycerol kinase. Found by blasting GlpK from Bsub to LC300 genome. Not a
% gene found in the original automatic annotation from Cordova et al. Lies
% between IB49_17040 and IB49_17050, so suggested new name is IB49_17045.

ReactionID='R00847';
ReactionName='glycerol kinase';    
ReactionFormula='C00002[c] + C00116[c] <=> C00008[c] + C00093[c]';
Reversibility=false;
lb=0;
ub=1000;
grRules='IB49_17045';
model=addReaction(model, ReactionID, 'reactionName',ReactionName,...
    'reactionFormula', ReactionFormula, 'reversible', Reversibility,...
    'lowerBound', lb, 'upperBound', ub, 'geneRule', grRules);


%%
% Mannitol import and exchange
metID='C00392[e]';
metName='Mannitol[e]';
Formula='C6H14O6';
KEGGID='C00392';
PubChemID='6251';
InChi='InChI=1S/C6H14O6/c7-1-3(9)5(11)6(12)4(10)2-8/h3-12H,1-2H2/t3-,4-,5-,6-/m1/s1';
Charge=0;
model=addMetabolite(model, metID, 'metName', metName, 'metFormula', Formula,...
    'KEGGId', KEGGID, 'PubChemID', PubChemID, 'InChi', InChi, 'Charge', Charge)
clear metID metName Formula KEGGID PubChemID InChi Charge
model = addExchangeRxn(model, 'C00392[e]', -1000, 1000);

ReactionID='Mannitol_PTS';
ReactionName='Mannitol import PTS';    
ReactionFormula='C00644[c] + C00022[c] => C00392[e] + C00074[c]';
Reversibility=false;
lb=-1000;
ub=0;
grRules='IB49_01725 and IB49_01740';
model=addReaction(model, ReactionID, 'reactionName',ReactionName,...
    'reactionFormula', ReactionFormula, 'reversible', Reversibility,...
    'lowerBound', lb, 'upperBound', ub, 'geneRule', grRules);

ReactionID='R00758';
ReactionName='D-mannitol-1-phosphate:NAD+ 5-oxidoreductase';    
ReactionFormula='C00644[c] + C00003[c] <=> C00085[c] + C00004[c] + C00080[c]';
Reversibility=true;
lb=-1000;
ub=1000;
grRules='IB49_01720';
model=addReaction(model, ReactionID, 'reactionName',ReactionName,...
    'reactionFormula', ReactionFormula, 'reversible', Reversibility,...
    'lowerBound', lb, 'upperBound', ub, 'geneRule', grRules);


%% Cellobiose exchange
model.rxns(find(contains(model.rxns, 'C00185[e]'))) = {'EX_Cellobiose'};

%%
% Remove 3HB and Ethanol bullshit stuff
model = removeRxns(model, {'EX_3HB', '3HB_thioesterase', 'EX_Ethanol'});

%%
model.lb(find(contains(model.rxns, 'EX_Glucose'))) = -30.06;
model.lb(find(contains(model.rxns, 'EX_C00185[e]'))) = 0;
model.lb(find(contains(model.rxns, 'EX_C00124[e]'))) = 0;
model.lb(find(contains(model.rxns, 'EX_C00392[e]'))) = 0;
model.lb(find(contains(model.rxns, 'EX_C00208[e]'))) = 0;

sol=solveLP(model, 1);
printFluxes(model, sol.x, false);
