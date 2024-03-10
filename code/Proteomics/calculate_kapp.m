clear

load('../../../Models/eciGEL629.mat');
Proteomics = loadProtData(4,'../../../Databases/proteomics.tsv');
C13fluxes = readtable('../../../Databases/Cordova_fluxes_rxns.xlsx', 'ReadVariableNames', true);

% Rescale Cordovas fluxes to their measured qS
C13fluxes.Flux = C13fluxes.Flux.*(30.6/100);

% Find CCM reactions in ecModel
CCM.Reactions = ecModel.ec.rxns(find(contains(ecModel.ec.rxns, C13fluxes.RxnID)));
CCM.Reactions(contains(CCM.Reactions, 'REV')) = '';
CCM.IDXs = find(ismember(ecModel.ec.rxns, CCM.Reactions));
idx2remove=[];
% Find protein abundances
for i=1:length(CCM.Reactions)
    if sum(ecModel.ec.rxnEnzMat(CCM.IDXs(i), :))>1 % multi-subunit enzymes, remove.
        idx2remove = [idx2remove; i];
        % Subunits = ecModel.ec.enzymes(find(ecModel.ec.rxnEnzMat(CCM.IDXs(i), :)));
        % CCM.abundance(i) = sum(Proteomics.abundances(contains(Proteomics.uniprotIDs, Subunits)));
        % CCM.mw(i) = sum(ecModel.ec.mw(contains(ecModel.ec.enzymes, Subunits)));
    else
        Protein = ecModel.ec.enzymes(find(ecModel.ec.rxnEnzMat(CCM.IDXs(i), :)));
        CCM.abundance(i) = Proteomics.abundances(contains(Proteomics.uniprotIDs, Protein));
        CCM.mw(i) = ecModel.ec.mw(contains(ecModel.ec.enzymes, Protein));
    end
end
CCM.Reactions(idx2remove) = '';
CCM.IDXs(idx2remove) = '';
CCM.abundance(idx2remove) = '';
CCM.mw(idx2remove) = '';

CCM.abundance = CCM.abundance';
CCM.mw = CCM.mw';

% Add fluxes
for i=1:length(C13fluxes.RxnID)
    Structidxs = find(contains(CCM.Reactions, C13fluxes.RxnID(i)));
    CCM.flux(Structidxs) = C13fluxes.Flux(i);
end
CCM.flux = CCM.flux';

% Convert enzyme abundance from mg/gDW to mmol/gDW
CCM.abundance = CCM.abundance./CCM.mw;

% For isoenzymes, pool abundances of each isoenzyme
IsoRxns = CCM.Reactions(find(contains(CCM.Reactions, 'EXP')));
IsoRxns = extractBefore(IsoRxns, 7);
IsoRxns = unique(IsoRxns);
for i=1:length(IsoRxns)
    Idxs = find(contains(CCM.Reactions, IsoRxns(i)));
    CCM.abundance(Idxs) = sum(CCM.abundance(Idxs));
end

% Calculate kapp
CCM.kapp = CCM.flux./CCM.abundance;

% Rescale kapp to s-1 from h-1
CCM.kapp = CCM.kapp/3600;

kapptable = table(CCM.Reactions, CCM.kapp, 'VariableNames', {'Reactions', 'kapp'});
writetable(kapptable, '../../../Databases/kapp.tsv', 'filetype','text', 'delimiter','\t');


