Proteomics_data = readtable('../../../Databases/Proteome_allocation_rep_original.tsv',"FileType","text",'Delimiter', '\t');
%load('../../../Models/eciGEL627_posterior.mat');

%model_proteins = ecModel.ec.enzymes;

ReplicateA = Proteomics_data(find(strcmp(Proteomics_data.rep,'Geo_0_1')),:);
ReplicateB = Proteomics_data(find(strcmp(Proteomics_data.rep,'Geo_0_2')),:);
ReplicateC = Proteomics_data(find(strcmp(Proteomics_data.rep,'Geo_0_3')),:);
ReplicateD = Proteomics_data(find(strcmp(Proteomics_data.rep,'Geo_0_4')),:);

prot_ints_replicates_idxs = ["ProteinID", "ReplicateA", "ReplicateB", "ReplicateC", "ReplicateD"];


%% Extract Relative abundance

% Extract relative abundance for each protein. Gecko wants to calculate
% stddev etc. at a later stage, so need to keep data for each replicate.

for i=1:length(ReplicateA.protein)

    repA_abundance(i) = ReplicateA.fraction(i);

    if any(strcmp(ReplicateB.protein, ReplicateA.protein(i)))
        repB_abundance(i) = ReplicateB.fraction(find(strcmp(ReplicateB.protein, ReplicateA.protein(i))));
    else
        repB_abundance(i) = NaN;
    end

    if any(strcmp(ReplicateC.protein, ReplicateA.protein(i)))
        repC_abundance(i) = ReplicateC.fraction(find(strcmp(ReplicateC.protein, ReplicateA.protein(i))));
    else
        repC_abundance(i) = NaN;
    end

    if any(strcmp(ReplicateD.protein, ReplicateA.protein(i)))
        repD_abundance(i) = ReplicateD.fraction(find(strcmp(ReplicateD.protein, ReplicateA.protein(i))));
    else
        repD_abundance(i) = NaN;
    end
end

% Convert Relative abundances to Absolute abundances (mg prot / g CDW).
% Protein is 51% of CDW, or 510 mg / g.
repA_abundance = repA_abundance*510;
repB_abundance = repB_abundance*510;
repC_abundance = repC_abundance*510;
repD_abundance = repD_abundance*510;

All_proteins = table(ReplicateA.protein, repA_abundance', repB_abundance', repC_abundance', repD_abundance');
All_proteins.Properties.VariableNames = prot_ints_replicates_idxs;
writetable(All_proteins, '../../../Databases/proteomics.tsv', 'filetype','text', 'delimiter','\t');


