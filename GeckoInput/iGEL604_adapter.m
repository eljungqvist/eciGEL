classdef iGEL604_adapter < ModelAdapter
    methods
        function obj = iGEL604_adapter()
            % Set initial values of the obj.params - they can be changed by the user
            
            % Directory where all model-specific files and scripts are kept.
            % Is assumed to follow the GECKO-defined folder structure. The
            % code below refers to userData/yourModel in the GECKO path.
            geckoPath = findGECKOroot;
            obj.params.path = fullfile(geckoPath,'eciGEL');

			% Path to the conventional GEM that this ecModel will be based on.
			%obj.params.convGEM = fullfile(obj.params.path,'models','iGEL604_before_ec.mat');
            obj.params.convGEM = fullfile(obj.params.path,'models','iGEL629.xml');


			% Average enzyme saturation factor
			obj.params.sigma = 0.5;

			% Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.51;

			% Fraction of enzymes in the model [g enzyme/g protein]
			obj.params.f = 0.4261;
            
            % Growth rate the model should be able to reach when not
            % constraint by nutrient uptake (e.g. max growth rate) [1/h]
			obj.params.gR_exp = 2.24;

			% Provide your organism scientific name
			obj.params.org_name = 'Geobacillus sp. LC300';
            
            % Matching name for Complex Portal
            obj.params.complex_org_name = '';

			% Provide your organism KEGG ID, selected at
			% https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.kegg.ID = 'gel';
            % Field for KEGG gene identifier; should match the gene
            % identifiers used in the model. With 'kegg', it takes the
            % default KEGG Entry identifier (for example YER023W here:
            % https://www.genome.jp/dbget-bin/www_bget?sce:YER023W).
            % Alternatively, gene identifiers from the "Other DBs" section
            % of the KEGG page can be selected. For example "NCBI-GeneID",
            % "UniProt", or "Ensembl". Not all DB entries are available for
            % all organisms and/or genes.
            obj.params.kegg.geneID = 'kegg';

			% Provide what identifier should be used to query UniProt.
            % Select proteome IDs at https://www.uniprot.org/proteomes/
            % or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprot.type = 'taxonomy_id'; % 'proteome' or 'taxonomy_id'
			obj.params.uniprot.ID = '1519377'; % This is now for E coli K12 generally
            % should match the ID type
            % Field for Uniprot gene ID - should match the gene ids used in the 
            % model. It should be one of the "Returned Field" entries under
            % "Names & Taxonomy" at this page: https://www.uniprot.org/help/return_fields
            obj.params.uniprot.geneIDfield = 'gene_oln';
            % Whether only reviewed data from UniProt should be considered.
            % Reviewed data has highest confidence, but coverage might be (very)
            % low for non-model organisms
            obj.params.uniprot.reviewed = true;

			% Reaction ID for glucose exchange reaction (or other preferred carbon source)
			obj.params.c_source = 'EX_Glucose'; 

			% Reaction ID for biomass pseudoreaction
			obj.params.bioRxn = 'Biomass';

			% Reaction ID for non-growth associated maitenance pseudoreaction
			obj.params.NGAM = 'R00086';

			% Compartment name in which the added enzymes should be located
			obj.params.enzyme_comp = 'Cytosol';
        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            % Indicates how spontaneous reactions are identified. Here it
            % is done by the reaction have 'spontaneous' in its name.
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end
