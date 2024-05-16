% Create uniprot database of all predicted proteins in G. LC300 genomes.
% Requires matlab bioinformatics addon

[Header, Sequence] = fastaread('../../../Databases/FASTA_GLC300.fasta');

% Load Martin's python-generated uniprot file as it contains essential
% pseduoproteins
old_uniprot = readtable('../../../../GECKO-3.0.1/userData/iGEL604/data/uniprot_old.tsv','FileType','text','Delimiter','\t');

% Goal is to get a file with ProteinID, GeneID, MW, Sequence

for i=1:length(Header)
    
    ProteinID_idx = strfind(Header(i), 'protein_id=');
    ProteinID(i) = extractBetween(Header(i),ProteinID_idx{1}+11,ProteinID_idx{1}+20);
    
    locustag_idx = strfind(Header(i), 'locus_tag=');
    gene(i) = extractBetween(Header(i),locustag_idx{1}+10,locustag_idx{1}+19);

    MW(i) = molweight(Sequence{i});
    
    EC_number(i) = {'N/A'};
end

 
TableHeader = ["Entry", "Gene names", "EC Number", "Mass", "Sequence"];
ProteinDb = table(ProteinID', gene', EC_number', MW', Sequence');
ProteinDb.Properties.VariableNames = TableHeader;
old_uniprot.Properties.VariableNames = TableHeader;

% Append a few genes from the old uniprot database that martin made.
ProteinDb = [ProteinDb; old_uniprot(find(~ismember(old_uniprot.("Gene names"), ProteinDb.("Gene names"))),:)];

writetable(ProteinDb, '../../../Databases/uniprot.tsv', 'filetype','text', 'delimiter','\t');

