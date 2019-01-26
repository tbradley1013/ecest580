%% ECES T580 Lab 3 - Tyler Bradley
clc;close all; clear;

%% Lab 3.1.1:
% 1. load the sequence in dinoDNA.txt in Matlab.
dino_dna = textread("labs/lab3_files/dinoDNA.txt", "%q");

dino_seq = strcat(dino_dna{2:length(dino_dna)});

if isfile("labs/lab3_files/dinoDNA.fasta")
    delete labs/lab3_files/dinoDNA.fasta
end

fastawrite("labs/lab3_files/dinoDNA.fasta", "Dino DNA", dino_seq);

dino_seq = fastaread("labs/lab3_files/dinoDNA.fasta");

%% Lab 3.2.1:
% 1. BLAST mystery_sequence.txt online.

% The best match was to "Francisella tularensis subsp. novicida strain
% AL97-2214, complete genome" with a 68% identity and a E value = 3e-07


%% Lab 3.2.2:
% 1. BLAST dinoDNA.txt in Matlab.
[dino_blastn, dino_ROTE] = blastncbi(dino_seq, "blastn", "database", "nr");
dino_bn = getblast(dino_blastn, "ToFile", "labs/lab3_files/dino_blast.rpt");

% The best blast match
dino_bn.Hits(1).Definition

%% Lab 3.2.3:
% Run the code in Ex 3.2.2 and answer the following questions:
% 1. What was the difference between the blastn and tblastx searches?
% 2. Why were the shorter lengths chosen?
% 3. For the full length BLAST, what organisms were closely related? What taxa do these belong to?
% 4. Did this sequence share homology (similarity) with any known and functionally annotated genes?
% 5. What was the difference between the tblastx searches performed by the full-length sequence vs.
% 300 bp sequence vs. 40 bp sequence?

% Perform a tblastx and retrieve the report
mys_dna = textread("labs/lab3_files/mystery_sequence.txt", "%q");

mys_seq = strcat(mys_dna{2:length(mys_dna)});

if isfile("labs/lab3_files/mys_seq.fasta")
    delete labs/lab3_files/mys_seq.fasta
end

fastawrite("labs/lab3_files/mys_seq.fasta", "Mystery Sequence", mys_seq);

mys_seq = fastaread("labs/lab3_files/mys_seq.fasta");
[Data_tblastx, RTOE_tx] = blastncbi("labs/lab3_files/mys_seq.fasta", 'tblastx', 'database', 'nr');
tblastx = getblast(Data_tblastx);
% Perform a tblastx for the first 40 base pairs
Seq_40 = mys_seq.Sequence(1:40);
if isfile("labs/lab3_files/seq_40.fasta")
    delete labs/lab3_files/seq_40.fasta
end
fastawrite('labs/lab3_files/seq_40.fasta', 'Sequence', Seq_40);
Seq_40 = fastaread('labs/lab3_files/seq_40.fasta');
[Data_tblastx_40, RTOE_40] = blastncbi("labs/lab3_files/seq_40.fasta", 'tblastx', 'database', 'nr');
blast_40 = getblast(Data_tblastx_40);
% Perform a tblastx for the first 300 base pairs
Seq_300 = mys_seq.Sequence(1:300);
if isfile("labs/lab3_files/seq_300.fasta")
    delete labs/lab3_files/seq_300.fasta
end
fastawrite('labs/lab3_files/seq_300.fasta', 'Sequence', Seq_300);
Seq_300 = fastaread('labs/lab3_files/seq_300.fasta');
[Data_tblastx_300, RTOE_300] = blastncbi("labs/lab3_files/seq_300.fasta", 'tblastx', 'database', 'nr');
blast_300 = getblast(Data_tblastx_300);

% Answers:
% 1. The difference between the two blast methods is that blastn compares
% the sequence to reference gene's nucleotides while tblastx is translating
% the nucleotide sequence into amino acid and comparing the results to
% proteins in ncbi's database
%
% 2. Shorter lengths were chosen because they were comparing amino acids
% rather than nucleotide and every amino acid represents a combination of 3
% nucleotides
%
% 3. The top five closest matches are shown here
full_close_species = repelem("a", 5);
for i = 1:5
    full_close_species(i) = tblastx.Hits(i).Definition;
end
full_close_species

% 4. 