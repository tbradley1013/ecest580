%% ECES T580 Homework 1 - Tyler Bradley
clc;close all;clear;
%% Question 1
% a) 
% Find the protein sequence of the hypothetical protein AF1226 precursor of Archaeoglobus
% fulgidus. What is the sequence of amino acids (in single letter representation) from
% positions 141-147?
af_1226 = getgenpept('O29042.1');

seq = af_1226.Sequence;

seq_sub = seq(141:147);

% b)
% How many nucleotide sequences can give rise to an arbitrary amino acid sequence (hint:
% revgeneticcode)? Write a function that returns the possible sequences. The input to
% the function should be a string (amino acid sequence) and the output should be the
% potential nucleotide sequences.
revcode = revgeneticcode();

temp_len = length(seq_sub);
possible_seqs = [""];
len_out = 1;
for i = 1:temp_len
    prot = seq_sub(i);
    
    prot_seq = revcode.(upper(prot));
    
    [x, y] = meshgrid(possible_seqs, prot_seq);
    pairs = [x(:) y(:)];
    
    temp_out = repelem("a", length(pairs));
    for j = 1:length(pairs)
      temp_out(j) = pairs(j, 1) + pairs(j, 2);
    end
    possible_seqs = temp_out;
end



function output = rev_seq(seq)
revcode = revgeneticcode();

temp_len = length(seq);
possible_seqs = [""];
len_out = 1;
for i = 1:temp_len
    prot = seq(i);
    
    prot_seq = revcode.(upper(prot));
    
    [x, y] = meshgrid(possible_seqs, prot_seq);
    pairs = [x(:) y(:)];
    
    temp_out = repelem("a", length(pairs));
    for j = 1:length(pairs)
      temp_out(j) = pairs(j, 1) + pairs(j, 2);
    end
    possible_seqs = temp_out;
end

end