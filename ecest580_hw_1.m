%% ECES T580 Homework 1 - Tyler Bradley

%% Question 1
% a) 
% Find the protein sequence of the hypothetical protein AF1226 precursor of Archaeoglobus
% fulgidus. What is the sequence of amino acids (in single letter representation) from
% positions 141-147?
af_1226 = getgenpept('O29042.1');

seq = af_1226.Sequence;

seq_sub = seq(141:147)

% b)
% How many nucleotide sequences can give rise to an arbitrary amino acid sequence (hint:
% revgeneticcode)? Write a function that returns the possible sequences. The input to
% the function should be a string (amino acid sequence) and the output should be the
% potential nucleotide sequences.
function output = rev_seq(seq)


end