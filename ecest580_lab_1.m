%% Lab 1 - Tyler Bradley
clc;close all;clear;

%% Lab 1.1.1 
% Download sequence with accession number nm_000520 from GenBank 
% and load it in Matlab.

% Download the sequence information from GenBank and save to file
getgenbank("nm_000520", "ToFile", "NM000520.txt")

% read the saved file
s=genbankread("NM000520.txt")

% Download the sequence only
seq=getgenbank("nm_000520", "SequenceOnly", true)

%% Lab 1.2.1
% 1. Repeat Ex 1.2.1 for sequence nm_000520.
% 2. Default window size for function ntdensity is length(Seq)/20. Try different windows. What is the
% advantage or disadvantage of longer windows?

% Format the long sequence output for easy viewing
seqdisp(s.Sequence)

% Count the nucleotides in sequence
[seq_counts]=basecount(s.Sequence)

% Plot density of nucleotides along sequence
figure(1)
seq_density = ntdensity(s.Sequence)

% Count dimers in nucleotide sequence
figure(2)
[Dimers, Percent] = dimercount(s.Sequence, "chart", "pie")

% Count 3-mer in nucleotide sequence
trimer = nmercount(s.Sequence, 3)

% Trying different window sizes for ntdensity
% doubling the defualt window size
figure(3)
ntdensity(s.Sequence, "Window", round(length(s.Sequence)/10))

% halving the default window size
figure(4)
ntdensity(s.Sequence, "Window", round(length(s.Sequence)/40))

% The advantages and disadvatages to longer window sizes go hand in hand.
% There may be trends in nucleotide density that is missed or has its
% effect dampened if the window is too large and similarly if the window is
% too small it may underestimate the size of the effect for a given trend
% in necleotide densities. 


%% Lab 1.3.1
% count the codons in each of the six reading frames, and plot the results
% in a heat map for sequence nm_000520

figure(1)
r1codons = codoncount(s.Sequence, "frame", 1, "figure", true)

figure(2)
r2codons = codoncount(s.Sequence, "frame", 2, "figure", true)

figure(3)
r3codons = codoncount(s.Sequence, "frame", 3, "figure", true)

% There are no recognized 4th, 5th, or 6th reading frames for nm_000520
% r4codons = codoncount(s.Sequence, "frame", 4, "figure", true)
% r5codons = codoncount(s.Sequence, "frame", 5, "figure", true)
% r6codons = codoncount(s.Sequence, "frame", 6, "figure", true)


%% Lab 1.3.2
% 1. Find the ORFs of length > 50 in Frame 1 for sequence nm_000520.
orf_great_50 = seqshoworfs(s.Sequence, "Frame", 1, "MinimumLength", 50)

% 2. Find the ORFs of length > 500 in Frame 1 for sequence nm_000520.
orf_great_500 = seqshoworfs(s.Sequence, "Frame", 1, "MinimumLength", 500)



%% Lab 1.3.3
% Estimate P(stop) from the sequence nm_000520 and determine the threshold given ? = 0.05
% P(k nonstops) = (1 ? P(stop))^k
% P(k nonstops) <= alpha
% get the count of each codon
all_codon_n = codoncount(s.Sequence)

% sum the total count of all stop codons in the sequence
stop_n = all_codon_n.TTA + all_codon_n.TAG + all_codon_n.TGA

% get the total number of codons in the sequence
total_n = sum(cell2mat(struct2cell(all_codon_n)))

% find the frequency of stop codons i.e. P(stop)
p_stop = stop_n/total_n

% from the example, k >= log(alpha)/log(1-p_stop)
k = log(0.05)/log(1-p_stop)

% add 2 codons to k for the start and stop codons in a ORF
k_final = k + 2



