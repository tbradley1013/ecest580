%% Homework 2
% Tyler Bradley

value1 = {'X', 'Y', 'Z', 'W'};
value2 = {'TACCCGAT', 'TAAACGAT', 'AAAACGCG', 'AAAACGAT'};

tree_seqs = struct('Header', value1, 'Sequence', value2);

% original tree
seq_dist = seqpdist(tree_seqs);
tree = seqlinkage(seq_dist, 'average', tree_seqs);
view(tree)


num_boot = 5;
seq_len = length(tree_seqs(1).Sequence);
num_seqs = length(tree_seqs);
boots = cell(num_boot, 1);
boots_dist = cell(num_boot, 1);
boots_trees = cell(num_boot, 1);


for i = 1:num_boot
    idx = randsample(seq_len, seq_len, 'true');
    for j = 1:num_seqs
        boot_seq(j).Header = tree_seqs(j).Header;
        boot_seq(j).Sequence = tree_seqs(j).Sequence(idx);
    end
    boots{i} = boot_seq;
    boot_dist = seqpdist(boot_seq);
    boots_dist{i} = boot_dist;
    boot_tree = seqlinkage(boot_dist, 'average', boot_seq);
    boots_trees{i} = boot_tree;
end



%% Question 3
% Go to a database of E. Coli binding sites: http://arep.med.harvard.edu/ecoli_matrices/.
% Click on the lexA link. Then, the alignment link contains the list of binding sites.
% Unfortunately between each sequence, there are notes about the sequence (you will need
% to strip these when importing the sequence into Matlab, hint: fastaread).
lexA = fastaread("hw2-files/lexA.fasta");

% Find R_sequence(l) for this alignment
len_seq = length(lexA(1).Sequence);
num_seqs = length(lexA);

% create empty vectors for output values
prob_A = repelem(0, len_seq);
prob_C = repelem(0, len_seq);
prob_G = repelem(0, len_seq);
prob_T = repelem(0, len_seq);
H = repelem(0, len_seq);
b_table = zeros(len_seq, 4);

for i = 1:len_seq
  count_A = 0;
  count_C = 0;
  count_G = 0;
  count_T = 0;
  
  for j = 1:num_seqs
      if lexA(j).Sequence(i) == "a"
          count_A = count_A + 1;
      elseif lexA(j).Sequence(i) == "c"
          count_C = count_C + 1;
      elseif lexA(j).Sequence(i) == "g"
          count_G = count_G + 1;
      elseif lexA(j).Sequence(i) == "t"
          count_T = count_T + 1;
      end
  end
  
  prob_A(i) = count_A/num_seqs + 0.0000001;
  prob_C(i) = count_C/num_seqs + 0.0000001;
  prob_G(i) = count_G/num_seqs + 0.0000001;
  prob_T(i) = count_T/num_seqs + 0.0000001;
  H(i) = -1*(prob_A(i)*log2(prob_A(i)) + prob_C(i)*log2(prob_C(i)) + ...
             prob_G(i)*log2(prob_G(i)) + prob_T(i)*log2(prob_T(i)));
         
  R_seq(i) = 2 - H(i);
         
  b_table(i, 1) = prob_A(i)*R_seq(i);
  b_table(i, 2) = prob_C(i)*R_seq(i);
  b_table(i, 3) = prob_G(i)*R_seq(i);
  b_table(i, 4) = prob_T(i)*R_seq(i);
end

%R_seq = 2 - H;

bar(R_seq)
ylabel("R(l) (bits)")
title("Sequence Logo Plot")
xlabel("Nucleotide position")

rownames = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', ...
    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'};
array2table(b_table, "VariableNames", {'A', 'C', 'G', 'T'}, "RowNames", rownames)

% The plot created above is fairly similar to the one created using the
% online tool. Bases 3, 4, 5 and 16, 17, 18 are the highest groupings 
% of R values in both of the plots. The e(n) value in the online tool 
% (set to zero here) likely corrects for potential bias that may be
% introduced when only a small amount of data is used to determine the
% entropy in a sequence alignment.
