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





