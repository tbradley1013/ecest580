
%% Lab 5.2.1 Matlab Portion
seqs = fastaread("SARS_Data_Part_2.fasta");

seq_dist = seqpdist(seqs);

tree = seqlinkage(seq_dist, "average", seqs);

view(tree)

% Overall, the two tree building methods make the same major groupings of 
% sequences. The two December 2002 sequences are separated from the rest of
% the sequences indicating that they are likely they starting point and the
% norovirus evolved from that point. However, way in which they are
% separated is slightly different. In the Cipres tree, it appears that the
% 12/16/02 sequence is more closely related to the future sequences, but in
% the matlab tree, they appear to be equally distant from the rest of the
% sequences. The rest of the groupings appear to be fairly similar. There
% are some apparent discrepencies that arise due to CIPRES use of actual
% branch lengths while matlab appears to normalize branch length in favor
% of more clear sub groupings. However, when looking closely, the shorted
% branch lengths on the sub nodes near the tree tips correspond with more
% closely related seqeunces in the matlab tree. So overall, there is good
% agreement between the two methods. 