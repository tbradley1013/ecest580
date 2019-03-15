%% In Class 8 - Tyler Bradley
clc;clear

% read in fasta files
limno = fastaread('Limnohabitans.fasta');
thermo = fastaread('s_thermotolerans.fasta');
vulcan = fastaread('t_vulcanus_rbcl.fasta');
unk = fastaread('uncultured.fasta');

% thermo should be highest

limno_probs = tri_prob(limno.Sequence);
thermo_probs = tri_prob(thermo.Sequence);
vulcan_probs = tri_prob(vulcan.Sequence);

probs_mat = [limno_probs, thermo_probs, vulcan_probs];
seq_names = ["limno", "thermo", "vulcan"];

% get the liklihood of the unknown sequence being related to 
% each of the reference sequences
likeli = zeros(3,1);
for i = 1:3
    n_probs = length(unk.Sequence)-3+1;
    probs_unk = zeros(n_probs,1);
    
    for j = i:n_probs
        tri = unk.Sequence(j:j+2);
        hash = hashtag(tri);
        probs_unk = probs_mat(hash,i);
    end
    
    likeli(i) = sum(probs_unk);
end

% likeli =

%   -6.203884135059216
%   -5.119323599846757
%   -5.227089857259669

% Get the argmax value of the liklihoods
[~,like_max] = max(likeli);

% get the sequence name of the most likely match
most_prob = seq_names(like_max)
% most_prob = 
% 
%    "thermo"


% Get the probability of each trimer for a given sequence
function output = tri_prob(seq)
  len_seq = length(seq)-3+1;
  tri_count = trimer_count(seq);
  len_tri = length(tri_count);
  output = zeros(len_tri, 1);
  n_feat = 4^3;
  for i = 1:len_tri
      output(i) = log((tri_count(i)+1)/(len_seq+n_feat));
  end


end


% get the count of the trimers in the sequence
function output = trimer_count(seq)
  ncounts = nmercount(seq,3);
  output = zeros(length(ncounts),1);
  for i = 1:length(ncounts)
      hash_num = hashtag(ncounts{i,1});
      output(hash_num) = ncounts{i,2};
  end

end

% get the numeric value of the trimer
function output = hashtag(seq)
  temp_vec = zeros(3, 1);
  
  idx = 1;
  for i = 1:length(seq)
      if seq(i) == 'A'
         temp_vec(idx) = 0;
      elseif seq(i) == "C"
         temp_vec(idx) = 1;
      elseif seq(i) == 'G'
          temp_vec(idx) = 2;
      else
          temp_vec(idx) = 3;
      end
      idx = idx+1;
  end
  
  output = temp_vec(1)*(4^2)+temp_vec(2)*(4)+temp_vec(3)+1;
end