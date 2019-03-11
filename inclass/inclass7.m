%% In Class 7 - Tyler Bradley
A = [0.5, 0.1, 0.4; 0.2, 0.7, 0.1; 0.5, 0.3, 0.2];
B = [0.5, 0.3, 0.2; 0.3, 0.2, 0.5; 0.1, 0.6, 0.3];

I = [1/3, 1/3, 1/3];
O = [1, 3, 2, 3, 1];

ans = viterbi(A,B,I,O)
% ans =
%
%     1     3     1     3     1


function output = viterbi(A, B, I, O)
  len_o = length(O);
  nrow_a = length(A(:,1));
  ncol_a = length(A(1,:));
  temp_like = transpose(I);
  hidden_states = [;];
  
  temp_like = temp_like.*B(:,O(1));
  
  idx = 1;
  
  for j = 2:len_o
      trace_back = [];
      
      for i = 1:nrow_a
         state_temp_like = temp_like.*A(:,i);
         [max_state, argmax_state] = max(state_temp_like);
         trace_back(i) = argmax_state;
         
         next_temp_like = max_state .* B(:, O(j));
      end
      hidden_states(idx, :) = trace_back;
      idx = idx+1;
      temp_like = next_temp_like;
  end

  idx = len_o;
  output = [];
  [~, out_argmax] = max(temp_like);
  
  output(idx) = out_argmax;
  for i = 1:(len_o-1)
      step = len_o-i;
      output(idx-1) = hidden_states(step, output(idx));
      idx = idx-1;
  end
end