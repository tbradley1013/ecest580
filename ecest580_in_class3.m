%% ECES T580 In class Lecture 3 - Tyler Bradley

% 1. Create a function to compute fibonacci sequence
fib(13)


% 2. Calculate levensthein distance between two strings
lev_dist("kitten", "sitting")


function output = fib(n)
  out_vec = repelem(0, n);
  
  next = 1;
  
  out_vec(1) = 0;
  for i = 2:n
      out_vec(i) = next;
      next = out_vec(i-1) + next;
  end
  output = out_vec;

end

function output = lev_dist(a, b)
len_a = length(a);
len_b = length(b);
my_mat = zeros(0, len_b+1, len_a+1);

for i = 2:len_b+1
    for j = 2:len_a+1
      if a(j-1) == b(i-1)
          cost = 0;
      else
          cost = 1;
      end   
      hor = my_mat(i-1, j) + 1;
      ver = my_mat(i, j-1) + 1;
      diag = my_mat(i-1, j-1) + cost;
          
      my_mat(i, j) = min([hor, ver, diag]);
    end   
    
end
output = my_mat;
end