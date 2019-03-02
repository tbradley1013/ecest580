%% In Class 6 - Tyler Bradley

%% Part 1
A = [1/3, 2/3, 0; 1/3, 0, 2/3; 0, 2/3, 1/3];

% Prob of state = 1 when initial state is 3 after 2 steps
hmm_one(A, 3, 1, 2)
% ans = 0.2222

% Prob of state = 1 when initial state is 3 after 2 steps
hmm_one(A, 3, 1, 10)
% and = 0.2023

% Prob of state = 1 when initial state is 3 after 2 steps
hmm_one(A, 3, 1, 50)
% and = 0.2000

% Prob of state = 1 when initial state = 3 after 100 steps
hmm_one(A, 3, 1, 100)
% and = 0.2000

%% Part 2
%Define initial conditions and givens
A = [1/3, 2/3, 0; 1/3, 0, 2/3; 0, 2/3, 1/3];
B = [9/10, 1/10; 1/10, 9/10; 1/2, 1/2];
pi0 = [1, 0, 0];
o = [0, 1, 1];

probs_out = hmm_two(A, B, pi0, o);
probs_out
%probs_out =
%
%    0.3667
%    0.2144
%    0.0908

%% Functions
% Function to calculate probability for part 1
function prob = hmm_one(mat, init, final, steps)
    final_mat = mat^steps;
    
    prob = final_mat(init, final);
end


function probs = hmm_two(A_mat, B_mat, I, O)
    % Define length of input sequence
    O_len = length(O);
    % Create empty output vector that will hold the probability at each
    % position in O
    probs = zeros(O_len,1);
    %Define initial pi value for for loop
    pi = I;
    % loop over each element in O and calculate probability of getting to
    % that point with the given sequence
    for i = 1:O_len
        % get the column index (1 or 2) from o (0 or 1, respectively)
        b_i = O(i) + 1;
        % Multiply pi by the A matrix and then do element wise
        % mulitplication of the transposed B - B is transposed to get
        % dimensions correct in the matrix multiplication
        new_pi = pi*A_mat.*transpose(B_mat(:,b_i));
        % Sum the new pi value to get the probability for that element in O
        probs(i) = sum(new_pi);
        % Assign the new_pi value to the looping value of pi
        pi = new_pi;
    end
    
end
