%% In-class Assignment 2 - Tyler Bradley

%% a) 
% Given x=[1, 2, 5, 4, 3, 4, 10, 2, 1] of N=9 length, make a 5-length 
% window (k=5), that moves by 1 point at a time (the window overlap is 4), 
% that will compute the average value of the window for x. 
% (output is: [3.0, 3.6, 5.2, 4.6, 4.0])

% define x (used in all three parts)
x=[1, 2, 5, 4, 3, 4, 10, 2, 1];

% define empty vector to assign elements to
moving_mean_a = [];

% write for loop to calculate the mean for each window
for i = 1:5
    moving_mean_a(i) = mean(x(i:i+4));
end

moving_mean_a

% moving_mean_a =
%
%    3.0000    3.6000    5.2000    4.6000    4.0000

%% b) 
% By exploiting the N-1 overlap of values, can you make an algorithm 
% that is faster than O(n*k)?

% create a vector that has a predefined length to speed up computation
window_mean_b = repelem(0, 5);

% calculate the sum of the first window
my_sum = sum(x(1:5))
for i = 1:5
    % add the mean of the current window to the output vector
    window_mean_b(i) = my_sum/5;
    
    % exclude the 5th loop because there are no more elements to add
    if i < 5
      % calculate the new sum by subtracting the first element and 
      % adding the next element to my_sum
      my_sum = my_sum - x(i) + x(i + 5);
    end
end

window_mean_b

%window_mean_b =
%
%    3.0000    3.6000    5.2000    4.6000    4.0000

%% c) 
% Can you make an algorithm that will efficiently compute the maximum of 
% each N-length window moving with N-1 overlap?

% create predefine vector
window_max_c = repelem(0, 5);

% calculate the starting maximum
current_max = max(x(1:5))

% This method compares the current max value with the newest element 
% to determine if it is greater than the current max
for i = 1:5
    window_max_c(i) = current_max;
    if i < 5 & x(i+5) > current_max
        current_max = x(i+5);
    end
end

window_max_c

%window_max_c =
%
%     5     5    10    10    10
