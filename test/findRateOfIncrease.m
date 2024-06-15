% Example usage
S = 100; % Example sum
a = 2;   % Example first term
n = 5;   % Example number of terms

r = findRateOfIncrease(S, a, n);
disp(['The rate of increase r is: ', num2str(r)]);

function r = findRateOfIncrease(S, a, n)
    % Define the function representing the sum of the geometric series
    sumGP = @(r) a * (r^n - 1) / (r - 1) - S;
    
    % Initial guess for r
    r_initial_guess = 1.1; % This can be adjusted based on the problem context
    
    % Use fsolve to find the rate of increase r
    options = optimoptions('fsolve', 'Display', 'off');
    r = fsolve(sumGP, r_initial_guess, options);
end


