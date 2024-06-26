% Example usage
S = 10; % Example sum
a = 2;   % Example first term
n = 5;   % Example number of terms

r = findRateOfIncrease(S, a, n);
disp(['The rate of increase r is: ', num2str(r)]);

function r = findRateOfIncrease(S, a, n)
    % Define the function representing the sum of the geometric series
    sumGP = @(r) a * (r^n - 1) / (r - 1) - S;
    
    % Initial guess for r, ensuring it's not 1
    r_initial_guess = 1.1; % This can be adjusted based on the problem context
    
    % Define a function that avoids the value r = 1
    function F = sumGP_nonUnity(r)
        if abs(r - 1) ==0
            F = inf; % Arbitrarily large number to avoid r = 1
        else
            F = sumGP(r);
        end
    end
    
    % Use fsolve to find the rate of increase r
    options = optimoptions('fsolve', 'Display', 'off');
    r = fsolve(@sumGP_nonUnity, r_initial_guess, options);
    
    % Check if the solution is valid and not close to 1
    if abs(r - 1) < 1e-6
        error('The solver converged to a value close to 1, which is not allowed.');
    end
end


