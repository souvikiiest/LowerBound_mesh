% Constants
L = 40;  % Total extent of boundary (constant)
D = 10;  % Depth below fan mesh (constant)

% Initial values
Df = 0.5; % Initial guess for Df
B = 1;   % Initial guess for B

% Target Lf
target_Lf_multiplier = 2;
tolerance = 0.01; % Tolerance for how close Lf should be to target value

% Iterate to find appropriate Df and B
found = false;
while ~found
    % Calculate Lf
    Lf = Df * (L - B) / (Df + D) + B;
    
    % Check if Lf is within the desired tolerance
    if abs(Lf - target_Lf_multiplier * B) < tolerance
        found = true;
    else
        % Adjust Df and B for next iteration
        if Lf > target_Lf_multiplier * B
            Df = Df * 0.99;  % Decrease Df
        else
            B = B * 1.01;  % Increase B
        end
    end
end

% Display results
fprintf('Found values:\n');
fprintf('Df = %.4f\n', Df);
fprintf('B = %.4f\n', B);
fprintf('Lf = %.4f\n', Lf);
fprintf('Target Lf = %.4f\n', target_Lf_multiplier * B);
