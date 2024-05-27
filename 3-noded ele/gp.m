clear all;

% Equilibrium
E1 = (1/12.5)*[5, 0, 0, -5, 0, 5, 0, 0, -5; 0, 0, 5, 0, 5, -5, 0, -5, 0];
E2 = (1/34)*[1.8, 0, 11.8, -6.8, 0, -6.8, 5, 0, -5; 0, 11.8, 1.8, 0, -6.8, -6.8, 0, -5, 5];
E3 = (1/46.58)*[6.8, 0, 6.8, -6.8, 0, 6.9, 0, 0, -13.7; 0, 6.8, 6.8, 0, 6.9, -6.8, 0, -13.7, 0];
B_E = repmat([0; 18], 3, 1);
A_E = zeros(6, 27);
A_E(1:2, 1:9) = E1;
A_E(3:4, 10:18) = E2;
A_E(5:6, 19:27) = E3;

% Discontinuity
D11 = [0, 1, 0; 0, 0.5, 1];
D22 = [0.5, 0.5, 1; 0.5, -0.5, 0];
A_D = zeros(16, 27);
A_D(1:2, 1:3) = D11;
A_D(1:2, 10:12) = -D11;
A_D(3:4, 7:9) = D11;
A_D(3:4, 13:15) = -D11;
A_D(5:6, 10:12) = D22;
A_D(5:6, 19:21) = -D22;
A_D(7:8, 16:18) = D22;
A_D(7:8, 22:24) = -D22;
B_D = zeros(16, 1);

% Boundary
B11 = [0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 1];
B22 = [0, 0, -1, 0, 0, 0; 0, 0, 0, 0, 0, -1];
B7 = [0, 1, 0; 0, 0, 1];
B9 = [0, 1, 0; 0, 0, 1];

A_B = zeros(8, 27);
A_B(1:2, 1:6) = B11;
A_B(3:4, 4:9) = B22;
A_B(5:6, 19:21) = B7;
A_B(7:8, 25:27) = B9;
B_B = zeros(8, 1);

% Yield Condition
% p=3 c=20 gamma=18
Y1 = [-0.5, 0.5, 1.732; -0.5, 0.5, -1.732; 1, -1, 0];
A_Y = zeros(27, 27);
A_Y = kron(eye(9), Y1);
B_Y = repmat(20, 27, 1);

% Objective Function
C_z = zeros(21, 1);
C_1 = 2.5 * [0; 1; 0; 0; 1; 0];
C = [C_1; C_z];

% Global matrix
G_A1 = A_Y;
% G_A2_dash = zeros(5, 27);
% G_B_dash = zeros(5, 1);
G_A2 = [A_E; A_B; A_D];
G_B = [B_E; B_B; B_D];

% Check if sizes match
disp('Sizes of matrices and vectors:')
disp(['A_Y: ', num2str(size(A_Y))])
disp(['B_Y: ', num2str(size(B_Y))])
disp(['G_A2: ', num2str(size(G_A2))])
disp(['G_B: ', num2str(size(G_B))])
disp(['f: ', num2str(size(C))])

% Solution
f = -C;


% Start with yield condition constraints
[solution, fval, exitflag, output] = linprog(f, A_E,B_E);



