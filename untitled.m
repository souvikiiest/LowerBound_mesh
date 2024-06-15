% Parameters
B = 1; % You can set this to any desired value
numDivisions = 13; % Number of divisions for the geometric progression
L = 40 * B;
H = 40 * B;
a = B;
r = 2; % Common ratio for the geometric progression

% Generate the geometric progression for the y-axis
y = a * r .^ (0:numDivisions-1);

% Generate x-coordinates
x = linspace(0, L, numDivisions+1);

% Plotting the mesh
figure;
hold on;

% Draw the horizontal lines
for i = 1:numDivisions
    plot([0, L], [y(i), y(i)], 'k');
end

% Draw the vertical lines
for i = 1:length(x)
    plot([x(i), x(i)], [0, max(y)], 'k');
end

% Draw the boundary lines in red
plot([0, L], [0, 0], 'r'); % Bottom boundary
plot([0, L], [H, H], 'r'); % Top boundary
plot([0, 0], [0, H], 'r'); % Left boundary
plot([L, L], [0, H], 'r'); % Right boundary

% Set labels and view
xlabel('X');
ylabel('Y');
title('2D Mesh with Geometric Progression in Y-axis');
grid on;
axis equal;
ylim([0, H]);

hold off;
