% Define the coordinates of nodes
nodes = [0, 0; 1, 0; 0.5, sqrt(3)/2];

% Define the connectivity of elements (triangles)
elements = [1, 2, 3];

% Plot the nodes
scatter(nodes(:, 1), nodes(:, 2), 'filled');
hold on;

% Plot the elements (connectivity)
for i = 1:size(elements, 1)
    element_nodes = nodes(elements(i, :), :);
    element_nodes = [element_nodes; element_nodes(1, :)]; % Closing the loop
    plot(element_nodes(:, 1), element_nodes(:, 2), 'b');
end

xlabel('X');
ylabel('Y');
title('Triangular Mesh');
axis equal;
grid on;
