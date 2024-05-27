% Lower Bound Limit Analysis for Simple Footing (Axisymmetric) with Finite Element Mesh
% Sloan (1988) Method

% Define material properties and footing geometry
phi = deg2rad(30); % Friction angle of soil (in radians)
c = 0; % Cohesion of soil
gamma = 18; % Unit weight of soil (kN/m^3)
q = 100; % Applied load (kN)
R = 3; % Radius of footing (m)

% Define mesh parameters
num_rings = 5; % Number of concentric rings
num_points_per_ring = 10; % Number of points per ring

% Initialize variables
node_coords = []; % Node coordinates
element_connectivity = []; % Element connectivity
num_nodes = 0;

% Generate mesh
for i = 1:num_rings
    radius = R * (i - 1) / (num_rings - 1);
    theta = linspace(0, pi/2, num_points_per_ring);
    for j = 1:num_points_per_ring
        x = radius * cos(theta(j));
        y = radius * sin(theta(j));
        node_coords = [node_coords; x, y];
        num_nodes = num_nodes + 1;
    end
end

% Generate element connectivity
for i = 1:num_rings - 1
    for j = 1:num_points_per_ring - 1
        node1 = (i - 1) * num_points_per_ring + j;
        node2 = node1 + 1;
        node3 = node1 + num_points_per_ring;
        node4 = node3 + 1;
        element_connectivity = [element_connectivity; node1, node2, node4, node3];
    end
end

% Initialize stresses and safety factor arrays
gamma_v = zeros(num_nodes, 1);
gamma_h = zeros(num_nodes, 1);
gamma_z = zeros(num_nodes, 1);
safety_factor = zeros(num_nodes, 1);

% Calculate stresses and safety factors for each node
for i = 1:num_nodes
    x = node_coords(i, 1);
    y = node_coords(i, 2);
    r = sqrt(x^2 + y^2);
    theta = atan2(y, x);
    
    % Calculate stresses at each node
    gamma_v(i) = gamma * r * cos(theta); % Vertical stress
    gamma_h(i) = gamma * r * sin(theta); % Horizontal stress
    gamma_z(i) = sqrt(gamma_v(i)^2 + gamma_h(i)^2); % Total stress
    
    % Calculate cohesion and friction components of yield surface
    phi_c = c * cos(phi) * r * cos(theta);
    phi_f = (q * R^2 / (2 * pi)) * (sin(phi) * cos(theta) + cos(phi) * sin(theta));
    
    % Calculate total yield surface
    phi_total = sqrt(phi_c^2 + phi_f^2);
    
    % Calculate safety factor
    safety_factor(i) = phi_total / gamma_z(i);
end

% Find minimum safety factor and corresponding node
[min_safety_factor, min_index] = min(safety_factor);
min_node_coords = node_coords(min_index, :);

% Output results
fprintf('Minimum Safety Factor: %.4f\n', min_safety_factor);
fprintf('Location of Minimum Safety Factor (x, y): %.2f m, %.2f m\n', min_node_coords(1), min_node_coords(2));

% Plot mesh with safety factor color-coded
figure;
trisurf(element_connectivity, node_coords(:,1), node_coords(:,2), zeros(num_nodes,1), safety_factor);
colorbar;
xlabel('X (m)');
ylabel('Y (m)');
title('Finite Element Mesh with Safety Factor');
