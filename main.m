clear;
clc;

% Load the load distribution data
load_data = xlsread('load_dist.xlsx');

x1 = load_data(:,1); % Starting point
x2 = load_data(:,2); % Ending point
l1 = load_data(:,3); % Load at x1
l2 = load_data(:,4); % Load at x2

% Define ship length and discretize x
x_ship = linspace(min(x1), max(x2), 1000); 
load_distribution = zeros(size(x_ship));

% Compute load distribution
for i = 1:length(x1)
    idx = (x_ship >= x1(i)) & (x_ship <= x2(i));
    load_distribution(idx) = l1(i) + (l2(i)-l1(i)) .* (x_ship(idx)-x1(i)) / (x2(i)-x1(i));
end

% Calculate total weight and LCG
W_total = trapz(x_ship, load_distribution); % Total weight
LCG = trapz(x_ship, load_distribution .* x_ship) / W_total; % Center of gravity

% Initial mean draft assumption
z_i = ones(107,1) * 5; % Initial draft guess

% Iteration parameters
alpha = 0.1; % Learning rate
max_iterations = 100;
tolerance_w = 0.005; % Weight tolerance
tolerance_lcg = 0.005; % LCG-LCB tolerance

for iter = 1:max_iterations
    % Get hydrostatic properties at current draft
    [VCB, LCB, LCF, Disp] = hydrostat_properties(z_i);
    
    % Compute moment and check convergence
    Moment = Disp * (LCG - LCB);
    if abs((W_total - Disp) / W_total) < tolerance_w && abs((LCG - LCB) / LCG) < tolerance_lcg
        break;
    end
    
    % Update trim using alpha factor
    z_i = z_i + alpha * abs(LCB - LCG) * sign(LCB - LCG);
    
    % Adjust displacement and draft
    z_i = z_i * (W_total / Disp)^(1/3);
end

% Final trim values
z_aft = z_i(1);
z_bow = z_i(end);

fprintf('Final Aft Trim: %.3f\n', z_aft);
fprintf('Final Bow Trim: %.3f\n', z_bow);
