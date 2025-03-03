clear;
clc;
%%
% Load the load distribution data
load('load_data.mat');
L = 160.8;
x1 = load_data(:,1); % Starting point
x2 = load_data(:,2); % Ending point
l1 = load_data(:,3); % Load at x1
l2 = load_data(:,4); % Load at x2

% Load sections data
load('sections.mat');
ns = sections(1,1); % Number of sections

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
%%
% Initial mean draft assumption by equating total weight to buoyant force
z_i = ones(ns,1) * 5; % Initial guess for z_i

% Iteratively solve for mean draft
for iter = 1:100
    [~, ~, ~, Disp] = hydrostat_properties(z_i);
    if abs((W_total - Disp) / W_total) < 1e-3
        break;
    end
    z_i = z_i * (W_total / Disp)^(1/3);
end
%%
% Iteration parameters
alpha = 0.1; % Learning rate
max_iterations = 100;
tolerance_w = 0.005; % Weight tolerance
tolerance_lcg = 0.005; % LCG-LCB tolerance

for iter = 1:max_iterations
    % Get hydrostatic properties at current draft
    [VCB, LCB, LCF, Disp, I] = hydrostat_properties(z_i);
    KB = VCB;
    BM = I/Disp;
    KG = 5; %vertical postion of CG
    MCT = (KB + BM - KG)*W_total/L; 
    % Compute moment and check convergence
    Moment = Disp * (LCG - LCB);
    if MCT == 0
        warning('MCT is zero, causing division error.');
        break;
    end

    fprintf('Final Aft Trim: %.3f\n', z_i(1));
    fprintf('Final Bow Trim: %.3f\n', z_i(end));

    if isnan(VCB) || isnan(LCB) || isnan(LCF) || isnan(Disp) || isnan(MCT)
        Disp(z_i);
        error('Hydrostat properties returned NaN values.');        
    end
    if abs((W_total - Disp) / W_total) < tolerance_w && abs((LCG - LCB) / LCG) < tolerance_lcg
        break;
    end
    % Compute change in trim using MCT
    d_trim = Moment / MCT;
    
    % Normalize so LCB remains on the same side of LCG
    if LCB > LCG
        d_trim = -abs(d_trim);
    else
        d_trim = abs(d_trim);
    end
    % Ensure array compatibility for z_i update
    x_ship_interp = linspace(min(x_ship), max(x_ship), ns)';
    z_i = z_i + d_trim * (x_ship_interp - LCF) / L;

    % Adjust displacement and draft
    z_i = z_i * (W_total / Disp)^(1/3);
    iter


end

% % Final trim values
% z_aft = z_i(1);
% z_bow = z_i(end);
% 
% fprintf('Final Aft Trim: %.3f\n', z_aft);
% fprintf('Final Bow Trim: %.3f\n', z_bow);