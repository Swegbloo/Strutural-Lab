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
% for i = 1:length(x1)
%     idx = (x_ship >= x1(i)) & (x_ship <= x2(i));
%     load_distribution(idx) = l1(i) + (l2(i)-l1(i)) .* (x_ship(idx)-x1(i)) / (x2(i)-x1(i));
% end

for j = 1:length(x1)
    for i = 1:1000
        if (x_ship(i) >= x1(j)) && (x_ship(i) <= x2(j))
            load_distribution(i) = load_distribution(i) + l1(j) + (l2(j)-l1(j)) * (x_ship(i)-x1(j)) / (x2(j)-x1(j));
        end
    end
end
%%
% Calculate total weight and LCG
W_total = trapz(x_ship, load_distribution); % Total weight
LCG = trapz(x_ship, load_distribution .* x_ship) / W_total; % Center of gravity
%%
% Initial mean draft assumption by equating total weight to buoyant force
z_i = ones(ns,1) * 5; % Initial guess for z_i

% Iteratively solve for mean draft
for iter = 1:100
    [~, ~, ~, Disp,~,x_ship_interp,~] = hydrostat_properties(z_i);
    if abs((W_total - Disp) / W_total) < 1e-3
        break;
    end
    z_i = z_i * (W_total / Disp)^(1/3);
end
zz=ones(ns,1);
% %%
% % Iteration parameters
% alpha = 0.1; % Learning rate
% max_iterations = 20;
% tolerance_w = 0.01; % Weight tolerance
% tolerance_lcg = 0.01; % LCG-LCB tolerance
% 
% for iter = 1:max_iterations
%     % Get hydrostatic properties at current draft
%     [VCB, LCB, LCF, Disp, I,~,AWP] = hydrostat_properties(z_i);
%     KB = VCB;
%     BM = I/Disp;
%     KG = 9.393; %vertical postion of CG
%     MCT = (KB + BM - KG)*W_total/(100*L); 
% 
%     % Compute moment and check convergence
%     Moment = Disp * (LCG - LCB);
%     if MCT == 0
%         warning('MCT is zero, causing division error.');
%         break;
%     end
%     % fprintf('Mean draft change: %.3f\n', (Disp-W_total)/AWP);
%     % fprintf('trim angle: %.3f\n', 180/pi*atan((z_i(end)-z_i(1))/L));
%     % fprintf('Final Aft Trim: %.3f\n', z_i(1));
%     % fprintf('Final Bow Trim: %.3f\n', z_i(end));
% disp(z_i);
%     if isnan(VCB) || isnan(LCB) || isnan(LCF) || isnan(Disp) || isnan(MCT)
%         Disp(z_i);
%         break;      
%     end
%     if abs((W_total - Disp) / W_total) < tolerance_w && abs((LCG - LCB) / LCG) < tolerance_lcg
%         disp("sx");
%         break;
%     end
%     % Compute change in trim using MCT
%     d_trim = Moment / MCT;
% 
%     % Normalize so LCB remains on the same side of LCG
%     if LCB > LCG
%         d_trim = -abs(d_trim);
%     else
%         d_trim = abs(d_trim);
%     end
% 
%     % Ensure array compatibility for z_i update
%     %x_ship_interp = linspace(min(x_ship), max(x_ship), ns)';
%     z_i = z_i + d_trim * (x_ship_interp - LCF - L*(Disp-W_total)/(AWP*d_trim)) / L;
% 
%     % Adjust displacement and draft
%     %z_i = z_i * (W_total / Disp)^(1/3);
% 
% 
% 
% end
% 
% % % Final trim values
% % z_aft = z_i(1);
% % z_bow = z_i(end);
% % 
% % fprintf('Final Aft Trim: %.3f\n', z_aft);
% % fprintf('Final Bow Trim: %.3f\n', z_bow);

%%
% Iteration parameters
trim = 0.01; % Learning rate
max_iterations = 100000;
tolerance_w = 0.01; % Weight tolerance
tolerance_lcg = 0.01; % LCG-LCB tolerance

for iter = 1:max_iterations
    % Get hydrostatic properties at current draft
    [VCB, LCB, LCF, Disp, I,~,~] = hydrostat_properties(z_i);
    KB = VCB;
    BM = I/Disp;
    KG = 9.393; %vertical postion of CG
    
    % Normalize so LCB remains on the same side of LCG
    if LCB > LCG
        trim = -abs(trim);
    else
        trim = abs(trim);
    end
    z_i = z_i + trim * (x_ship_interp - LCF) / L;

    % fprintf('Mean draft change: %.3f\n', (Disp-W_total)/AWP);
    % fprintf('trim angle: %.3f\n', 180/pi*atan((z_i(end)-z_i(1))/L));
    % fprintf('Final Aft Trim: %.3f\n', z_i(1));
    % fprintf('Final Bow Trim: %.3f\n', z_i(end));

    % if isnan(VCB) || isnan(LCB) || isnan(LCF) || isnan(Disp) || isnan(MCT)
    %     Disp(z_i);
    %     break;      
    % end

    
    LCB_old = LCB;

    [VCB, LCB, LCF, Disp, I,~,AWP] = hydrostat_properties(z_i);

    % Ensure array compatibility for z_i update
    %x_ship_interp = linspace(min(x_ship), max(x_ship), ns)';

    z_i = z_i - ones(ns,1)*(Disp-W_total)/AWP;


    % Adjust displacement and draft
    %z_i = z_i * (W_total / Disp)^(1/3);
    trim = trim*(LCB-LCG)/(LCB_old-LCG);

    % if abs((W_total - Disp) / W_total) < tolerance_w && abs((LCG - LCB) / LCG) < tolerance_lcg
    %     disp("sx");
    %     break;
    % end
end

% % Final trim values
% z_aft = z_i(1);
% z_bow = z_i(end);
% 
% fprintf('Final Aft Trim: %.3f\n', z_aft);
% fprintf('Final Bow Trim: %.3f\n', z_bow);

