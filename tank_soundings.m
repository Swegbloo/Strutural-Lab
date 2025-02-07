clc;

function [volume, lcg, tcg, vcg] = calculate_tank_properties(file_path, draft)
    % Reads section data and calculates volume, LCG, TCG, and VCG for a given draft.

    %% Step 1: Read and parse the .txt file
    data = read_section_data(file_path);

    for i = 1:length(data)
        section = data(i);
        points = section.points;
        flag = 0;
        for j = 1:length(points)
            if points(j,2)<draft
                flag = flag+1;
                disp(flag);
            end
        end
    end

    if flag == 0
      volume = 0;
      lcg = 0;
      tcg = 0;
      vcg = 0;
    else
    %% Step 2: Filter data based on draft
    filtered_data = filter_data_by_draft(data, draft);

    %% Step 3: Calculate volume and centers of gravity
    [volume, lcg, tcg, vcg] = compute_properties(filtered_data);

    %% Step 4: Display the results
    fprintf('Calculated Properties for Draft %.2f:\n', draft);
    fprintf('Volume: %.2f cubic units\n', volume);
    fprintf('LCG: %.2f units\n', lcg);
    fprintf('TCG: %.2f units\n', tcg);
    fprintf('VCG: %.2f units\n', vcg);
    end
end



function data = read_section_data(file_path)
    % Reads the tank data from a .txt file.
    % Each section's data starts with an x-coordinate, followed by z-y pairs.
    
    fileID = fopen(file_path, 'r');
    raw_data = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    
    raw_data = raw_data{1};
    data = struct();
    
    idx = 1;
    for i = 1:length(raw_data)
        line = strtrim(raw_data{i});
        if isempty(line)
            continue;
        end
        
        nums = sscanf(line, '%f,');
        if length(nums) == 1 % New section starting with x-coordinate
            data(idx).x = nums(1);
            data(idx).points = [];
            idx = idx + 1;
        else % z-y coordinate pair
            data(idx-1).points = [data(idx-1).points; nums(1), nums(2)];
        end
    end
end


function filtered_data = filter_data_by_draft(data, draft)
    % Filters and interpolates section data up to the given draft.
    % section = data;
    % points = section.points;
    % for i = 1:length(data)
    %     section = data(i);
    %     points = section.points;
    %     flag = 0;
    %     for j = 1:length(points)
    %         if points(j)<draft
    %             flag = flag+1;
    %             disp(flag);
    %         end
    %     end
    % end
 
    filtered_data = struct();
    for i = 1:length(data)
        section = data(i);
        points = section.points;
        % Keep points below or at the draft level
        below_draft = points(points(:,2) <= draft, :);
        % Interpolate the waterline if necessary
        if any(points(:,2) > draft)
            above_draft = points(points(:,2) > draft, :);
            below_draft = [below_draft; interpolate_waterline(points, draft)];
        end
        filtered_data(i).x = section.x;
        filtered_data(i).points = below_draft;
    end
end

function interpolated_point = interpolate_waterline(points, draft)
    % Interpolates the intersection of the tank edge with the waterline (draft).
    for j = 2:size(points, 1)
        if points(j-1, 2) <= draft && points(j, 2) > draft
            z1 = points(j-1, 2); y1 = points(j-1, 1);
            z2 = points(j, 2); y2 = points(j, 1);
            
            % Linear interpolation
            t = (draft - z1) / (z2 - z1);
            interpolated_y = y1 + t * (y2 - y1);
            interpolated_point = [interpolated_y, draft];
            return;
        end
    end
end


function [volume, lcg, tcg, vcg] = compute_properties(data)
    % Calculates the tank properties (volume, LCG, TCG, VCG).
    total_volume = 0;
    total_lcg = 0;
    total_tcg = 0;
    total_vcg = 0;
    Xmom = [];
    Ymom = [];
    Zmom = [];
    
    sectional_data = compute_section_properties(data);

    %Calculate Xmom
    for i = 1:length([sectional_data.x])
        Xm = sectional_data(i).x*sectional_data(i).sectional_area;
            Xmom = [Xmom; Xm];
    end
   
    %Calculate Ymom
    for i = 1:length([sectional_data.x])
        Ym = sectional_data(i).centroid_tcg*sectional_data(i).sectional_area;
            Ymom = [Ymom; Ym];
    end

     %Calculate Zmom
    for i = 1:length([sectional_data.x])
        Zm = sectional_data(i).centroid_vcg*sectional_data(i).sectional_area;
            Zmom = [Zmom; Zm];
    end
    
    %Total_volume
    for i = 1:length([sectional_data.x])-1
        h = sectional_data(i+1).x- sectional_data(i).x;
        total_volume = total_volume + h*0.5*(sectional_data(i).sectional_area+sectional_data(i+1).sectional_area);
    end
  
    % LCG
    for i = 1:length([sectional_data.x])-1
        h = sectional_data(i+1).x- sectional_data(i).x;
        total_lcg = total_lcg + h*0.5*(Xmom(i)+Xmom(i+1));
    end

    % TCG
    for i = 1:length([sectional_data.x])-1
        h = sectional_data(i+1).x- sectional_data(i).x;
        total_tcg = total_tcg + h*0.5*(Ymom(i)+Ymom(i+1));
    end

     % VCG
    for i = 1:length([sectional_data.x])-1
        h = sectional_data(i+1).x- sectional_data(i).x;
        total_vcg = total_vcg + h*0.5*(Zmom(i)+Zmom(i+1));
    end
    
    

    % Final properties
    volume = -1*total_volume;
    lcg = total_lcg / total_volume;
    tcg = total_tcg / total_volume;
    vcg = total_vcg / total_volume;

end


function [Sectional_data] = compute_section_properties(filtered_data)
    % Initializes structure to store sectional properties for each X
    Sectional_data = struct();

    for i = 1:length(filtered_data)
        points1 = filtered_data(i).points;
        x_value = filtered_data(i).x;
        
        % Initialize variables
        sectional_area = 0;
        y_over_2 = [0];
        y_squared_2 = [0];
        y_times_z = [0];
        f_z_yz = 0;
        f_z_y2 = 0;
        
        % Calculate sectional area and related values for points
        for j = 1:length(points1)-1
            h = points1(j+1,2) - points1(j,2);
            sectional_area = sectional_area + h * 0.5 * (points1(j,1) + points1(j+1,1));
        end

        for j = 2:length(points1)
            % Y/2 values
            y_2 = points1(j,1) / 2;
            y_over_2 = [y_over_2; y_2];
            
            % Y^2/2 values
            y_squared_2_val = points1(j,1)^2 / 2;
            y_squared_2 = [y_squared_2; y_squared_2_val];
            
            % Y*Z values
            y_times_z_val = points1(j,1) * points1(j,2);
            y_times_z = [y_times_z; y_times_z_val];
        end

        % Compute f(Z,Y*Z)
        for j = 1:length(y_times_z)-1
            h = points1(j+1,2) - points1(j,2);
            f_z_yz = f_z_yz + h * 0.5 * (y_times_z(j) + y_times_z(j+1));
        end
        
        % Compute f(Z,Y^2/2)
        for j = 1:length(y_squared_2)-1
            h = points1(j+1,2) - points1(j,2);
            f_z_y2 = f_z_y2 + h * 0.5 * (y_squared_2(j) + y_squared_2(j+1));
        end
        % Compute area and centroids
        area = sectional_area;
        if sectional_area == 0
            centroid_tcg = 0;
            centroid_vcg = 0;
        else
            centroid_tcg = f_z_y2 / sectional_area; % TCG (breadth center)
            centroid_vcg = f_z_yz / sectional_area; % VCG (vertical center)
        end
        
        % Store data in structure
        Sectional_data(i).x = x_value;
        Sectional_data(i).sectional_area = area;
        Sectional_data(i).centroid_tcg = centroid_tcg;
        Sectional_data(i).centroid_vcg = centroid_vcg;
    end
end

% Define the range for the increment
increments = 0.1:0.1:13;

% Initialize empty arrays to store the results
volumes = [];
lcgs = [];
tcgs = [];
vcgs = [];

% Loop through the increments
for i = 1:length(increments)
    % Call the function and store the results
    [volume, lcg, tcg, vcg] = calculate_tank_properties('F.O.T.S.txt', increments(i));
    
    % Append the results to the arrays
    volumes = [volumes; volume];
    lcgs = [lcgs; lcg];
    tcgs = [tcgs; tcg];
    vcgs = [vcgs; vcg];
end

% Create a table for the data
T = table(increments', volumes, lcgs, tcgs, vcgs, 'VariableNames', {'Increment', 'Volume', 'LCG', 'TCG', 'VCG'});

% Write the table to an Excel file
writetable(T, 'tank_properties.xlsx');