%% section 1
load('ElecPosXYZ');
load('Interictal');

ModelParams.R = [8 8.5 9.2]; 
layers = [8 8.5 9.2]; 
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3]; 
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];
Resolution = 1;
[SourceLocations, LeadFieldMatrix] = ForwardModel_3shell(Resolution, ModelParams);
scatter3(SourceLocations(1,:), SourceLocations(2,:), SourceLocations(3,:));

%% section 2
electrode_x = cellfun(@(elec) elec.XYZ(1) * layers(3), ElecPos);
electrode_y = cellfun(@(elec) elec.XYZ(2) * layers(3), ElecPos);
electrode_z = cellfun(@(elec) elec.XYZ(3) * layers(3), ElecPos);
electrode_labels = cellfun(@(elec) elec.Name, ElecPos, 'UniformOutput', false);

scatter3(SourceLocations(1,:), SourceLocations(2,:), SourceLocations(3,:));
hold on
scatter3(electrode_x, electrode_y, electrode_z, 'filled');
text(electrode_x, electrode_y, electrode_z, electrode_labels);

%% section 3
% rand_index = int32(1316 * rand(1,1) + 1);

scatter3(SourceLocations(1, rand_index), SourceLocations(2, rand_index), SourceLocations(3, rand_index), 'g', 'filled');

line_x = SourceLocations(1, rand_index) + (0:0.01:1);
line_y = SourceLocations(2, rand_index) + (0:0.01:1);
line_z = SourceLocations(3, rand_index) + (0:0.01:1);

scatter3(SourceLocations(1,:), SourceLocations(2,:), SourceLocations(3,:));
hold on
scatter3(electrode_x, electrode_y, electrode_z, 'filled');
text(electrode_x, electrode_y, electrode_z, electrode_labels);
hold on
plot3(line_x, line_y, line_z, 'LineWidth', 2);

%% section 4
interictal_signal = Interictal(5,:);
dipole_direction = [SourceLocations(1,rand_index), SourceLocations(2,rand_index), SourceLocations(3,rand_index)]./norm([SourceLocations(1,rand_index), SourceLocations(2,rand_index), SourceLocations(3,rand_index)]);
M = LeadFieldMatrix(:, (3 * rand_index - 2:3 * rand_index)) * dipole_direction' * interictal_signal;

disp_eeg(M, max(abs(M(:))), [], electrode_labels, 'Electrode Potentials');
xlim("tight");

%% section 5
peak_times_all = cell(1, 21);
avg_potential = zeros(1, 21);

for i = 1:21
    [peak_values, peak_times] = findpeaks(M(i,:), 'MinPeakHeight', 0.3 * max(M(i,:)));
    peak_times_all{i} = peak_times;
end

for i = 1:21
    peak_times = peak_times_all{i};
    window = peak_times + (-3:1:3)'; 
    valid_window = window(window > 0 & window <= size(M, 2));
    avg_potential(i) = mean(M(i, valid_window), 'all');
end

Display_Potential_3D(ModelParams.R(3), avg_potential);

%% section 6
DipoleLocationEstimate = LeadFieldMatrix' * ((LeadFieldMatrix * LeadFieldMatrix' + 0.6 * eye(21)) \ M);

%% section 7
DipoleLocationAvg = max(DipoleLocationEstimate, [], 2);
DipoleLocationAvg_reshaped = reshape(DipoleLocationAvg, 3, []);
DipoleLocationNorm = vecnorm(DipoleLocationAvg_reshaped, 2, 1);

[NormMax, LocationIndex] = max(DipoleLocationNorm);
DipoleDirection = DipoleLocationAvg_reshaped(:, LocationIndex) / NormMax;

%% section 8
position_error = norm(SourceLocations(:, rand_index) - SourceLocations(:, LocationIndex));
direction_error = norm(dipole_direction - DipoleDirection);
disp("position error is:" + position_error)
disp("direction error is:" + direction_error)
%% section 9.1
[max_value, index_1] = max(vecnorm(SourceLocations, 2, 1));

interictal_signal = Interictal(5,:);
dipole_direction = SourceLocations(:, index_1) / norm(SourceLocations(:, index_1));

M = LeadFieldMatrix(:, (3 * index_1 - 2):(3 * index_1)) * (dipole_direction .* interictal_signal);

disp_eeg(M, max(abs(M(:))), [], electrode_labels, 'Electrode Potentials');
xlim("tight");

avg_potential = zeros(1, 21);
for i = 1:21
    [peak_values, peak_times] = findpeaks(M(i,:), 'MinPeakHeight', 0.3 * max(M(i,:)));
    window_indices = peak_times + (-3:3)';
    window_indices = window_indices(window_indices > 0 & window_indices <= size(M, 2));
    avg_potential(i) = mean(M(i, window_indices), 'all');
end

figure;
Display_Potential_3D(ModelParams.R(3), avg_potential);

DipoleLocationEstimate = LeadFieldMatrix' * ((LeadFieldMatrix * LeadFieldMatrix' + 0.6 * eye(21)) \ M);

DipoleLocationAvg = max(DipoleLocationEstimate, [], 2);
DipoleLocationAvg_reshaped = reshape(DipoleLocationAvg, 3, []);
DipoleLocationNorm = vecnorm(DipoleLocationAvg_reshaped, 2, 1);

[NormMax, LocationIndex] = max(DipoleLocationNorm);
DipoleDirection = DipoleLocationAvg_reshaped(:, LocationIndex) / NormMax;

position_error_1 = norm(SourceLocations(:, index_1) - SourceLocations(:, LocationIndex));
direction_error_1 = abs(norm(dipole_direction) - norm(DipoleDirection));

%% section 9.2
index_2 = 104;
interictal_signal = Interictal(5,:);
dipole_direction = SourceLocations(:, index_2) / norm(SourceLocations(:, index_2));
M = LeadFieldMatrix(:, (3 * index_2 - 2):(3 * index_2)) * (dipole_direction .* interictal_signal);

disp_eeg(M, max(abs(M(:))), [], electrode_labels, 'Electrode Potentials');
xlim("tight");

avg_potential = zeros(1, 21);
for i = 1:21
    [peak_values, peak_times] = findpeaks(M(i,:), 'MinPeakHeight', 0.3 * max(M(i,:)));
    window = peak_times + (-3:3)';
    window = window(window > 0 & window <= size(M, 2));
    avg_potential(i) = mean(M(i, window), 'all');
end

figure;
Display_Potential_3D(ModelParams.R(3), avg_potential);

DipoleLocationEstimate = LeadFieldMatrix' * ((LeadFieldMatrix * LeadFieldMatrix' + 0.6 * eye(21)) \ M);

DipoleLocationAvg = max(DipoleLocationEstimate, [], 2);
DipoleLocationAvg_reshaped = reshape(DipoleLocationAvg, 3, []);
DipoleLocationNorm = vecnorm(DipoleLocationAvg_reshaped, 2, 1);

[NormMax, LocationIndex] = max(DipoleLocationNorm);
DipoleDirection = DipoleLocationAvg_reshaped(:, LocationIndex) / NormMax;

position_error_2 = norm(SourceLocations(:, index_2) - SourceLocations(:, LocationIndex));
direction_error_2 = abs(norm(dipole_direction) - norm(DipoleDirection));

%% section 9.3
[min_value, index_3] = min(vecnorm(SourceLocations, 2, 1));

interictal_signal = Interictal(5,:);
dipole_direction = SourceLocations(:, index_3) / norm(SourceLocations(:, index_3));
M = LeadFieldMatrix(:, (3 * index_3 - 2):(3 * index_3)) * (dipole_direction .* interictal_signal);

disp_eeg(M, max(abs(M(:))), [], electrode_labels, 'Electrode Potentials');
xlim("tight");

avg_potential = zeros(1, 21);
for i = 1:21
    [peak_values, peak_times] = findpeaks(M(i,:), 'MinPeakHeight', 0.3 * max(M(i,:)));
    window = peak_times + (-3:3)';
    window = window(window > 0 & window <= size(M, 2));
    avg_potential(i) = mean(M(i, window), 'all');
end

figure;
Display_Potential_3D(ModelParams.R(3), avg_potential);

DipoleLocationEstimate = LeadFieldMatrix' * ((LeadFieldMatrix * LeadFieldMatrix' + 0.6 * eye(21)) \ M);

DipoleLocationAvg = max(DipoleLocationEstimate, [], 2);
DipoleLocationAvg_reshaped = reshape(DipoleLocationAvg, 3, []);
DipoleLocationNorm = vecnorm(DipoleLocationAvg_reshaped, 2, 1);

[NormMax, LocationIndex] = max(DipoleLocationNorm);
DipoleDirection = DipoleLocationAvg_reshaped(:, LocationIndex) / NormMax;

position_error_3 = norm(SourceLocations(:, index_3) - SourceLocations(:, LocationIndex));
direction_error_3 = abs(norm(dipole_direction) - norm(DipoleDirection));
%% section 10: Simulated Annealing for Dipole Location and Direction
costFunction = @(params) dipoleCostFunction(params, LeadFieldMatrix, Interictal, SourceLocations, rand_index);

initialGuess = [rand*max(SourceLocations(1,:)); rand*max(SourceLocations(2,:)); rand*max(SourceLocations(3,:))];

initialDirection = rand(3,1);
initialDirection = initialDirection / norm(initialDirection); % Normalize

initialParams = [initialGuess; initialDirection];

options = optimoptions('simulannealbnd', 'Display', 'iter', 'MaxIterations', 2000, 'FunctionTolerance', 1e-6);

[optimalParams, optimalCost] = simulannealbnd(costFunction, initialParams, [], [], options);

optimalLocation = optimalParams(1:3);
optimalDirection = optimalParams(4:6);

disp('Optimal Dipole Location:');
disp(optimalLocation);

disp('Optimal Dipole Direction:');
disp(optimalDirection);

% Calculate the position and direction errors
position_error_10 = norm(optimalLocation - SourceLocations(:, rand_index));
direction_error_10 = norm(optimalDirection - DipoleDirection);

disp("Position Error (Simulated Annealing): " + position_error_10);
disp("Direction Error (Simulated Annealing): " + direction_error_10);

%% Function for the cost function (dipoleCostFunction)
function cost = dipoleCostFunction(params, LeadFieldMatrix, Interictal, SourceLocations, rand_index)
    location = params(1:3);
    direction = params(4:6);
    
    direction = direction / norm(direction);
    
    source_columns = (double(rand_index)-1)*3 + (1:3);  
    
    lead_field = LeadFieldMatrix(:, source_columns); 
    M = lead_field * direction;  
    
    M = M(1:20 , :);
    

    cost = norm(M - Interictal); 
end


