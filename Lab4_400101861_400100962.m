%% section 1.1

EEG_ERP = load('ERP_EEG.mat');
EEG = EEG_ERP.ERP_EEG;
EEG = EEG.';

N_values = 100:100:2500;
fs = 240;
t = 0:1/fs:(size(EEG, 2)-1)/fs;

average_ERPs = zeros(length(N_values), length(t));

figure;
hold on;
for i = 1:length(N_values)
    N = N_values(i);
    
    avg_erp = mean(EEG(1:N, :), 1);
    average_ERPs(i, :) = avg_erp;
    
    plot(t, avg_erp, 'DisplayName', sprintf('N = %d', N));
end

xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Average ERP response for different values of N');
legend;
hold off;

%% section 1.2

N_values = 1:2550;
max_amplitudes = zeros(size(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    
    avg_erp = mean(EEG(1:N, :), 1);
    
    max_amplitudes(i) = max(abs(avg_erp));
end

figure;
plot(N_values, max_amplitudes, 'LineWidth', 1.5);
xlabel('Number of Trials (N)');
ylabel('Maximum Absolute Amplitude (\muV)');
title('Maximum Absolute Amplitude of Averaged ERP vs. Number of Trials (N)');
grid on;

%% section 1.3

N_values = 1:2550;
rms_errors = zeros(size(N_values));

previous_avg_erp = zeros(1, size(EEG, 2));
for i = 1:length(N_values)
    N = N_values(i);
    
    avg_erp = mean(EEG(1:N, :), 1);
    
    if i > 1
        rms_errors(i) = sqrt(mean((avg_erp - previous_avg_erp).^2));
    end
    
    previous_avg_erp = avg_erp;
end

figure;
plot(N_values, rms_errors, 'LineWidth', 1.5);
xlabel('Number of Trials (N)');
ylabel('RMS Error');
title('RMS Error between consecutive averaged ERPs vs. Number of Trials (N)');
grid on;

%% section 1.5

fs = 240;
t = 0:1/fs:(size(EEG, 2)-1)/fs;

N_0 = 500;
N_values = [2550, N_0, round(N_0/3)];

average_ERPs = zeros(length(N_values) + 2, length(t));

average_ERPs(1, :) = mean(EEG(1:2550, :), 1);

average_ERPs(2, :) = mean(EEG(1:N_0, :), 1);

average_ERPs(3, :) = mean(EEG(1:round(N_0/3), :), 1);

rng(1);
random_indices_N0 = randperm(2550, N_0);
average_ERPs(4, :) = mean(EEG(random_indices_N0, :), 1);

random_indices_N0_div3 = randperm(2550, round(N_0/3));
average_ERPs(5, :) = mean(EEG(random_indices_N0_div3, :), 1);

figure;
plot(t, average_ERPs(1, :), 'DisplayName', 'N = 2550');
hold on;
plot(t, average_ERPs(2, :), 'DisplayName', sprintf('N = %d', N_0));
plot(t, average_ERPs(3, :), 'DisplayName', sprintf('N = %d (N_0 / 3)', round(N_0/3)));
plot(t, average_ERPs(4, :), 'DisplayName', sprintf('Random N = %d', N_0));
plot(t, average_ERPs(5, :), 'DisplayName', sprintf('Random N = %d (N_0 / 3)', round(N_0/3)));
hold off;

xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Average ERP Responses for Different Numbers of Trials');
legend;
grid on;

%% section 2.1

SSVEP_data = load('SSVEP_EEG.mat');
SSVEP_Signal = SSVEP_data.SSVEP_Signal;
fs = 250;
channels = {'O1', 'O2', 'P8', 'P7', 'Oz', 'Pz'};

low_cutoff = 1;
high_cutoff = 40;
[b, a] = butter(4, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

figure('Name', 'Time Domain - Original and Filtered Signals');
for i = 1:size(SSVEP_Signal, 1)
    original_signal = SSVEP_Signal(i, :);
    
    filtered_signal = filtfilt(b, a, original_signal);
    
    subplot(size(SSVEP_Signal, 1), 1, i);
    plot((1:length(original_signal)) / fs, original_signal, 'b', 'DisplayName', 'Original Signal');
    hold on;
    plot((1:length(filtered_signal)) / fs, filtered_signal, 'r', 'DisplayName', 'Filtered Signal (1-40 Hz)');
    hold off;
    title(['Channel ' channels{i} ' - Time Domain']);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    legend;
end
sgtitle('Original and Filtered EEG Signals in Time Domain for Each Channel');

figure('Name', 'Frequency Domain - Original and Filtered Signals (Zoomed)');
for i = 1:size(SSVEP_Signal, 1)
    original_signal = SSVEP_Signal(i, :);
    
    filtered_signal = filtfilt(b, a, original_signal);
    
    L = length(original_signal);
    f = (0:L-1) * (fs / L);
    
    original_fft = abs(fft(original_signal));
    filtered_fft = abs(fft(filtered_signal));
    
    subplot(size(SSVEP_Signal, 1), 1, i);
    plot(f, original_fft, 'b', 'DisplayName', 'Original Signal');
    hold on;
    plot(f, filtered_fft, 'r', 'DisplayName', 'Filtered Signal (1-40 Hz)');
    hold off;
    xlim([0 0.5]);
    title(['Channel ' channels{i} ' - Frequency Domain (0-10 Hz)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    legend;
end
sgtitle('Original and Filtered EEG Signals in Frequency Domain (0-10 Hz) for Each Channel');


%% section 2.2

SSVEP_data = load('SSVEP_EEG.mat');
SSVEP_Signal = SSVEP_data.SSVEP_Signal;
Event_Samples = SSVEP_data.Event_samples;
fs = 250;
window_size = 5 * fs;

num_trials = 15;
trials = cell(num_trials, 1);

for i = 1:num_trials
    start_sample = Event_Samples(i);
    end_sample = start_sample + window_size - 1;
    
    trials{i} = SSVEP_Signal(:, start_sample:end_sample);
end

for i = 1:num_trials
    fprintf('Size of trial %d: %s\n', i, mat2str(size(trials{i})));
end

%% section 2.3

fs = 250;
window_size = 5 * fs;
channels = {'O1', 'O2', 'P8', 'P7', 'Oz', 'Pz'};

num_trials = 15;
trials = cell(num_trials, 1);

for i = 1:num_trials
    start_sample = Event_Samples(i);
    end_sample = start_sample + window_size - 1;
    
    trials{i} = SSVEP_Signal(:, start_sample:end_sample);
end

figure('Name', 'Power Spectral Density (PSD) for Each Channel Across 15 Trials');

for trial = 1:num_trials
    subplot(3, 5, trial);
    
    hold on;
    for ch = 1:length(channels)
        [pxx, f] = pwelch(trials{trial}(ch, :), [], [], [], fs);
        
        plot(f, 10*log10(pxx), 'DisplayName', channels{ch});
    end
    hold off;
    
    title(['Trial ' num2str(trial)]);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    legend;
    xlim([0 150]);
end

sgtitle('Power Spectral Density (PSD) for Each Channel Across 15 Trials');

%% section 3.1
load('FiveClass_EEG.mat');
fs = 256;

bands = struct('Delta', [1 4], 'Theta', [4 8], 'Alpha', [8 13], 'Beta', [13 30]);
band_names = fieldnames(bands);
filtered_signals = struct();

for i = 1:numel(band_names)
    band_name = band_names{i};
    Fpass1 = bands.(band_name)(1);
    Fpass2 = bands.(band_name)(2);
    
    N = 4;
    Apass = 1;
    
    h = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, fs);
    Hd = design(h, 'cheby1');
    
    filtered_signals.(band_name) = zeros(size(X));
    for c = 1:30
        filtered_signals.(band_name)(:,c) = filter(Hd, X(:,c));
    end
end

figure;
subplot(5,1,1);
plot(0:1/fs:5-1/fs, X(1:5*fs, 1));
xlabel('time', 'Interpreter', 'latex');
title('Original Signal', 'Interpreter', 'latex');

subplot(5,1,2);
plot(0:1/fs:5-1/fs, filtered_signals.Delta(1:5*fs, 1));
xlabel('time', 'Interpreter', 'latex');
title('Delta', 'Interpreter', 'latex');

subplot(5,1,3);
plot(0:1/fs:5-1/fs, filtered_signals.Theta(1:5*fs, 1));
xlabel('time', 'Interpreter', 'latex');
title('Theta', 'Interpreter', 'latex');

subplot(5,1,4);
plot(0:1/fs:5-1/fs, filtered_signals.Alpha(1:5*fs, 1));
xlabel('time', 'Interpreter', 'latex');
title('Alpha', 'Interpreter', 'latex');

subplot(5,1,5);
plot(0:1/fs:5-1/fs, filtered_signals.Beta(1:5*fs, 1));
xlabel('time', 'Interpreter', 'latex');
title('Beta', 'Interpreter', 'latex');

%% section 3.2
num_samples = 10 * fs; 

Alpha_X = filtered_signals.Alpha;
Alpha_Trials = zeros(num_samples, size(Alpha_X, 2), length(trial));
for i = 1:length(trial)
    Alpha_Trials(:, :, i) = Alpha_X(trial(i) : trial(i) + num_samples - 1, :);
end

Beta_X = filtered_signals.Beta;
Beta_Trials = zeros(num_samples, size(Beta_X, 2), length(trial));
for i = 1:length(trial)
    Beta_Trials(:, :, i) = Beta_X(trial(i) : trial(i) + num_samples - 1, :);
end

Theta_X = filtered_signals.Theta;
Theta_Trials = zeros(num_samples, size(Theta_X, 2), length(trial));
for i = 1:length(trial)
    Theta_Trials(:, :, i) = Theta_X(trial(i) : trial(i) + num_samples - 1, :);
end

Delta_X = filtered_signals.Delta;
Delta_Trials = zeros(num_samples, size(Delta_X, 2), length(trial));
for i = 1:length(trial)
    Delta_Trials(:, :, i) = Delta_X(trial(i) : trial(i) + num_samples - 1, :);
end

trial_index = 1; 

figure;
for channel = 1:3
    subplot(3, 1, channel);
    plot((0:num_samples-1)/fs, Alpha_Trials(:, channel, trial_index));
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel(['Channel ' num2str(channel)], 'Interpreter', 'latex');
    title(['Alpha - Trial ' num2str(trial_index) ', Channel ' num2str(channel)], 'Interpreter', 'latex');
end
sgtitle('10-Second EEG Segments for Selected Channels - Alpha Band', 'Interpreter', 'latex');

figure;
for channel = 1:3
    subplot(3, 1, channel);
    plot((0:num_samples-1)/fs, Beta_Trials(:, channel, trial_index));
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel(['Channel ' num2str(channel)], 'Interpreter', 'latex');
    title(['Beta - Trial ' num2str(trial_index) ', Channel ' num2str(channel)], 'Interpreter', 'latex');
end
sgtitle('10-Second EEG Segments for Selected Channels - Beta Band', 'Interpreter', 'latex');

figure;
for channel = 1:3
    subplot(3, 1, channel);
    plot((0:num_samples-1)/fs, Theta_Trials(:, channel, trial_index));
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel(['Channel ' num2str(channel)], 'Interpreter', 'latex');
    title(['Theta - Trial ' num2str(trial_index) ', Channel ' num2str(channel)], 'Interpreter', 'latex');
end
sgtitle('10-Second EEG Segments for Selected Channels - Theta Band', 'Interpreter', 'latex');

figure;
for channel = 1:3
    subplot(3, 1, channel);
    plot((0:num_samples-1)/fs, Delta_Trials(:, channel, trial_index));
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel(['Channel ' num2str(channel)], 'Interpreter', 'latex');
    title(['Delta - Trial ' num2str(trial_index) ', Channel ' num2str(channel)], 'Interpreter', 'latex');
end
sgtitle('10-Second EEG Segments for Selected Channels - Delta Band', 'Interpreter', 'latex');
%% section 3.3 and 3.4
bands = {'Delta', 'Theta', 'Alpha', 'Beta'};
trial_data = struct('Delta', Delta_Trials, 'Theta', Theta_Trials, 'Alpha', Alpha_Trials, 'Beta', Beta_Trials);
averaged_data = struct();

for b = 1:numel(bands)
    band_name = bands{b};
    current_data = trial_data.(band_name);
    
    averaged_data.(band_name) = zeros(fs * 10, size(current_data, 2), 5);
    
    for i = 1:5
        label_index = find(y == i);
        averaged_data.(band_name)(:, :, i) = mean(current_data(:, :, label_index).^2, 3);
        
        fprintf('Band: %s, Class: %d\n', band_name, i);
        fprintf('Mean Power Values (First 5 Samples, Channel 1):\n');
        disp(averaged_data.(band_name)(1:5, 1, i)); % Display first 5 samples for Channel 1
    end
end

%% section 3.5

window_size = 200;
newWin = ones(1, window_size) / sqrt(window_size);

smoothed_data = struct();
for b = 1:numel(bands)
    band_name = bands{b};
    current_band_data = averaged_data.(band_name);
    
    smoothed_data.(band_name) = zeros(size(current_band_data));
    
    for i = 1:5 
        for ch = 1:size(current_band_data, 2) 
            smoothed_data.(band_name)(:, ch, i) = conv(current_band_data(:, ch, i), newWin, 'same');
        end
        
        fprintf('Smoothed Data - Band: %s, Class: %d\n', band_name, i);
        fprintf('Smoothed Power Values (First 5 Samples, Channel 1):\n');
        disp(smoothed_data.(band_name)(1:5, 1, i)); 
    end
end

%% section 3.6

bands = {'Delta', 'Theta', 'Alpha', 'Beta'};
selected_channel = 15; 
time_vector = (0:size(smoothed_data.Delta, 1)-1) / fs;
colors = lines(5); 

figure;

for b = 1:numel(bands)
    band_name = bands{b};
    subplot(2, 2, b); 
    hold on;
    
    for i = 1:5
        plot(time_vector, smoothed_data.(band_name)(:, selected_channel, i), ...
            'DisplayName', ['Class ' num2str(i)], 'Color', colors(i,:));
    end
    
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Power', 'Interpreter', 'latex');
    title([band_name ' Band - Channel CPz'], 'Interpreter', 'latex');
    
    
    legend('show', 'Location', 'best', 'Interpreter', 'latex');
    
    hold off;
end

sgtitle('Smoothed Power for Each Frequency Band - Channel CPz', 'Interpreter', 'latex');

%% section 3.7

