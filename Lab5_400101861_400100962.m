%% Section 1
clc; clear; close all;

load('normal.mat');

Fs = 250; 
t = normal(:, 1); 
sig = normal(:, 2); 
sig2 = randn(1, 1000); 


t_total = length(sig) / Fs;

clean_range = 5*Fs+1 : 15*Fs; 
noisy_range = 250*Fs+1 : 260*Fs;

sig_clean = sig(clean_range);
sig_noisy = sig(noisy_range);

window = gausswin(128); 
noverlap = 64;
nfft = 128;

figure;

% Clean signal plot
subplot(4 , 1, 1);
plot(t(clean_range), sig_clean, 'Color', '#49d9f2', 'LineWidth', 1.5);
title('Clean Signal Segment', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude (V)', 'Interpreter', 'latex');

% Noisy signal plot
subplot(4 , 1 , 2);
plot(t(noisy_range), sig_noisy, 'Color', '#49d9f2', 'LineWidth', 1.5);
title('Noisy Signal Segment', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude (V)', 'Interpreter', 'latex');

% PSD of clean signal
subplot(4 , 1 , 3);
[pxx_clean, f_clean] = pwelch(sig_clean, window, noverlap, nfft, Fs);
plot(f_clean, db(pxx_clean), 'Color', '#49d9f2', 'LineWidth', 1.5);
title('Power Spectrum of Clean Signal', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Power/Frequency (dB/Hz)', 'Interpreter', 'latex');
xlim([0, 120]);

% PSD of noisy signal
subplot(4 , 1, 4);
[pxx_noisy, f_noisy] = pwelch(sig_noisy, window, noverlap, nfft, Fs);
plot(f_noisy, db(pxx_noisy), 'Color', '#49d9f2', 'LineWidth', 1.5);
title('Power Spectrum of Noisy Signal', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Power/Frequency (dB/Hz)', 'Interpreter', 'latex');
xlim([0, 120]);

%% section 1.2

signal = normal(:, 2);
n = length(signal);
freq_axis = (0:n/2-1) * (Fs / n);

signal_fft = abs(fft(signal));
signal_half = signal_fft(1:n/2);

low_idx = round(1 * n / Fs);
high_idx = round(90 * n / Fs);

signal_band_pwr = sum(signal_half(low_idx:high_idx).^2);

filter_order = 3;
cutoff = [1, 90] / (Fs / 2);
[b, a] = butter(filter_order, cutoff, 'stop');

signal_filtered = filter(b, a, signal);

filtered_fft = abs(fft(signal_filtered));
filtered_half = filtered_fft(1:n/2);

filtered_band_pwr = sum(filtered_half(low_idx:high_idx).^2);

figure;
freqz(b, a, 512, Fs);
title('Frequency Response of Bandstop Filter', 'Interpreter', 'latex');
grid on;

figure;
impz(b, a, 50);
title('Impulse Response of Bandstop Filter', 'Interpreter', 'latex');
grid on;
%% Section 1.3: Filter and Analyze Clean and Noisy Signals

% Filter the clean and noisy segments
clean_filtered = filter(b, a, sig_clean);
noisy_filtered = filter(b, a, sig_noisy);

% Parameters for PSD calculation
window = gausswin(128); 
noverlap = 64; 
nfft = 128;

figure;

% PSD of filtered clean signal
subplot(2, 1, 1);
[pxx_clean_filtered, f_clean_filtered] = pwelch(clean_filtered, window, noverlap, nfft, Fs);
plot(f_clean_filtered, db(pxx_clean_filtered), 'Color', '#49d9f2', 'LineWidth', 1.5);
title('Filtered Clean Signal - Power Spectrum', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Power/Frequency (dB/Hz)', 'Interpreter', 'latex');
xlim([0, 120]);
grid on;

% PSD of filtered noisy signal
subplot(2, 1, 2);
[pxx_noisy_filtered, f_noisy_filtered] = pwelch(noisy_filtered, window, noverlap, nfft, Fs);
plot(f_noisy_filtered, db(pxx_noisy_filtered), 'Color', '#49d9f2', 'LineWidth', 1.5);
title('Filtered Noisy Signal - Power Spectrum', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Power/Frequency (dB/Hz)', 'Interpreter', 'latex');
xlim([0, 120]);
grid on;

%% Part 2 - Ventricular Arrhythmia Detection

clc;
clear;
close all;
fs = 250;

%% Section 1: Power Spectral Density (PSD) Analysis

dataPath = 'Lab 5_data\n_422.mat';
load(dataPath, 'n_422');

segments = struct(...
    'normal1', n_422(1:10710), ...
    'abnormal1', n_422(10711:11210), ...
    'normal2', n_422(11211:11442), ...
    'abnormal2', n_422(11442:59710), ...
    'abnormal3', n_422(61288:end));

windowLength = 100;
noverlap = 90;
nfft = 100;
window = hamming(windowLength);

[PSD_normal1, f_normal] = pwelch(segments.normal1(1:10*fs), window, noverlap, nfft, fs);
[PSD_abnormal2, ~] = pwelch(segments.abnormal2(1:10*fs), window, noverlap, nfft, fs);
[PSD_abnormal3, f_abnormal] = pwelch(segments.abnormal3(1:10*fs), window, noverlap, nfft, fs);

figure('Name', 'Power Spectral Density Analysis', 'NumberTitle', 'off');
subplot(3,1,1);
plot(f_normal, 20*log10(PSD_normal1), 'LineWidth', 1.5);
title('PSD of Normal Signal 1');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([0 fs/2]);

subplot(3,1,2);
plot(f_abnormal, 20*log10(PSD_abnormal3), 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
title('PSD of Abnormal Signal (VFIB)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([0 fs/2]);

subplot(3,1,3);
plot(f_normal, 20*log10(PSD_abnormal2), 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
title('PSD of Abnormal Signal (VT)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([0 fs/2]);

%% Section 2: Time-Domain Signal Visualization

timeVec = (0:(10*fs-1))/fs;

figure('Name', 'Time-Domain Signals', 'NumberTitle', 'off');
subplot(3,1,1);
plot(timeVec, segments.normal1(1:10*fs), 'LineWidth', 1.2);
title('Normal Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(timeVec, segments.abnormal3(1:10*fs), 'LineWidth', 1.2, 'Color', [0.4660, 0.6740, 0.1880]);
title('Abnormal Signal (VFIB)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(timeVec, segments.abnormal2(1:10*fs), 'LineWidth', 1.2, 'Color', [0.9290, 0.6940, 0.1250]);
title('Abnormal Signal (VT)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Section 3: Frame Segmentation and Labeling

frameDuration = 10;
overlapRatio = 0.5;

frameLength = frameDuration * fs;
frameStep = round(frameLength * (1 - overlapRatio));
ecgLength = length(n_422);
numFrames = floor((ecgLength - (frameLength - frameStep)) / frameStep);

labels_n422 = zeros(1, numFrames);
endpoints_n422 = zeros(1, numFrames);

for i = 1:numFrames
    startIdx = (i-1)*frameStep + 1;
    endIdx = (i-1)*frameStep + frameLength;
    
    if endIdx <= 10710
        labels_n422(i) = 1;
    elseif startIdx > 10710 && endIdx <= 11210
        labels_n422(i) = 3;
    elseif startIdx > 11210 && endIdx <= 11442
        labels_n422(i) = 1;
    elseif startIdx > 11442 && endIdx <= 59710
        labels_n422(i) = 3;
    elseif startIdx > 59710 && endIdx <= 61288
        labels_n422(i) = 4;
    elseif startIdx >= 61288
        labels_n422(i) = 2;
    else
        labels_n422(i) = 0;
    end
    endpoints_n422(i) = endIdx;
end

%% Section 4: Feature Extraction

features_n422 = zeros(4, numFrames);

for i = 1:numFrames
    frameData = n_422((i-1)*frameStep + 1 : (i-1)*frameStep + frameLength);
    features_n422(:,i) = [
        bandpower(frameData, fs, [0 60]);
        bandpower(frameData, fs, [60 120]);
        meanfreq(frameData, fs);
        medfreq(frameData, fs)
    ];
end

%% Section 5: Feature Histograms

isNormal = labels_n422 == 1;
isVFIB = labels_n422 == 2;

featuresNormal = features_n422(:, isNormal);
featuresVFIB = features_n422(:, isVFIB);

featureNames = {'Bandpower [0-60 Hz]', 'Bandpower [60-120 Hz]', 'Mean Frequency', 'Median Frequency'};

figure('Name', 'Feature Histograms', 'NumberTitle', 'off');
for i = 1:4
    subplot(2,2,i);
    histogram(featuresNormal(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    histogram(featuresVFIB(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold off;
    title(['Feature ', num2str(i), ': ', featureNames{i}]);
    xlabel('Value');
    ylabel('Probability');
    legend('Normal', 'Abnormal');
    grid on;
end

%% Section 6: Alarm Detection Using Bandpower and Median Frequency

[alarm_bandpower, t_bandpower] = va_detect_bandpower(n_422, fs);
[alarm_medfreq, t_medfreq] = va_detect_medfreq(n_422, fs);

%% Section 7: Confusion Matrix and Performance Metrics

indices = find(isNormal | isVFIB);
groundTruth = double(isVFIB(indices));
alarm_pred = alarm_bandpower(indices);

groundTruth_cat = categorical(groundTruth, [0 1], {'Normal', 'VFIB'});
alarm_pred_cat = categorical(alarm_pred, [0 1], {'Normal', 'VFIB'});

[cm_bandpower, order] = confusionmat(groundTruth_cat, alarm_pred_cat);

disp('Confusion Matrix for Bandpower Detection:');
disp(cm_bandpower);

accuracy_bp = (cm_bandpower(1,1) + cm_bandpower(2,2)) / sum(cm_bandpower(:));
sensitivity_bp = cm_bandpower(2,2) / (cm_bandpower(2,2) + cm_bandpower(1,2));
specificity_bp = cm_bandpower(1,1) / (cm_bandpower(1,1) + cm_bandpower(2,1));

fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n\n', ...
    accuracy_bp*100, sensitivity_bp*100, specificity_bp*100);

alarm_medfreq_pred = alarm_medfreq(indices);
alarm_medfreq_cat = categorical(alarm_medfreq_pred, [0 1], {'Normal', 'VFIB'});

[cm_medfreq, order_medfreq] = confusionmat(groundTruth_cat, alarm_medfreq_cat);

disp('Confusion Matrix for Median Frequency Detection:');
disp(cm_medfreq);

accuracy_med = (cm_medfreq(1,1) + cm_medfreq(2,2)) / sum(cm_medfreq(:));
sensitivity_med = cm_medfreq(2,2) / (cm_medfreq(2,2) + cm_medfreq(1,2));
specificity_med = cm_medfreq(1,1) / (cm_medfreq(1,1) + cm_medfreq(2,1));

fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n', ...
    accuracy_med*100, sensitivity_med*100, specificity_med*100);

%% Section 8: Additional Feature Extraction

additionalFeatures_n422 = zeros(6, numFrames);

for i = 1:numFrames
    frameData = n_422((i-1)*frameStep + 1 : (i-1)*frameStep + frameLength);
    peakData = findpeaks(frameData);
    additionalFeatures_n422(:,i) = [
        max(frameData);
        min(frameData);
        max(frameData) - min(frameData);
        mean(peakData);
        sum(frameData == 0);
        var(frameData)
    ];
end

%% Section 9: Histograms for Additional Features

featuresNormal_add = additionalFeatures_n422(:, isNormal);
featuresVFIB_add = additionalFeatures_n422(:, isVFIB);

additionalFeatureNames = {
    'Max Amplitude', 
    'Min Amplitude', 
    'Peak-to-Peak', 
    'Mean Peak', 
    'Zero Crossings', 
    'Variance'
};

figure('Name', 'Additional Feature Histograms', 'NumberTitle', 'off');
for i = 1:6
    subplot(2,3,i);
    histogram(featuresNormal_add(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    histogram(featuresVFIB_add(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold off;
    title(['Feature ', num2str(i), ': ', additionalFeatureNames{i}]);
    xlabel('Value');
    ylabel('Probability');
    legend('Normal', 'Abnormal');
    grid on;
end


%% Section 10: Alarm Detection Using Additional Features (Mean Peak and Max Amplitude)

[alarm_meanpeak, t_meanpeak] = va_detect_meanpeak(n_422, fs);
[alarm_maxamp, t_maxamp] = va_detect_maxamp(n_422, fs);

%% Section 11: Confusion Matrix and Performance Metrics for Additional Alarms

indices_jk = find(isNormal | isVFIB);
groundTruth_jk = double(isVFIB(indices_jk));

alarm_meanpeak_pred = double(alarm_meanpeak(indices_jk) > 0);
alarm_maxamp_pred = double(alarm_maxamp(indices_jk) > 0);

groundTruth_cat = categorical(groundTruth_jk, [0 1], {'Normal', 'VFIB'});
alarm_meanpeak_cat = categorical(alarm_meanpeak_pred, [0 1], {'Normal', 'VFIB'});
alarm_maxamp_cat = categorical(alarm_maxamp_pred, [0 1], {'Normal', 'VFIB'});

[cm_meanpeak, ~] = confusionmat(groundTruth_cat, alarm_meanpeak_cat, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Mean Peak Detection:');
disp(cm_meanpeak);

if size(cm_meanpeak,1) < 2 || size(cm_meanpeak,2) < 2
    cm_meanpeak = padarray(cm_meanpeak, [2 - size(cm_meanpeak,1) 2 - size(cm_meanpeak,2)], 0, 'post');
end

accuracy_meanpeak = (cm_meanpeak(1,1) + cm_meanpeak(2,2)) / sum(cm_meanpeak(:));
sensitivity_meanpeak = cm_meanpeak(2,2) / (cm_meanpeak(2,2) + cm_meanpeak(1,2));
specificity_meanpeak = cm_meanpeak(1,1) / (cm_meanpeak(1,1) + cm_meanpeak(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n\n', ...
    accuracy_meanpeak*100, sensitivity_meanpeak*100, specificity_meanpeak*100);

[cm_maxamp, ~] = confusionmat(groundTruth_cat, alarm_maxamp_cat, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Max Amplitude Detection:');
disp(cm_maxamp);

if size(cm_maxamp,1) < 2 || size(cm_maxamp,2) < 2
    cm_maxamp = padarray(cm_maxamp, [2 - size(cm_maxamp,1) 2 - size(cm_maxamp,2)], 0, 'post');
end

accuracy_maxamp = (cm_maxamp(1,1) + cm_maxamp(2,2)) / sum(cm_maxamp(:));
sensitivity_maxamp = cm_maxamp(2,2) / (cm_maxamp(2,2) + cm_maxamp(1,2));
specificity_maxamp = cm_maxamp(1,1) / (cm_maxamp(1,1) + cm_maxamp(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n', ... 
    accuracy_maxamp*100, sensitivity_maxamp*100, specificity_maxamp*100);



%% Section 12: Processing Second Dataset (n_424)

dataPath_n424 = 'Lab 5_data\n_424.mat';
load(dataPath_n424, 'n_424');

%% Section 12.1: PSD Analysis for n_424

segments_n424 = struct(...
    'normal1', n_424(1:10710), ...
    'abnormal1', n_424(10711:11210), ...
    'normal2', n_424(11211:11442), ...
    'abnormal2', n_424(11442:59710), ...
    'abnormal3', n_424(61288:end));

[PSD_normal1_n424, f_normal_n424] = pwelch(segments_n424.normal1(1:10*fs), window, noverlap, nfft, fs);
[PSD_abnormal2_n424, ~] = pwelch(segments_n424.abnormal2(1:10*fs), window, noverlap, nfft, fs);
[PSD_abnormal3_n424, f_abnormal_n424] = pwelch(segments_n424.abnormal3(1:10*fs), window, noverlap, nfft, fs);

figure('Name', 'PSD Analysis for n\_424', 'NumberTitle', 'off');
subplot(3,1,1);
plot(f_normal_n424, 20*log10(PSD_normal1_n424), 'LineWidth', 1.5);
title('PSD of Normal Signal 1 (n\_424)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([0 fs/2]);

subplot(3,1,2);
plot(f_abnormal_n424, 20*log10(PSD_abnormal3_n424), 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
title('PSD of Abnormal Signal (VFIB) (n\_424)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([0 fs/2]);

subplot(3,1,3);
plot(f_normal_n424, 20*log10(PSD_abnormal2_n424), 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
title('PSD of Abnormal Signal (VT) (n\_424)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
xlim([0 fs/2]);

%% Section 12.2: Time-Domain Signal Visualization for n_424

timeVec_n424 = (0:(10*fs-1))/fs;

figure('Name', 'Time-Domain Signals for n\_424', 'NumberTitle', 'off');
subplot(3,1,1);
plot(timeVec_n424, segments_n424.normal1(1:10*fs), 'LineWidth', 1.2);
title('Normal Signal (n\_424)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(timeVec_n424, segments_n424.abnormal3(1:10*fs), 'LineWidth', 1.2, 'Color', [0.4660, 0.6740, 0.1880]);
title('Abnormal Signal (VFIB) (n\_424)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(timeVec_n424, segments_n424.abnormal2(1:10*fs), 'LineWidth', 1.2, 'Color', [0.9290, 0.6940, 0.1250]);
title('Abnormal Signal (VT) (n\_424)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Section 12.3: Frame Segmentation and Labeling for n_424

frameDuration_n424 = 10;
overlapRatio_n424 = 0.5;

frameLength_n424 = frameDuration_n424 * fs;
frameStep_n424 = round(frameLength_n424 * (1 - overlapRatio_n424));
ecgLength_n424 = length(n_424);
numFrames_n424 = floor((ecgLength_n424 - (frameLength_n424 - frameStep_n424)) / frameStep_n424);

labels_n424 = zeros(1, numFrames_n424);
endpoints_n424 = zeros(1, numFrames_n424);

for i = 1:numFrames_n424
    startIdx_n424 = (i-1)*frameStep_n424 + 1;
    endIdx_n424 = (i-1)*frameStep_n424 + frameLength_n424;
    
    if endIdx_n424 <= 27249
        labels_n424(i) = 1;
    elseif startIdx_n424 > 27249 && endIdx_n424 <= 53673
        labels_n424(i) = 2;
    elseif startIdx_n424 > 53673 && endIdx_n424 <= 55134
        labels_n424(i) = 4;
    elseif startIdx_n424 > 55134 && endIdx_n424 <= 55134
        labels_n424(i) = 3;
    elseif startIdx_n424 >= 55134
        labels_n424(i) = 5;
    else
        labels_n424(i) = 0;
    end
    endpoints_n424(i) = endIdx_n424;
end

%% Section 12.4: Feature Extraction for n_424

features_n424 = zeros(4, numFrames_n424);

for i = 1:numFrames_n424
    frameData_n424 = n_424((i-1)*frameStep_n424 + 1 : (i-1)*frameStep_n424 + frameLength_n424);
    features_n424(:,i) = [
        bandpower(frameData_n424, fs, [0 40]);
        bandpower(frameData_n424, fs, [40 120]);
        meanfreq(frameData_n424, fs);
        medfreq(frameData_n424, fs)
    ];
end

%% Section 12.5: Feature Histograms for n_424

isNormal_n424 = labels_n424 == 1;
isVFIB_n424 = labels_n424 == 2;

featuresNormal_n424 = features_n424(:, isNormal_n424);
featuresVFIB_n424 = features_n424(:, isVFIB_n424);

featureNames_n424 = {'Bandpower [0-40 Hz]', 'Bandpower [40-120 Hz]', 'Mean Frequency', 'Median Frequency'};

figure('Name', 'Feature Histograms for n\_424', 'NumberTitle', 'off');
for i = 1:4
    subplot(2,2,i);
    histogram(featuresNormal_n424(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    histogram(featuresVFIB_n424(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold off;
    title(['Feature ', num2str(i), ': ', featureNames_n424{i}]);
    xlabel('Value');
    ylabel('Probability');
    legend('Normal', 'Abnormal');
    grid on;
end

%% Section 12.6: Alarm Detection Using Bandpower and Mean Frequency for n_424

[alarm_bandpower_n424, t_bandpower_n424] = va_detect_bandpower_424(n_424, fs);
[alarm_meanfreq_n424, t_meanfreq_n424] = va_detect_meanfreq_424(n_424, fs);

%% Section 12.7: Confusion Matrix and Performance Metrics for n_424 Alarms

indices_lg = find(isNormal_n424 | isVFIB_n424);
groundTruth_lg = double(isVFIB_n424(indices_lg));

alarm_bandpower_pred_n424 = double(alarm_bandpower_n424(indices_lg) > 0);
alarm_meanfreq_pred_n424 = double(alarm_meanfreq_n424(indices_lg) > 0);

groundTruth_cat_lg = categorical(groundTruth_lg, [0 1], {'Normal', 'VFIB'});
alarm_bandpower_cat_n424 = categorical(alarm_bandpower_pred_n424, [0 1], {'Normal', 'VFIB'});
alarm_meanfreq_cat_n424 = categorical(alarm_meanfreq_pred_n424, [0 1], {'Normal', 'VFIB'});

[cm_bp_n424, ~] = confusionmat(groundTruth_cat_lg, alarm_bandpower_cat_n424, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Bandpower Detection (n\_424):');
disp(cm_bp_n424);

if size(cm_bp_n424,1) < 2 || size(cm_bp_n424,2) < 2
    cm_bp_n424 = padarray(cm_bp_n424, [2 - size(cm_bp_n424,1), 2 - size(cm_bp_n424,2)], 0, 'post');
end

accuracy_bp_n424 = (cm_bp_n424(1,1) + cm_bp_n424(2,2)) / sum(cm_bp_n424(:));
sensitivity_bp_n424 = cm_bp_n424(2,2) / (cm_bp_n424(2,2) + cm_bp_n424(1,2));
specificity_bp_n424 = cm_bp_n424(1,1) / (cm_bp_n424(1,1) + cm_bp_n424(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n\n', ...
    accuracy_bp_n424*100, sensitivity_bp_n424*100, specificity_bp_n424*100);

[cm_mf_n424, ~] = confusionmat(groundTruth_cat_lg, alarm_meanfreq_cat_n424, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Mean Frequency Detection (n\_424):');
disp(cm_mf_n424);

if size(cm_mf_n424,1) < 2 || size(cm_mf_n424,2) < 2
    cm_mf_n424 = padarray(cm_mf_n424, [2 - size(cm_mf_n424,1), 2 - size(cm_mf_n424,2)], 0, 'post');
end

accuracy_mf_n424 = (cm_mf_n424(1,1) + cm_mf_n424(2,2)) / sum(cm_mf_n424(:));
sensitivity_mf_n424 = cm_mf_n424(2,2) / (cm_mf_n424(2,2) + cm_mf_n424(1,2));
specificity_mf_n424 = cm_mf_n424(1,1) / (cm_mf_n424(1,1) + cm_mf_n424(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n', ...
    accuracy_mf_n424*100, sensitivity_mf_n424*100, specificity_mf_n424*100);


%% Section 12.8: Additional Feature Extraction for n_424

additionalFeatures_n424 = zeros(6, numFrames_n424);

for i = 1:numFrames_n424
    frameData_n424 = n_424((i-1)*frameStep_n424 + 1 : (i-1)*frameStep_n424 + frameLength_n424);
    peakData_n424 = findpeaks(frameData_n424);
    additionalFeatures_n424(:,i) = [
        max(frameData_n424);
        min(frameData_n424);
        max(frameData_n424) - min(frameData_n424);
        mean(peakData_n424);
        sum(frameData_n424 == 0);
        var(frameData_n424)
    ];
end

%% Section 12.9: Histograms for Additional Features in n_424

featuresNormal_add_n424 = additionalFeatures_n424(:, isNormal_n424);
featuresVFIB_add_n424 = additionalFeatures_n424(:, isVFIB_n424);

additionalFeatureNames_n424 = {
    'Max Amplitude', 
    'Min Amplitude', 
    'Peak-to-Peak', 
    'Mean Peak', 
    'Zero Crossings', 
    'Variance'
};

figure('Name', 'Additional Feature Histograms for n\_424', 'NumberTitle', 'off');
for i = 1:6
    subplot(2,3,i);
    histogram(featuresNormal_add_n424(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    histogram(featuresVFIB_add_n424(i,:), 10, 'Normalization', 'probability', 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold off;
    title(['Feature ', num2str(i), ': ', additionalFeatureNames_n424{i}]);
    xlabel('Value');
    ylabel('Probability');
    legend('Normal', 'Abnormal');
    grid on;
end

%% Section 12.10: Alarm Detection Using Zero Crossings and Peak-to-Peak for n_424

[alarm_zero_n424, t_zero_n424] = va_detect_zero_424(n_424, fs);
[alarm_peaktopeak_n424, t_peaktopeak_n424] = va_detect_peaktopeak_424(n_424, fs);

%% Section 12.11: Confusion Matrix and Performance Metrics for Zero Crossings and Peak-to-Peak Alarms

indices_lk = find(isNormal_n424 | isVFIB_n424);
groundTruth_lk = double(isVFIB_n424(indices_lk));

alarm_zero_pred_n424 = double(alarm_zero_n424(indices_lk) > 0);
alarm_peaktopeak_pred_n424 = double(alarm_peaktopeak_n424(indices_lk) > 0);

disp('Unique values in groundTruth_lk:');
disp(unique(groundTruth_lk));

disp('Unique values in alarm_zero_pred_n424:');
disp(unique(alarm_zero_pred_n424));

disp('Unique values in alarm_peaktopeak_pred_n424:');
disp(unique(alarm_peaktopeak_pred_n424));

groundTruth_cat_lk = categorical(groundTruth_lk, [0 1], {'Normal', 'VFIB'});
alarm_zero_cat_n424 = categorical(alarm_zero_pred_n424, [0 1], {'Normal', 'VFIB'});
alarm_peaktopeak_cat_n424 = categorical(alarm_peaktopeak_pred_n424, [0 1], {'Normal', 'VFIB'});

[c_zero_n424, cm_zero_n424] = confusionmat(groundTruth_cat_lk, alarm_zero_cat_n424, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Zero Crossings Detection (n\_424):');
disp(cm_zero_n424);

if size(cm_zero_n424,1) < 2 || size(cm_zero_n424,2) < 2
    cm_zero_padded = zeros(2,2);
    cm_zero_padded(1:size(cm_zero_n424,1), 1:size(cm_zero_n424,2)) = cm_zero_n424;
    cm_zero_n424 = cm_zero_padded;
end

accuracy_zero_n424 = (cm_zero_n424(1,1) + cm_zero_n424(2,2)) / sum(cm_zero_n424(:));
sensitivity_zero_n424 = cm_zero_n424(2,2) / (cm_zero_n424(2,2) + cm_zero_n424(1,2));
specificity_zero_n424 = cm_zero_n424(1,1) / (cm_zero_n424(1,1) + cm_zero_n424(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n\n', ...
    accuracy_zero_n424*100, sensitivity_zero_n424*100, specificity_zero_n424*100);

[c_peaktopeak_n424, cm_peaktopeak_n424] = confusionmat(groundTruth_cat_lk, alarm_peaktopeak_cat_n424, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Peak-to-Peak Detection (n\_424):');
disp(cm_peaktopeak_n424);

if size(cm_peaktopeak_n424,1) < 2 || size(cm_peaktopeak_n424,2) < 2
    cm_peaktopeak_padded = zeros(2,2);
    cm_peaktopeak_padded(1:size(cm_peaktopeak_n424,1), 1:size(cm_peaktopeak_n424,2)) = cm_peaktopeak_n424;
    cm_peaktopeak_n424 = cm_peaktopeak_padded;
end

accuracy_peaktopeak_n424 = (cm_peaktopeak_n424(1,1) + cm_peaktopeak_n424(2,2)) / sum(cm_peaktopeak_n424(:));
sensitivity_peaktopeak_n424 = cm_peaktopeak_n424(2,2) / (cm_peaktopeak_n424(2,2) + cm_peaktopeak_n424(1,2));
specificity_peaktopeak_n424 = cm_peaktopeak_n424(1,1) / (cm_peaktopeak_n424(1,1) + cm_peaktopeak_n424(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n', ...
    accuracy_peaktopeak_n424*100, sensitivity_peaktopeak_n424*100, specificity_peaktopeak_n424*100);


%% Section 12.11: Confusion Matrix and Performance Metrics for Zero Crossings and Peak-to-Peak Alarms

indices_lk = find(isNormal_n424 | isVFIB_n424);
groundTruth_lk = double(isVFIB_n424(indices_lk));

alarm_zero_pred_n424 = double(alarm_zero_n424(indices_lk) > 0);
alarm_peaktopeak_pred_n424 = double(alarm_peaktopeak_n424(indices_lk) > 0);

disp('Unique values in groundTruth_lk:');
disp(unique(groundTruth_lk));

disp('Unique values in alarm_zero_pred_n424:');
disp(unique(alarm_zero_pred_n424));

disp('Unique values in alarm_peaktopeak_pred_n424:');
disp(unique(alarm_peaktopeak_pred_n424));

groundTruth_cat_lk = categorical(groundTruth_lk, [0 1], {'Normal', 'VFIB'});
alarm_zero_cat_n424 = categorical(alarm_zero_pred_n424, [0 1], {'Normal', 'VFIB'});
alarm_peaktopeak_cat_n424 = categorical(alarm_peaktopeak_pred_n424, [0 1], {'Normal', 'VFIB'});

[c_zero_n424, cm_zero_n424] = confusionmat(groundTruth_cat_lk, alarm_zero_cat_n424, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Zero Crossings Detection (n\_424):');
disp(cm_zero_n424);

if size(cm_zero_n424,1) < 2 || size(cm_zero_n424,2) < 2
    cm_zero_padded = zeros(2,2);
    cm_zero_padded(1:size(cm_zero_n424,1), 1:size(cm_zero_n424,2)) = cm_zero_n424;
    cm_zero_n424 = cm_zero_padded;
end

accuracy_zero_n424 = (cm_zero_n424(1,1) + cm_zero_n424(2,2)) / sum(cm_zero_n424(:));
sensitivity_zero_n424 = cm_zero_n424(2,2) / (cm_zero_n424(2,2) + cm_zero_n424(1,2));
specificity_zero_n424 = cm_zero_n424(1,1) / (cm_zero_n424(1,1) + cm_zero_n424(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n\n', ...
    accuracy_zero_n424*100, sensitivity_zero_n424*100, specificity_zero_n424*100);

[c_peaktopeak_n424, cm_peaktopeak_n424] = confusionmat(groundTruth_cat_lk, alarm_peaktopeak_cat_n424, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Peak-to-Peak Detection (n\_424):');
disp(cm_peaktopeak_n424);

if size(cm_peaktopeak_n424,1) < 2 || size(cm_peaktopeak_n424,2) < 2
    cm_peaktopeak_padded = zeros(2,2);
    cm_peaktopeak_padded(1:size(cm_peaktopeak_n424,1), 1:size(cm_peaktopeak_n424,2)) = cm_peaktopeak_n424;
    cm_peaktopeak_n424 = cm_peaktopeak_padded;
end

accuracy_peaktopeak_n424 = (cm_peaktopeak_n424(1,1) + cm_peaktopeak_n424(2,2)) / sum(cm_peaktopeak_n424(:));
sensitivity_peaktopeak_n424 = cm_peaktopeak_n424(2,2) / (cm_peaktopeak_n424(2,2) + cm_peaktopeak_n424(1,2));
specificity_peaktopeak_n424 = cm_peaktopeak_n424(1,1) / (cm_peaktopeak_n424(1,1) + cm_peaktopeak_n424(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n', ...
    accuracy_peaktopeak_n424*100, sensitivity_peaktopeak_n424*100, specificity_peaktopeak_n424*100);


%% Section 13: Cross-Dataset Alarm Application

% === Applying Best Detector from Data1 (n_424) to Data2 (n_422) ===

[alarm_zero_from1_to2, t_zero_from1_to2] = va_detect_zero_424(n_422, fs);

indices_cross1 = find(labels_n424 == 1 | labels_n424 == 2);
groundTruth_cross1 = double(labels_n424(indices_cross1) == 2);

alarm_zero_pred_cross1 = double(alarm_zero_from1_to2(indices_cross1) > 0);

disp('--- Cross-Application: Data1 Detector on Data2 ---');
disp('Unique values in groundTruth_cross1:');
disp(unique(groundTruth_cross1));

disp('Unique values in alarm_zero_pred_cross1:');
disp(unique(alarm_zero_pred_cross1));

groundTruth_cat_cross1 = categorical(groundTruth_cross1, [0 1], {'Normal', 'VFIB'});
alarm_zero_cat_cross1 = categorical(alarm_zero_pred_cross1, [0 1], {'Normal', 'VFIB'});

[cm_zero_cross1, ~] = confusionmat(groundTruth_cat_cross1, alarm_zero_cat_cross1, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Zero Crossings Alarm Applied to n\_424:');
disp(cm_zero_cross1);

if size(cm_zero_cross1,1) < 2 || size(cm_zero_cross1,2) < 2
    cm_zero_padded_cross1 = zeros(2,2);
    cm_zero_padded_cross1(1:size(cm_zero_cross1,1), 1:size(cm_zero_cross1,2)) = cm_zero_cross1;
    cm_zero_cross1 = cm_zero_padded_cross1;
end

accuracy_zero_cross1 = (cm_zero_cross1(1,1) + cm_zero_cross1(2,2)) / sum(cm_zero_cross1(:));
sensitivity_zero_cross1 = cm_zero_cross1(2,2) / (cm_zero_cross1(2,2) + cm_zero_cross1(1,2));
specificity_zero_cross1 = cm_zero_cross1(1,1) / (cm_zero_cross1(1,1) + cm_zero_cross1(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n\n', ...
    accuracy_zero_cross1*100, sensitivity_zero_cross1*100, specificity_zero_cross1*100);

% === Applying Best Detector from Data2 (n_424) to Data1 (n_422) ===

[alarm_bandpower_from2_to1, t_bandpower_from2_to1] = va_detect_bandpower(n_424, fs);

indices_cross2 = find(labels_n422 == 1 | labels_n422 == 2);
groundTruth_cross2 = double(labels_n422(indices_cross2) == 2);

alarm_bandpower_pred_cross2 = double(alarm_bandpower_from2_to1(indices_cross2) > 0);

disp('--- Cross-Application: Data2 Detector on Data1 ---');
disp('Unique values in groundTruth_cross2:');
disp(unique(groundTruth_cross2));

disp('Unique values in alarm_bandpower_pred_cross2:');
disp(unique(alarm_bandpower_pred_cross2));

groundTruth_cat_cross2 = categorical(groundTruth_cross2, [0 1], {'Normal', 'VFIB'});
alarm_bandpower_cat_cross2 = categorical(alarm_bandpower_pred_cross2, [0 1], {'Normal', 'VFIB'});

[cm_bandpower_cross2, ~] = confusionmat(groundTruth_cat_cross2, alarm_bandpower_cat_cross2, 'Order', {'Normal', 'VFIB'});
disp('Confusion Matrix for Bandpower Alarm Applied to n\_422:');
disp(cm_bandpower_cross2);

if size(cm_bandpower_cross2,1) < 2 || size(cm_bandpower_cross2,2) < 2
    cm_bandpower_padded_cross2 = zeros(2,2);
    cm_bandpower_padded_cross2(1:size(cm_bandpower_cross2,1), 1:size(cm_bandpower_cross2,2)) = cm_bandpower_cross2;
    cm_bandpower_cross2 = cm_bandpower_padded_cross2;
end

accuracy_bandpower_cross2 = (cm_bandpower_cross2(1,1) + cm_bandpower_cross2(2,2)) / sum(cm_bandpower_cross2(:));
sensitivity_bandpower_cross2 = cm_bandpower_cross2(2,2) / (cm_bandpower_cross2(2,2) + cm_bandpower_cross2(1,2));
specificity_bandpower_cross2 = cm_bandpower_cross2(1,1) / (cm_bandpower_cross2(1,1) + cm_bandpower_cross2(2,1));
fprintf('Accuracy: %.2f%%\nSensitivity: %.2f%%\nSpecificity: %.2f%%\n', ...
    accuracy_bandpower_cross2*100, sensitivity_bandpower_cross2*100, specificity_bandpower_cross2*100);

%% Section 15: Best Detector Application and Performance Evaluation

load("Lab 5_data/n_426.mat");
fs = 250;

labels = zeros(59, 11);

for i = 1:59
    t_start = (i-1) * 5 * 250;
    t_end = ((i-1) * 5 + 10) * 250;
    labels(i, 2) = t_end;
    if t_start >= 1 && t_end < 26432
        labels(i, 1) = 1;
    elseif t_start >= 26432
        labels(i, 1) = 2;
    else
        labels(i, 1) = 0;
    end
end

d = designfilt('bandpassiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 10, 'HalfPowerFrequency2', 30, 'DesignMethod', 'butter', 'SampleRate', fs);
filtered_signal = filtfilt(d, n_426);

for i = 1:59
    if labels(i, 1) == 1 || labels(i, 1) == 2
        signal = n_426(labels(i, 2) - 10 * 250 + 1 : labels(i, 2));
        filtered_signal = filtfilt(d, signal);
        labels(i, 3) = meanfreq(signal);
        labels(i, 4) = medfreq(signal);
        labels(i, 5) = sum(filtered_signal.^2, 'all');
    end
end

normal_features = labels(2:20, :);
VFIB_features = labels(23:59, :);

for i = 1:59
    if labels(i, 1) == 1 || labels(i, 1) == 2
        signal = n_426(labels(i, 2) - 10 * 250 + 1 : labels(i, 2));
        labels(i, 6) = min(signal);         
        labels(i, 7) = max(signal);        
        labels(i, 8) = labels(i, 7) - labels(i, 6);
        labels(i, 9) = mean(findpeaks(signal), 'all');
        labels(i, 10) = sum(signal == 0);
        labels(i, 11) = var(signal);
    end
end

normal_features = labels(2:20, [6:11]);
VFIB_features = labels(23:59, [6:11]);

[alarm_zero_crossing, t_zero_crossing] = va_detect(n_426, fs, 'zero-crossing');

targets = labels((labels(:, 1) == 1 | labels(:, 1) == 2), 1) == 2;

outputs_zero = alarm_zero_crossing((labels(:, 1) == 1 | labels(:, 1) == 2)) == 1;

cm1 = confusionmat(targets, outputs_zero);

disp('Confusion Matrix for Normal and VFIB Detection (Zero-Crossing Detector):');
disp(cm1);

accuracy_zero = (cm1(1,1) + cm1(2,2)) / sum(cm1(:));
sensitivity_zero = cm1(2,2) / (cm1(2,2) + cm1(1,2));
specificity_zero = cm1(1,1) / (cm1(1,1) + cm1(2,1));

fprintf('Zero-Crossing Detector Performance:\n');
fprintf('Accuracy: %.2f%%\n', accuracy_zero * 100);
fprintf('Sensitivity: %.2f%%\n', sensitivity_zero * 100);
fprintf('Specificity: %.2f%%\n', specificity_zero * 100);

disp('False Positives:');
disp(cm1(1, 2));
disp('Missed Detections:');
disp(cm1(2, 1));

%% Function Definitions (Extended)


function [alarm, t] = va_detect_meanpeak(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        peaks = findpeaks(seg);
        if isempty(peaks)
            mean_peak = 0;
        else
            mean_peak = mean(peaks);
        end
        if mean_peak < 326
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_maxamp(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        max_amp = max(seg);
        if max_amp < 326
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_zero_424(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        zero_crossings = sum(seg == 0);
        if zero_crossings > 10
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_peaktopeak_424(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        peak_to_peak = max(seg) - min(seg);
        if peak_to_peak < 317
            alarm(i) = 1;
        end
    end
end

%% Additional Function: Confusion Matrix Display

function display_confusion_matrix(cm, titleStr)
    figure('Name', titleStr, 'NumberTitle', 'off');
    confusionchart(cm, {'Normal', 'Abnormal'});
    title(titleStr);
end


%% Function Definitions

function [alarm, t] = va_detect_bandpower(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        feature = bandpower(seg, Fs, [60 120]);
        if feature < 2.3
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_medfreq(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        feature = medfreq(seg, Fs);
        if feature > 3
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_zero_crossings(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        zero_crossings = sum(seg == 0);
        if zero_crossings > 10
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_peak_to_peak(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        peak_to_peak = max(seg) - min(seg);
        if peak_to_peak < 317
            alarm(i) = 1;
        end
    end
end

function [alarm, t] = va_detect_bandpower_424(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        feature = bandpower(seg, Fs, [0 40]);
        threshold = 3167.2323;
        
        if feature < threshold
            alarm(i) = 1;
        end
    end
end


function [alarm, t] = va_detect_meanfreq_424(ecg_data, Fs)
    frame_sec = 10;
    overlap = 0.5;
    frame_length = round(frame_sec * Fs);
    frame_step = round(frame_length * (1 - overlap));
    num_frames = floor((length(ecg_data) - (frame_length - frame_step)) / frame_step);
    alarm = zeros(num_frames, 1);
    t = ((0:num_frames-1)' * frame_step + frame_length) / Fs;
    
    for i = 1:num_frames
        seg = ecg_data((i-1)*frame_step + 1 : (i-1)*frame_step + frame_length);
        feature = meanfreq(seg, Fs);
        threshold = 0.7714; 
        
        if feature > threshold
            alarm(i) = 1;
        end
    end
end


