%% section 1.1

load('E:\uni\term 7\signalProccessingLab\Lab1\Lab 1_data\EEG_sig.mat');

Fs = 256;
channel_5 = Z(5, :);
time = (0:length(channel_5)-1) / Fs;

figure;
plot(time, channel_5 , 'Color' , "#2596be");
xlabel('Time (s)', Interpreter='latex');
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');
title('EEG Signal - Channel 5' , Interpreter='latex');
grid on;

%% section 1.2
intervals = [0 15; 18 40; 45 50; 50 length(channel_5)/Fs]; 

figure;

for i = 1:size(intervals, 1)
    start_idx = round(intervals(i, 1) * Fs) + 1;
    end_idx = round(intervals(i, 2) * Fs);

    signal_part = channel_5(start_idx:end_idx);
    time_part = time(start_idx:end_idx);

    subplot(4, 1, i);
    plot(time_part, signal_part, 'Color', "#2596be"); 
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
    title(sprintf('EEG Signal - Channel 5 (%.1f s to %.1f s)', intervals(i, 1), intervals(i, 2)), 'Interpreter', 'latex');
    xlim("tight")
    grid on;

end

%% section 1.3
Fs = 256;
channel_8 = Z(8, :);
time = (0:length(channel_8)-1) / Fs;

figure;
plot(time, channel_8 , 'Color' , "#2596be");
xlabel('Time (s)', Interpreter='latex');
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');
title('EEG Signal - Channel 8' , Interpreter='latex');
grid on;

%% section 1.4
offset = max(max(abs(Z)))/3 ;
feq = 256 ;
ElecName = des.channelnames ;
disp_eeg(Z,offset,feq,ElecName) ;
xlim("tight")
%% section 1.5
%% section 1.6
Fs = 256;
channel_C3 = Z(5, :);
time = (0:length(channel_C3)-1) / Fs;

intervals = [2 7; 30 35; 42 47; 50 55];

tiledlayout(4, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

for i = 1:size(intervals, 1)
    start_idx = round(intervals(i, 1) * Fs) + 1;
    end_idx = round(intervals(i, 2) * Fs);
    
    signal_part = channel_C3(start_idx:end_idx);
    time_part = time(start_idx:end_idx);
    
    L = length(signal_part);
    Y = fft(signal_part);
    P2 = abs(Y / L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f = Fs * (0:(L/2)) / L;
    
    nexttile(2*i-1);
    plot(time_part, signal_part, 'Color', [37 150 190]/255);
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Amplitude ($\mu$V)', 'Interpreter', 'latex');
    title(sprintf('EEG Signal - Channel C3 (%.1f s to %.1f s)', intervals(i, 1), intervals(i, 2)), 'Interpreter', 'latex');
    grid on;
    
    nexttile(2*i);
    plot(f, P1, 'Color', [37 150 190]/255);
    xlabel('Frequency (Hz)', 'Interpreter', 'latex');
    ylabel('Amplitude', 'Interpreter', 'latex');
    title(sprintf('Frequency Spectrum - Channel C3 (%.1f s to %.1f s)', intervals(i, 1), intervals(i, 2)), 'Interpreter', 'latex');
    grid on;
    
    xlim([0 128]);
end
%% section 1.7
intervals = [2 7; 30 35; 42 47; 50 55];

window_length = 256;
overlap = 128;
nfft = 512; 

tiledlayout(4, 1, 'TileSpacing', 'loose', 'Padding', 'loose');

for i = 1:size(intervals, 1)
    start_idx = round(intervals(i, 1) * Fs) + 1;
    end_idx = round(intervals(i, 2) * Fs);
    
    signal_part = channel_C3(start_idx:end_idx);
    
    [pxx, f] = pwelch(signal_part, window_length, overlap, nfft, Fs);
    
    nexttile;
    plot(f, 10*log10(pxx), 'Color', [37 150 190]/255);
    xlabel('Frequency (Hz)', 'Interpreter', 'latex');
    ylabel('Power/Frequency (dB/Hz)', 'Interpreter', 'latex');
    title(sprintf('Power Spectral Density - Channel C3 (%.1f s to %.1f s)', intervals(i, 1), intervals(i, 2)), 'Interpreter', 'latex');
    grid on;
    
    xlim([0 128]);
end

%% section 1.8
intervals = [2 7; 30 35; 42 47; 50 55];
window_length = 128;
overlap = 64;
nfft = 128;

tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

for i = 1:size(intervals, 1)
    start_idx = round(intervals(i, 1) * Fs) + 1;
    end_idx = round(intervals(i, 2) * Fs);
    
    signal_part = channel_C3(start_idx:end_idx);
    
    nexttile;
    spectrogram(signal_part, hamming(window_length), overlap, nfft, Fs, 'yaxis');
    xlabel('Time (seconds)', 'Interpreter', 'latex');
    ylabel('Frequency (Hz)', 'Interpreter', 'latex');
    title(sprintf('Spectrogram - Channel C3 (%.1f s to %.1f s)', intervals(i, 1), intervals(i, 2)), 'Interpreter', 'latex');
    ylim([0 128]);
end

%% section 1.9
start_idx = round(30 * Fs) + 1;
end_idx = round(35 * Fs);
signal_part = channel_C3(start_idx:end_idx);

Fc = 64;
[b, a] = butter(4, Fc / (Fs / 2));
signal_filtered = filtfilt(b, a, signal_part);

Fs_new = 128;
signal_downsampled = downsample(signal_filtered, 2);

figure;
subplot(3,1,1);
plot((0:length(signal_part)-1)/Fs, signal_part);
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Original Signal', 'Interpreter', 'latex');

L = length(signal_downsampled);
Y = fft(signal_downsampled);
P2 = abs(Y / L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1);
f = Fs_new * (0:(L/2)) / L;

subplot(3,1,2);
plot(f, P1);
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
title('Frequency Spectrum (Simple FFT) - Downsampled', 'Interpreter', 'latex');
xlim([0 64]);

window_length = 128;
overlap = 64;
nfft = 128;
subplot(3,1,3);

spectrogram(signal_downsampled, hamming(window_length), overlap, nfft, Fs_new, 'yaxis');
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('Frequency (Hz)', 'Interpreter', 'latex');
title('Spectrogram (STFT) - Downsampled', 'Interpreter', 'latex');
ylim([0 64]);





%% section 2.1
ECG = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');
disp(ECG);
ecg_signal = ECG.Sig;
Fs = 256;
n = length(ecg_signal);
t = (0:n-1)/Fs;

figure;
plot(t, ecg_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal');

%% section 2.1 zoom
figure;
plot(t, ecg_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Signal');
xlim([0,50])

%% section 2.1 divided channels
ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');
ecg_channel1 = ECG_data.Sig(:,1);
ecg_channel2 = ECG_data.Sig(:,2);
Fs = ECG_data.sfreq;
n = length(ecg_channel1);
t = (0:n-1)/Fs;

figure;
subplot(2,1,1);
plot(t, ecg_channel1);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 1');
xlim([0 max(t)]);

subplot(2,1,2);
plot(t, ecg_channel2);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 2');
xlim([0 max(t)]);

zoom_time_window = [0 5];

figure;
subplot(2,1,1);
plot(t, ecg_channel1);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 1 (Zoomed)');
xlim(zoom_time_window);

subplot(2,1,2);
plot(t, ecg_channel2);
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 2 (Zoomed)');
xlim(zoom_time_window);


%% section 2.1 ++
ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');
ecg_channel1 = ECG_data.Sig(:,1);
ecg_channel2 = ECG_data.Sig(:,2);
Fs = ECG_data.sfreq;
n = length(ecg_channel1);
t = (0:n-1)/Fs;

[R_peaks1, locs1] = findpeaks(ecg_channel1, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);
[R_peaks2, locs2] = findpeaks(ecg_channel2, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);

P_waves1 = zeros(size(locs1)); Q_waves1 = zeros(size(locs1));
S_waves1 = zeros(size(locs1)); T_waves1 = zeros(size(locs1));
P_waves2 = zeros(size(locs2)); Q_waves2 = zeros(size(locs2));
S_waves2 = zeros(size(locs2)); T_waves2 = zeros(size(locs2));

P_offset = round(0.05 * Fs);
Q_offset = round(0.01 * Fs);
S_offset = round(0.01 * Fs);
T_offset = round(0.05 * Fs);

P_waves1 = locs1 - P_offset;
Q_waves1 = locs1 - Q_offset;
S_waves1 = locs1 + S_offset;
T_waves1 = locs1 + T_offset;

P_waves2 = locs2 - P_offset;
Q_waves2 = locs2 - Q_offset;
S_waves2 = locs2 + S_offset;
T_waves2 = locs2 + T_offset;

figure;
subplot(2,1,1);
plot(t, ecg_channel1);
hold on;
plot(t(locs1), ecg_channel1(locs1), 'ro', 'MarkerSize', 6);
plot(t(P_waves1), ecg_channel1(P_waves1), 'go', 'MarkerSize', 6);
plot(t(Q_waves1), ecg_channel1(Q_waves1), 'bo', 'MarkerSize', 6);
plot(t(S_waves1), ecg_channel1(S_waves1), 'mo', 'MarkerSize', 6);
plot(t(T_waves1), ecg_channel1(T_waves1), 'co', 'MarkerSize', 6);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 1');
xlim([0,2])
legend('ECG', 'R-peak', 'P', 'Q', 'S', 'T');

subplot(2,1,2);
plot(t, ecg_channel2);
hold on;
plot(t(locs2), ecg_channel2(locs2), 'ro', 'MarkerSize', 6);
plot(t(P_waves2), ecg_channel2(P_waves2), 'go', 'MarkerSize', 6);
plot(t(Q_waves2), ecg_channel2(Q_waves2), 'bo', 'MarkerSize', 6);
plot(t(S_waves2), ecg_channel2(S_waves2), 'mo', 'MarkerSize', 6);
plot(t(T_waves2), ecg_channel2(T_waves2), 'co', 'MarkerSize', 6);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 2');
xlim([0,2])
legend('ECG', 'R-peak', 'P', 'Q', 'S', 'T');


%% section 2.2
ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');
ecg_signal = ECG_data.Sig;
r_times = ECG_data.ATRTIMED;
annotations = ECG_data.ANNOTD;
Fs = 256;
n = length(ecg_signal);
t = (0:n-1)/Fs;

beat_types = {'NOTQRS', 'NORMAL', 'LBBB', 'RBBB', 'ABERR', 'PVC', ...
              'FUSION', 'NPC', 'APC', 'SVPB', 'VESC', 'NESC', 'PACE', ...
              'UNKNOWN', 'NOISE', 'ARFCT', 'STCH', 'TCH', 'SYSTOLE', ...
              'DIASTOLE', 'NOTE', 'MEASURE', 'PWAVE', 'BBB', 'PACESP', ...
              'TWAVE', 'RHYTHM', 'UWAVE', 'LEARN', 'FLWAV', 'VFON', ...
              'VFOFF', 'AESC', 'SVESC', 'LINK', 'NAPC', 'PFUS', ...
              'WFON', 'WFOFF', 'RONT'};

figure;
plot(t, ecg_signal);
hold on;

for i = 1:length(r_times)
    r_index = round(r_times(i) * Fs);
    
    if annotations(i) <= length(beat_types)
        beat_type = beat_types{annotations(i)};
    else
        beat_type = 'UNKNOWN';
    end
    
    plot(t(r_index), ecg_signal(r_index), 'ro', 'MarkerSize', 8);
    
    text(t(r_index), ecg_signal(r_index) - 0.05, beat_type, 'FontSize', 8, 'Color', 'blue', 'HorizontalAlignment', 'center');
end

xlabel('Time (s)');
ylabel('ECG Amplitude');
xlim([0,1845])
title('ECG Signal with R-Peak Annotations and Beat Types');
hold off;

%% section 2.1 last
ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');
ecg_channel1 = ECG_data.Sig(:,1);
ecg_channel2 = ECG_data.Sig(:,2);
Fs = ECG_data.sfreq;
n = length(ecg_channel1);
t = (0:n-1)/Fs;

[R_peaks1, locs1] = findpeaks(ecg_channel1, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);
[R_peaks2, locs2] = findpeaks(ecg_channel2, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);

P_waves1 = zeros(size(locs1)); Q_waves1 = zeros(size(locs1));
S_waves1 = zeros(size(locs1)); T_waves1 = zeros(size(locs1));
P_waves2 = zeros(size(locs2)); Q_waves2 = zeros(size(locs2));
S_waves2 = zeros(size(locs2)); T_waves2 = zeros(size(locs2));

P_offset = round(0.05 * Fs);
Q_offset = round(0.01 * Fs);
S_offset = round(0.01 * Fs);
T_offset = round(0.05 * Fs);

for i = 1:length(locs1)
    P_waves1(i) = max(1, locs1(i) - P_offset);
    Q_waves1(i) = max(1, locs1(i) - Q_offset);
    S_waves1(i) = min(n, locs1(i) + S_offset);
    T_waves1(i) = min(n, locs1(i) + T_offset);
end

for i = 1:length(locs2)
    P_waves2(i) = max(1, locs2(i) - P_offset);
    Q_waves2(i) = max(1, locs2(i) - Q_offset);
    S_waves2(i) = min(n, locs2(i) + S_offset);
    T_waves2(i) = min(n, locs2(i) + T_offset);
end

figure;
subplot(2,1,1);
plot(t, ecg_channel1);
hold on;
plot(t(locs1), ecg_channel1(locs1), 'ro', 'MarkerSize', 6);
plot(t(P_waves1), ecg_channel1(P_waves1), 'go', 'MarkerSize', 6);
plot(t(Q_waves1), ecg_channel1(Q_waves1), 'bo', 'MarkerSize', 6);
plot(t(S_waves1), ecg_channel1(S_waves1), 'mo', 'MarkerSize', 6);
plot(t(T_waves1), ecg_channel1(T_waves1), 'co', 'MarkerSize', 6);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 1');
legend('ECG', 'R-peak', 'P', 'Q', 'S', 'T');

subplot(2,1,2);
plot(t, ecg_channel2);
hold on;
plot(t(locs2), ecg_channel2(locs2), 'ro', 'MarkerSize', 6);
plot(t(P_waves2), ecg_channel2(P_waves2), 'go', 'MarkerSize', 6);
plot(t(Q_waves2), ecg_channel2(Q_waves2), 'bo', 'MarkerSize', 6);
plot(t(S_waves2), ecg_channel2(S_waves2), 'mo', 'MarkerSize', 6);
plot(t(T_waves2), ecg_channel2(T_waves2), 'co', 'MarkerSize', 6);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('ECG Channel 2');
legend('ECG', 'R-peak', 'P', 'Q', 'S', 'T');

%% section 2.3

ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');

ecg_channel1 = ECG_data.Sig(:,1); 

Fs = ECG_data.sfreq;

[R_peaks, locs] = findpeaks(ecg_channel1, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);

normal_beat_idx = find(ECG_data.ANNOTD == 1, 1); 
normal_beat_loc = locs(normal_beat_idx);

abnormal_beat_idx = find(ECG_data.ANNOTD == 5, 1); 
abnormal_beat_loc = locs(abnormal_beat_idx);

if isempty(abnormal_beat_idx)
    error('No PVC (Premature Ventricular Contraction) detected in the dataset.');
end

beat_window = round(0.2 * Fs); 

normal_beat = ecg_channel1(normal_beat_loc-beat_window:normal_beat_loc+beat_window);
abnormal_beat = ecg_channel1(abnormal_beat_loc-beat_window:abnormal_beat_loc+beat_window);

t_normal = (-beat_window:beat_window)/Fs;
t_abnormal = (-beat_window:beat_window)/Fs;

figure;
subplot(2,1,1);
plot(t_normal, normal_beat);
xlabel('Time (s)');
ylabel('Amplitude');
title('Normal Beat');
grid on;

subplot(2,1,2);
plot(t_abnormal, abnormal_beat);
xlabel('Time (s)');
ylabel('Amplitude');
title('Abnormal Beat (PVC)');
grid on;

%% section 2.4

ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');

ecg_channel1 = ECG_data.Sig(:,1); 
ecg_channel2 = ECG_data.Sig(:,2); 

Fs = ECG_data.sfreq; 

[R_peaks, locs] = findpeaks(ecg_channel1, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);

normal_beats_idx = find(ECG_data.ANNOTD == 1, 3); 
normal_beats_loc = locs(normal_beats_idx);

abnormal_beats_idx = find(ECG_data.ANNOTD ~= 1, 3); 
abnormal_beats_loc = locs(abnormal_beats_idx);

beat_window = round(1.5 * Fs); 

n = length(ecg_channel1);

normal_start_idx = max(1, normal_beats_loc(1) - beat_window); 
normal_end_idx = min(n, normal_beats_loc(3) + beat_window);  

abnormal_start_idx = max(1, abnormal_beats_loc(1) - beat_window); 
abnormal_end_idx = min(n, abnormal_beats_loc(3) + beat_window); 

normal_segment = ecg_channel1(normal_start_idx:normal_end_idx);
abnormal_segment = ecg_channel1(abnormal_start_idx:abnormal_end_idx);

t_normal = (0:length(normal_segment)-1) / Fs;
t_abnormal = (0:length(abnormal_segment)-1) / Fs;

figure;
subplot(2,1,1);
plot(t_normal, normal_segment);
xlabel('Time (s)');
ylabel('Amplitude');
title('Normal Segment (3 Consecutive Beats)');
grid on;

subplot(2,1,2);
plot(t_abnormal, abnormal_segment);
xlabel('Time (s)');
ylabel('Amplitude');
title('Abnormal Segment (3 Consecutive Beats)');
grid on;
%xlim([0,4.5])

figure;

subplot(2,1,1);
spectrogram(normal_segment, round(0.2*Fs), round(0.15*Fs), 512, Fs, 'yaxis');
title('Normal Segment - Time-Frequency Spectrum');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(2,1,2);
spectrogram(abnormal_segment, round(0.2*Fs), round(0.15*Fs), 512, Fs, 'yaxis');
title('Abnormal Segment - Time-Frequency Spectrum');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

%% section 2.4 + twoside fft

ECG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\ECG_sig.mat');

ecg_channel1 = ECG_data.Sig(:,1); 
ecg_channel2 = ECG_data.Sig(:,2); 

Fs = ECG_data.sfreq; 

[R_peaks, locs] = findpeaks(ecg_channel1, 'MinPeakHeight', 0.5, 'MinPeakDistance', Fs*0.6);

normal_beats_idx = find(ECG_data.ANNOTD == 1, 3); 
normal_beats_loc = locs(normal_beats_idx);

abnormal_beats_idx = find(ECG_data.ANNOTD ~= 1, 3); 
abnormal_beats_loc = locs(abnormal_beats_idx);

beat_window = round(1.5 * Fs); 

n = length(ecg_channel1);

normal_start_idx = max(1, normal_beats_loc(1) - beat_window); 
normal_end_idx = min(n, normal_beats_loc(3) + beat_window);  

abnormal_start_idx = max(1, abnormal_beats_loc(1) - beat_window); 
abnormal_end_idx = min(n, abnormal_beats_loc(3) + beat_window); 

normal_segment = ecg_channel1(normal_start_idx:normal_end_idx);
abnormal_segment = ecg_channel1(abnormal_start_idx:abnormal_end_idx);

t_normal = (0:length(normal_segment)-1) / Fs;
t_abnormal = (0:length(abnormal_segment)-1) / Fs;

figure;
subplot(2,1,1);
plot(t_normal, normal_segment);
xlabel('Time (s)');
ylabel('Amplitude');
title('Normal Segment (3 Consecutive Beats)');
grid on;

subplot(2,1,2);
plot(t_abnormal, abnormal_segment);
xlabel('Time (s)');
ylabel('Amplitude');
title('Abnormal Segment (3 Consecutive Beats)');
xlim([0,4.5])
grid on;

figure;
subplot(2,1,1);
spectrogram(normal_segment, round(0.2*Fs), round(0.15*Fs), 512, Fs, 'yaxis');
title('Normal Segment - Time-Frequency Spectrum');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(2,1,2);
spectrogram(abnormal_segment, round(0.2*Fs), round(0.15*Fs), 512, Fs, 'yaxis');
title('Abnormal Segment - Time-Frequency Spectrum');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

fft_normal = fftshift(fft(normal_segment));  
fft_abnormal = fftshift(fft(abnormal_segment));

n_normal = length(normal_segment);
n_abnormal = length(abnormal_segment);

f_normal = (-n_normal/2:n_normal/2-1)*(Fs/n_normal);  
f_abnormal = (-n_abnormal/2:n_abnormal/2-1)*(Fs/n_abnormal);  

figure;
subplot(2,1,1);
plot(f_normal, abs(fft_normal));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Normal Segment');
xlim([-500 500]);  

subplot(2,1,2);
plot(f_abnormal, abs(fft_abnormal));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Abnormal Segment');
xlim([-500 500]);  

%% section 3.1

EOG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\EOG_sig.mat');

eog_left = EOG_data.Sig(1,:); 
eog_right = EOG_data.Sig(2,:); 

Fs = EOG_data.fs;  

n = length(eog_left); 
t = (0:n-1)/Fs;

figure;
plot(t, eog_left,Color='b',DisplayName='Left Eye');
xlabel('Time (s)');
ylabel('Amplitude');
legend
hold on;
plot(t, eog_right,Color='r',DisplayName='Right Eye');
title('EOG Signal');
legend;
grid on;

horizontal_movement = eog_left - eog_right;

vertical_movement = eog_left + eog_right;

figure;
subplot(2,1,1);
plot(t, horizontal_movement);
xlabel('Time (s)');
ylabel('Amplitude');
title('Horizontal Eye Movement');
grid on;

subplot(2,1,2);
plot(t, vertical_movement);
xlabel('Time (s)');
ylabel('Amplitude');
title('Vertical Eye Movement');
grid on;

%% section 3.2

n_left = length(eog_left);
n_right = length(eog_right);

f_left = (-n_left/2:n_left/2-1)*(Fs/n_left);
f_right = (-n_right/2:n_right/2-1)*(Fs/n_right);

fft_left = fftshift(fft(eog_left));
fft_right = fftshift(fft(eog_right));

figure;
subplot(2,1,1);
plot(f_left, abs(fft_left));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Left Eye');
xlim([-50 50]);  

subplot(2,1,2);
plot(f_right, abs(fft_right));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Right Eye');
xlim([-50 50]);

window_length = min(64, n_left); 
noverlap = round(0.5 * window_length); 

figure;

subplot(2,1,1);
spectrogram(eog_left, window_length, noverlap, 128, Fs, 'yaxis');
title('Time-Frequency Spectrum - Left Eye');

subplot(2,1,2);
spectrogram(eog_right, window_length, noverlap, 128, Fs, 'yaxis');
title('Time-Frequency Spectrum - Right Eye');

%% section 4.1

EMG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\EMG_sig.mat');

emg_healthy = EMG_data.emg_healthym;   
emg_myopathy = EMG_data.emg_myopathym; 
emg_neuropathy = EMG_data.emg_neuropathym; 

Fs = EMG_data.fs;  

t_healthy = (0:length(emg_healthy)-1)/Fs;
t_myopathy = (0:length(emg_myopathy)-1)/Fs;
t_neuropathy = (0:length(emg_neuropathy)-1)/Fs;

figure;
subplot(3,1,1);
plot(t_healthy, emg_healthy);
xlabel('Time (s)');
ylabel('Amplitude');
title('EMG Signal - Healthy Individual');
grid on;

subplot(3,1,2);
plot(t_myopathy, emg_myopathy);
xlabel('Time (s)');
ylabel('Amplitude');
title('EMG Signal - Individual with Myopathy');
grid on;

subplot(3,1,3);
plot(t_neuropathy, emg_neuropathy);
xlabel('Time (s)');
ylabel('Amplitude');
title('EMG Signal - Individual with Neuropathy');
grid on;

zoom_window = 0.5; 

figure;
subplot(3,1,1);
plot(t_healthy, emg_healthy);
xlim([0 zoom_window]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Zoomed EMG Signal - Healthy Individual');
grid on;

subplot(3,1,2);
plot(t_myopathy, emg_myopathy);
xlim([0 zoom_window]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Zoomed EMG Signal - Individual with Myopathy');
grid on;

subplot(3,1,3);
plot(t_neuropathy, emg_neuropathy);
xlim([0 zoom_window]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Zoomed EMG Signal - Individual with Neuropathy');
grid on;

%% section 4.2

EMG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\EMG_sig.mat');

emg_healthy = EMG_data.emg_healthym;   
emg_myopathy = EMG_data.emg_myopathym; 
emg_neuropathy = EMG_data.emg_neuropathym; 

Fs = EMG_data.fs;  

n_healthy = length(emg_healthy);
n_myopathy = length(emg_myopathy);
n_neuropathy = length(emg_neuropathy);

f_healthy = (0:n_healthy-1)*(Fs/n_healthy);
f_myopathy = (0:n_myopathy-1)*(Fs/n_myopathy);
f_neuropathy = (0:n_neuropathy-1)*(Fs/n_neuropathy);

fft_healthy = abs(fft(emg_healthy));
fft_myopathy = abs(fft(emg_myopathy));
fft_neuropathy = abs(fft(emg_neuropathy));

figure;
subplot(3,1,1);
plot(f_healthy, fft_healthy);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum - Healthy Individual');
xlim([0 500]);  

subplot(3,1,2);
plot(f_myopathy, fft_myopathy);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum - Individual with Myopathy');
xlim([0 500]);

subplot(3,1,3);
plot(f_neuropathy, fft_neuropathy);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum - Individual with Neuropathy');
xlim([0 500]);

figure;

subplot(3,1,1);
spectrogram(emg_healthy, 128, 120, 128, Fs, 'yaxis');
title('Time-Frequency Spectrum - Healthy Individual');

subplot(3,1,2);
spectrogram(emg_myopathy, 128, 120, 128, Fs, 'yaxis');
title('Time-Frequency Spectrum - Individual with Myopathy');

subplot(3,1,3);
spectrogram(emg_neuropathy, 128, 120, 128, Fs, 'yaxis');
title('Time-Frequency Spectrum - Individual with Neuropathy');

%% section 4.2 fft twoside

EMG_data = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab1\Lab 1_data\EMG_sig.mat');

emg_healthy = EMG_data.emg_healthym;   
emg_myopathy = EMG_data.emg_myopathym; 
emg_neuropathy = EMG_data.emg_neuropathym; 

Fs = EMG_data.fs;  

n_healthy = length(emg_healthy);
n_myopathy = length(emg_myopathy);
n_neuropathy = length(emg_neuropathy);

fft_healthy = fftshift(fft(emg_healthy));  
fft_myopathy = fftshift(fft(emg_myopathy));
fft_neuropathy = fftshift(fft(emg_neuropathy));

f_healthy = (-n_healthy/2:n_healthy/2-1)*(Fs/n_healthy);  
f_myopathy = (-n_myopathy/2:n_myopathy/2-1)*(Fs/n_myopathy);
f_neuropathy = (-n_neuropathy/2:n_neuropathy/2-1)*(Fs/n_neuropathy);

figure;
subplot(3,1,1);
plot(f_healthy, abs(fft_healthy));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Healthy Individual');
xlim([-500 500]);  

subplot(3,1,2);
plot(f_myopathy, abs(fft_myopathy));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Individual with Myopathy');
xlim([-500 500]);

subplot(3,1,3);
plot(f_neuropathy, abs(fft_neuropathy));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Full Frequency Spectrum - Individual with Neuropathy');
xlim([-500 500]);

