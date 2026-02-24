%% section 1.1
load('X_org.mat')
load('X_noise.mat')
load('Electrodes') ;
offset = max(abs(X_org(:))) ;
feq = 250 ; 
ElecName = Electrodes.labels ;
disp_eeg(X_org,offset,feq,ElecName);
xlim("tight")
%% section 1.2
offset = max(abs(X_noise(:))) ;
feq = 250 ; 
ElecName = Electrodes.labels ;
disp_eeg(X_noise,offset,feq,ElecName);
xlim("tight")
%% section 1.3
calculate_energy = @(signal) sum(signal(:).^2);
Psignal = calculate_energy(X_org);
Pnoise = calculate_energy(X_noise);

SNR1 = -5;
SNR2 = -15;

sigma_1 = sqrt(Psignal / Pnoise * 10^(-SNR1 / 10));
sigma_2 = sqrt(Psignal / Pnoise * 10^(-SNR2 / 10));

noisy_sig_1 = X_org + sigma_1 * X_noise;
noisy_sig_2 = X_org + sigma_2 * X_noise;

offset = max(max(abs(noisy_sig_1)))/2;
feq = 250;
ElecName = Electrodes.labels;

disp_eeg(noisy_sig_1, offset, feq, ElecName);
title('Noisy EEG Signal with SNR = -5 dB' , 'Interpreter','latex');
offset = max(max(abs(noisy_sig_2)))/2;
disp_eeg(noisy_sig_2, offset, feq, ElecName);
title('Noisy EEG Signal with SNR = -15 dB' , 'Interpreter','latex');
%% section 1.4

Pest = 32;

[F, W, K] = COM2R(noisy_sig_1, Pest);

sources = W * noisy_sig_1;

offset = max(abs(noisy_sig_1(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(noisy_sig_1,offset,feq,ElecName);
title('noisy Signal with SNR = -5db' , 'Interpreter','latex');

offset = max(abs(sources(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(sources,offset,feq,ElecName);
title('Extracted Independent Components using COM2R' , 'Interpreter','latex');


Pest = 32;

[F2, W2, K2] = COM2R(noisy_sig_2, Pest);

sources2 = W2 * noisy_sig_2;

offset = max(abs(noisy_sig_2(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(noisy_sig_2,offset,feq,[]);
title('noisy Signal with SNR = -15db' , 'Interpreter','latex');

offset = max(abs(sources2(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(sources2,offset,feq,[]);
title('Extracted Independent Components using COM2R' , 'Interpreter','latex');
%% section 1.5 & 1.6
selected_components = [2,5];

mask = zeros(size(sources));
mask(selected_components, :) = sources(selected_components, :);

X_den = F * mask;

offset = max(abs(X_den(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X_den,offset,feq,ElecName);
title('Reconstructed Signal after Removing Unwanted Components for SNR = -5 db', 'Interpreter','latex');

selected_components2 = [15,18,19];

mask2 = zeros(size(sources2));
mask2(selected_components2, :) = sources2(selected_components2, :);

X_den2 = F2 * mask2;

offset2 = max(abs(X_den2(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X_den2,offset2,feq,ElecName);
title('Reconstructed Signal after Removing Unwanted Components for SNR = -15 db', 'Interpreter','latex');
%% section 1.7
figure;
[N , K] = size(X_org);
t = (1:K)/feq;
subplot(3,1,1)
plot(t, X_org(13 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("original signal in channel 13" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(noisy_sig_1);
t = (1:K)/feq;
subplot(3,1,2)
plot(t, noisy_sig_1(13 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("noisy signal with SNR = -5 db in channel 13" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(X_den);
t = (1:K)/feq;
subplot(3,1,3)
plot(t, X_den(13 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("reconstructed signal with SNR = -5 db in channel 13" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

figure
[N , K] = size(X_org);
t = (1:K)/feq;
subplot(3,1,1)
plot(t, X_org(24 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("original signal in channel 24" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(noisy_sig_1);
t = (1:K)/feq;
subplot(3,1,2)
plot(t, noisy_sig_1(24 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("noisy signal with SNR = -5 db in channel 24" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(X_den);
t = (1:K)/feq;
subplot(3,1,3)
plot(t, X_den(24 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("reconstructed signal with SNR = -5 db in channel 24" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

figure;
[N , K] = size(X_org);
t = (1:K)/feq;
subplot(3,1,1)
plot(t, X_org(13 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("original signal in channel 13" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(noisy_sig_2);
t = (1:K)/feq;
subplot(3,1,2)
plot(t, noisy_sig_2(13 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("noisy signal with SNR = -15 db in channel 13" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(X_den2);
t = (1:K)/feq;
subplot(3,1,3)
plot(t, X_den2(13 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("reconstructed signal with SNR = -15 db in channel 13" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

figure
[N , K] = size(X_org);
t = (1:K)/feq;
subplot(3,1,1)
plot(t, X_org(24 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("original signal in channel 24" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(noisy_sig_2);
t = (1:K)/feq;
subplot(3,1,2)
plot(t, noisy_sig_2(24 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("noisy signal with SNR = -15 db in channel 24" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

[N , K] = size(X_den2);
t = (1:K)/feq;
subplot(3,1,3)
plot(t, X_den2(24 , :) , "Color","#2596be")
xlim("tight")
ylim("tight")
title("reconstructed signal with SNR = -15 db in channel 24" , 'Interpreter','latex')
xlabel("Time(s)" , 'Interpreter','latex')
ylabel('Amplitude ($\mu$V)' , Interpreter='latex');

%% section 1.8
numerator = sum(sum((X_org - X_den).^2));
denominator = sum(sum(X_org.^2));

RRMSE = sqrt(numerator / denominator);

disp(['RRMSE for SNR = -5 db : ', num2str(RRMSE)]);

numerator = sum(sum((X_org - X_den2).^2));
denominator = sum(sum(X_org.^2));

RRMSE = sqrt(numerator / denominator);

disp(['RRMSE for SNR = -15 db : ', num2str(RRMSE)]);
%% section 2.1
clc; clear;
data1 = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Lab2_2\NewData1.mat');
data2 = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Lab2_2\NewData2.mat');
data3 = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Lab2_2\NewData3.mat');
data4 = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Lab2_2\NewData4.mat');

fs = 250;
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Lab2_2\Electrodes.mat');
electrode_labels = Electrodes.labels;  

EEG_Sig_1 = data1.EEG_Sig;
EEG_Sig_2 = data2.EEG_Sig;
EEG_Sig_3 = data3.EEG_Sig;
EEG_Sig_4 = data4.EEG_Sig;

load('Electrodes.mat');  
ElecName = Electrodes.labels; 

X = EEG_Sig_1;
offset = max(abs(X(:))); 
feq = 250; 
titre = 'EEG Signal 1';
disp_eeg(X, offset, feq, ElecName, titre);

X = EEG_Sig_2;
offset = max(abs(X(:)));  
titre = 'EEG Signal 2';
disp_eeg(X, offset, feq, ElecName, titre);

X = EEG_Sig_3;
offset = max(abs(X(:)));  
titre = 'EEG Signal 3';
disp_eeg(X, offset, feq, ElecName, titre);

X = EEG_Sig_4;
offset = max(abs(X(:)));  
titre = 'EEG Signal 4';
disp_eeg(X, offset, feq, ElecName, titre);

%% section 2.2

% no need for any code in this part

%% section 2.3 for signal 1

EEG_Sig = data1.EEG_Sig; 

Pest = size(EEG_Sig, 1);  

[F, W, K] = COM2R(EEG_Sig, Pest);

offset = max(abs(EEG_Sig(:))); 
feq = 250;  
titre = 'Original EEG Signal 1';
disp_eeg(EEG_Sig, offset, feq, electrode_labels, titre);

disp('Mixing matrix (W):');
disp(W);

disp('Number of sweeps (K):');
disp(K);

EEG_Sig_ICA=W*EEG_Sig;
offset = max(abs(EEG_Sig_ICA(:)));
titre ="After apply ICA Signal 1";
disp_eeg(EEG_Sig_ICA, offset, feq, electrode_labels, titre);

%% section 2.3 for signal 2

EEG_Sig2 = data2.EEG_Sig; 

Pest2 = size(EEG_Sig2, 1);  

[F2, W2, K2] = COM2R(EEG_Sig2, Pest2);

offset2 = max(abs(EEG_Sig2(:))); 
feq = 250;  
titre = 'Original EEG Signal 2';
disp_eeg(EEG_Sig2, offset2, feq, electrode_labels, titre);

disp('Mixing matrix (W):');
disp(W2);

disp('Number of sweeps (K):');
disp(K2);

EEG_Sig_ICA2=W2*EEG_Sig2;
offset2 = max(abs(EEG_Sig_ICA2(:)));
titre ="After apply ICA Signal 2";
disp_eeg(EEG_Sig_ICA2, offset2, feq, electrode_labels, titre);

%% section 2.4

figure;
for i = 1:size(EEG_Sig_ICA, 1)
    subplot(ceil(size(EEG_Sig_ICA, 1) / 4), 4, i);  % Adjust subplot grid based on the number of components
    plot(EEG_Sig_ICA(i, :));
    title(['Time domain of Component ', num2str(i)]);
end
%titre="Time Domain of ICA Components";

figure;
for i = 1:size(EEG_Sig_ICA, 1)
    subplot(ceil(size(EEG_Sig_ICA, 1) / 4), 4, i);
    [Pxx, Fn] = pwelch(EEG_Sig_ICA(i, :), [], [], [], fs);  % Welch's method to compute power spectrum
    plot(Fn, 10 * log10(Pxx));
    title(['Frequency Spectrum of Component ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
end
%titre= "Frequency Domain of ICA Components";

% Spatial Domain: Correcting the F matrix size mismatch
Electrodes = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Electrodes.mat');


% topography
elocsX = Electrodes.Electrodes.X;
elocsY = Electrodes.Electrodes.Y;
elabels = Electrodes.Electrodes.labels;

figure;
for i = 1:size(EEG_Sig_ICA, 1)
    subplot(6, 4, i);
    plottopomap(elocsX,elocsY,elabels(:), F(:,i));
    title(['Spatial Map of Component ', num2str(i)]);
end

%% section 2.4 for signal 2

figure;
for i = 1:size(EEG_Sig_ICA2, 1)
    subplot(ceil(size(EEG_Sig_ICA2, 1) / 4), 4, i);
    plot(EEG_Sig_ICA2(i, :));
    title(['Time domain of Component ', num2str(i)]);
end

figure;
for i = 1:size(EEG_Sig_ICA2, 1)
    subplot(ceil(size(EEG_Sig_ICA2, 1) / 4), 4, i);
    [Pxx, Fn] = pwelch(EEG_Sig_ICA2(i, :), [], [], [], fs);
    plot(Fn, 10 * log10(Pxx));
    title(['Frequency Spectrum of Component ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
end

Electrodes = load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab2\Electrodes.mat');

elocsX = Electrodes.Electrodes.X;
elocsY = Electrodes.Electrodes.Y;
elabels = Electrodes.Electrodes.labels;

figure;
for i = 1:size(EEG_Sig_ICA2, 1)
    subplot(6, 4, i);
    plottopomap(elocsX,elocsY,elabels(:), F2(:,i));
    title(['Spatial Map of Component ', num2str(i)]);
end
%% section 2.5

SelSources = [4,10];

EEG_Sig_ICA_denoised = EEG_Sig_ICA;
EEG_Sig_ICA_denoised(SelSources, :) = 0;

X_denoised = F * EEG_Sig_ICA_denoised;

%% section 2.5 for signal 2

SelSources2 = [4,9,12,18];

EEG_Sig_ICA_denoised2 = EEG_Sig_ICA2;
EEG_Sig_ICA_denoised2(SelSources2, :) = 0;

X_denoised2 = F2 * EEG_Sig_ICA_denoised2;

%% section 2.6

figure;
subplot(2, 1, 1);
plot(EEG_Sig(1, :));
title('Original EEG Signal (Channel 1)');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(X_denoised(1, :));
title('Denoised EEG Signal (Channel 1)');
xlabel('Time');
ylabel('Amplitude');

original_signal = EEG_Sig(1, :);
denoised_signal = X_denoised(1, :);

snr_original = snr(original_signal);
snr_denoised = snr(denoised_signal);

disp(['SNR of Original Signal (Channel 1): ', num2str(snr_original), ' dB']);
disp(['SNR of Denoised Signal (Channel 1): ', num2str(snr_denoised), ' dB']);

[Pxx_original, f_original] = pwelch(original_signal, [], [], [], fs);
[Pxx_denoised, f_denoised] = pwelch(denoised_signal, [], [], [], fs);

figure;
subplot(2, 1, 1);
plot(f_original, 10 * log10(Pxx_original));
title('Power Spectral Density of Original Signal (Channel 1)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

subplot(2, 1, 2);
plot(f_denoised, 10 * log10(Pxx_denoised));
title('Power Spectral Density of Denoised Signal (Channel 1)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%% section 2.6 for signal 2

figure;
subplot(2, 1, 1);
plot(EEG_Sig2(1, :));
title('Original EEG Signal 2(Channel 1)');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(X_denoised2(1, :));
title('Denoised EEG Signal 2(Channel 1)');
xlabel('Time');
ylabel('Amplitude');

original_signal2 = EEG_Sig2(1, :);
denoised_signal2 = X_denoised2(1, :);

snr_original = snr(original_signal2);
snr_denoised = snr(denoised_signal2);

disp(['SNR of Original Signal (Channel 1): ', num2str(snr_original), ' dB']);
disp(['SNR of Denoised Signal (Channel 1): ', num2str(snr_denoised), ' dB']);

[Pxx_original, f_original] = pwelch(original_signal2, [], [], [], fs);
[Pxx_denoised, f_denoised] = pwelch(denoised_signal2, [], [], [], fs);

figure;
subplot(2, 1, 1);
plot(f_original, 10 * log10(Pxx_original));
title('Power Spectral Density of Original Signal (Channel 1)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

subplot(2, 1, 2);
plot(f_denoised, 10 * log10(Pxx_denoised));
title('Power Spectral Density of Denoised Signal (Channel 1)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

%% section 2.6 example for bad choices

SelSources = [3,5,8];

EEG_Sig_ICA_denoised = EEG_Sig_ICA;
EEG_Sig_ICA_denoised(SelSources, :) = 0;

X_denoised = F * EEG_Sig_ICA_denoised;

figure;
subplot(2, 1, 1);
plot(EEG_Sig(1, :));
title('Original EEG Signal (Channel 1)');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(X_denoised(1, :));
title('Denoised EEG Signal (Channel 1)');
xlabel('Time');
ylabel('Amplitude');

original_signal = EEG_Sig(1, :);
denoised_signal = X_denoised(1, :);

snr_original = snr(original_signal);
snr_denoised = snr(denoised_signal);

disp(['SNR of Original Signal (Channel 1): ', num2str(snr_original), ' dB']);
disp(['SNR of Denoised Signal (Channel 1): ', num2str(snr_denoised), ' dB']);

[Pxx_original, f_original] = pwelch(original_signal, [], [], [], fs);
[Pxx_denoised, f_denoised] = pwelch(denoised_signal, [], [], [], fs);

figure;
subplot(2, 1, 1);
plot(f_original, 10 * log10(Pxx_original));
title('Power Spectral Density of Original Signal (Channel 1)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

subplot(2, 1, 2);
plot(f_denoised, 10 * log10(Pxx_denoised));
title('Power Spectral Density of Denoised Signal (Channel 1)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');