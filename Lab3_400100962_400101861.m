%% section 1.1


mother_signal = mecg1;
fetus_signal = fecg1;    
noise_signal = noise1;

Fs = 256;

t = (0:length(mother_signal)-1) / Fs;

figure;
subplot(4, 1, 1);
plot(t, mother_signal);
title('Mother Signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude (mV)', 'Interpreter', 'latex');
grid on;

subplot(4, 1, 2);
plot(t, fetus_signal);
title('Fetus Signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude (mV)', 'Interpreter', 'latex');
grid on;

subplot(4, 1, 3);
plot(t, noise_signal);
title('Noise Signal', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude (mV)', 'Interpreter', 'latex');
grid on;

combined_signal = mother_signal + fetus_signal + noise_signal;

subplot(4, 1, 4);
plot(t, combined_signal);
title('Combined Signal (Mother + Fetus + Noise)', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Amplitude (mV)', 'Interpreter', 'latex');
grid on;

%% section 1.2

Fs = 256;

figure;
subplot(3, 1, 1);
pwelch(mother_signal, [], [], [], Fs);
title('Power Spectrum of Mother Signal', 'Interpreter', 'latex');

subplot(3, 1, 2);
pwelch(fetus_signal, [], [], [], Fs);
title('Power Spectrum of Fetus Signal', 'Interpreter', 'latex');

subplot(3, 1, 3);
pwelch(noise_signal, [], [], [], Fs);
title('Power Spectrum of Noise Signal', 'Interpreter', 'latex');

%% section 1.3

mean_mother = mean(mother_signal);
mean_fetus = mean(fetus_signal);
mean_noise = mean(noise_signal);

var_mother = var(mother_signal);
var_fetus = var(fetus_signal);
var_noise = var(noise_signal);

fprintf('Mean of Mother Signal: %.4f\n', mean_mother);
fprintf('Variance of Mother Signal: %.4f\n', var_mother);

fprintf('Mean of Fetus Signal: %.4f\n', mean_fetus);
fprintf('Variance of Fetus Signal: %.4f\n', var_fetus);

fprintf('Mean of Noise Signal: %.4f\n', mean_noise);
fprintf('Variance of Noise Signal: %.4f\n', var_noise);

%% section 1.4

figure;
subplot(3, 1, 1);
histogram(mother_signal, 'Normalization', 'pdf');
title('Histogram of Mother Signal (PDF Approximation)', 'Interpreter', 'latex');
xlabel('Amplitude', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');

subplot(3, 1, 2);
histogram(fetus_signal, 'Normalization', 'pdf');
title('Histogram of Fetus Signal (PDF Approximation)', 'Interpreter', 'latex');
xlabel('Amplitude', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');

subplot(3, 1, 3);
histogram(noise_signal, 'Normalization', 'pdf');
title('Histogram of Noise Signal (PDF Approximation)', 'Interpreter', 'latex');
xlabel('Amplitude', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');

kurtosis_mother = kurtosis(mother_signal);
kurtosis_fetus = kurtosis(fetus_signal);
kurtosis_noise = kurtosis(noise_signal);

fprintf('Kurtosis of Mother Signal: %.4f\n', kurtosis_mother);
fprintf('Kurtosis of Fetus Signal: %.4f\n', kurtosis_fetus);
fprintf('Kurtosis of Noise Signal: %.4f\n', kurtosis_noise);

%% section 2.1

data_double = X;

plot3ch(data_double , Fs , "3d plot of data")

[U, S, V] = svd(data_double);

singular_values = diag(S);

%% section 2.2

num_vectors = size(V, 2);

for i = 1:num_vectors
    singular_value = S(i,i);
    eigenvector = V(:,i);
    plot3dv(eigenvector, singular_value, 'r'); 
end
title('3D Plot of Eigenvectors (V) with Singular Values', 'Interpreter', 'latex');
xlabel('X-axis', 'Interpreter', 'latex');
ylabel('Y-axis', 'Interpreter', 'latex');
zlabel('Z-axis', 'Interpreter', 'latex');
grid on;

saveas(gcf, 'svd_eigenvectors_plot.fig');

save('SVD_matrices.mat', 'U', 'S', 'V');

%% section 2.3

figure;
subplot(3, 1, 1);
plot(U(:,1));
title('First Column of U', 'Interpreter', 'latex');
xlabel('Index', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
xlim("tight")
subplot(3, 1, 2);
plot(U(:,2));
title('Second Column of U', 'Interpreter', 'latex');
xlabel('Index', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
xlim("tight")

subplot(3, 1, 3);
plot(U(:,3));
title('Third Column of U', 'Interpreter', 'latex');
xlabel('Index', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
xlim("tight")

singular_values = diag(S);
eigenspectrum = singular_values;

figure;
stem(eigenspectrum);
title('Eigenspectrum (Variance Explained)', 'Interpreter', 'latex');
xlabel('Component', 'Interpreter', 'latex');
ylabel('Variance Explained (%)', 'Interpreter', 'latex');
grid on;
xlim([0 4])
fprintf('Largest variance is explained by the first component: %.2f%%\n', eigenspectrum(1));
%% section 2.4
num_components = 2;

S_reduced = S;
S_reduced(num_components+1:end, num_components+1:end) = 0;
S_reduced(1,1)=0;
data_reconstructed = U * S_reduced * V';

figure;
subplot(3, 1, 1);
plot(data_reconstructed(:,1));
title('Reconstructed Sensor Channel 1', 'Interpreter', 'latex');
xlabel('Index', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
xlim("tight")
subplot(3, 1, 2);
plot(data_reconstructed(:,2));
title('Reconstructed Sensor Channel 2', 'Interpreter', 'latex');
xlabel('Index', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
xlim("tight")

subplot(3, 1, 3);
plot(data_reconstructed(:,3));
title('Reconstructed Sensor Channel 3', 'Interpreter', 'latex');
xlabel('Index', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');
grid on;
xlim("tight")

reconstruction_error = norm(data_double - data_reconstructed, 'fro');
fprintf('Reconstruction Error: %.4f\n', reconstruction_error);

if reconstruction_error < 0.01
    disp('The fetal signal reconstruction was successful.');
else
    disp('The fetal signal reconstruction was not successful.');
end

%% Section 3.1: Plotting the signals

clear;
close all;

load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\mecg1.dat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\fecg1.dat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\noise1.dat');

fs = 256;

X_r = mecg1 + fecg1 + noise1;

t = 0:1/fs:(length(X_r)-1)/fs;

subplot(4,1,1)
plot(t, mecg1);
xlabel('time (sec)');
ylabel('amp (mV)');
title('Mother ECG');

subplot(4,1,2)
plot(t, fecg1);
xlabel('time (sec)');
ylabel('amp (mV)');
title('Fetal ECG');

subplot(4,1,3)
plot(t, noise1);
xlabel('time (sec)');
ylabel('amp (mV)');
title('Noise');

subplot(4,1,4)
plot(t, X_r);
xlabel('time (sec)');
ylabel('amp (mV)');
title('Recorded Mixed Signal');

%% Section 3.2: Scatter plot of original data using plot3ch and plot of W^(-1) using plot3dv

clc;
close all;
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\X.dat');

figure;
plot3ch(X, fs, 'Scatter Plot of Original Mixed Data');

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\scatter_plot_original_data.fig');

[W, Z_hat] = ica(X');

A = inv(W);

figure;
plot3dv(A(:,1), [], 'r');
hold on;
plot3dv(A(:,2), [], 'g');
plot3dv(A(:,3), [], 'b');
title('Visualization of W^{-1} Columns Using plot3dv');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\W_inv_visualization.fig');


%% Section 3.3: Visualization of A columns using plot3dv

figure();
plot3dv(A(:,1), [], 'r');
hold on;
plot3dv(A(:,2), [], 'g');
plot3dv(A(:,3), [], 'b');
title('Visualization of A Vectors');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\A_vectors_visualization.fig');

%% Section 3.4: Plot Z_hat (3-column data)

fs = 256;
t = 0:1/fs:(size(Z_hat, 2) - 1) / fs;

figure();
plot3ch(Z_hat.', fs, 'Separated Sources (Z_hat)');

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\Z_hat_visualization.fig');

%% Section 3.5: Reconstruct fECG Signal

Z_des = [zeros(2, size(Z_hat, 2)); Z_hat(3,:)];
X_des = A * Z_des; 

plot3ch(X_des.', fs, 'Reconstructed fECG Signal Using ICA');
title('Reconstructed Signal Using ICA');

save('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\X_des_reconstructed.mat', 'X_des');

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\reconstructed_signal_visualization.fig');

%% Section 4.1: Comparison of Scatter Plots, Calculation of Angles, and Observation Matrices

clc;
close all;

load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\X.dat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\fecg2.dat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\W_matrix.mat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\Z_hat_matrix.mat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\A_matrix.mat');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\X_des_reconstructed.mat');
load('SVD_matrices.mat', 'U', 'S', 'V');

figure;
plot3ch(X, fs, 'Scatter Plot of Original Observation Matrix (X)');
savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\scatter_plot_observation_matrix.fig');

figure;
plot3ch(X_des.', fs, 'Scatter Plot of Reconstructed Signal (X_{\text{reconstructed}}) from ICA');
savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\scatter_plot_reconstructed_signal_ica.fig');

S_reduced = S;
data_reconstructed_SVD = U * S_reduced * V';
figure;
plot3ch(data_reconstructed_SVD, fs, 'Scatter Plot of Reconstructed Signal from SVD');
savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\scatter_plot_reconstructed_signal_svd.fig');

figure;
plot3dv(A(:,1), [], 'r'); 
hold on;
plot3dv(A(:,2), [], 'g');
plot3dv(A(:,3), [], 'b');
title('Visualization of W^{-1} Columns');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\W_inv_visualization.fig');

figure;
plot3dv(V(:,1), [], 'r');
hold on;
plot3dv(V(:,2), [], 'g');
plot3dv(V(:,3), [], 'b');
title('Visualization of V Columns (SVD)');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\V_columns_visualization.fig');

v1 = A(:,1); v2 = A(:,2); v3 = A(:,3); 
w1 = W(:,1); w2 = W(:,2); w3 = W(:,3);

v1 = v1 / norm(v1);
v2 = v2 / norm(v2);
v3 = v3 / norm(v3);
w1 = w1 / norm(w1);
w2 = w2 / norm(w2);
w3 = w3 / norm(w3);

angle_v1_w1 = acosd(dot(v1, w1));
angle_v2_w2 = acosd(dot(v2, w2));
angle_v3_w3 = acosd(dot(v3, w3));  

disp(['Angle between v1 and w1: ', num2str(angle_v1_w1), ' degrees']);
disp(['Angle between v2 and w2: ', num2str(angle_v2_w2), ' degrees']);
disp(['Angle between v3 and w3: ', num2str(angle_v3_w3), ' degrees']);

figure;
plot(fecg2); 
hold on;
plot(X_des(3,:), '--'); 
title('Comparison of Ground Truth fECG and Reconstructed fECG (ICA)');
xlabel('Samples');
ylabel('Amplitude');
legend('Ground Truth fECG', 'ICA Reconstructed fECG');

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\comparison_fECG_ICA.fig');

figure;
plot(fecg2);
hold on;
plot(data_reconstructed_SVD(:,3), '--');
title('Comparison of Ground Truth fECG and Reconstructed fECG (SVD)');
xlabel('Samples');
ylabel('Amplitude');
legend('Ground Truth fECG', 'SVD Reconstructed fECG');

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\comparison_fECG_SVD.fig');



%% Section 4.2: Plot and Compare ICA and SVD Reconstructed Signals in Separate Subplots

load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\fecg2.dat');
load('SVD_matrices.mat', 'U', 'S', 'V');
load('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\X_des_reconstructed.mat');

fecg_ground_truth = fecg2 / max(abs(fecg2));

num_components = 3;
S_reduced = S;
S_reduced(num_components+1:end, num_components+1:end) = 0;
data_reconstructed_SVD = U * S_reduced * V';

fECG_reconstructed_SVD = data_reconstructed_SVD(:, 3) / max(abs(data_reconstructed_SVD(:, 3)));

fECG_reconstructed_ICA = X_des(3,:) / max(abs(X_des(3,:)));

min_len = min([length(fecg_ground_truth), length(fECG_reconstructed_SVD), length(fECG_reconstructed_ICA)]);
fecg_ground_truth = fecg_ground_truth(1:min_len);
fECG_reconstructed_SVD = fECG_reconstructed_SVD(1:min_len);
fECG_reconstructed_ICA = fECG_reconstructed_ICA(1:min_len);

figure;

subplot(2, 1, 1);
plot(fecg_ground_truth, 'LineWidth', 1.5);
hold on;
plot(fECG_reconstructed_ICA, '--', 'LineWidth', 1.5);
title('Comparison of Ground Truth fECG and ICA Reconstructed fECG');
xlabel('Samples');
ylabel('Normalized Amplitude');
legend('Ground Truth fECG', 'ICA Reconstructed fECG');
grid on;

subplot(2, 1, 2);
plot(fecg_ground_truth, 'LineWidth', 1.5);
hold on;
plot(fECG_reconstructed_SVD, '--', 'LineWidth', 1.5);
title('Comparison of Ground Truth fECG and SVD Reconstructed fECG');
xlabel('Samples');
ylabel('Normalized Amplitude');
legend('Ground Truth fECG', 'SVD Reconstructed fECG');
grid on;

savefig('C:\Users\ASUS\Desktop\uni\term 7\Medical Signal Proccesing Lab\Lab3\Lab3_data\comparison_SVD_ICA_subplots.fig');

%% Calculate and Print Reconstruction Errors for ICA and SVD

ica_reconstruction_error = norm(fecg_ground_truth - fECG_reconstructed_ICA', 'fro');

svd_reconstruction_error = norm(fecg_ground_truth - fECG_reconstructed_SVD', 'fro');

fprintf('Reconstruction Error for ICA: %.4f\n', ica_reconstruction_error);
fprintf('Reconstruction Error for SVD: %.4f\n', svd_reconstruction_error);


%% Section 4.3: Calculate Correlation Coefficients for ICA and SVD Reconstructed Signals

min_len = min([length(fecg_ground_truth), length(fECG_reconstructed_SVD), length(fECG_reconstructed_ICA)]);
fecg_ground_truth = fecg_ground_truth(1:min_len);
fECG_reconstructed_SVD = fECG_reconstructed_SVD(1:min_len);
fECG_reconstructed_ICA = fECG_reconstructed_ICA(1:min_len);

corr_coef_ICA = corrcoef(fecg_ground_truth, fECG_reconstructed_ICA');
corr_ica = corr_coef_ICA(1,2);

corr_coef_SVD = corrcoef(fecg_ground_truth, fECG_reconstructed_SVD);
corr_svd = corr_coef_SVD(1,2);

fprintf('Correlation coefficient between ground truth and ICA reconstructed fECG: %.4f\n', corr_ica);
fprintf('Correlation coefficient between ground truth and SVD reconstructed fECG: %.4f\n', corr_svd);
