%% section 4.1
close all
ct = imread("S1_Q4_utils\ct.jpg");
img = ct(:, :, 1);

FT_img = fft2(img);
FT_img_shifted = fftshift(FT_img);

x_shift = 20;
y_shift = 40;

[rows, cols] = size(img);
[xF, yF] = meshgrid(-cols/2:cols/2-1, -rows/2:rows/2-1);

shift_kernel = exp(-1i * 2 * pi * (xF * x_shift / cols + yF * y_shift / rows));
FT_shifted_img = FT_img_shifted .* shift_kernel;

shifted_img = ifft2(ifftshift(FT_shifted_img));
shifted_img_abs = abs(shifted_img);

figure;
subplot(1, 2, 1);
imshow(img, []);
title('Original Image');

subplot(1, 2, 2);
imshow(shifted_img_abs, []);
title('Shifted Image (Spatial Domain)');

shift_kernel_abs = abs(shift_kernel);
figure;
plot(shift_kernel_abs);
title('Shift Kernel (Real Part)');
%% section 4.2
close all;

rotation_angle = -30;

rotated_image_spatial = imrotate(img, rotation_angle);

figure;
subplot(1, 2, 1);
imshow(img, []);
title('Original Image', 'Interpreter', 'latex');

subplot(1, 2, 2);
imshow(rotated_image_spatial, []);
title('Rotated Image (Spatial Domain)', 'Interpreter', 'latex');

FT_original = fftshift(fft2(img));
FT_rotated_spatial = fftshift(fft2(rotated_image_spatial));

figure;
subplot(1, 2, 1);
imshow(log(abs(FT_original) + 1), []);
title('Original Image Fourier Transform', 'Interpreter', 'latex');

subplot(1, 2, 2);
imshow(log(abs(FT_rotated_spatial) + 1), []);
title('Fourier Transform of Spatially Rotated Image', 'Interpreter', 'latex');

FT_original_shifted = fftshift(fft2(ifftshift(img)));
FT_rotated_FT = imrotate(FT_original_shifted, rotation_angle);

rotated_image_fft = abs(fftshift(ifft2(ifftshift(FT_rotated_FT))));

figure;
subplot(1, 2, 1);
imshow(img, []);
title('Original Image', 'Interpreter', 'latex');

subplot(1, 2, 2);
imshow(rotated_image_fft, []);
title('Rotated Image (via FFT)', 'Interpreter', 'latex');

figure;
subplot(1, 2, 1);
imshow(log(abs(FT_original) + 1), []);
title('Original Image Fourier Transform', 'Interpreter', 'latex');

subplot(1, 2, 2);
imshow(log(abs(FT_rotated_FT) + 1), []);
title('Fourier Transform of Rotated Image (via FFT)', 'Interpreter', 'latex');


%% section 5

clc;
clear;

img = imread("S1_Q5_utils\t1.jpg");
img = im2double(img(:, :, 1));

horizontal_derivative = (circshift(img, 1, 1) - circshift(img, -1, 1)) / 2;
vertical_derivative = (circshift(img, 1, 2) - circshift(img, -1, 2)) / 2;
gradient_magnitude = sqrt(horizontal_derivative.^2 + vertical_derivative.^2);

normalize = @(x) (x - min(x(:))) / (max(x(:)) - min(x(:)));
normalized_img = normalize(img);
normalized_h_derivative = normalize(horizontal_derivative);
normalized_v_derivative = normalize(vertical_derivative);
normalized_gradient_magnitude = normalize(gradient_magnitude);

figure;
subplot(2, 2, 1);
imshow(img, []);
title('Original Image', 'Interpreter', 'latex');

subplot(2, 2, 2);
imshow(horizontal_derivative, []);
title('Horizontal Derivative', 'Interpreter', 'latex');

subplot(2, 2, 3);
imshow(vertical_derivative, []);
title('Vertical Derivative', 'Interpreter', 'latex');

subplot(2, 2, 4);
imshow(gradient_magnitude, []);
title('Gradient Magnitude', 'Interpreter', 'latex');

figure;
subplot(2, 2, 1);
imshow(normalized_img, []);
title('Normalized Original Image', 'Interpreter', 'latex');

subplot(2, 2, 2);
imshow(normalized_h_derivative, []);
title('Normalized Horizontal Derivative', 'Interpreter', 'latex');

subplot(2, 2, 3);
imshow(normalized_v_derivative, []);
title('Normalized Vertical Derivative', 'Interpreter', 'latex');

subplot(2, 2, 4);
imshow(normalized_gradient_magnitude, []);
title('Normalized Gradient Magnitude', 'Interpreter', 'latex');

%% section 6
clc;
close all;

sobel_edges = edge(img, 'sobel');
canny_edges = edge(img, 'canny');

figure;
t = tiledlayout(1, 3);

nexttile;
imshow(img, []);
title('Original Image', 'Interpreter', 'latex');

nexttile;
imshow(sobel_edges, []);
title('Sobel Edge Detection', 'Interpreter', 'latex');

nexttile;
imshow(canny_edges, []);
title('Canny Edge Detection', 'Interpreter', 'latex');


