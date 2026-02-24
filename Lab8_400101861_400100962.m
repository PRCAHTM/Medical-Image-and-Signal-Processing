%% section 1
img = imread("S2_Q1_utils/t2.jpg");
figure;
imshow(img);
title('S2-Q1-utils/t2.jpg', 'Interpreter', 'latex');

img_d = img(:,:,1);

noiseVariance = 0.015;
noiseMean = 0;
img_d_noisy = imnoise(img_d, 'gaussian', noiseMean, noiseVariance);

k_size = 4;
[X, Y] = ndgrid(-size(img,1)/2:size(img,1)/2-1, -size(img,2)/2:size(img,2)/2-1);
kernel = abs(X) < k_size & abs(Y) < k_size;
imshow(kernel);
title(['Binary Kernel, $k$ size = ', num2str(k_size)], 'Interpreter', 'latex');

kernel = kernel / sum(kernel(:));

kernel_fft = fftshift(fft2(kernel));
img_d_noisy_fft = fftshift(fft2(img_d_noisy));

filtered_fft = kernel_fft .* img_d_noisy_fft;
img_filtered = fftshift(abs(ifft2(filtered_fft)));
img_filtered = img_filtered / max(img_filtered(:));

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
imshow(img_d);
title('Original Image - First Slice', 'Interpreter', 'latex');

subplot(2,2,2);
imshow(img_d_noisy);
title('Image with Added Gaussian Noise', 'Interpreter', 'latex');

subplot(2,2,3);
imshow(img_filtered);
title('Filtered Using Rectangular Kernel', 'Interpreter', 'latex');

img_filtered_2 = imgaussfilt(img_d, 1);
subplot(2,2,4);
imshow(img_filtered_2);
title('Filtered Using Gaussian Kernel on Original Image', 'Interpreter', 'latex');

%% section 2 part 1
clc;clear;close all;
image = imread("S2_Q2_utils/t2.jpg");
image=image(:,:,1);
f = double(image);
h = Gaussian(1, [256 256]);
g=conv2(f,h,'same');
figure()
subplot(1,2,1)
imshow(f/255)
title('original')
subplot(1,2,2)
imshow(g/max(g,[],'all'))
title('blurred')

G_2=fftshift(fft2(ifftshift((g))));
H=fftshift(fft2(ifftshift((h))));
F=G_2./H;
f_hat=fftshift(ifft2(ifftshift(F)));
figure()
subplot(1,3,1)
imshow(f/255)
title('original')

subplot(1,3,2)
imshow(g/max(g,[],'all'))
title('blurred')

subplot(1,3,3)
imshow(abs(f_hat)/max(f_hat,[],'all'));
title('deblur')
g_noise=g+randn(size(g))*0.001;
figure()
subplot(1,2,1)
imshow(f/255)
title('original')
subplot(1,2,2)
imshow(g_noise/max(g_noise,[],'all'))
title('blurred noisy')
G_2=fftshift(fft2(ifftshift((g_noise))));
H=fftshift(fft2(ifftshift((h))));
F=G_2./H;
f_hat=fftshift(ifft2(ifftshift(F)));
f_hat=f_hat-min(f_hat,[],'all');
figure()
subplot(1,2,1)
imshow(f/255)
title('original')

subplot(1,2,2)
imshow(f_hat/max(f_hat,[],'all'));
title('deblur noisy')
%% section 3
clc; clear; close all;

image = imread("S2_Q2_utils/t2.jpg");
img_resized = imresize(image(:,:,1), [64, 64], 'bilinear');
f = double(img_resized);

h = [0 1 0; 1 2 1; 0 1 0];
K = zeros(64);
D = zeros(64*64);

K(1:3, 1:3) = h;
for count = 1:64*64
    [r, c] = ind2sub([64, 64], count);
    temp = circshift(K, [r-1, c-1]);
    D(count, :) = temp(:)';
end

figure();
spy(D);

g = D * f(:);
g_reshape = reshape(g, 64, 64);

figure();
subplot(1, 3, 1);
imshow(g_reshape, []);
title('Blurred', 'Interpreter', 'latex');

g_noise = g_reshape + randn(size(g_reshape)) * 0.05;

subplot(1, 3, 2);
imshow(g_noise, []);
title('Noisy', 'Interpreter', 'latex');

f_hat = pinv(D) * g_noise(:);
f_reconstructed = reshape(f_hat, 64, 64);

subplot(1, 3, 3);
imshow(f_reconstructed, []);
title('Reconstructed', 'Interpreter', 'latex');

%% section 4
beta = 0.01;
f_k = zeros(64*64, 1);

for i = 1:50
    f_k = f_k + beta * D' * (g - D * f_k);
end

figure;
imshow(reshape(f_k, 64, 64), []);
title('Reconstructed with Gradient', 'Interpreter', 'latex');

%% functions
function g = Gaussian(sigma, dims)

	rows = dims(1);
	cols = dims(2);
    slices = 1;
    D = 2;
    if length(dims)>2
        slices = dims(3);
        D = 3;
    end
    
	cr = ceil( (rows-1)/2 ) + 1;
	cc = ceil( (cols-1)/2 ) + 1;
    cs = ceil( (slices-1)/2) + 1;
    
    a = 1 / (2*sigma^2);
	g = zeros(rows,cols,slices);

    for s = 1:slices
        for c = 1:cols
            for r = 1:rows
                r_sh = r - cr;
                c_sh = c - cc;
                s_sh = s - cs;
                g(r,c,s) = exp( -a * (r_sh^2 + c_sh^2 + s_sh^2) );
            end
        end
    end
    
    g = g / (sqrt(2*pi)*sigma)^D;
    
end