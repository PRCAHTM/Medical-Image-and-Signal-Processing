%% section Q1 part1

image_path = "S3_Q1_utils/thorax_t1.jpg";
image = imread(image_path);

image = image(:, :, 1);

I = double(image);

imshow(I, []);
title('Input Image');

center1 = [90, 95];
center2 = [90, 175];
radius = 30;

flag = zeros(size(I));

[row, col] = size(I);
for i = 1:row
    for j = 1:col
        dist1 = sqrt((i - center1(1))^2 + (j - center1(2))^2);
        dist2 = sqrt((i - center2(1))^2 + (j - center2(2))^2);
        if (dist1 < radius || dist2 < radius) && I(i, j) <= 50
            flag(i, j) = 1;
        end
    end
end

Overlay(I, flag);

%% section Q1 part2

image_path = "S3_Q1_utils/thorax_t1.jpg";
image = imread(image_path);

image = image(:, :, 1);

I = double(image);

imshow(I, []);
title('Input Image');

center1 = [180, 88];
center2 = [140, 88];
center3 = [150, 130];

radius = 32;

flag = zeros(size(I));

[row, col] = size(I);
for i = 1:row
    for j = 1:col
        dist1 = sqrt((i - center1(1))^2 + (j - center1(2))^2);
        dist2 = sqrt((i - center2(1))^2 + (j - center2(2))^2);
        dist3 = sqrt((i - center3(1))^2 + (j - center3(2))^2);

        if (dist1 < radius) && (I(i, j) >= 85) && (I(i, j) <= 135)
            flag(i, j) = 1;
        end
        if (dist2 < radius) && (I(i, j) >= 70) && (I(i, j) <= 200)
            flag(i, j) = 1;
        end
        if (dist3 < radius) && (I(i, j) >= 70) && (I(i, j) <= 120)
            flag(i, j) = 1;
        end
    end
end

Overlay(I, flag);


%% section Q2

image_t1 = imread("S3_Q2_utils/t1.jpg");
image_t2 = imread("S3_Q2_utils/t2.jpg");
image_pd = imread("S3_Q2_utils/pd.jpg");

if size(image_t1, 3) == 3
    image_t1 = rgb2gray(image_t1);
end
if size(image_t2, 3) == 3
    image_t2 = rgb2gray(image_t2);
end
if size(image_pd, 3) == 3
    image_pd = rgb2gray(image_pd);
end

image_t1 = double(image_t1);
image_t2 = double(image_t2);
image_pd = double(image_pd);

[rows, cols] = size(image_t1);
features = [image_t1(:), image_t2(:), image_pd(:)];

num_clusters = 6;
[cluster_idx, cluster_centers] = kmeans(features, num_clusters, 'MaxIter', 1000);

clustered_image = reshape(cluster_idx, rows, cols);

figure;
imagesc(clustered_image);
colormap('jet');
colorbar;
title('K-means Clustering (6 Clusters)');

imwrite(uint8(clustered_image * (255 / num_clusters)), 'clustered_output.jpg');

figure;
for k = 1:num_clusters
    cluster_mask = (clustered_image == k);
    
    subplot(2, 3, k);
    imagesc(cluster_mask);
    colormap('gray');
    title(['Cluster ', num2str(k)]);
    axis off;
end

sgtitle('Clusters from K-means Clustering');
