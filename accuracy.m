% function metrics = accuracy()
data_dir = '/Users/kabir/Desktop/Research/MTP/Figures/finalGLCM/';
img_dir = '/Users/kabir/Desktop/Research/MTP/MTP_data/';
mask_dir = '/Users/kabir/Desktop/Research/MTP/';
feature_dir = '/Users/kabir/Desktop/Research/MTP/Figures/finalGLCM/'

% img_file = strcat(feature_dir,'Linear_4_7_1_1.jpg');
img_file = strcat(img_dir,'DWI_b800.mat');
input_img = extract_img(img_file);


mask_file = strcat(mask_dir,'mask1.nii');
mask_img = extract_img(mask_file);

input_img = input_img/(max(max(input_img)));
otsu_mask = regiongrow(input_img,102,93,0.15);
% 	refined_img  = input_img.*otsu_mask ;
dice_metric = dice(otsu_mask,mask_img);
jaccard_metric = jaccard(otsu_mask, mask_img);
dice_metric
jaccard_metric
metrics = dice_metric;
    
%     imshow(mat2gray(refined_img));
%     [x,y] = getpts;

% function img = extract_image()