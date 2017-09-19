% function Labels = cluster()

img_dir = '/Users/kabir/Desktop/Research/MTP/MTP_data/';
img_file = strcat(img_dir,'DWI_b800.mat');
mask_dir = '/Users/kabir/Desktop/Research/MTP/';
feature_dir = '/Users/kabir/Desktop/Research/MTP/Figures/finalGLCM/';
feature_fil = strcat(feature_dir,'Linear_4_7_1_1.jpg');
predict_dir = '/Users/kabir/Desktop/OSData/01-/';
predict_fil = strcat(predict_dir,'DWI.nii');
predict_img = extract_img(predict_fil);
predict_mask_fil = strcat(predict_dir,'tumorbaseline.nii');
predict_mask = extract_img(predict_mask_fil);

mask_file = strcat(mask_dir,'mask1.nii');
mask_img = extract_img(mask_file);




input_img = extract_img(img_file);
mask_img = extract_img(mask_file);
sz = size(input_img);
ImageInColumnFormat = reshape(input_img,[],1);
MaskInColumnFormat = reshape(mask_img,[],1);

ImageInColumnFormat2 = reshape(predict_img,[],1);

Num_clusters = 4;
[~,U] = fcm(double(ImageInColumnFormat),Num_clusters);
[Values,Labels] = max(U,[],1);
LabelsInImageFormat = reshape(Labels, sz(1),sz(2));

[~,U2] = fcm(double(ImageInColumnFormat2),Num_clusters);
[Values2,Labels2] = max(U2,[],1);
LabelsInImageFormat2 = reshape(Labels2, sz(1),sz(2));


SVMModel = fitcsvm([ImageInColumnFormat Labels.'],MaskInColumnFormat,'KernelFunction','rbf','ClassNames',[0,1]);
if not(SVMModel.ConvergenceInfo.Converged)
    error('SVM training did not reach convergence')
end

[prlabels,~] = predict(SVMModel,[ImageInColumnFormat2 Labels2.']);

final_predict = reshape(prlabels, sz(1),sz(2));
imshow(mat2gray(final_predict));




%     figure; imagesc(LabelsInImageFormat)
% end
