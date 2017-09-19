%%%%%%%%%%%%  Parametric maps texture analysis  %%%%%%%%%%%%

close all;
clear all;

b=[0,10,20,30,40,50,80,100,200,400,800];

offsets0 = [zeros(10,1) (1:10)'];

% offsets1 = [ 0 1; -1 1; -1 0; -1 -1];
% offsets2 = [ 0 2; -2 2; -2 0; -2 -2]; 
% offsets3 = [ 0 3; -3 3; -3 0; -3 -3]; 
% ffsets4 = [ 0 4; -4 4; -4 0; -4 -4];
offsets5 = [ 0 5; -5 5; -5 0; -5 -5];

offsets = [offsets5];

total_feature = 28;

total_direc = 4;      % total direction -0, 45,  90,  135 degree
total_d = length(offsets)/total_direc;     %total distance - 1,2,3,4,5 in pixcel


% read  patient data in baseline & follow-op
baseline_data=load('D:\OS result\25-prasu\25-prasu latest\prasu_1.mat', 'allslice', 'ADCmap200_lsqall', 'Dmap_tvall','Dpmap_tvall', 'fmap_tvall', 'tumor','normal', 'start_slice', 'end_slice');
folup_data1=load('D:\OS result\25-prasu\25-prasu latest\prasu_2.mat', 'allslice', 'ADCmap200_lsqall', 'Dmap_tvall','Dpmap_tvall', 'fmap_tvall', 'tumor','start_slice', 'end_slice');
folup_data2=load('D:\OS result\25-prasu\25-prasu latest\prasu_3.mat', 'allslice', 'ADCmap200_lsqall', 'Dmap_tvall','Dpmap_tvall', 'fmap_tvall', 'tumor','start_slice', 'end_slice');




%%%%%  change for baseline proessing
% baseline_data=load('C:\Users\dalim\Documents\MATLAB\OS result\divyansh\divyansh_1.mat', 'allslice',  'Dmap_tvall','Dpmap_tvall', 'fmap_tvall', 'ADCmap200_lsqall');
% tumor=load_untouch_niigz('C:\Users\DELL\Documents\MATLAB\OS\15-Chotu-UHID-102081675\1-chotu-17-10-2016\tumor 1.nii.gz');
%   normal=load_untouch_niigz('D:\OS\03-Pappu-UHID-101646593\0-Pappu_29-04-2016\normal.nii.gz');
%  


%%%%%%%%%%%%%%   Process baseline data 

allslice = baseline_data.allslice;
soft_tss_tumor = baseline_data.tumor;
% % normal = baseline_data.normal;
% tslice= baseline_data.start_slice : baseline_data.end_slice;
% nslice=baseline_data.nslice;
Dmap_0 = baseline_data.Dmap_tvall;
Dpmap_0 = baseline_data.Dpmap_tvall;
fmap_0 = baseline_data.fmap_tvall;
adcmap_0 = baseline_data.ADCmap200_lsqall;

%%% find the slice no. on which ROI mask is present
count =1;
for i=1:size(soft_tss_tumor.img,3)
  if(max(max(soft_tss_tumor.img(:,:,i)))==1)
      tslice(count) = i;
      count=count+1;
  end
end  
% 
% 
% 
% 
% % count =1;
% % for i=1:size(normal.img,3)
% %   if(max(max(normal.img(:,:,i)))==1)
% %       nslice(count) = i;
% %       count=count+1;
% %   end
% % end  

%get image dimention
 [row,col] = size(soft_tss_tumor.img(:,:,1));

  D_feature_tumor = zeros(total_feature, length(tslice));
   Dp_feature_tumor  = zeros(total_feature, length(tslice));
    f_feature_tumor  = zeros(total_feature, length(tslice));
      adc_feature_tumor  = zeros(total_feature, length(tslice));
      b800_feature_tumor  = zeros(total_feature, length(tslice));
      
% %         D_feature_normal = zeros(total_feature, length(nslice));
% %    Dp_feature_normal = zeros(total_feature, length(nslice));
% %     f_feature_normal = zeros(total_feature, length(nslice));
% %       adc_feature_normal = zeros(total_feature, length(nslice));
% %       b800_feature_normal = zeros(total_feature, length(nslice));



 %%%%%%%%%%  Analisis on tumour vilume %%%%%%%%%%%%%%%%%

for i=1:length(tslice)
    
D_img= double(Dmap_0(:,:,tslice(i)));
Dp_img= double(Dpmap_0(:,:,tslice(i)));
f_img= double(fmap_0(:,:,tslice(i)));
adc_img= double(adcmap_0(:,:,tslice(i)));
b800_img= double(allslice.img(:,:,tslice(i),11));

% figure, imagesc(img);

tmask= double(soft_tss_tumor.img(:,:,tslice(i)));
% figure, imagesc(tmask);

D_t=immultiply(D_img,tmask);
Dp_t=immultiply(Dp_img,tmask);
f_t=immultiply(f_img,tmask);
adc_t=immultiply(adc_img,tmask);
b800_t=immultiply(b800_img,tmask);

% figure, imagesc(t);

%%%%%  normalize masked image

maxtroi=max(D_t(:));
  mintroi = min(D_t(:));
  D_norm_troi= ((D_t-mintroi)/(maxtroi-mintroi+0.0001)).*255;
  
 maxtroi=max(Dp_t(:));
  mintroi = min(Dp_t(:));
  Dp_norm_troi= ((Dp_t-mintroi)/(maxtroi-mintroi+0.0001)).*255; 
  
  maxtroi=max(f_t(:));
  mintroi = min(f_t(:));
  f_norm_troi= ((f_t-mintroi)/(maxtroi-mintroi+0.0001)).*255;
  
  maxtroi=max(adc_t(:));
  mintroi = min(adc_t(:));
  adc_norm_troi= ((adc_t-mintroi)/(maxtroi-mintroi+0.0001)).*255;
  
 maxtroi=max(b800_t(:));
  mintroi = min(b800_t(:));
 b800_norm_troi= ((b800_t-mintroi)/(maxtroi-mintroi+0.0001)).*255;


D_norm_troi(D_norm_troi==0) = NaN;   %%  norm_troi(find(~norm_troi)) = NaN;
Dp_norm_troi(Dp_norm_troi==0) = NaN;
f_norm_troi(f_norm_troi==0) = NaN;
adc_norm_troi(adc_norm_troi==0) = NaN;
b800_norm_troi(b800_norm_troi==0) = NaN;
%  figure, imagesc(norm_troi);

D_norm_troi_all(:,:,i) = D_norm_troi;
Dp_norm_troi_all(:,:,i) = Dp_norm_troi;
f_norm_troi_all(:,:,i) = f_norm_troi;
adc_norm_troi_all(:,:,i) = adc_norm_troi;
b800_norm_troi_all(:,:,i) =b800_norm_troi;
 
 %%%%  GLCM 

  warning('off','Images:graycomatrix:scaledImageContainsNan');

  [D_glcm_tumor, It] = graycomatrix(D_norm_troi,'NumLevels',256,'Offset',offsets,'G',[]);
D_stats_t = graycoprops((D_glcm_tumor),'Contrast Correlation Energy Homogeneity');

[Dp_glcm_tumor, It] = graycomatrix(Dp_norm_troi,'NumLevels',256,'Offset',offsets,'G',[]);
Dp_stats_t = graycoprops((Dp_glcm_tumor),'Contrast Correlation Energy Homogeneity');

[f_glcm_tumor, It] = graycomatrix(f_norm_troi,'NumLevels',256,'Offset',offsets,'G',[]);
f_stats_t = graycoprops((f_glcm_tumor),'Contrast Correlation Energy Homogeneity');

[adc_glcm_tumor, It] = graycomatrix(adc_norm_troi,'NumLevels',256,'Offset',offsets,'G',[]);
adc_stats_t = graycoprops((adc_glcm_tumor),'Contrast Correlation Energy Homogeneity');

[b800_glcm_tumor, It] = graycomatrix(b800_norm_troi,'NumLevels',256,'Offset',offsets,'G',[]);
b800_stats_t = graycoprops((b800_glcm_tumor),'Contrast Correlation Energy Homogeneity');

 warning('on','Images:graycomatrix:scaledImageContainsNan');
 
% figure, imagesc(It);
% texture(i,1)= mean(stats_t.Contrast);
% texture(i,2)= mean(stats_t.Correlation);
% texture(i,3) = mean(stats_t.Energy);
% texture(i,4) = mean(stats_t.Homogeneity);

%  stat_tumor = GLCM_Features4(glcm_tumor,0);

% feature(1:20, i) = stats_t.Contrast';
% feature(21:40, i) = stats_t.Correlation';
% feature(41:60, i) = stats_t.Energy';
% feature(61:80, i) = stats_t.Homogeneity';


%% store the GLCM features after averaging values in all direction(0,45,90,135) at
%% each  distance (1,2,3,4,5 pixcel)


for x=1: total_d  

D_feature_tumor (0*total_d+x , i) = mean(D_stats_t.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
D_feature_tumor (1*total_d+x , i) = mean(D_stats_t.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
D_feature_tumor (2*total_d+x , i) = mean(D_stats_t.Energy(((x-1)*total_direc)+1 : x*total_direc ));
D_feature_tumor (3*total_d+x , i) = mean(D_stats_t.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));

Dp_feature_tumor (0*total_d+x , i) = mean(Dp_stats_t.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
Dp_feature_tumor (1*total_d+x , i) = mean(Dp_stats_t.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
Dp_feature_tumor (2*total_d+x , i) = mean(Dp_stats_t.Energy(((x-1)*total_direc)+1 : x*total_direc ));
Dp_feature_tumor (3*total_d+x , i) = mean(Dp_stats_t.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));

f_feature_tumor (0*total_d+x , i) = mean(f_stats_t.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
f_feature_tumor (1*total_d+x , i) = mean(f_stats_t.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
f_feature_tumor (2*total_d+x , i) = mean(f_stats_t.Energy(((x-1)*total_direc)+1 : x*total_direc ));
f_feature_tumor (3*total_d+x , i) = mean(f_stats_t.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));

adc_feature_tumor (0*total_d+x , i) = mean(adc_stats_t.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
adc_feature_tumor (1*total_d+x , i) = mean(adc_stats_t.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
adc_feature_tumor (2*total_d+x , i) = mean(adc_stats_t.Energy(((x-1)*total_direc)+1 : x*total_direc ));
adc_feature_tumor (3*total_d+x , i) = mean(adc_stats_t.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));

b800_feature_tumor (0*total_d+x , i) = mean(b800_stats_t.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
b800_feature_tumor (1*total_d+x , i) = mean(b800_stats_t.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
b800_feature_tumor (2*total_d+x , i) = mean(b800_stats_t.Energy(((x-1)*total_direc)+1 : x*total_direc ));
b800_feature_tumor (3*total_d+x , i) = mean(b800_stats_t.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));

end

%%%%  Run Length Matrix

  [D_GLRLMS_t,SI_t] = grayrlmatrix(D_norm_troi,'NumLevels',256,'G',[]);
  D_stats_rlm_t = grayrlprops(D_GLRLMS_t);
  
    [Dp_GLRLMS_t,SI_t] = grayrlmatrix(Dp_norm_troi,'NumLevels',256,'G',[]);
  Dp_stats_rlm_t = grayrlprops(Dp_GLRLMS_t);
  
    [f_GLRLMS_t,SI_t] = grayrlmatrix(f_norm_troi,'NumLevels',256,'G',[]);
  f_stats_rlm_t = grayrlprops(f_GLRLMS_t);
  
  [adc_GLRLMS_t,SI_t] = grayrlmatrix(adc_norm_troi,'NumLevels',256,'G',[]);
  adc_stats_rlm_t = grayrlprops(adc_GLRLMS_t);
  
  [b800_GLRLMS_t,SI_t] = grayrlmatrix(b800_norm_troi,'NumLevels',256,'G',[]);
  b800_stats_rlm_t = grayrlprops(b800_GLRLMS_t);
%  figure, imagesc(SI_t);

% D_feature(21, i) = mean(D_stats_rlm_t(:,1));      % Short Run Emphasis (SRE)
% D_feature(22, i) = mean(D_stats_rlm_t(:,2));      % Long Run Emphasis (LRE)
% D_feature(23, i) = mean(D_stats_rlm_t(:,3));      % Gray-Level Nonuniformity (GLN)
% D_feature(24, i) = mean(D_stats_rlm_t(:,4));      % Run Length Nonuniformity (RLN)
% D_feature(25, i) = mean(D_stats_rlm_t(:,5));     % Run Percentage (RP)
% D_feature(26, i) = mean(D_stats_rlm_t(:,6));    % Low Gray-Level Run Emphasis (LGRE) 
% D_feature(27, i) = mean(D_stats_rlm_t(:,7));    % High Gray-Level Run Emphasis (HGRE)
% D_feature(28, i) = mean(D_stats_rlm_t(:,8));    % Short Run Low Gray-Level Emphasis (SRLGE)
% D_feature(29, i) = mean(D_stats_rlm_t(:,9));    % Short Run High Gray-Level Emphasis (SRHGE)
% D_feature(30, i) = mean(D_stats_rlm_t(:,10));   % Long Run Low Gray-Level Emphasis (LRLGE)
% D_feature(31, i) = mean(D_stats_rlm_t(:,11));   % Long Run High Gray-Level Emphasis (LRHGE)

w=total_d*total_direc;
for j=1:11

D_feature_tumor (w+j, i) = mean(D_stats_rlm_t(:,j));     
Dp_feature_tumor (w+j, i) = mean(Dp_stats_rlm_t(:,j));     
f_feature_tumor (w+j, i) = mean(f_stats_rlm_t(:,j));      
adc_feature_tumor (w+j, i) = mean(adc_stats_rlm_t(:,j)); 
b800_feature_tumor (w+j, i) = mean(b800_stats_rlm_t(:,j)); 

end


% histogra analysis at each slice

D_histo_result_tumor=histogram_analysis(D_norm_troi);
D_feature_tumor(16:22,i) = D_histo_result_tumor;

Dp_histo_result_tumor=histogram_analysis(Dp_norm_troi);
Dp_feature_tumor(16:22,i) = Dp_histo_result_tumor;

f_histo_result_tumor=histogram_analysis(f_norm_troi);
f_feature_tumor(16:22,i) = f_histo_result_tumor;

adc_histo_result_tumor=histogram_analysis(adc_norm_troi);
adc_feature_tumor(16:22,i) = adc_histo_result_tumor;

b800_histo_result_tumor=histogram_analysis(b800_norm_troi);
b800_feature_tumor(16:22,i) = b800_histo_result_tumor;

%gradient image analysis at each slice

D_grdimg_result_tumor=gradientimg_analysis(D_norm_troi);
D_feature_tumor(23:28,i) = D_grdimg_result_tumor;

Dp_grdimg_result_tumor=gradientimg_analysis(Dp_norm_troi);
Dp_feature_tumor(23:28,i) = Dp_grdimg_result_tumor;

f_grdimg_result_tumor=gradientimg_analysis(f_norm_troi);
f_feature_tumor(23:28,i) = f_grdimg_result_tumor;

adc_grdimg_result_tumor=gradientimg_analysis(adc_norm_troi);
adc_feature_tumor(23:28,i) = adc_grdimg_result_tumor;

b800_grdimg_result_tumor=gradientimg_analysis(b800_norm_troi);
b800_feature_tumor(23:28,i) = b800_grdimg_result_tumor;


end

 %%%%%%%%%%  Analisis on normal volume %%%%%%%%%%%%%%%%%

% % for i=1:length(nslice)
% %     
% % D_img= double(Dmap_0(:,:,nslice(i)));
% % Dp_img= double(Dpmap_0(:,:,nslice(i)));
% % f_img= double(fmap_0(:,:,nslice(i)));
% % adc_img= double(adcmap_0(:,:,nslice(i)));
% % b800_img= double(allslice.img(:,:,nslice(i),11));
% % 
% % % figure, imagesc(img);
% % 
% % nmask= double(normal.img(:,:,nslice(i)));
% % % figure, imagesc(tmask);
% % 
% % D_n=immultiply(D_img,nmask);
% % Dp_n=immultiply(Dp_img,nmask);
% % f_n=immultiply(f_img,nmask);
% % adc_n=immultiply(adc_img,nmask);
% % b800_n=immultiply(b800_img,nmask);
% % 
% % % figure, imagesc(t);
% % 
% % %%%%%  normalize masked image
% % 
% % maxnroi=max(D_n(:));
% %   minnroi = min(D_n(:));
% %   D_norm_nroi= ((D_n-minnroi)/(maxnroi-minnroi+1e-10)).*255;
% %   
% %  maxnroi=max(Dp_n(:));
% %   minnroi = min(Dp_n(:));
% %   Dp_norm_nroi= ((Dp_n-minnroi)/(maxnroi-minnroi+1e-10)).*255; 
% %   
% %   maxnroi=max(f_n(:));
% %   minnroi = min(f_n(:));
% %   f_norm_nroi= ((f_n-minnroi)/(maxnroi-minnroi+1e-10)).*255;
% %   
% %   maxnroi=max(adc_n(:));
% %   minnroi = min(adc_n(:));
% %   adc_norm_nroi= ((adc_n-minnroi)/(maxnroi-minnroi+1e-10)).*255;
% %   
% %  maxnroi=max(b800_n(:));
% %   minnroi = min(b800_n(:));
% %  b800_norm_nroi= ((b800_n-minnroi)/(maxnroi-minnroi+1e-10)).*255;
% % 
% % 
% % D_norm_nroi(D_norm_nroi==0) = NaN;   %%  norm_troi(find(~norm_troi)) = NaN;
% % Dp_norm_nroi(Dp_norm_nroi==0) = NaN;
% % f_norm_nroi(f_norm_nroi==0) = NaN;
% % adc_norm_nroi(adc_norm_nroi==0) = NaN;
% % b800_norm_nroi(b800_norm_nroi==0) = NaN;
% % %  figure, imagesc(norm_troi);
% %  
% %  %%%%  GLCM 
% % 
% %   warning('off','Images:graycomatrix:scaledImageContainsNan');
% % 
% %   [D_glcm_normal, It] = graycomatrix(D_norm_nroi,'NumLevels',256,'Offset',offsets,'G',[]);
% % D_stats_n = graycoprops((D_glcm_normal),'Contrast Correlation Energy Homogeneity');
% % 
% % [Dp_glcm_normal, It] = graycomatrix(Dp_norm_nroi,'NumLevels',256,'Offset',offsets,'G',[]);
% % Dp_stats_n = graycoprops((Dp_glcm_normal),'Contrast Correlation Energy Homogeneity');
% % 
% % [f_glcm_normal, It] = graycomatrix(f_norm_nroi,'NumLevels',256,'Offset',offsets,'G',[]);
% % f_stats_n = graycoprops((f_glcm_normal),'Contrast Correlation Energy Homogeneity');
% % 
% % [adc_glcm_normal, It] = graycomatrix(adc_norm_nroi,'NumLevels',256,'Offset',offsets,'G',[]);
% % adc_stats_n = graycoprops((adc_glcm_normal),'Contrast Correlation Energy Homogeneity');
% % 
% % [b800_glcm_normal, It] = graycomatrix(b800_norm_nroi,'NumLevels',256,'Offset',offsets,'G',[]);
% % b800_stats_n = graycoprops((b800_glcm_normal),'Contrast Correlation Energy Homogeneity');
% % 
% %  warning('on','Images:graycomatrix:scaledImageContainsNan');
% %  
% % % figure, imagesc(It);
% % % texture(i,1)= mean(stats_t.Contrast);
% % % texture(i,2)= mean(stats_t.Correlation);
% % % texture(i,3) = mean(stats_t.Energy);
% % % texture(i,4) = mean(stats_t.Homogeneity);
% % 
% % %  stat_tumor = GLCM_Features4(glcm_tumor,0);
% % 
% % % feature(1:20, i) = stats_t.Contrast';
% % % feature(21:40, i) = stats_t.Correlation';
% % % feature(41:60, i) = stats_t.Energy';
% % % feature(61:80, i) = stats_t.Homogeneity';
% % 
% % 
% % %% store the GLCM features after averaging values in all direction(0,45,90,135) at
% % %% each  distance (1,2,3,4,5 pixcel)
% % 
% % 
% % for x=1: total_d  
% % 
% % D_feature_normal(0*total_d+x , i) = mean(D_stats_n.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
% % D_feature_normal(1*total_d+x , i) = mean(D_stats_n.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
% % D_feature_normal(2*total_d+x , i) = mean(D_stats_n.Energy(((x-1)*total_direc)+1 : x*total_direc ));
% % D_feature_normal(3*total_d+x , i) = mean(D_stats_n.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));
% % 
% % Dp_feature_normal(0*total_d+x , i) = mean(Dp_stats_n.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
% % Dp_feature_normal(1*total_d+x , i) = mean(Dp_stats_n.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
% % Dp_feature_normal(2*total_d+x , i) = mean(Dp_stats_n.Energy(((x-1)*total_direc)+1 : x*total_direc ));
% % Dp_feature_normal(3*total_d+x , i) = mean(Dp_stats_n.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));
% % 
% % f_feature_normal(0*total_d+x , i) = mean(f_stats_n.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
% % f_feature_normal(1*total_d+x , i) = mean(f_stats_n.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
% % f_feature_normal(2*total_d+x , i) = mean(f_stats_n.Energy(((x-1)*total_direc)+1 : x*total_direc ));
% % f_feature_normal(3*total_d+x , i) = mean(f_stats_n.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));
% % 
% % adc_feature_normal(0*total_d+x , i) = mean(adc_stats_n.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
% % adc_feature_normal(1*total_d+x , i) = mean(adc_stats_n.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
% % adc_feature_normal(2*total_d+x , i) = mean(adc_stats_n.Energy(((x-1)*total_direc)+1 : x*total_direc ));
% % adc_feature_normal(3*total_d+x , i) = mean(adc_stats_n.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));
% % 
% % b800_feature_normal(0*total_d+x , i) = mean(b800_stats_n.Contrast(((x-1)*total_direc)+1 : x*total_direc ));
% % b800_feature_normal(1*total_d+x , i) = mean(b800_stats_n.Correlation(((x-1)*total_direc)+1 : x*total_direc ));
% % b800_feature_normal(2*total_d+x , i) = mean(b800_stats_n.Energy(((x-1)*total_direc)+1 : x*total_direc ));
% % b800_feature_normal(3*total_d+x , i) = mean(b800_stats_n.Homogeneity(((x-1)*total_direc)+1 : x*total_direc ));
% % 
% % end
% % 
% % %%%%  Run Length Matrix
% % 
% %   [D_GLRLMS_n,SI_n] = grayrlmatrix(D_norm_nroi,'NumLevels',256,'G',[]);
% %   D_stats_rlm_n = grayrlprops(D_GLRLMS_n);
% %   
% %     [Dp_GLRLMS_n,SI_n] = grayrlmatrix(Dp_norm_nroi,'NumLevels',256,'G',[]);
% %   Dp_stats_rlm_n = grayrlprops(Dp_GLRLMS_n);
% %   
% %     [f_GLRLMS_n,SI_n] = grayrlmatrix(f_norm_nroi,'NumLevels',256,'G',[]);
% %   f_stats_rlm_n = grayrlprops(f_GLRLMS_n);
% %   
% %   [adc_GLRLMS_n,SI_n] = grayrlmatrix(adc_norm_nroi,'NumLevels',256,'G',[]);
% %   adc_stats_rlm_n = grayrlprops(adc_GLRLMS_n);
% %   
% %   [b800_GLRLMS_n,SI_n] = grayrlmatrix(b800_norm_nroi,'NumLevels',256,'G',[]);
% %   b800_stats_rlm_n = grayrlprops(b800_GLRLMS_n);
% % %  figure, imagesc(SI_t);
% % 
% % % D_feature(21, i) = mean(D_stats_rlm_t(:,1));      % Short Run Emphasis (SRE)
% % % D_feature(22, i) = mean(D_stats_rlm_t(:,2));      % Long Run Emphasis (LRE)
% % % D_feature(23, i) = mean(D_stats_rlm_t(:,3));      % Gray-Level Nonuniformity (GLN)
% % % D_feature(24, i) = mean(D_stats_rlm_t(:,4));      % Run Length Nonuniformity (RLN)
% % % D_feature(25, i) = mean(D_stats_rlm_t(:,5));     % Run Percentage (RP)
% % % D_feature(26, i) = mean(D_stats_rlm_t(:,6));    % Low Gray-Level Run Emphasis (LGRE) 
% % % D_feature(27, i) = mean(D_stats_rlm_t(:,7));    % High Gray-Level Run Emphasis (HGRE)
% % % D_feature(28, i) = mean(D_stats_rlm_t(:,8));    % Short Run Low Gray-Level Emphasis (SRLGE)
% % % D_feature(29, i) = mean(D_stats_rlm_t(:,9));    % Short Run High Gray-Level Emphasis (SRHGE)
% % % D_feature(30, i) = mean(D_stats_rlm_t(:,10));   % Long Run Low Gray-Level Emphasis (LRLGE)
% % % D_feature(31, i) = mean(D_stats_rlm_t(:,11));   % Long Run High Gray-Level Emphasis (LRHGE)
% % 
% % w=total_d*total_direc;
% % for j=1:11
% % 
% % D_feature_normal(w+j, i) = mean(D_stats_rlm_n(:,j));     
% % Dp_feature_normal(w+j, i) = mean(Dp_stats_rlm_n(:,j));     
% % f_feature_normal(w+j, i) = mean(f_stats_rlm_n(:,j));      
% % adc_feature_normal(w+j, i) = mean(adc_stats_rlm_n(:,j)); 
% % b800_feature_normal(w+j, i) = mean(b800_stats_rlm_n(:,j)); 
% % 
% % end
% % 
% % % histogra analysis at each slice
% % 
% % D_histo_result_normal=histogram_analysis(D_norm_nroi);
% % D_feature_normal(16:22,i) = D_histo_result_normal;
% % 
% % Dp_histo_result_normal=histogram_analysis(Dp_norm_nroi);
% % Dp_feature_normal(16:22,i) = Dp_histo_result_normal;
% % 
% % f_histo_result_normal=histogram_analysis(f_norm_nroi);
% % f_feature_normal(16:22,i) = f_histo_result_normal;
% % 
% % adc_histo_result_normal=histogram_analysis(adc_norm_nroi);
% % adc_feature_normal(16:22,i) = adc_histo_result_normal;
% % 
% % b800_histo_result_normal=histogram_analysis(b800_norm_nroi);
% % b800_feature_normal(16:22,i) = b800_histo_result_normal;
% % 
% % 
% % %gradient image analysis at each slice
% % 
% % 
% % D_grdimg_result_normal=gradientimg_analysis(D_norm_nroi);
% % D_feature_normal(23:28,i) = D_grdimg_result_normal;
% % 
% % Dp_grdimg_result_normal=gradientimg_analysis(Dp_norm_nroi);
% % Dp_feature_normal(23:28,i) = Dp_grdimg_result_normal;
% % 
% % f_grdimg_result_normal=gradientimg_analysis(f_norm_nroi);
% % f_feature_normal(23:28,i) = f_grdimg_result_normal;
% % 
% % adc_grdimg_result_normal=gradientimg_analysis(adc_norm_nroi);
% % adc_feature_normal(23:28,i) = adc_grdimg_result_normal;
% % 
% % b800_grdimg_result_normal=gradientimg_analysis(b800_norm_nroi);
% % b800_feature_normal(23:28,i) = b800_grdimg_result_normal;
% % 
% % end

%%--------only for base-lie data texture analysis

%%%%%%% claculate mean value for all feature for all slice -tumour volume %%%%%%%%%%%%%%

D_mean_feature_tumor = zeros(total_feature,1);
Dp_mean_feature_tumor = zeros(total_feature,1);
f_mean_feature_tumor = zeros(total_feature,1);
 adc_mean_feature_tumor = zeros(total_feature,1);
b800_mean_feature_tumor = zeros(total_feature,1);

for i= 1:total_feature
    D_mean_feature_tumor(i,1) = nanmean(D_feature_tumor(i,:));
    
    Dp_mean_feature_tumor(i,1) = nanmean(Dp_feature_tumor(i,:));
        
    f_mean_feature_tumor(i,1) = nanmean(f_feature_tumor(i,:));
 
    adc_mean_feature_tumor(i,1) = nanmean(adc_feature_tumor(i,:));
  
    b800_mean_feature_tumor(i,1) = nanmean(b800_feature_tumor(i,:));
    
end    