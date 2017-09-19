%load image volume
data_dir = '/Users/kabir/Desktop/Research/MTP/OSData/01-/'

DWI = load_untouch_niigz(strcat(data_dir,'DWI.nii.gz'));
T1=load_untouch_niigz(strcat(data_dir,'T1.nii.gz'));
T2=load_untouch_niigz(strcat(data_dir,'T2.nii.gz'));

fixed_volume = DWI.img(:,:,:,11);
moving_T1 = T1.img(:,:,:);
moving_T2_nonfatsat = T2.img(:,:,:,1);
moving_T2_fatsat = T2.img(:,:,:,2);

%spatial referencing information about the resolution and/or location
Rfixed  = imref3d(size(fixed_volume),DWI.hdr.dime.pixdim(1,2),DWI.hdr.dime.pixdim(1,3),DWI.hdr.dime.pixdim(1,4));
Rmoving = imref3d(size(moving_T1),T1.hdr.dime.pixdim(1,2),T1.hdr.dime.pixdim(1,3),T1.hdr.dime.pixdim(1,4));

%optimizer and metric configuration
[optimizer,metric] = imregconfig('multimodal');

%Specify a non-default setting for the InitialRadius property of the optimizer
optimizer.InitialRadius = 0.004;

%Use a similarity transformation to register the images.
movingRegisteredVolume = imregister(moving_T1,Rmoving, fixed_volume,Rfixed, 'rigid', optimizer, metric);

for i=10:50
	w = waitforbuttonpress;
	figure
	imshowpair(movingRegisteredVolume(:,:,i), fixed_volume(:,:,i));
	title('Axial slice of registered volume.')
end

% % view registered volume
% helperVolumeRegistration(fixed_volume,movingRegisteredVolume);

%3-D Geometric Transformation That Aligns Moving With Fixed
% geomtform = imregtform(moving_T1,Rmoving, fixed_volume,Rfixed, 'rigid', optimizer, metric);

% %determine the transformed location of the center of the moving image in the world coordinate system
% centerXWorld = mean(Rmoving.XWorldLimits);
% centerYWorld = mean(Rmoving.YWorldLimits);
% centerZWorld = mean(Rmoving.ZWorldLimits);
% [xWorld,yWorld,zWorld] = transformPointsForward(geomtform,centerXWorld,centerYWorld,centerZWorld);

% % use the worldToSubscript method of Rfixed to determine the element of the fixed volume that aligns with the center of the moving volume
% [r,c,p] = worldToSubscript(Rfixed,xWorld,yWorld,zWorld);

% movingRegisteredVolume1 = imwarp(moving_T1,Rmoving,geomtform,'bicubic','OutputView',Rfixed);

% figure
% imshowpair(movingRegisteredVolume1(:,:,40), fixed_volume(:,:,40));
% title('Axial slice of registered volume.');

% %check difference of two methods
% a=movingRegisteredVolume(:,:,40) - movingRegisteredVolume1(:,:,40);
% figure, imagesc(a);
