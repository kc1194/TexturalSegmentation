function loaded_img = extract_img(img_file)
	extension = img_file(find(img_file=='.',1,'last'):end);
	if (strcmp(extension,'.nii'))
		Image1 = load_untouch_nii(img_file);
		Image1_data_comp = Image1.img;
		loaded_img = Image1_data_comp(:,:,32,1);
	elseif (strcmp(extension,'.jpg'))
		imread(img_file);
		loaded_img = ans;
	elseif (strcmp(extension,'.mat'))
		load(img_file);
		loaded_img = original_img_all(:,:,32);
	else
		disp('Error');
	end
	loaded_img = double(loaded_img);
end