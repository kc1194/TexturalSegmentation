function generated =  generate(task,modality,window_size,saver,limits,save_name,save_path)
	clear all;

	%comments kaun daalega?? :P
	%suppress all warnings
	warning('off','all');
	% Saving  : task limits Modality Patient no 

	%set default values
	if ~exist('task', 'var') || isempty(task); task='glcm'; end

	if strcmp(task,'glcm')
		if ~exist('saver', 'var') || isempty(saver); saver=true; end
		if ~exist('save_name', 'var') || isempty(save_name); save_name='Linear'; end
		if ~exist('save_path', 'var') || isempty(save_path); save_path='/Users/kabir/Desktop/Research/MTP/Figures/'; end
		if ~exist('window_size', 'var') || isempty(window_size); window_size=7; end
		if ~exist('limits', 'var') || isempty(limits); limits='lin'; end
		if ~exist('modality', 'var') || isempty(modality); modality=1; end
	elseif strcmp(task,'glrlm')
		if ~exist('saver', 'var') || isempty(saver); saver=true; end
		if ~exist('save_name', 'var') || isempty(save_name); save_name='Linear'; end
		if ~exist('save_path', 'var') || isempty(save_path); save_path='/Users/kabir/Desktop/Research/MTP/Figures/'; end
		if ~exist('window_size', 'var') || isempty(window_size); window_size=25; end
		if ~exist('limits', 'var') || isempty(limits); limits='lin'; end
		if ~exist('modality', 'var') || isempty(modality); modality=1; end
	else
		disp('Wrong Task!');
		exit
	end

	%loading and writing directories
	data_dir = '/Users/kabir/Desktop/Research/MTP/MTP_data/';
	files = {}; 
	images = {};
	generated = [];
	contents = dir(data_dir);
	index=0;
	for i=1:size(contents,1)
		if strcmp(contents(i).name(1),'.')
			continue
		end
		index = index+1;
		files{index} = (strcat(data_dir,contents(i).name));
	end
	
	files
	

	% if strcmp(task,'glcm')
	% 	generated = glcm(img,window_size,limits);
	% else
	% 	generated = glrlm(img,window_size,limits);
	% end

% Saving format
% feature/modality_window_featureno_patientno

	if (saver)
		% for i=1:size(files,2)
		% 	disp('Displaying file no');
		% 	disp(i);
		% 	cur_img = load_img(files{i});
		% 	cur_gen = glcm(cur_img,7,limits);
		% 	for j=1:size(cur_gen,1)
		% 		disp('Displaying feature no');
		% 		disp(j);
		% 		fil_to_write = strcat(save_path,'finalGLCM/',save_name,'_',num2str(i),'_',num2str(7),'_',num2str(j),'_1.jpg');
		% 		% if (exist(dir_to_write,'dir')~=7)
		% 		% 	mkdir(dir_to_write);
		% 		% end
		% 		imwrite(mat2gray(cur_gen{j}),fil_to_write);
		% 	end
		% end
		for i=1:size(files,2)
			disp('Displaying file no');
			disp(i);
			cur_img = extract_img(files{i});
			cur_gen = glrlm(cur_img,25,limits);
			for j=1:size(cur_gen,1)
				disp('Displaying feature no');
				disp(j);
				fil_to_write = strcat(save_path,'finalGLRLM/',save_name,'_',num2str(i),'_',num2str(25),'_',num2str(j),'_1.jpg');
				% if (exist(dir_to_write,'dir')~=7)
				% 	mkdir(dir_to_write);
				% end
				imwrite(mat2gray(cur_gen{j}),fil_to_write);
			end
		end
	end
	exit
end

