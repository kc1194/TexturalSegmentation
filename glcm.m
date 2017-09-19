function glcm_images = glcm(img,window_size,limits)
	[w h] = size(img);
	upp_lim = ceil(window_size/2);
	low_lim = floor(window_size/2);
	glcm_images = cell([3 1]);
	max_quant = max(max(img));
	for k=1:3
		glcm_images{k} = zeros(w,h);
    end
    iter=0;
	for i=1:w
		for j=1:h
            
            if (mod(iter,10000)==0)
            	disp(iter)
            end
            iter = iter+1;
			if ((i<upp_lim) || (j<upp_lim))
				continue;
			end
			if ((i>(w-low_lim)) || (j>(h-low_lim)))
				continue;
			end
			glcm_window = img(i-low_lim:i+low_lim,j-low_lim:j+low_lim);
			if strcmp(limits,'lin')
				[glcms,si] = graycomatrix(glcm_window,'Offset',[0 1],'Symmetric',true,'GrayLimits',[min(min(glcm_window)) max(max(glcm_window))],'NumLevels',max(max(glcm_window))-min(min(glcm_window))+1);
			else
				[glcms,si] = graycomatrix(glcm_window,'Offset',[0 1],'Symmetric',true,'GrayLimits',[min(min(glcm_window)) max(max(glcm_window))],'NumLevels',8);
			end
			stats = graycoprops(glcms);
			glcm_images{1}(i,j)=stats.Contrast;
            glcm_images{2}(i,j)=stats.Energy;
            glcm_images{3}(i,j)=stats.Homogeneity;
		end
	end
end
