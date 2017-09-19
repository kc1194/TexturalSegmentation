function glrlm_images = glrlm(img,window_size,limits)
	[w h] = size(img)
	max_quant =max(max(img));
	glrlm_images = cell([11 1]);
	for i=1:11
		glrlm_images{i} =(zeros(w,h));
	end

	iteration = 0;
	for i = 1:w
		for j = 1:h
			if ((i<13) | (j<13))
				continue;
			end
			if ((i>w-12) | (j>h-12))
				continue;
			end
			iteration = iteration + 1;
			if (mod(iteration,10000)==0)
				disp(iteration)
			end
			glrlm_window = img(i-12:i+12,j-12:j+12);
			if strcmp(limits,'lin')
				[glrlms,si] = grayrlmatrix(glrlm_window,'Offset',[1], 'Numlevels',max(max(glrlm_window))-min(min(glrlm_window))+1,'G',[]);
			else
				[glrlms,si] = grayrlmatrix(glrlm_window,'Offset',[1], 'Numlevels',8,'G',[]);
			end	

			stats = grayrlprops(glrlms);
			for k=1:11
				glrlm_images{k}(i,j) = stats(k);
			end
		end
	end
end


