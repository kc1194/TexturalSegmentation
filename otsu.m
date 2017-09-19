function otsu_mask = otsu(img)
	[level,EM] = graythresh(img);
	otsu_mask = im2bw(img,level);
	