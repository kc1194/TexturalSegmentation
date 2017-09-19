function index = dice(mask,gold)
	inter_img = mask & gold;
	index = 2*sum(inter_img(:))/(sum(mask(:))+sum(gold(:)));
end