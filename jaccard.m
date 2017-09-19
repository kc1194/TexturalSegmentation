function index = jaccard(mask,gold)
	inter_img = mask & gold;
	union_img = mask | gold;
	index = sum(inter_img(:))/sum(union_img(:));
end