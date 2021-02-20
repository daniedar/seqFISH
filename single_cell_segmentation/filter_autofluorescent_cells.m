function [output_mask] = filter_autofluorescent_cells(bg_im,bg_im_crop,cell_mask,percentile) 
	
	fprintf(' - thickening mask for calcs\n')
    cell_mask_thick = logical(Thicken_Mask(cell_mask,3));

	fprintf(' - labeling mask\n')
    cell_mask_lab = bwlabel(cell_mask_thick);
    cell_mask_filt = zeros(size(cell_mask));

    %measure mean intensity per cell and remove outliers
    percentile_cutoff = prctile(bg_im(:),percentile); %calculate the cutoff
    max_num_of_high_bg_pixels = 0;
	pixel_int_threshold = 500;
    avg_bg_int_list = zeros(1,max(max(cell_mask_lab))); %
    for cell_idx = 1:max(max(cell_mask_lab))
		if rem(cell_idx,1000) == 0
			disp(cell_idx)
		end
		
		specific_cell_mask = cell_mask_lab == cell_idx;
        number_of_high_bg_pixels = sum(bg_im_crop(specific_cell_mask) >= pixel_int_threshold); 
        cell_mean_bg = mean(bg_im_crop(specific_cell_mask));
        if cell_mean_bg <= 230 & number_of_high_bg_pixels <= max_num_of_high_bg_pixels
           cell_mask_filt = cell_mask_filt + specific_cell_mask;
        end
    end
	fprintf(' -intersect thick and thin masks to get output\n')
	output_mask = cell_mask_filt & cell_mask; 
end
