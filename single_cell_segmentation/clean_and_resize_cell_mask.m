function [clean_img] = clean_and_resize_cell_mask(main_path,fov,cell_mask,min_cell_pixel_count)

    cell_mask = bwareaopen(cell_mask, min_cell_pixel_count);
    dapi_im = imread(sprintf('%s/zproj_aligned/fov_%i_hyb_%i.maxInt.DAPI.aligned.tif',main_path,fov,0)); %read dapi

    %clear border cells - any object touching the border is eliminated
    %adjust border in cases where image registration caused a shift
    im_size = size(dapi_im);
    col_sum = sum(dapi_im);
    col_sum_zeros = col_sum == 0;
    cols_on_right = sum(col_sum_zeros(im_size(2)/2:im_size));
    cols_on_left = sum(col_sum_zeros(1:im_size(2)/2));

    row_sum = sum(dapi_im,2);
    row_sum_zeros = row_sum == 0;
    rows_on_top = sum(row_sum_zeros(1:im_size(1)/2));
    rows_on_bottom = sum(row_sum_zeros(im_size(1)/2:im_size));

    %clear border cells - any object touching the border is eliminated 
    cell_mask_bord = cell_mask;
    cell_mask_bord(1:(3+rows_on_top),:) = 1;
    cell_mask_bord(end-(2+rows_on_bottom):end,:) = 1;
    cell_mask_bord(:,1:(3+cols_on_left)) = 1;
    cell_mask_bord(:,end-(2+cols_on_right):end) = 1;

    cell_mask_bord_clean = imclearborder(cell_mask_bord);

    % clear irregular shapes 
    BW = cell_mask_bord_clean;
    stats = struct2table(regionprops(BW,{'Area','Solidity','PixelIdxList','MinorAxisLength','MajorAxisLength'}));

    idx = stats.Solidity < 0.80 | stats.Area < min_cell_pixel_count | stats.MinorAxisLength < 5.5 | stats.MinorAxisLength > 10.5  | stats.MajorAxisLength > 35;
    for kk = find(idx)'
      BW(stats.PixelIdxList{kk}) = false;
    end

    clean_img = BW;
end
