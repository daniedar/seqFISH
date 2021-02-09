function [im_seg_stack] = load_and_process_ilastik_mask(main_path,fov,thin_mask_by,min_object_size)
   
    fprintf('Collecting simple segmentation Ilastik mask for fov %i\n',fov)
    simple_segmentation = h5read(sprintf('%s\\segmentation\\fov_%i\\fov_%i_simple_segmentation.h5',main_path,fov,fov),sprintf('/fov_%i',fov));

    num_of_rows = size(simple_segmentation,3);
    num_of_columns = size(simple_segmentation,2);
    num_of_z = size(simple_segmentation,4);

    im_seg_stack = zeros(num_of_rows,num_of_columns,num_of_z);

    for z = 1:num_of_z
        fprintf(' - z = %i\n',z)
        bw = transpose(squeeze(simple_segmentation(1,:,:,z)));
        bw_original = bw == 1;
        bw = bw_original;

        % thin only the large segments only 
        bw_large_objects = bwareaopen(bw, 30);
        bw_smaller_objects = logical(logical(bw - bw_large_objects));
        
        bw_large_objects_thin = bwmorph(bw_large_objects, 'thin',thin_mask_by);
        
        bw_thin = logical(bw_large_objects_thin + bw_smaller_objects);

        bw_thin_fill = imfill(bw_thin, 'holes'); %fill in holes

        bw_thin_fill_spur = bwmorph(bw_thin_fill, 'spur'); 
        bw_thin_fill_spur_diag = bwmorph(bw_thin_fill_spur, 'diag');
        
        %save
        im_seg_stack(:,:,z) = bwareafilt(bw_thin_fill_spur_diag,[min_object_size 150]);
    end %z



end

        