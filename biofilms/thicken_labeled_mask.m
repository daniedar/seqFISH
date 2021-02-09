function [output_mask] = thicken_labeled_mask(im_seg_stack_lab,ncomp_3d_min,thicken_by)

    %thicken a labled mask while keeping the same labeling...
    num_of_z = size(im_seg_stack_lab,3);
    
    output_mask = zeros(size(im_seg_stack_lab));
    
    for z = 1:num_of_z
        fprintf(' - z = %i\n',z)
        im_seg_stack_lab_z = im_seg_stack_lab(:,:,z);        
        mask_z_bw = logical(im_seg_stack_lab_z);
        try
            %re-label in 2D + change the pixel 
            bw_lab =  bwlabeln(im_seg_stack_lab_z,4); %2D labeled comps
            bw_lab_thick =  bwlabeln(bwmorph(mask_z_bw, 'thick',thicken_by),4); %2D labeled comps - thick
            
            s = regionprops(bw_lab_thick,{'PixelIdxList'});
            NumObjects = size(s,1);
            
            bw_lab_thick_output = zeros(size(bw_lab_thick));
            for thick_label = 1:NumObjects
               
                object_pixel_idx = s(thick_label).PixelIdxList;
                labels_in_original_im = im_seg_stack_lab_z(object_pixel_idx);
                labels_in_original_im = unique(labels_in_original_im(labels_in_original_im > 0));
                
                if length(labels_in_original_im) > 1
                    fprintf('multi labels for %i...\n',thick_label)
                    continue
                end
                bw_lab_thick_output(object_pixel_idx) = labels_in_original_im; 
            end
            output_mask(:,:,z) = bw_lab_thick_output;

        end
    end % z
       
end
   