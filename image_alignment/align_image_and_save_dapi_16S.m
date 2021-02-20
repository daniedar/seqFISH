function [] = align_image_and_save_dapi_16S(fixed_dapi,ch_tforms,exp_dir,fov,hyb)
    
    zproj_dir = sprintf('%s%s',exp_dir,'zproj\');
    zproj_output_dir = sprintf('%s%s',exp_dir,'zproj_aligned\');

    hyb = num2str(hyb);
    moving_dapi = imread(sprintf('%sfov_%i_hyb_%s.maxInt.DAPI.tif',zproj_dir,fov,hyb));
    Rfixed = imref2d(size(fixed_dapi));
    
    fprintf('corr\n')
    tform_corr = imregcorr(moving_dapi,fixed_dapi,'translation');
    moving_dapi = imwarp(moving_dapi,tform_corr,'OutputView',Rfixed);
    
    fprintf('reg\n')
    tform_reg = image_registration(fixed_dapi,moving_dapi);
    movingReg_dapi = imwarp(moving_dapi,tform_reg,'OutputView',Rfixed);

    imwrite(convert_to_8bits(movingReg_dapi),sprintf('%sfov_%i_hyb_%s.maxInt.DAPI.aligned.tif',zproj_output_dir,fov,hyb));
    
    ch_list = ["A647","cy3B","A488"]; 
    for i = 1:length(ch_list) 
        ch_name = ch_list(1,i);
        moving_ch_im = imread(sprintf('%sfov_%i_hyb_%s.maxInt.%s.tif',zproj_dir,fov,hyb,ch_name));
        movingReg_ch_cor = imwarp(moving_ch_im,tform_corr,'OutputView',Rfixed);
        movingReg_ch = imwarp(movingReg_ch_cor,tform_reg,'OutputView',Rfixed);      
        ch_tform = ch_tforms.(ch_name);
        movingReg_ch = imwarp(movingReg_ch,ch_tform,'OutputView',Rfixed);
        imwrite(movingReg_ch,sprintf('%sfov_%i_hyb_%s.maxInt.%s.aligned.tif',zproj_output_dir,fov,hyb,ch_name)); %write
    end
    
    zstacks_dir = sprintf('%s%s',exp_dir,'zstacks\');
    zstacks_output_dir = sprintf('%s%s',exp_dir,'images_aligned\');
    for i = 1:length(ch_list) 
        ch_name = ch_list(1,i);
        zstack_path = sprintf('%sfov_%i_hyb_%s.%s.zstack.tif',zstacks_dir,fov,hyb,ch_name); 
        num_zslices = numel(imfinfo(zstack_path));
        for j = 1:num_zslices
            moving_zstack_slice = imread(zstack_path,j);
            moving_zstack_slice_reg_cor = imwarp(moving_zstack_slice,tform_corr,'OutputView',Rfixed); %correct using corr
            movingReg_slice = imwarp(moving_zstack_slice_reg_cor,tform_reg,'OutputView',Rfixed); %correct using reg
            ch_tform = ch_tforms.(ch_name); %second correction for abbreration
            movingReg_slice = imwarp(movingReg_slice,ch_tform,'OutputView',Rfixed);
            slice_out_path = sprintf('%sfov_%i_hyb_%s.%s.z%i.tif',zstacks_output_dir,fov,hyb,ch_name,j);
            imwrite(movingReg_slice,slice_out_path);
        end
    end

end