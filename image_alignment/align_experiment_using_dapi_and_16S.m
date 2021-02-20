function [] = align_experiment_using_dapi_and_16S(main_path,fov,num_of_hybs,alignment_round,correction_16S_round)

    %correct mechanical shifts between cycles using DAPI (mech_tform)
    %correct channel differences using 16S and DAPI (channel_tform)

    cd(main_path);
    mkdir zproj_aligned;
    mkdir images_aligned;
    mkdir segmentation;
    
    zproj_dir = sprintf('%s%s',main_path,'zproj\');

    fixed_dapi = imread(sprintf('%sfov_%i_hyb_%i.maxInt.DAPI.tif',zproj_dir,fov,alignment_round));
    
    ch_fixed_dapi = imread(sprintf('%sfov_%i_hyb_%i.maxInt.DAPI.tif',zproj_dir,fov,correction_16S_round));
    moving_A488 =  imread(sprintf('%sfov_%i_hyb_%i.maxInt.A488.tif',zproj_dir,fov,correction_16S_round));
    movin_A647 =  imread(sprintf('%sfov_%i_hyb_%i.maxInt.A647.tif',zproj_dir,fov,correction_16S_round));
    movin_cy3B =  imread(sprintf('%sfov_%i_hyb_%i.maxInt.A647.tif',zproj_dir,fov,correction_16S_round));
     
    disp('calc tform channel A488')
    [A488_tform,Rfixed] = calc_tform_channel_to_dapi(ch_fixed_dapi,moving_A488);
    disp('calc tform channel A647')
    [A647_tform,Rfixed] = calc_tform_channel_to_dapi(ch_fixed_dapi,movin_A647);
    disp('calc tform channel cy3B')
    [cy3B_tform,Rfixed] = calc_tform_channel_to_dapi(ch_fixed_dapi,movin_cy3B);
    
    ch_tforms = struct();
    ch_tforms.A488 = A488_tform;
    ch_tforms.A647 = A647_tform;
    ch_tforms.cy3B = cy3B_tform;
    ch_tforms.A647 = A647_tform;
    ch_tforms.cy3B = cy3B_tform;
     
    zproj_output_dir = sprintf('%s%s',main_path,'zproj_aligned\');
    align_image_and_save_dapi_16S(fixed_dapi,ch_tforms,main_path,fov,alignment_round);
    for hyb = 0:num_of_hybs-1
        disp(hyb)
        disp(sprintf('aligning hyb %i',hyb))
        align_image_and_save_dapi_16S(fixed_dapi,ch_tforms,main_path,fov,hyb); 
     end
end