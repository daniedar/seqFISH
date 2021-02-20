function [cell_masks_output,size_filtered_mask] = get_superSegger_cell_mask(main_path,seg_dir,fov,bg_hyb)

    main_path = char(main_path);
    seg_dir = char(seg_dir);
    min_cell_pixel_count = 30;

    %% get supersegger data
    all_fovs_dir = dir([sprintf('%s',seg_dir),'/xy*']);
    path_to_chosen_fov_dir = sprintf('%s%s%s%s',all_fovs_dir(1).folder,'/',all_fovs_dir(1).name,'/');

    seg_data_dir = dir([sprintf('%s%s',path_to_chosen_fov_dir,'seg/'),'*_seg.mat']); %get 
    seg_data = load([sprintf('%s%s',path_to_chosen_fov_dir,'seg/'),seg_data_dir(1).name]);

    cell_mask = seg_data.mask_cell; %get the cell mask
    phase = seg_data.phaseNorm; %get phase image (DAPI)
    cell_mask_unfiltered = cell_mask;

    %% Clean and discard odd shapes - some poor segmentations
    disp('Clean and discard small and odd shapes')
    [cleaned_cell_mask] = clean_and_resize_cell_mask(main_path,fov,cell_mask,min_cell_pixel_count);

    %% filter no DNA cells and small cells
    %some cases of regular shaped and sized cells with not enough dapi.
    disp('filter out small segments and those with low DAPI stain'); 
    min_cell_size_cutoff = 30;
    min_dapi_signal = 350; %a little arbitrary

    dapi_im = imread(sprintf('%s/zproj_aligned/fov_%i_hyb_%i.maxInt.DAPI.aligned.tif',main_path,fov,0)); %read dapi

    bw_lab_dapi = bwlabel(cleaned_cell_mask); %label cells
    dapi_filtered_cell_mask = zeros(size(cleaned_cell_mask));
    for i = 1:max(max(bw_lab_dapi))
        if rem(i,1000) == 0
            fprintf(' - %i\n',i)
        end
        specific_cell_mask = bw_lab_dapi == i; %filter mask for cell i
        cell_size = sum(sum(specific_cell_mask));
        mean_dapi = mean(dapi_im(specific_cell_mask)); %mean dapi in this mask
        num_of_dapi_positive_pixels = sum(dapi_im(specific_cell_mask) >= min_dapi_signal);

        if cell_size >= min_cell_size_cutoff & mean_dapi >= min_dapi_signal & num_of_dapi_positive_pixels/cell_size >= 0.8
            dapi_filtered_cell_mask = dapi_filtered_cell_mask + specific_cell_mask;
        end
    end


    %% remove autofluoresnct cells 
    disp('Remove highly autofluor cells')
    fluor_percentile_cutoff = 99.5;
    bg_im = imread(sprintf('%s/zproj_aligned/fov_%i_hyb_%i.maxInt.A488.aligned.tif',main_path,fov,bg_hyb)); %read dapi
    [autofluo_filtered_mask] = filter_autofluorescent_cells(bg_im,bg_im,dapi_filtered_cell_mask,fluor_percentile_cutoff); 

    %% final mask calculations, save and plot
    disp('Thicken cell mask and label');
    
    output_cell_mask = logical(autofluo_filtered_mask); %bw mask
    output_cell_mask_thick = logical(Thicken_Mask(output_cell_mask,3)); %thicken the mask without overlapping cells

    output_cell_mask_thick_lab = bwlabel(output_cell_mask_thick);
    output_cell_mask_lab = output_cell_mask_thick_lab .* imcomplement(logical(output_cell_mask_thick_lab) - output_cell_mask); %save exact numbering as in thick but keep normal segment size (labels comes out different otherwise)

    disp('removing small_cells cut outs after thick to thin mask adjustments')
    size_filtered_output_cell_mask_lab = zeros(size(output_cell_mask_lab));
    for i = 1:max(max(output_cell_mask_lab))
        specific_cell_mask = output_cell_mask_lab == i; %filter mask for cell i
        cell_size = sum(sum(specific_cell_mask));

        if cell_size >= min_cell_size_cutoff
            size_filtered_output_cell_mask_lab = size_filtered_output_cell_mask_lab + specific_cell_mask .* i;
        end
    end
    output_cell_mask_lab = size_filtered_output_cell_mask_lab;
    output_cell_mask = logical(output_cell_mask_lab);

    %% save options and plot results
    disp('Saving masks');
    total_num_of_segments = max(max(output_cell_mask_lab));
    fprintf('total number of segments: cell_mask = %i\n',max(max(bwlabel(output_cell_mask))));
    fprintf('total number of segments: cell_mask_lab = %i\n',max(max(bwlabel(output_cell_mask_lab))));
    fprintf('total number of segments: cell_mask_thick = %i\n',max(max(bwlabel(output_cell_mask_thick))));
    fprintf('total number of segments: cell_mask_thick_lab = %i\n',max(max(bwlabel(output_cell_mask_thick_lab))));
    cell_masks_output = struct('cell_mask',output_cell_mask,...
                               'cell_mask_lab',output_cell_mask_lab,...
                               'cell_mask_thick',output_cell_mask_thick,...
                               'cell_mask_thick_lab',output_cell_mask_thick_lab,...
                               'cell_mask_unfiltered',cell_mask_unfiltered,...
                               'phase',phase);

    save(sprintf('%s/fov_%i_cell_masks_db.mat',seg_dir,fov),'-struct','cell_masks_output'); %a = load(sprintf('%sgene_spot_data.mat',main_path)); 
end
