function [] = demultiplex_seqFISH_standalone_hpc(main_path,fov,total_bg_hyb_num)

    main_path = char(main_path);
    cd(main_path);
    mkdir demultiplexing;

    output_path = sprintf('%s//demultiplexing//fov_%i_demult_table.txt',main_path,fov);

    fprintf(sprintf('%s//demult_hybs_table.txt\n',main_path))
    demult_hybs_table = readtable(sprintf('%s//demult_hybs_table.txt',main_path)); %read sample sheet;

    fprintf('loading segments\n')

    mask_db_path = sprintf('%s//segmentation//fov_%i_seg//fov_%i_cell_masks_db.mat',main_path,fov,fov);

    cell_mask_db = load(mask_db_path);
    cell_mask = cell_mask_db.cell_mask;
    cell_mask_lab = cell_mask_db.cell_mask_lab;
    cell_mask_thick = cell_mask_db.cell_mask_thick;
    cell_mask_thick_lab = cell_mask_db.cell_mask_thick_lab; 

    num_of_cells = size(unique(cell_mask_thick_lab(cell_mask_thick_lab>0)),1);
    cell_id_list = unique(cell_mask_thick_lab(cell_mask_thick_lab>0));

    fprintf('%i cells in field\n',num_of_cells)

    %% dapi total-DNA im
    disp('loading_DAPI_image')
    total_DNA_path = sprintf('%s//zproj_aligned//fov_%i_hyb_0.maxInt.DAPI.aligned.tif',main_path,fov);
    total_DNA_im = imread(total_DNA_path);

    %% calculate background signal cutoffs from distributions
    disp('loading_backgrounds_images')
    bg_images = struct();
    ch_img_A488_path = sprintf('%s//zproj_aligned//fov_%i_hyb_%i.maxInt.A488.aligned.tif',main_path,fov,total_bg_hyb_num);
    ch_img_A647_path = sprintf('%s//zproj_aligned//fov_%i_hyb_%i.maxInt.A647.aligned.tif',main_path,fov,total_bg_hyb_num);
    ch_img_cy3B_path = sprintf('%s//zproj_aligned//fov_%i_hyb_%i.maxInt.cy3B.aligned.tif',main_path,fov,total_bg_hyb_num);

    bg_images.A647 = imread(ch_img_A647_path);
    bg_images.A488 = imread(ch_img_A488_path);
    bg_images.cy3B = imread(ch_img_cy3B_path);


    %% get raw data for classification (16S hyb images)
    disp('loading_channel_images')
    demult_images = struct();
    for ro_idx = 1:size(demult_hybs_table,1) 
        ro_num = char(demult_hybs_table.ro_num(ro_idx));
        ro_channel = char(demult_hybs_table.channel(ro_idx));
        ro_hyb = demult_hybs_table.hyb(ro_idx);
        im_path = sprintf('%s//zproj_aligned//fov_%i_hyb_%i.maxInt.%s.aligned.tif',main_path,fov,ro_hyb,ro_channel);
        ro_im = imread(im_path);
        demult_images.(ro_num) = ro_im;
    end


    %% measure background and reference channels per cell
    fprintf('Measuring background and reference signals per cell (%i cells)\n',num_of_cells)
    demult_16S_table = demult_hybs_table(strcmp(demult_hybs_table.desc,'demult'),:);
    ref_16S_table = demult_hybs_table(strcmp(demult_hybs_table.desc,'ref'),:);

    cell_class_db = struct('cell_id',transpose(1:num_of_cells),'fov',fov.*ones(num_of_cells,1),'cell_size',zeros(num_of_cells,1),...
                           'total_DNA_median',zeros(num_of_cells,1),...
                           'total_DNA_sum',zeros(num_of_cells,1),...
                           'bg_A647',zeros(num_of_cells,1),...
                           'bg_A488',zeros(num_of_cells,1),...
                           'bg_cy3B',zeros(num_of_cells,1),...
                           'ref_A647',zeros(num_of_cells,1),...
                           'ref_A488',zeros(num_of_cells,1),...
                           'ref_cy3B',zeros(num_of_cells,1));

    % measure background and reference signals per cell
    channel_list = {'A647','A488','cy3B'};

    for ch_idx = 1:size(channel_list,2)
        channel = channel_list{ch_idx}; 

        bg_im = bg_images.(channel);%background image
        ref_im = demult_images.(ref_16S_table(strcmp(ref_16S_table.channel,channel),:).ro_num{1});%ref image

        bg_key = sprintf('bg_%s',channel);%struct keys 
        ref_key = sprintf('ref_%s',channel);

        for cell_idx = 1:num_of_cells
            cell_idx_corrected = cell_id_list(cell_idx);
            cell_specific_mask = cell_mask_lab == cell_idx_corrected;
            cell_class_db.cell_size(cell_idx) = sum(sum(cell_specific_mask)); 

            cell_specific_mask = imerode(cell_specific_mask,strel('disk',1)); 

            cell_class_db.total_DNA_median(cell_idx) = median(total_DNA_im(cell_specific_mask)); 
            cell_class_db.total_DNA_sum(cell_idx) = sum(sum(total_DNA_im(cell_specific_mask)));

            cell_class_db.(bg_key)(cell_idx) = median(bg_im(cell_specific_mask));
            cell_class_db.(ref_key)(cell_idx) = median(ref_im(cell_specific_mask));
        end
    end
    fprintf(' - done\n')


    %% run calculations for individual 16S ROs
    fprintf('Measuring demultiplexing 16S ROs per cell (%i cells)\n',num_of_cells)

    for ro_idx = 1:size(demult_16S_table,1)  
        ro_num = char(demult_16S_table.ro_num(ro_idx));
        ro_channel = char(demult_16S_table.channel(ro_idx));

        ch_key = sprintf('%s_%s',ro_num,ro_channel);
        demult_im = demult_images.(ro_num);
        cell_class_db.(ch_key) = zeros(num_of_cells,1);

        for cell_idx = 1:num_of_cells
            cell_idx_corrected = cell_id_list(cell_idx);
            cell_specific_mask = cell_mask_lab == cell_idx_corrected;    
            cell_specific_mask = imerode(cell_specific_mask,strel('disk',1));
            cell_class_db.(ch_key)(cell_idx) = median(demult_im(cell_specific_mask));
        end
    end
    fprintf(' - done\n')

    %% go over cells and calculate the ratios
    fprintf('Demultiplexing and classification (%i cells)\n',num_of_cells)
    for ro_idx = 1:size(demult_16S_table,1) 
        ro_num = char(demult_16S_table.ro_num(ro_idx));
        ro_channel = char(demult_16S_table.channel(ro_idx));
        cell_class_db.(sprintf('%s_%s_sig2bg',ro_num,ro_channel)) = zeros(num_of_cells,1);
        cell_class_db.(sprintf('%s_ref2bg',ro_channel)) = zeros(num_of_cells,1);
        cell_class_db.(sprintf('%s_%s_sig2ref',ro_num,ro_channel)) = zeros(num_of_cells,1);
    end

    for cell_idx = 1:num_of_cells
        for ro_idx = 1:size(demult_16S_table,1)
            ro_num = char(demult_16S_table.ro_num(ro_idx));
            ro_channel = char(demult_16S_table.channel(ro_idx));
            ch_key = sprintf('%s_%s',ro_num,ro_channel);
            bg_key = sprintf('bg_%s',ro_channel);%struct keys 
            ref_key = sprintf('ref_%s',ro_channel);

            %get cell specific signals in bg, ref and demult 16S channel
            bg_signal = cell_class_db.(bg_key)(cell_idx);
            ref_signal = cell_class_db.(ref_key)(cell_idx);
            ch_signal = cell_class_db.(ch_key)(cell_idx);

            %calculate the ratios and save 
            ref_to_bg_ratio = ref_signal / bg_signal;
            signal_to_bg_ratio = ch_signal / bg_signal;
            signal_to_ref_ratio = ch_signal / ref_signal;

            cell_class_db.(sprintf('%s_%s_sig2bg',ro_num,ro_channel))(cell_idx) = signal_to_bg_ratio;
            cell_class_db.(sprintf('%s_ref2bg',ro_channel))(cell_idx) = ref_to_bg_ratio;
            cell_class_db.(sprintf('%s_%s_sig2ref',ro_num,ro_channel))(cell_idx) = signal_to_ref_ratio;
        end
    end   
    fprintf(' - done\n')

    %% write results file 
    cell_class_db_table = struct2table(cell_class_db); %save as table
    cell_class_db_table.cell_id = cell_id_list;

    %output
    fprintf('writing table to %s\n\n',output_path)
    writetable(cell_class_db_table,output_path,'Delimiter','\t') %writing full table to file 

end