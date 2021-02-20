%% Collect masks, cleanup and filter + label 3D connected components (CCs)

% Parameters
main_path = "";
thin_mask_by = 1;
min_object_size = 4;

% Fields to analyze
fov_list = [5];

% Collect
all_cell_masks_db = struct();
for fov = fov_list
    fprintf('Extracting initial mask fov = %i\n',fov)
    
    [im_seg_stack] = load_and_process_ilastik_mask(main_path,fov,thin_mask_by,min_object_size); %get the results simple_seg from Ilastik and process it
    [projected_mask] = prepare_projected_mask(main_path,fov,im_seg_stack); % prepare a global mask for mapping spots
    
    [im_seg_stack_lab,ncomp_3d_min] = bwlabeln(im_seg_stack,conndef(3,'minimal'));
    fov_regprops_stats = regionprops3(im_seg_stack_lab,{'Volume','VoxelList','Centroid'});
    
    fprintf('Found %i CCs in fov = %i\n',ncomp_3d_min,fov)    
    all_cell_masks_db.(sprintf('fov_%i',fov)) = struct();
    all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_lab = im_seg_stack_lab;
    all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_stats = fov_regprops_stats;
    
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_mask_initial.mat',main_path,fov,fov),'im_seg_stack');
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_mask_initial_lab.mat',main_path,fov,fov),'im_seg_stack_lab');
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_mask_initial_regprops_stats.mat',main_path,fov,fov),'fov_regprops_stats');
end

%% Identidy over-connected CCs and view their locations

% Parameters
main_path = "";
max_num_of_z = 5; %CCs with more than this many will be examined; z-slice = 0.5um

options = struct();
options.min_blob_size = 8; 
options.min_ecc = 0.65; 
options.min_orientation_change = 20; 
options.max_area = 25; 
options.curve_cutoff = 0.78; 
options.max_z_per_cc = 7; 
options.min_num_of_z = 2; 

fov_list = 0:6;
for fov = fov_list
    
    %get data 
    im_seg_stack_lab = all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_lab;
    fov_regprops_stats = all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_stats;
    
    % Extract CC (cells) information into a table
    [cell_by_gene_frame] = generate_cell_by_gene_frame(fov_regprops_stats,im_seg_stack_lab,fov);

    % find objects within reasonable z-num range
    good_comps =  cell_by_gene_frame.cell_id(cell_by_gene_frame.num_of_z <= max_num_of_z);
    over_connected_comps =  cell_by_gene_frame.cell_id(cell_by_gene_frame.num_of_z > max_num_of_z);
    
    %divide components into different labeled masks
    ok_connected_im_seg = im_seg_stack_lab;
    over_connected_im_seg = im_seg_stack_lab;
    over_connected_im_seg(~ismember(over_connected_im_seg,over_connected_comps)) = 0;
    ok_connected_im_seg(ismember(over_connected_im_seg,over_connected_comps)) = 0;
    
    fprintf('Fov=%i\n - Identified %i good CCs and %i over-connected CCs\n',fov,length(good_comps),length(over_connected_comps))
    
    %save masks
    all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_lab_ok_connected = ok_connected_im_seg;
    all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_lab_overconnected = over_connected_im_seg;
    all_cell_masks_db.(sprintf('fov_%i',fov)).initial_cell_by_gene_frame = cell_by_gene_frame;
    
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_mask_initial_lab_ok_connected.mat',main_path,fov,fov),'ok_connected_im_seg');
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_mask_initial_lab_overconnected.mat',main_path,fov,fov),'over_connected_im_seg');

end


%% Do the disconnecting on all of the data and save

% Parameters 
main_path = "";
thicken_by = 4;

fov_list = 0:6;
for fov = fov_list
    fprintf('\n\nDisconnecting components Fov = %i\n\n',fov)
    cell_by_gene_frame = all_cell_masks_db.(sprintf('fov_%i',fov)).initial_cell_by_gene_frame;

    ok_connected_im_seg = all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_lab_ok_connected;
    over_connected_im_seg = all_cell_masks_db.(sprintf('fov_%i',fov)).im_seg_stack_lab_overconnected;
    
    [disc_im_seg_lab,total_recovered,unslaveged_components_list] = disconnect_components(options,ok_connected_im_seg,over_connected_im_seg,cell_by_gene_frame);
    number_of_cells = length(unique(disc_im_seg_lab(:))) - 1;
    fov_regprops_stats = regionprops3(disc_im_seg_lab,{'Volume','VoxelList','VoxelIdxList','Centroid'}); %get the spatial info
    
    %save masks and regprop stats    
    all_cell_masks_db.(sprintf('fov_%i',fov)).disc_im_seg_lab = disc_im_seg_lab; %new labeled 
    all_cell_masks_db.(sprintf('fov_%i',fov)).disc_im_seg_lab_stats = fov_regprops_stats; %new lab stats
    
    %thicken masks
    [disc_im_seg_lab_thick] = thicken_labeled_mask(disc_im_seg_lab,number_of_cells,thicken_by);
    fov_regprops_stats_thick = regionprops3(disc_im_seg_lab_thick,{'Volume','VoxelList','VoxelIdxList','Centroid'});

    %get the unsalvaged CCs
    unslaveged_components = over_connected_im_seg;
    unslaveged_components(~ismember(unslaveged_components,unslaveged_components_list)) = 0;
    
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_disc_im_seg_lab.mat',main_path,fov,fov),'disc_im_seg_lab'); %save
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_disc_im_seg_lab_thick.mat',main_path,fov,fov),'disc_im_seg_lab_thick'); %save
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_lab_mask_unsalvaged.mat',main_path,fov,fov),'unslaveged_components'); %save

    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_disc_im_seg_lab_region_stats.mat',main_path,fov,fov),'fov_regprops_stats'); %save
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_disc_im_seg_lab_thick_region_stats.mat',main_path,fov,fov),'fov_regprops_stats_thick'); %save
    
    [cell_by_gene_frame] = generate_cell_by_gene_frame(fov_regprops_stats,disc_im_seg_lab,fov); %get component info 
    
    fov_mask_db = all_cell_masks_db.(sprintf('fov_%i',fov));
    
    save(sprintf('%s\\segmentation\\fov_%i\\fov_%i_all_cell_masks_db.mat',main_path,fov,fov),'fov_mask_db'); %save

end
