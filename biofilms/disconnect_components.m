function [output_im_seg,total_recovered,components_we_could_not_repair] = disconnect_components(options,ok_connected_im_seg,over_connected_im_seg,cell_by_gene_frame)
   
    %general params
    total_recovered = 0;
    
    %re-label the filtered matrix 
    [output_im_seg,num_of_ok_connected] = bwlabeln(ok_connected_im_seg,conndef(3,'minimal'));
    
    %list of bad comp ids
    list_of_components_to_repair = unique(over_connected_im_seg(:));
    list_of_components_to_repair = list_of_components_to_repair(list_of_components_to_repair > 0);
    
    components_we_could_not_repair = [];
    
    %disconnect algorithm params
    min_blob_size = options.min_blob_size; % tougher filter in over-connected (remove small segs)
    min_ecc = options.min_ecc; %how line like the seg has to be for orientation analysis
    min_orientation_change = options.min_orientation_change; %min change to disconnect
    max_area = options.max_area; %max num of pixels in skeletal seg (detects overly long ones)
    curve_cutoff = options.curve_cutoff; %threhold for detecting curved segs
    max_z_per_cc = options.max_z_per_cc; % allow up to this many z in component if the oreintation is right (ex, standing cell)
    min_num_of_z = options.min_num_of_z; % discard any component that appears in less than 2 z
    
    %do the disconnecting per each 3d-CC
    for bad_comp_idx = 1:length(list_of_components_to_repair)
        bad_comp_num = list_of_components_to_repair(bad_comp_idx);
                
        fprintf('Component idx= %i, num=%i\n',bad_comp_idx,bad_comp_num)
        comp_specific_mask = logical(over_connected_im_seg == bad_comp_num);        

        %find the z-range for the bad component
        cell_by_gene_frame_sp = cell_by_gene_frame(cell_by_gene_frame.cell_id==bad_comp_num,:); %find 
        num_of_z_sections_in_comp = cell_by_gene_frame_sp.num_of_z(1);
        min_z = cell_by_gene_frame_sp.min_z_list(1);
        max_z = min_z + num_of_z_sections_in_comp - 1;
        
        %compnent dbs
        comp_index_list = [0]; 
        cc_database = struct();
        disconnected_mat = zeros(size(comp_specific_mask)); % keeps the new comps as we gos

        % run over Z-sections and re-connect 
        for z = min_z:max_z

            seg_bw_full = comp_specific_mask(:,:,z);

            %remove smaller blobs
            seg_bw_full = bwareaopen(seg_bw_full,min_blob_size); %remove tiny blobs

            [seg_bw_full_lab,n_comps] = bwlabeln(seg_bw_full,conndef(2,'minimal')); %label componenets 

            for comp_idx = 1:n_comps

                seg_bw_comp = seg_bw_full_lab == comp_idx; %CC specific mask
                comp_bw_s = regionprops(seg_bw_comp,{'PixelIdxList','Centroid'}); %to fill in the disconnected mat

                centroid_xy = comp_bw_s.Centroid;
                size_of_cropped_img = 50;
                xmin = int16(centroid_xy(1)) - size_of_cropped_img/2;
                ymin = int16(centroid_xy(2)) - size_of_cropped_img/2;

                %for simplicity use a crop
                seg_bw = imcrop(seg_bw_comp,[xmin ymin size_of_cropped_img size_of_cropped_img]);
                bw_s = regionprops(seg_bw,{'PixelList','PixelIdxList','Orientation','Eccentricity','Area','MajorAxisLength','MinorAxisLength','Centroid'});

                %find out whether blob is curved or not
                [is_curved,curve_factor,area] = is_blob_curved(seg_bw,curve_cutoff);

                if is_curved == 1 | area > max_area 
                    continue
                end

                %is this a new comp or does it overlap one in the last z-section?
                seg_bw_full_z_minus = disconnected_mat(:,:,max(z-1,1));

                intersection_with_previous_z = unique(seg_bw_full_z_minus((seg_bw_comp & seg_bw_full_z_minus)));
                intersection_with_previous_z = intersection_with_previous_z(intersection_with_previous_z > 0);
                overlapping_comp = comp_index_list(find(ismember(comp_index_list,intersection_with_previous_z)));

                if length(overlapping_comp) > 2 %throw away anything that connects more than two segments...
                    continue
                end

                is_new_comp = 0;

                if ~isempty( overlapping_comp ) % if connected, set the componenet number to that one
                    current_comp = overlapping_comp; %until proven otherwise, we assume this will be the same component
                else % no overlap -> set a new comp
                    is_new_comp = 1;
                    current_comp = max(comp_index_list) + 1;
                    comp_index_list(end + 1) = current_comp; 
                end

                current_comp_to_save = [];
                comp_orientation = bw_s.Orientation;
                comp_eccentricity = bw_s.Eccentricity;
                %if new component
                if is_new_comp == 1 %only if there is no overlap with any other CCs
                    cc_database.(sprintf('cc_%i',current_comp)) = struct('orientation',comp_orientation,'num_of_z',1);
                    current_comp_to_save(end+1) = current_comp;
                else %there is an overlap - is there a shift?
                    for current_comp_idx = 1:length(current_comp) %cases where there are multi comps
                        current_comp_i = current_comp(current_comp_idx);

                        current_cc_orientation = cc_database.(sprintf('cc_%i',current_comp_i)).orientation;
                        current_cc_num_of_z = cc_database.(sprintf('cc_%i',current_comp_i)).num_of_z;
                        [is_shift,current_cc_orientation,orientation_diff] = is_cc_orientation_change(current_comp_i,min_ecc,min_orientation_change,max_z_per_cc,current_cc_orientation,current_cc_num_of_z,comp_orientation,comp_eccentricity);

                        if is_shift == 0
                            cc_database.(sprintf('cc_%i',current_comp_i)).num_of_z = cc_database.(sprintf('cc_%i',current_comp_i)).num_of_z + 1;
                            current_comp_to_save(end+1) = current_comp_i;
                        else %significant shift, this is possibly a new component
                            is_new_comp = is_new_comp + 1;
                            current_comp_i = max(comp_index_list) + 1;
                            current_comp_to_save(end+1) = current_comp_i;
                            comp_index_list(end + 1) = current_comp_i;
                            cc_database.(sprintf('cc_%i',current_comp_i)) = struct('orientation',comp_orientation,'num_of_z',1);
                        end
                    end 
                end

                if is_new_comp == 2 & length(current_comp_to_save) == 2%overlapped 2 but belongs to none
                    new_comp = min(current_comp_to_save);
                    for cc_num = current_comp_to_save
                        cc_database = rmfield(cc_database,sprintf('cc_%i',cc_num));
                        comp_index_list(comp_index_list==cc_num) = [];
                    end

                    cc_database.(sprintf('cc_%i',new_comp)) = struct('orientation',comp_orientation,'num_of_z',1);
                    current_comp_to_save = new_comp;
                    comp_index_list(end+1) = new_comp;
                end

                if length(current_comp_to_save) > 1 %convert the ones connected to one number, including this one
                    comp_for_all = min(current_comp_to_save);
                    disconnected_mat_z_minus = disconnected_mat(:,:,max(z-1,1));
                    disconnected_mat_z = disconnected_mat(:,:,z);

                    for comp_to_change_idx = 1:length(current_comp_to_save)
                        comp_2_change = current_comp_to_save(comp_to_change_idx);
                        disconnected_mat_z_minus(disconnected_mat_z_minus == comp_2_change) = comp_for_all;
                        disconnected_mat_z_minus(disconnected_mat_z == comp_2_change) = comp_for_all;
                    end

                    disconnected_mat_z(comp_bw_s.PixelIdxList) = comp_for_all; %fill it in

                    disconnected_mat(:,:,max(z-1,1)) = disconnected_mat_z_minus; %save it 
                    disconnected_mat(:,:,z) = disconnected_mat_z; %save it 

                else
                    disconnected_mat_z = disconnected_mat(:,:,z);
                    for current_comp_idx = 1:length(current_comp_to_save)
                        current_comp_i = current_comp_to_save(current_comp_idx);
                        disconnected_mat_z(comp_bw_s.PixelIdxList) = current_comp_i; %fill it in
                    end

                    disconnected_mat(:,:,z) = disconnected_mat_z; %save it 
                end
            end %comps in Z

        end
        
        fnames = fieldnames(cc_database);
        for i = 1:length(fnames)
            cc_fn = fnames{i};
            cc_database.(cc_fn).num_of_z = 0; %zero them 
        end

        %count again
        for z = min_z:max_z
            seg_bw_full = disconnected_mat(:,:,z);
            cc_numbers = unique(seg_bw_full(seg_bw_full > 0));
            for cc_idx = 1:length(cc_numbers)
                cc_num = cc_numbers(cc_idx);
                cc_database.(sprintf('cc_%i',cc_num)).num_of_z = cc_database.(sprintf('cc_%i',cc_num)).num_of_z + 1;

            end
        end

        %filter
        num_of_pass_filter_comp = 0;
        for i = 1:length(fnames)
            cc_fn = fnames{i};
            num_of_z = cc_database.(cc_fn).num_of_z;
            cc_num = split(cc_fn,'_');
            cc_num = str2num(cc_num{2});
            if num_of_z < min_num_of_z | num_of_z > max_z_per_cc 
                disconnected_mat(disconnected_mat == cc_num) = 0;
            else
                num_of_pass_filter_comp = num_of_pass_filter_comp + 1;
            end
        end

        if num_of_pass_filter_comp > 0
            s = regionprops3(disconnected_mat,{'VoxelIdxList'}); %apply to new relabeled matrix
            for i = 1:size(s,1)
               si = s.VoxelIdxList{i};
               if ~isempty( si )
                  new_id = num_of_ok_connected + 1; %get new 
                  output_im_seg(si) = new_id; %use voxel index to add labels to output matrix
                  num_of_ok_connected = new_id; %update total count
                  fprintf('  - adding new label: %i\n',new_id)
                  total_recovered = total_recovered + 1;
               end
            end
        else
            components_we_could_not_repair(end+1) = bad_comp_num;
        end
   
    end 
    
end