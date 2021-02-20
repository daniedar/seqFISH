function [cell_by_gene_frame] = generate_cell_by_gene_frame(regprops_stats,im_seg_stack_lab,fov)


    cell_by_gene_frame = struct();
    num_of_cells = size(regprops_stats,1);


    cell_by_gene_frame.cell_id = [1:num_of_cells]';
    cell_by_gene_frame.fov = repelem(fov,num_of_cells)';

    num_of_z_list = zeros(num_of_cells,1);
    min_z_list = zeros(num_of_cells,1);
    for i = 1:num_of_cells
        voxels = regprops_stats.VoxelList(i);
        voxels = voxels{1};
        num_of_z_list(i) = length(unique(voxels(:,3)));
        min_z_list(i) = min(unique(voxels(:,3)));

    end
    
    centroids = regprops_stats.Centroid;
    centroid_mat_ind = sub2ind(size(im_seg_stack_lab),int16(centroids(:,2)),int16(centroids(:,1)),int16(centroids(:,3))); %convert to mat index

    
    cell_by_gene_frame.num_of_z = num_of_z_list;
    cell_by_gene_frame.min_z_list = min_z_list;
    
    cell_by_gene_frame.centroid_ind = centroid_mat_ind; %position in a 3D matrix to get xyz 
    
    cell_by_gene_frame.centroid_z = int16(centroids(:,3)); %position in a 3D matrix to get xyz 

    
    cell_by_gene_frame.Volume = regprops_stats.Volume;

    cell_by_gene_frame.total_rRNA = repelem(0,num_of_cells)';
    cell_by_gene_frame.total_mRNA = repelem(0,num_of_cells)';

    cell_by_gene_frame = struct2table(cell_by_gene_frame);


end