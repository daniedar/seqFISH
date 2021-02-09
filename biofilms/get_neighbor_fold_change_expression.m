function [output_vector] = get_neighbor_fold_change_expression(main_path,fov_list,max_neighbors_in_output_matrix,neighborhood_window_size,max_num_of_neighbors,center_cells_id_list,sample_sheet)
    
    gene_list = sample_sheet.gene_name;

    merged_cell_by_gene = [];
    merged_neighborhood_cell_by_gene = [];
    for fov = fov_list
        fprintf(' - calculating local cell neighborhoods: fov %i\n',fov)
        
        fov_center_cells = center_cells_id_list(center_cells_id_list.fov == fov,:); %Get the relevant center cells for this fov
        
        %get pre-calcualted neighborhoods of all cells (cell id and distances)
        neighbor_distances = load(sprintf('%s\\spatial_analysis\\fov_%i\\fov_%i_cell_distance_matrix_closest_%i_neighbors_distances.mat',main_path,fov,fov,max_neighbors_in_output_matrix));
        neighbor_cell_numbers = load(sprintf('%s\\spatial_analysis\\fov_%i\\fov_%i_cell_distance_matrix_closest_%i_neighbors_cell_ids.mat',main_path,fov,fov,max_neighbors_in_output_matrix));
        neighbor_distances = neighbor_distances.neighbor_distances;
        neighbor_cell_numbers = neighbor_cell_numbers.neighbor_cell_numbers;
        
        % slice by the center_cells ids (look only at those)
        neighbor_distances = neighbor_distances(fov_center_cells.cell_id,:); %select the center cell rows. These do not include the actual center cell
        neighbor_cell_numbers = neighbor_cell_numbers(fov_center_cells.cell_id,:);%select the center cell rows
        
        %get up to the maximal number of neighbors
        neighbor_distances = neighbor_distances(:,1:max_num_of_neighbors); % get only the X closest neighboring cells
        neighbor_cell_numbers = neighbor_cell_numbers(:,1:max_num_of_neighbors);
        
        %clear cells that are too far away
        neighbor_distances_above_cutoff = neighbor_distances > neighborhood_window_size; %identify cells too far away and erase their numbers from the list
        neighbor_distances(neighbor_distances_above_cutoff) = 0; %erasing
        neighbor_cell_numbers(neighbor_distances_above_cutoff) = 0; %erasing
        
        %Get the center cells and get the unique cell_ids for the neighborhoods in one list
        neighbor_cell_numbers = neighbor_cell_numbers(:); %reshape into a single vector
        neighbor_cell_numbers = sort(unique(neighbor_cell_numbers)); %discard duplicates/zeros and sort
               
        center_cells_in_neighborhood = neighbor_cell_numbers(ismember(neighbor_cell_numbers,fov_center_cells.cell_id)); %center cells in use
        center_cell_id_to_remove = fov_center_cells.cell_id;
        center_cell_id_to_remove = center_cell_id_to_remove(~ismember(center_cell_id_to_remove,center_cells_in_neighborhood));

        %Analyze the cell by genes
        fov_cell_by_gene = readtable(sprintf('%s\\count_data\\cell_by_gene\\fov_%i_cell_by_gene_table.txt',main_path,fov)); %all cellss
        
        %slice out the neighborhood only + save to merge
        neighborhood_cell_by_gene = fov_cell_by_gene(ismember(fov_cell_by_gene.cell_id,neighbor_cell_numbers),:); %all neighborhood cells
        merged_neighborhood_cell_by_gene = [merged_neighborhood_cell_by_gene;neighborhood_cell_by_gene]; %save into the merged fov db
        
        % remove the center cells from the total as they were not analyzed
        % in the neighborhoods (unless they are in the neighborhoods
        % themselves)
        fov_cell_by_gene_minus_center = fov_cell_by_gene(~ismember(fov_cell_by_gene.cell_id,center_cell_id_to_remove),:); %remove useed center cells         
        merged_cell_by_gene = [merged_cell_by_gene; fov_cell_by_gene_minus_center]; %save for the "all cell analysis minus center cells"

    end
    
    output_vector = zeros(length(gene_list),1);
    
    % calculate the fold change difference in table per gene
    for i = 1:length(gene_list)
        gene_name = gene_list{i};
        all_cells_mean = mean(merged_cell_by_gene.(gene_name));
        neighborhood_cells_mean = mean(merged_neighborhood_cell_by_gene.(gene_name));
%         ratio = double(neighborhood_cells_mean +0.01) / double(all_cells_mean+0.01);
        ratio_log2 = log2(double(neighborhood_cells_mean +0.01) / double(all_cells_mean+0.01));
        output_vector(i,1) = ratio_log2;
    end
    
end