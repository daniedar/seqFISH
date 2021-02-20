%% Measure the distances between cells and save neigbor ID and distances (um) for references

main_path = "D:\imaging_data\081520_10h_biofilm\Analysis\";
fov_list =[0 1 2 3 4 5 6];

max_neighbors_in_output_matrix = 2500; %save up to 2500 closest neighborhoods

for fov = fov_list
    
    [fov_cell_distance_matrix] = calculate_cell_distance_matrix(main_path,fov);
    
    [neighbor_distances,neighbor_cell_numbers] = mink(fov_cell_distance_matrix,max_neighbors_in_output_matrix + 1);
    neighbor_distances = neighbor_distances';
    neighbor_cell_numbers = neighbor_cell_numbers';
    neighbor_distances = neighbor_distances(:,2:end); %remove self
    neighbor_cell_numbers = neighbor_cell_numbers(:,2:end); 
   
    %write fov indepdent tables
    distance_output_name = sprintf('%s\\spatial_analysis\\fov_%i\\fov_%i_cell_distance_matrix_closest_%i_neighbors_distances.mat',main_path,fov,fov,max_neighbors_in_output_matrix);
    cell_id_output_name = sprintf('%s\\spatial_analysis\\fov_%i\\fov_%i_cell_distance_matrix_closest_%i_neighbors_cell_ids.mat',main_path,fov,fov,max_neighbors_in_output_matrix);
    save(distance_output_name,'neighbor_distances');
    save(cell_id_output_name,'neighbor_cell_numbers');
    
end

%% Local gene enrichment

main_path = "";
fov_list =[0:6];
max_neighbors_in_output_matrix = 2500;
gene_name_list = {'gpD','rpoA'};

max_size = 300;
go_by = 5;

for gene_idx = 1:length(gene_name_list)
    gene_name = gene_name_list(gene_idx);
    gene_name = gene_name{1};
    sample_sheet = readtable(sprintf('%s\\experiment_sample_sheet.txt',main_path)); %all gene list
    sample_sheet = sample_sheet(ismember(sample_sheet.gene_name,{gene_name}),:);

    neighborhood_window_size_list = [20];
    max_num_of_neighbors_list = 0:go_by:max_size;
    max_num_of_neighbors_list = max_num_of_neighbors_list(2:end);
    percentile_list = [99.5];

    all_output_matrix = zeros(length(max_num_of_neighbors_list),length(neighborhood_window_size_list));
    for neighborhood_size_idx = 1:length(neighborhood_window_size_list)
        combo_idx = 1;
        neighborhood_window_size = neighborhood_window_size_list(neighborhood_size_idx);
        for max_num_of_neighbors_idx = 1:length(max_num_of_neighbors_list)
            max_num_of_neighbors = max_num_of_neighbors_list(max_num_of_neighbors_idx);
            for percentile_cutoff_idx = 1:length(percentile_list)
                percentile_cutoff = percentile_list(percentile_cutoff_idx);
                gene_parameter_db = struct();
                for i = 1:length(sample_sheet.gene_name)
                    fprintf(' - gene %i of combo %i\n',i,combo_idx)
                    gene_name = sample_sheet.gene_name{i};

                    %identify center cells (expressing levels > cut_off)
                    [center_cells_id_list,exp_cutoff] = get_cells_expressing_gene_of_interest(main_path,fov_list,gene_name,percentile_cutoff);
                    gene_parameter_db.(gene_name) = struct('num_of_center_cells',size(center_cells_id_list,1),'exp_cutoff',exp_cutoff);

                    % get the list of cells that are neighbors 
                    [gene_neighborhood_mean_exp_fold_change] = get_neighbor_fold_change_expression(main_path,fov_list,max_neighbors_in_output_matrix,neighborhood_window_size,max_num_of_neighbors,center_cells_id_list,sample_sheet);

                    %save the fold change accordingly in the output matrix
                    all_output_matrix(combo_idx,neighborhood_size_idx) = 2^gene_neighborhood_mean_exp_fold_change; 
                end
                combo_idx = combo_idx + 1;
            end
        end
    end

    columns = {};
    for i = 1:length(neighborhood_window_size_list)
        columns{end+1} = sprintf('window_%0.1f',neighborhood_window_size_list(i));
    end

    rows = {};
    for i = 1:length(max_num_of_neighbors_list)
        rows{end+1} = sprintf('%i',max_num_of_neighbors_list(i));
    end

    gene_neighborhood_mean_exp_fold_change_table = array2table(all_output_matrix,'VariableNames',columns,'RowNames',rows);
    outpath = sprintf('%s\\spatial_analysis\\%s_local_fold_change_percentile_%0.1f_max_%i_by_%i.txt',main_path,gene_name,percentile_list,max_size,go_by);

    writetable(gene_neighborhood_mean_exp_fold_change_table,outpath,'Delimiter','\t','WriteRowNames',true) %writing full table to file 
end