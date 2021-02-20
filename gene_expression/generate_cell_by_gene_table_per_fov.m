function [cell_by_gene_table] = generate_cell_by_gene_table_per_fov(main_path,fov,spot_int_and_size_threshold_db)
    
    path_to_cell_by_gene_dir = sprintf('%s\\count_data\\cell_by_gene',main_path);
    if ~exist(path_to_cell_by_gene_dir, 'dir')
       mkdir(path_to_cell_by_gene_dir)
    end

    sample_sheet = readtable(sprintf('%s\\experiment_sample_sheet.txt',main_path)); %read sample sheet
    cell_class_table = readtable(sprintf('%s\\demultiplexing\\fov_%i_demult_table.decision.txt',main_path,fov));
    cells_to_discard = logical(strcmp(cell_class_table.classification,'only_one_signal') + strcmp(cell_class_table.classification,'undetermined') + strcmp(cell_class_table.classification,'false_signal'));
    cells_to_keep = logical(ones(length(cells_to_discard),1) - cells_to_discard);
    cell_class_table = cell_class_table(cells_to_keep,:); %filter out the unclassified cells and false positives
    
    fprintf(' - getting single-mol normalization factor table..\n')
    single_mol_norm_db = readtable(sprintf('%s\\count_data\\single_mol_normalization_factors.txt',main_path)); %single -mol norms
    
    A488_factor_correction = 1.5;
    
    cell_by_gene_table = table();
    cell_by_gene_table.cell_id = cell_class_table.cell_id;
    cell_by_gene_table.fov = cell_class_table.fov;
    cell_by_gene_table.classification = cell_class_table.classification;
    cell_by_gene_table.MajorAxisLength = cell_class_table.MajorAxisLength;
    cell_by_gene_table.MinorAxisLength = cell_class_table.MinorAxisLength;
    cell_by_gene_table.total_DNA_median = cell_class_table.total_DNA_median;
    cell_by_gene_table.total_rRNA = cell_class_table.ref_cy3B;
    cell_by_gene_table.total_mRNA = zeros(length(cell_class_table.cell_id),1); % to fill out after normalization
    
    num_of_non_gene_columns = size(cell_by_gene_table,2);
    
    for hyb_idx = 1:size(sample_sheet,1)

        t = sample_sheet(hyb_idx,:);
        hyb = t.hyb;
        readout = t.readout;
        channel = char(t.channel);
        gene_name = char(t.gene_name);
        num_of_probes = t.num_of_probes;

        fprintf('gene %s #%i\n',gene_name,hyb_idx)
        try
            fprintf(' - loading spot data for fov = %i\n',fov)
            min_peak_intensity = spot_int_and_size_threshold_db.(channel).min_peak_intensity;
            max_fit_size = spot_int_and_size_threshold_db.(channel).max_fit_size;
            min_fit_size = spot_int_and_size_threshold_db.(channel).min_fit_size;
            [gene_spot_data,gene_spots_discarded] = load_and_filter_gene_spots_data(main_path,fov,gene_name,min_peak_intensity,max_fit_size,min_fit_size);
        catch
            fprintf('could not open %s, exiting...\n\n',sprintf('%s\\spots_data\\fov_%i_hyb_%i_%s_spots.txt',main_path,fov,hyb,gene_name))
        end

        if ~ismember(gene_name,fieldnames(cell_by_gene_table))
           cell_by_gene_table.(gene_name) = zeros(size(cell_by_gene_table,1),1);
        end
        
        fprintf(' - counting spot intensities...\n')
        for spot_idx = 1:size(gene_spot_data,1)
           cell_idx = gene_spot_data(spot_idx,12); 
           spot_intensity = gene_spot_data(spot_idx,14);
           cell_by_gene_table{cell_by_gene_table.cell_id == cell_idx,gene_name}  = cell_by_gene_table{cell_by_gene_table.cell_id == cell_idx,gene_name} + spot_intensity;
        end

        norm_table = single_mol_norm_db(strcmp(single_mol_norm_db.channel,channel),:);
        norm_factor = norm_table.mean * num_of_probes;
        if strcmp(channel,'A488')
            norm_factor = norm_factor/A488_factor_correction;
        end
        cell_by_gene_table.(gene_name) = double(int16(cell_by_gene_table.(gene_name) / norm_factor));
        cell_by_gene_table.total_mRNA = cell_by_gene_table.total_mRNA + cell_by_gene_table.(gene_name);
    end

    fov_cell_by_gene_path = sprintf('%s\\fov_%i_cell_by_gene_table.txt',path_to_cell_by_gene_dir,fov);
    fprintf('Saving table:\n  - %s...\n\n',fov_cell_by_gene_path)
    writetable(cell_by_gene_table,fov_cell_by_gene_path,'Delimiter','\t') %writing full table to file 

    
end