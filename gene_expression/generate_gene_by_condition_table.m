function [norm_exp_gene_by_condition_table] = generate_gene_by_condition_table(main_path,fov_list,min_expression)

    sample_sheet = readtable(sprintf('%s\\experiment_sample_sheet.txt',main_path)); %read sample sheet
    demult_sample_sheet = readtable(sprintf('%s\\demultiplexing\\sample_sheets\\demult_sample_sheet.txt',main_path));
    demult_sample_sheet.Sample_name = matlab.lang.makeValidName(demult_sample_sheet.Sample_name); 
    
    fprintf(' - merging fov cell by gene tables...\n')
    merged_cell_by_gene_table = [];
    for fov = fov_list
        [fov_cell_by_gene_table] = load_cell_by_gene_table(main_path,fov,'cell_by_gene');
        merged_cell_by_gene_table = [merged_cell_by_gene_table;fov_cell_by_gene_table];
    end
    
    merged_cell_by_gene_table.classification = matlab.lang.makeValidName(merged_cell_by_gene_table.classification);
    
    fprintf(' - measuring gene means per condition...\n')
    mean_gene_by_condition_db = struct();
    for cond_idx = 1:length(demult_sample_sheet.Sample_name)
        condition = char(demult_sample_sheet.Sample_name{cond_idx});
        fprintf('  - %s\n',condition)
        mean_gene_by_condition_db.(condition) = zeros(size(sample_sheet,1),1);
        
        cells_in_conditon = merged_cell_by_gene_table(ismember(merged_cell_by_gene_table.classification,condition),:);
        
        for gene_idx = 1:size(sample_sheet,1)
            gene_name = char(sample_sheet.gene_name{gene_idx});
            gene_mean_expression = mean(cells_in_conditon.(gene_name));
            mean_gene_by_condition_db.(condition)(gene_idx) = gene_mean_expression;
        end
    end
    
    %produce a normalized expression (by max of each gene)
    norm_exp_gene_by_condition_table = struct2table(mean_gene_by_condition_db,'RowNames',sample_sheet.gene_name);
    for gene_idx = 1:size(sample_sheet,1)
        gene_name = char(sample_sheet.gene_name{gene_idx});
        
        max_exp = max(norm_exp_gene_by_condition_table(ismember(norm_exp_gene_by_condition_table.Row,gene_name),:).Variables);
        if max_exp >= min_expression 
            norm_exp_gene_by_condition_table(ismember(norm_exp_gene_by_condition_table.Row,gene_name),:).Variables = norm_exp_gene_by_condition_table(ismember(norm_exp_gene_by_condition_table.Row,gene_name),:).Variables / max_exp;
        else 
            norm_exp_gene_by_condition_table(ismember(norm_exp_gene_by_condition_table.Row,gene_name),:) = [];
        end
        
    end
    
    %raw means tables
    mean_gene_by_condition_table = struct2table(mean_gene_by_condition_db,'RowNames',sample_sheet.gene_name);
    output_path = sprintf('%s\\count_data\\cell_by_gene\\mean_gene_exp_by_condition_table.txt',main_path);
    fprintf('writing results to file:\n - %s\n\n',output_path)
    writetable(mean_gene_by_condition_table,output_path,'Delimiter','\t','WriteRowNames',true) %writing full table to file 
    
    %normalized by max tables
    output_path = sprintf('%s\\count_data\\cell_by_gene\\normalized_gene_exp_by_condition_table.txt',main_path);
    fprintf('writing results to file:\n - %s\n\n',output_path)
    writetable(norm_exp_gene_by_condition_table,output_path,'Delimiter','\t','WriteRowNames',true) %writing full table to file 

end