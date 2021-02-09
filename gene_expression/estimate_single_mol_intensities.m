function [single_mol_normalization_factors_table] = estimate_single_mol_intensities(main_path,fov_list,min_spots_for_fit,max_spots_per_total_cells)
     
    sample_sheet = readtable(sprintf('%s\\experiment_sample_sheet.txt',main_path)); %Experiment sample sheet
    single_mol_res_db = cell2table(cell(0,7),'VariableNames',{'gene_name','hyb','RO','channel','num_of_probes','num_of_spots','mean_signal_per_probe'});
    
    try         
        demultiplexing_results_statistics = readtable(sprintf('%s\\demultiplexing\\demultiplexing_results_statistics.txt',main_path)); %Experiment sample sheet
        total_num_of_cells = 0;
        for fov = fov_list 
            fov_demult_sub_table = demultiplexing_results_statistics(demultiplexing_results_statistics.fov ==fov,:);
            total_num_of_cells = total_num_of_cells  + fov_demult_sub_table.classified;
        end
    catch
        fprintf('could not open %s, exiting...\n\n',sprintf('%s\\demultiplexing\\demultiplexing_results_statistics.txt'))
    end

    %run over genes
    for hyb_idx = 1:size(sample_sheet,1)
        t = sample_sheet(hyb_idx,:);
        hyb = t.hyb;
        readout = t.readout;
        channel = char(t.channel);
        gene_name = char(t.gene_name);
        num_of_probes = t.num_of_probes;
        
        %get merged spots data per gene
        try
            fprintf('Gene = %s\n',gene_name)
            gene_spot_data = table2array(readtable(sprintf('%s\\count_data\\merged_spots_data\\%s_merged_spots.txt',main_path,gene_name)));
        catch
            fprintf('could not open %s, exiting...\n\n',sprintf('%s\\spots_data\\fov_%i_hyb_%i_%s_spots.txt',main_path,fov,hyb,gene_name))
        end
        
        num_of_spots = size(gene_spot_data,1);
        spots_per_cell = num_of_spots/total_num_of_cells;
        
        if num_of_spots >= min_spots_for_fit & spots_per_cell <= max_spots_per_total_cells % enough spots but not too much (low exp) ???
            options = statset('MaxIter',1000);
            f = fitgmdist(gene_spot_data(:,14),3,'RegularizationValue',0.1,'Options',options); %fit a 3 gaussian mixture model
            single_mRNA_intensity = min(f.mu);
            single_probe_int = single_mRNA_intensity / num_of_probes; %differences between genes in number of probes
            gene_res = {gene_name,hyb,readout,channel,num_of_probes,num_of_spots,single_probe_int};
            single_mol_res_db = [single_mol_res_db;gene_res];
        end
    end %genes
    
    %calculate the mean/median tables
    fprintf('Calculating normalization factors per channel...\n')
    channels = {'A647','A488','cy3B'};
    single_mol_res_db_sorted = [];
    channel_specific_estimates_table = [];
    for ch_idx = 1:length(channels)
       ch = char(channels(ch_idx));      
       ch_single_mol_res = single_mol_res_db(strcmp(single_mol_res_db.channel,ch),:);
       single_mol_res_db_sorted = [single_mol_res_db_sorted;ch_single_mol_res];
       
       total_genes_for_estimate = size(ch_single_mol_res,1);
       mean_ch = mean(ch_single_mol_res.mean_signal_per_probe);
       median_ch = median(ch_single_mol_res.mean_signal_per_probe);
       stdv_ch= std(ch_single_mol_res.mean_signal_per_probe);
       cv_ch = stdv_ch / mean_ch * 100;   
       
       line = {ch,mean_ch,median_ch,stdv_ch,cv_ch,total_genes_for_estimate};
       channel_specific_estimates_table = [channel_specific_estimates_table; line];
    end
    single_mol_normalization_factors_table = cell2table(channel_specific_estimates_table,'VariableNames', {'channel' 'mean' 'median' 'stdv','CV','num_of_genes'});
    
    %write raw fitting results
    single_mol_output_path = sprintf('%s\\count_data\\single_mol_intensity_estimates_raw.txt',main_path);
    writetable(single_mol_res_db_sorted,single_mol_output_path,'Delimiter','\t') %writing full table to file 
    
    %write normalization estiamtes 
    single_mol_normalization_output_path = sprintf('%s\\count_data\\single_mol_normalization_factors.txt',main_path);
    writetable(single_mol_normalization_factors_table,single_mol_normalization_output_path,'Delimiter','\t') %writing full table to file 
    
end