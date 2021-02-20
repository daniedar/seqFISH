function [] = merge_spot_data_per_gene(main_path,fov_list,spot_int_and_size_threshold_db)

    cd(main_path);
    mkdir count_data;
    cd(sprintf('%s\\count_data\\',main_path));
    mkdir merged_spots_data;
    
    sample_sheet = readtable(sprintf('%s\\experiment_sample_sheet.txt',main_path)); %Experiment sample sheet
    
    for hyb_idx = 1:size(sample_sheet,1)
        t = sample_sheet(hyb_idx,:);
        hyb = t.hyb;
        readout = t.readout;
        channel = char(t.channel);
        gene_name = char(t.gene_name);
        num_of_probes = t.num_of_probes;

        disp(gene_name)
        merged_spot_data = [];
        for fov = fov_list
            try
                fprintf(' - fov = %i\n',fov)
                fprintf(' - loading spot data for fov = %i\n',fov)
                min_peak_intensity = spot_int_and_size_threshold_db.(channel).min_peak_intensity;
                max_fit_size = spot_int_and_size_threshold_db.(channel).max_fit_size;
                min_fit_size = spot_int_and_size_threshold_db.(channel).min_fit_size;
                [fov_spot_data,gene_spots_discarded] = load_and_filter_gene_spots_data(main_path,fov,gene_name,min_peak_intensity,max_fit_size,min_fit_size);               
            catch
                fprintf('could not open %s, exiting...\n\n',sprintf('%s\\spots_data\\fov_%i_hyb_%i_%s_spots.txt',main_path,fov,hyb,gene_name))
                
            end     
            fov_spot_data(:,size(fov_spot_data,2)+1) = fov; %adding fov at column 22
            fov_spot_data(:,size(fov_spot_data,2)+1) = readout; %adding RO at column 24
            merged_spot_data = [merged_spot_data;fov_spot_data ];

        end 

        merged_spot_output_path = sprintf('%s\\count_data\\merged_spots_data\\%s_merged_spots.txt',main_path,gene_name);
        writetable(array2table(merged_spot_data),merged_spot_output_path,'Delimiter','\t')

    end 
    fprintf('done merging files...\n\n')
    
end