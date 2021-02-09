function [gene_spot_data_filtered] = load_and_filter_gene_spots_data(main_path,fov,gene,min_peak_intensity,max_fit_size,min_fit_size)

    sample_sheet = readtable(sprintf('%s\\experiment_sample_sheet.txt',main_path)); %read sample sheet
    sample_sheet = sample_sheet(strcmp(sample_sheet.gene_name,gene),:);
    hyb = sample_sheet.hyb;
    channel = char(sample_sheet.channel);
    num_of_probes = sample_sheet.num_of_probes;

    try
        gene_spot_data = table2array(readtable(sprintf('%s\\spots_data\\fov_%i_hyb_%i_%s_spots.txt',main_path,fov,hyb,gene)));
    catch
        fprintf('could not open %s, exiting...\n\n',sprintf('%s\\spots_data\\fov_%i_hyb_%i_%s_spots.txt',main_path,fov,hyb,gene_name))
    end
    
    %filter
    fit_size_list = gene_spot_data(:,15) .* gene_spot_data(:,16);
    fit_axis_ratio_list = gene_spot_data(:,16) ./ gene_spot_data(:,15);    
    gene_spot_data_output = gene_spot_data(fit_size_list <= max_fit_size,:);
    fit_size_list = gene_spot_data_output(:,15) .* gene_spot_data_output(:,16); %only large ones
    gene_spot_data_output = gene_spot_data_output(fit_size_list >= min_fit_size,:);
    gene_spot_data_filtered = gene_spot_data_output(gene_spot_data_output(:,1) >= min_peak_intensity ,:);
end