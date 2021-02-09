function [] =  classify_cells(main_path,fov_demult_table,fov_list,thresholds_db)

    RO_specific_sig2bg = thresholds_db.RO_specific_sig2bg;
    RO_specific_sig2ref = thresholds_db.RO_specific_sig2ref;

    demult_hybs_table = readtable(sprintf('%s\\demultiplexing\\sample_sheets\\demult_hybs_table.txt',main_path));
    demult_16S_table = demult_hybs_table(strcmp(demult_hybs_table.desc,'demult'),:);
    
    for fov_id = 1:length(fov_list)
        fov = fov_list(fov_id);
        fprintf('Demultiplexing and classification (fov %i)\n',fov)

        try 
            cell_class_db = readtable(sprintf('%s\\demultiplexing\\raw_demult_values\\fov_%i_demult_table.txt',main_path,fov));
            num_of_cells = size(cell_class_db,1);
        catch
            fprintf('could not load %s , existing...\n',sprintf('%s\\demultiplexing\\raw_demult_values\\fov_%i_demult_table.txt',main_path,fov))
            break
        end

        for ro_idx = 1:size(demult_16S_table,1)
           ro_num = char(demult_16S_table.ro_num(ro_idx));
           ro_channel = char(demult_16S_table.channel(ro_idx));
           cell_class_db.(sprintf('%s_%s_decision',ro_num,ro_channel)) = zeros(num_of_cells,1);
        end

        classification_db = {};
        for cell_idx = 1:num_of_cells 
            code_book = struct('sum',0);
            
            if rem(cell_idx,1000) == 0
                fprintf(' - cell_idx = %i\n',cell_idx)
            end
            
            for ro_idx = 1:size(demult_16S_table,1)
                ro_num = char(demult_16S_table.ro_num(ro_idx));
                ro_channel = char(demult_16S_table.channel(ro_idx));
                ch_key = sprintf('%s_%s',ro_num,ro_channel);
                bg_key = sprintf('bg_%s',ro_channel);
                ref_key = sprintf('ref_%s',ro_channel);

                code_book.(ro_num) = 0;

                bg_signal = cell_class_db.(bg_key)(cell_idx);
                ref_signal = cell_class_db.(ref_key)(cell_idx);
                ch_signal = cell_class_db.(ch_key)(cell_idx);

                %calculate the ratios
                signal_to_bg_ratio = ch_signal / bg_signal;
                signal_to_ref_ratio = ch_signal / ref_signal;
                cell_class_db.(sprintf('%s_%s_sig2bg',ro_num,ro_channel))(cell_idx) = signal_to_bg_ratio;
                cell_class_db.(sprintf('%s_%s_sig2ref',ro_num,ro_channel))(cell_idx) = signal_to_ref_ratio;

            end 

            cell_ro_intensities = struct();
            cell_ro_intensities.R112 = cell_class_db.R112_A647(cell_idx);
            cell_ro_intensities.R115 = cell_class_db.R115_A647(cell_idx);
            cell_ro_intensities.R113 = cell_class_db.R113_A488(cell_idx);
            cell_ro_intensities.R116 = cell_class_db.R116_A488(cell_idx);
            cell_ro_intensities.R114 = cell_class_db.R114_cy3B(cell_idx);
            cell_ro_intensities.R117 = cell_class_db.R117_cy3B(cell_idx);

            for ro_idx = 1:size(demult_16S_table,1)
                ro_num = char(demult_16S_table.ro_num(ro_idx));
                ro_channel = char(demult_16S_table.channel(ro_idx));
                sig2bg = cell_class_db.(sprintf('%s_%s_sig2bg',ro_num,ro_channel))(cell_idx);
                sig2ref = cell_class_db.(sprintf('%s_%s_sig2ref',ro_num,ro_channel))(cell_idx);            
                [is_ok] = compare_ro_intensities_in_cell(cell_ro_intensities,ro_num,thresholds_db);

                if sig2bg >= RO_specific_sig2bg.(ro_num) & sig2ref >= RO_specific_sig2ref.(ro_num) & is_ok == 1
                    cell_class_db.(sprintf('%s_%s_decision',ro_num,ro_channel))(cell_idx) = 1;
                    code_book.(ro_num) = 1;
                    code_book.sum = code_book.sum + 1; 
                end
            end
            
            c = 0; 
            if code_book.sum == 0
                classification_db{end+1} = 'undetermined';
                c = 1;
            elseif code_book.sum == 1
                classification_db{end+1} = 'only_one_signal';

            elseif code_book.sum > 2
                classification_db{end+1} = 'false_signal';
                c = 1;
            else 
                class_sample_name = 'false_signal';
                for i = 1:(size(fov_demult_table,1))
                    sample_name = char(fov_demult_table.Sample_name(i));
                    barcode_1 = char(fov_demult_table.barcode_1(i));
                    barcode_2 = char(fov_demult_table.barcode_2(i));

                    sample_code_score = code_book.(barcode_1) + code_book.(barcode_2); 
                    if sample_code_score == 2 
                        class_sample_name = sample_name;
                    end
                end
                classification_db{end+1} = class_sample_name; %save result

            end
        end
        fprintf(' - done -> writing results to file\n')
        cell_class_db.classification = classification_db';
        output_path = sprintf('%s\\demultiplexing\\fov_%i_demult_table.decision.txt',main_path,fov);
        fprintf('writing table to %s\n\n',output_path)
        writetable(cell_class_db,output_path,'Delimiter','\t')
    end
end