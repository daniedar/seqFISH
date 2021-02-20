main_path = "";

fov_demult_sheet = readtable(sprintf('%s\\demultiplexing\\sample_sheets\\demult_sample_sheet.txt',main_path)); %
fov_list = from_fov:1:to_fov;

%define the 16S signal thresholds for classification. done empirically
RO_specific_sig2bg = struct('R112',2.5,'R115',2.5,'R113',1.3,'R116',1.8,'R114',2,'R117',2); %fold- comparing to background
RO_specific_sig2ref = struct('R112',1.5,'R115',1.5,'R113',0.1,'R116',0.4,'R114',0.2,'R117',0.2); %fold-comparing to ref 16S
RO_within_channel_rel_int = struct('R112',1.75,'R115',2,'R113',1.30,'R116',2,'R114',1.8,'R117',1.8); % relative signal within same laser channel
RO_between_channel_rel_int = struct('A647',struct('A488',1,'Cy3B',1),'Cy3B',struct('A488',1)); %rel sign between different laser channels

thresholds_db = struct('RO_specific_sig2bg',RO_specific_sig2bg,'RO_specific_sig2ref',RO_specific_sig2ref, ... 
    'RO_within_channel_rel_int',RO_within_channel_rel_int,'RO_between_channel_rel_int',RO_between_channel_rel_int);

%running the classification
classify_cells(main_path,fov_demult_sheet,fov_list,thresholds_db);
