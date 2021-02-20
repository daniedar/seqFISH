function [distance_matrix] = calculate_cell_distance_matrix(main_path,fov)
    
    fprintf('Calculating distance matrix for fov = %i\n',fov)
    microns_per_pixel = 0.1;
    
    %get cell coordinates for fov
    fov_disc_im_seg_lab_stats = struct2array( load(sprintf('%s\\segmentation\\fov_%i\\fov_%i_disc_im_seg_lab_region_stats.mat',main_path,fov,fov)) );
    fov_cell_centroid_coordinates = fov_disc_im_seg_lab_stats.Centroid;
    
    %calculate distances
    distance_matrix = squareform(pdist(fov_cell_centroid_coordinates,@corrected_euclidean_dist_func)) * microns_per_pixel; %distance matrix
    
end