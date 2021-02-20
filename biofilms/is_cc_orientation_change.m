function [is_shift,current_cc_orientation,orientation_diff] = is_cc_orientation_change(current_comp,min_ecc,min_orientation_change,max_z_per_cc,current_cc_orientation,current_cc_num_of_z,comp_orientation,comp_eccentricity)
    
    is_shift = 0;
    
    orientation_diff = abs( max(comp_orientation,current_cc_orientation) - min(comp_orientation,current_cc_orientation) );
    if orientation_diff > 100
        orientation_diff = 180 - orientation_diff;
    end
    
    if comp_eccentricity >= min_ecc 
        if orientation_diff >= min_orientation_change
            is_shift = 1; %output
        end
    else current_cc_num_of_z <= max_z_per_cc
           is_shift = 0;
        
    end
    
    if is_shift == 1
        current_cc_orientation = comp_orientation;
    end

end