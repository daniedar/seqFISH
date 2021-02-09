function [is_curved,curve_factor,area] = is_blob_curved(seg_bw,curve_cutoff)
    
    is_curved = 0;
    seg_bw_skel = bwmorph(seg_bw,'skel');
    skel_s = regionprops(seg_bw_skel,{'PixelList','Area','MajorAxisLength'});
    area = skel_s.Area;
    
    x = skel_s.PixelList(:,1);
    y = skel_s.PixelList(:,2);
    distance = 0; %maximal distance between two pixels in this blob (end to end)
    for idx1 = 1:length(y)
        x1 = x(idx1);
        y1 = y(idx1);
        for idx2 = 1:length(x)
            x2 = x(idx2);
            y2 = y(idx2);
            d = sqrt((x1-x2).^2+(y1-y2).^2);
            if d > distance
                distance = d;
            end
        end
    end

    curve_factor = distance/skel_s.MajorAxisLength; % curve_factor = 1 is a straight line
    
    if curve_factor < curve_cutoff & area > 11
        is_curved = 1;
    end
    
end