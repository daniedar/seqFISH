function distance = corrected_euclidean_dist_func(XI,XJ)

z_pixels = 5; %correct for the z = ~5 pixels

diff = XI - XJ;
diff(3) = diff(3) * z_pixels;
sqdx = diff.^2;
sqdx_sum = sum(sqdx,2);

distance = sqrt(sqdx_sum);

end