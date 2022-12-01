function [dist_max, dist_last] = RadialDist(hand_x,hand_y,start_x,start_y,pix_per_mm)
reach_dist_fromStart = sqrt(((hand_x-start_x).^2)+((hand_y-start_y).^2));
dist_max = max(reach_dist_fromStart)/100*pix_per_mm;
dist_last = reach_dist_fromStart(end)/100*pix_per_mm;
end
