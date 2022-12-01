function [dist_time] = DistTime(hand_x,hand_y,hand_time,start_x,start_y,pix_per_mm,dist)
reach_dist_fromStart = (sqrt(((hand_x-start_x).^2)+((hand_y-start_y).^2)))/100*pix_per_mm;
wheredist=find(reach_dist_fromStart>dist);
if ~isempty(wheredist)
    dist_time=hand_time(wheredist(1));
else
    dist_time = [];
end
end
