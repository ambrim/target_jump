function velocity = RadialVel(hand_x,hand_y,hand_time)
%will do vel over whole range, input should be range of values wanted
x_vel = abs(diff(hand_x/1000)./diff(hand_time)); %cm/s time is in sec, dist is cm/1000            
y_vel = abs(diff(hand_y/1000)./diff(hand_time));
velocity = sqrt(x_vel.^2 + y_vel.^2);
end

