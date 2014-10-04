function [ traj_x, traj_y ] = calc_trial_trajectory( dx, dy, start_x, start_y )

traj_x = zeros( size(dx,2), 1 );
traj_y = zeros( size(dy,2), 1 );

if (nargin == 4)
    traj_x(1) = start_x;
    traj_y(1) = start_y;
else
    traj_x(1) = 0;
    traj_y(1) = 0;
end

i=2;
while i <= size(dx,2)
    traj_x(i) = traj_x(i-1) + dx(i);
    traj_y(i) = traj_y(i-1) + dy(i);
    i=i+1;
end

end

