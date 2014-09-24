%% Plot cummulative

[ traj_x, traj_y ] = calculate_traj(dx_all,dy_all);

figure;
plot(traj_x,traj_y)

%% Calculate velocity distributions

t_diff = diff(t_all);

v_x = dx_all(2:end) ./ t_diff;
v_y = dy_all(2:end) ./ t_diff;

v = sqrt(v_x.^2 + v_y.^2);

% hist(v,10000)
t_z = t_all-t_all(1);
figure;
plot(t_z(2:end),v)
axis tight;
xlabel('Time (s)');

%% 
figure;
hold on
plotyy([1:size(t_diff,2)],t_diff,[1:size(t_diff,2)],dy_all(2:end))

%%
%hist(t_diff,100000)
hist(v,100000)

%% Load all mat files from trials

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

STIM = 15.0;
PRE_STIM = 5.0;

figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);

max_vsize = -1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d = data{trial_idx, j};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_diff = diff(t);

        v_x = dx(2:end) ./ t_diff;
        v_y = dy(2:end) ./ t_diff;

        v = sqrt(v_x.^2 + v_y.^2);
        if( size(v,2) > max_vsize )
            max_vsize = size(v,2);
        end
        %all_data_v(trial_idx, j, 1:size(v,2)) = v;
    end
end

disp(['Max size: ' num2str(max_vsize)]);

clear all_data_v;
all_data = zeros(4*max(trial_type_cnt) * max_vsize,1 );

cur_idx = 1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d = data{trial_idx, j};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_diff = diff(t);

        v_x = dx(2:end) ./ t_diff;
        v_y = dy(2:end) ./ t_diff;

        v = sqrt(v_x.^2 + v_y.^2);
        
        %all_data(cur_idx:(cur_idx+size(v,2))-1) = v;
        all_data(cur_idx:(cur_idx+size(v,2))-1) = t_diff;
        cur_idx = cur_idx + size(v,2);
    end
end

clear all_data_trim;
all_data_trim = all_data(1:cur_idx);
hist(all_data_trim,100000)